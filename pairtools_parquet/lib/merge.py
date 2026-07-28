"""Merging several .pairs/.parquet files into one.

`pairtools merge` shells out to `sort --merge` over the decompressed bodies;
here the inputs are UNION ALL-ed in DuckDB and ordered by the same keys, so the
inputs may be a mix of formats and the output may be any of them.

The ordering matches `pairtools merge`: chrom1, chrom2, pos1, pos2, pair_type,
compared bytewise -- the same keys `pairtools sort` uses, which is what makes a
merged file still sorted.
"""

import glob
import os
import tempfile

from pairtools.lib import headerops

from . import arrowio, duckdb_utils, duckdbio
from .arrowio import PairsWriter
from .sort import DEFAULT_SORT_KEYS, quote_identifier

UTIL_NAME = "pairtools_parquet_merge"

#: Inputs merged in one pass, matching `pairtools merge`. Beyond this the merge
#: is staged through temporary files.
DEFAULT_MAX_NMERGE = 8


def expand_paths(patterns):
    """Expand shell wildcards, as `pairtools merge` does, preserving order.

    A pattern that matches nothing is kept as-is so that a missing file is
    reported by the reader, rather than silently vanishing from the merge.
    """
    paths = []
    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        paths.extend(matches if matches else [pattern])
    return paths


def merged_header(headers, util_name=UTIL_NAME):
    """Combine input headers the way `pairtools merge` does.

    Deliberately does not re-mark the header as sorted: the inputs already
    carry `#sorted`, merge_headers keeps it, and calling
    `headerops.mark_header_as_sorted` again would rewrite `#chromosomes: a b c`
    as `#chromosomes: : a b c` -- an upstream quirk that `pairtools merge` does
    not trigger, so neither do we.
    """
    if not headerops.all_same_columns(headers):
        raise ValueError("Input pairs cannot contain different columns")

    header = headerops.merge_headers(headers)
    return headerops.append_new_pg(header, ID=util_name, PN=util_name)


def merge_pairs(
    paths,
    output,
    concatenate=False,
    keep_first_header=False,
    max_nmerge=DEFAULT_MAX_NMERGE,
    nproc=8,
    tmpdir=None,
    memory=None,
    compress_program="auto",
    row_group_size=None,
    util_name=UTIL_NAME,
    **kwargs
):
    """Merge `paths` into `output`, maintaining the .pairs sort order.

    Parameters
    ----------
    paths : list of str
        Input paths or wildcards. Formats may be mixed.
    output : str
        Output path; its extension selects the format.
    concatenate : bool
        Concatenate rather than merging sorted inputs, leaving the output
        unsorted. Mirrors `pairtools merge --concatenate`.
    keep_first_header : bool
        Take the first input's header instead of merging all of them.
    max_nmerge : int
        How many inputs to merge in one pass. Beyond this the merge is staged
        through temporary Parquet files, as `pairtools merge` stages through
        temporary text files, so a fan-in of thousands does not need every
        input open at once. Pass 0 for no limit.
    """
    nproc_in = kwargs.get("nproc_in", 3)
    cmd_in = kwargs.get("cmd_in", None)

    paths = expand_paths(paths)
    if not paths:
        raise ValueError("No input paths")

    header_paths = paths[:1] if keep_first_header else paths
    headers = [
        arrowio.read_header(path, nproc_in=nproc_in, cmd_in=cmd_in)
        for path in header_paths
    ]

    if len(paths) == 1:
        # `pairtools merge` copies a lone input through untouched, header and
        # all, rather than recording itself in the @PG chain.
        out_header = headers[0]
    else:
        out_header = merged_header(headers, util_name=util_name)

    path_headers = _headers_for(paths, headers)

    if max_nmerge and len(paths) > max_nmerge:
        # The header is computed once, from the real inputs, and handed to every
        # stage -- so a staged merge writes exactly what a single-pass one would
        # rather than accumulating a @PG record per stage.
        _merge_staged(
            paths, path_headers, output, out_header, max_nmerge,
            concatenate=concatenate, nproc=nproc, tmpdir=tmpdir, memory=memory,
            compress_program=compress_program, row_group_size=row_group_size,
            **kwargs
        )
        return

    _merge_files(
        paths, path_headers, output, out_header,
        concatenate=concatenate, nproc=nproc, tmpdir=tmpdir, memory=memory,
        compress_program=compress_program, row_group_size=row_group_size,
        **kwargs
    )


def _merge_staged(
    paths, path_headers, output, out_header, max_nmerge, tmpdir=None, **kwargs
):
    """Merge in rounds, storing each round's result in a temporary file.

    Parquet for the intermediates: it round-trips the header and the column
    types without reparsing, which is the whole reason the staging is cheap
    here where `pairtools merge` pays for text on every round.
    """
    with tempfile.TemporaryDirectory(dir=tmpdir or None) as stage_dir:
        round_number = 0
        while len(paths) > max_nmerge:
            merged, merged_headers = [], []
            for i in range(0, len(paths), max_nmerge):
                group = paths[i : i + max_nmerge]
                if len(group) == 1:
                    merged.append(group[0])
                    merged_headers.append(path_headers[i])
                    continue
                part = os.path.join(
                    stage_dir, "round{}-{}.parquet".format(round_number, i)
                )
                _merge_files(
                    group,
                    path_headers[i : i + max_nmerge],
                    part,
                    out_header,
                    tmpdir=tmpdir,
                    **kwargs
                )
                merged.append(part)
                merged_headers.append(out_header)
            paths, path_headers = merged, merged_headers
            round_number += 1

        _merge_files(
            paths, path_headers, output, out_header, tmpdir=tmpdir, **kwargs
        )


def _merge_files(
    paths,
    path_headers,
    output,
    out_header,
    concatenate=False,
    nproc=8,
    tmpdir=None,
    memory=None,
    compress_program="auto",
    row_group_size=None,
    **kwargs
):
    """Merge `paths` into `output` under an already-decided header."""
    nproc_in = kwargs.get("nproc_in", 3)
    columns = headerops.extract_column_names(out_header)
    projection = ", ".join(quote_identifier(col) for col in columns)

    def run(use_chrom_enum):
        con = duckdb_utils.setup_duckdb_connection(
            temp_directory=tmpdir or None,
            memory_limit=memory or None,
            enable_progress_bar=False,
            enable_profiling="no_output",
            numb_threads=nproc,
        )
        try:
            chrom_type = (
                duckdbio.declare_chrom_enum(con, out_header)
                if use_chrom_enum
                else None
            )

            # Column names are projected explicitly: UNION ALL matches by
            # position, and all_same_columns only guarantees the same set.
            parts = []
            for path, path_header in zip(paths, path_headers):
                parts.append(
                    "SELECT {} FROM {}".format(
                        projection,
                        duckdbio.scan_sql(
                            path, path_header, nproc_in=nproc_in,
                            chrom_type=chrom_type,
                        ),
                    )
                )

            query = " UNION ALL ".join(parts)
            if not concatenate:
                query += " ORDER BY {}".format(
                    ", ".join(quote_identifier(key) for key in DEFAULT_SORT_KEYS)
                )

            if arrowio.is_parquet(output):
                duckdbio.copy_to_parquet(
                    con, query, output, out_header, row_group_size
                )
            else:
                result = con.execute(query)
                with PairsWriter(
                    output,
                    out_header,
                    compress_program=compress_program,
                    nproc_out=kwargs.get("nproc_out", 8),
                ) as writer:
                    for batch in duckdbio.result_batches(result):
                        writer.write(batch)
        finally:
            con.close()

    duckdbio.run_with_chrom_enum_fallback(run, output, description="merge")


def _headers_for(paths, headers):
    """Line up one header per path, re-reading if --keep-first-header trimmed them.

    Only the text reader needs a header (to know how many lines to skip and what
    the columns are), and each file's own header is the right one for that --
    the merged header is only for the output.
    """
    if len(headers) == len(paths):
        return headers
    return [arrowio.read_header(path) for path in paths]

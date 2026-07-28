"""Duplicate detection over Arrow batches.

The algorithm is pairtools': every chunk goes through
``pairtools.lib.dedup._dedup_chunk`` unchanged, with all of its backends and
options. What is reimplemented here is only the loop around it -- the ~30 lines
that chunk the input, carry the tail of each chunk into the next so duplicates
spanning a boundary are still found, and split the marked chunk into outputs.

That loop exists upstream as ``_dedup_stream``, but it hard-codes
``pd.read_table(in_stream, chunksize=...)`` for input and ``.to_csv()`` for
output, so it cannot be handed Arrow batches. UPSTREAM.md records the change
that would let this call into pairtools instead.

Chunks are rebuilt to exactly ``--chunksize`` rows before being marked, because
which pairs land in a chunk together determines what gets compared, and
matching pairtools row-for-row means matching its chunk boundaries.
"""

import numpy as np
import pandas as pd
import pyarrow as pa
from pairtools.lib import headerops
from pairtools.lib.dedup import _dedup_chunk
from pairtools.lib.stats import PairCounter

from . import arrowio
from .arrowio import PairsWriter, open_pairs
from .chunking import rechunk

UTIL_NAME = "pairtools_parquet_dedup"

#: Backends that work on DataFrames, and so can be driven from Arrow batches.
SUPPORTED_BACKENDS = ("scipy", "sklearn")


def mark_duplicate_chunks(
    frames,
    colnames,
    carryover,
    max_mismatch,
    method,
    mark_dups,
    keep_parent_id,
    extra_col_pairs,
    backend,
    n_proc,
    unmapped_chrom,
    c1="chrom1",
    c2="chrom2",
    p1="pos1",
    p2="pos2",
    s1="strand1",
    s2="strand2",
):
    """Yield each chunk with a `duplicate` column added.

    Mirrors ``pairtools.lib.dedup._dedup_stream``, differing only in taking
    DataFrames rather than reading them from a text stream.
    """
    df_prev_nodups = pd.DataFrame([])

    for df in frames:
        df = df.copy()
        df["carryover"] = False
        input_chunk = pd.concat(
            [df_prev_nodups, df], axis=0, ignore_index=True
        ).reset_index(drop=True)

        df_marked = _dedup_chunk(
            input_chunk,
            r=max_mismatch,
            method=method,
            keep_parent_id=keep_parent_id,
            extra_col_pairs=list(extra_col_pairs),
            backend=backend,
            n_proc=n_proc,
            c1=c1,
            c2=c2,
            p1=p1,
            p2=p2,
            s1=s1,
            s2=s2,
            unmapped_chrom=unmapped_chrom,
        )

        df_marked = (
            df_marked[~df_marked["carryover"]]
            .drop(columns=["carryover"])
            .reset_index(drop=True)
        )

        mask_duplicated = df_marked["duplicate"]
        if mark_dups:
            df_marked.loc[mask_duplicated, "pair_type"] = "DD"

        yield df_marked

        # The tail of the non-duplicates is prepended to the next chunk, so a
        # duplicate pair straddling the boundary is still seen as a pair.
        df_nodups = df_marked.loc[~mask_duplicated, colnames]
        df_prev_nodups = df_nodups.tail(carryover).reset_index(drop=True)
        df_prev_nodups["carryover"] = True


def _write(writer, df, columns):
    if writer is None or not len(df):
        return
    writer.write(pa.Table.from_pandas(df[columns], preserve_index=False))


def dedup_pairs(
    input_path,
    output,
    output_dups=None,
    output_unmapped=None,
    output_stats=None,
    max_mismatch=3,
    method="max",
    backend="scipy",
    chunksize=10_000,
    carryover=100,
    n_proc=1,
    mark_dups=True,
    keep_parent_id=False,
    extra_col_pair=(),
    unmapped_chrom="!",
    c1="chrom1",
    c2="chrom2",
    p1="pos1",
    p2="pos2",
    s1="strand1",
    s2="strand2",
    yaml=False,
    compress_program="auto",
    row_group_size=None,
    util_name=UTIL_NAME,
    **kwargs
):
    """Remove PCR duplicates from a sorted .pairs/.parquet file.

    Parameters
    ----------
    input_path : str
        Input path; must be sorted, as `pairtools dedup` requires.
    output : str
        Where the non-duplicate pairs go.
    output_dups, output_unmapped : str, optional
        Where duplicates and unmapped pairs go. Passing the same path as
        `output` for `output_dups` keeps duplicates in the main output, as in
        pairtools.
    output_stats : str, optional
        Where to write the duplication statistics.
    """
    if backend not in SUPPORTED_BACKENDS:
        raise NotImplementedError(
            "the {!r} backend is not available in pairtools_parquet: it works "
            "line-by-line on text streams rather than on dataframes. Use "
            "--backend scipy (the default) or --backend sklearn.".format(backend)
        )

    header, reader = open_pairs(
        input_path,
        batch_size=chunksize,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )
    columns = headerops.extract_column_names(header)

    out_header = headerops.append_new_pg(header, ID=util_name, PN=util_name)
    dups_header = out_header
    dups_columns = columns
    if keep_parent_id:
        dups_columns = columns + ["parent_readID"]
        # append_columns rewrites the #columns line of the list it is given and
        # returns that same list, so it has to be handed a copy.
        dups_header = headerops.append_columns(list(out_header), ["parent_readID"])

    stats = PairCounter() if output_stats else None

    writer_kwargs = dict(
        compress_program=compress_program,
        row_group_size=row_group_size,
        nproc_out=kwargs.get("nproc_out", 8),
    )

    frames = rechunk((batch.to_pandas() for batch in reader), chunksize)
    marked = mark_duplicate_chunks(
        frames,
        colnames=columns,
        carryover=carryover,
        max_mismatch=max_mismatch,
        method=method,
        mark_dups=mark_dups,
        keep_parent_id=keep_parent_id,
        extra_col_pairs=extra_col_pair,
        backend=backend,
        n_proc=n_proc,
        unmapped_chrom=unmapped_chrom,
        c1=c1,
        c2=c2,
        p1=p1,
        p2=p2,
        s1=s1,
        s2=s2,
    )

    dups_to_main_output = bool(output_dups) and str(output_dups) == str(output)

    writers = {}
    try:
        # pairtools drops parent_readID only just before writing the
        # non-duplicates, so the dups and unmapped streams still carry it.
        main_header = dups_header if dups_to_main_output else out_header
        main_columns = dups_columns if dups_to_main_output else columns

        writers["out"] = PairsWriter(output, main_header, **writer_kwargs)
        if output_dups and not dups_to_main_output:
            writers["dups"] = PairsWriter(output_dups, dups_header, **writer_kwargs)
        if output_unmapped:
            writers["unmapped"] = PairsWriter(
                output_unmapped, dups_header, **writer_kwargs
            )

        for df in marked:
            if stats is not None:
                stats.add_pairs_from_dataframe(df, unmapped_chrom=unmapped_chrom)

            mask_mapped = np.logical_and(
                (df[c1] != unmapped_chrom), (df[c2] != unmapped_chrom)
            )
            mask_duplicates = df["duplicate"].to_numpy(dtype=bool)
            df = df.drop(columns=["duplicate"])

            _write(writers.get("unmapped"), df.loc[~mask_mapped], dups_columns)

            if dups_to_main_output:
                _write(writers["out"], df.loc[mask_mapped], main_columns)
            else:
                _write(writers.get("dups"), df.loc[mask_duplicates], dups_columns)
                _write(
                    writers["out"], df.loc[mask_mapped & (~mask_duplicates)], columns
                )
    finally:
        for writer in writers.values():
            writer.close()

    if stats is not None:
        with open(output_stats, "w") as f:
            stats.save(f, yaml=yaml)

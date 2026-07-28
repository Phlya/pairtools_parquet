"""Removing pairs from high-coverage regions.

Coverage here means: for each end of each pair, how many other ends lie within
`max_dist` on the same chromosome. That is the same neighbour search `dedup`
does, and it gets the same treatment -- bucket the ends by position and the
search becomes an equi-join, since two ends within `max_dist` are either in the
same bucket or in an adjacent one.

pairtools computes it with ``_filterbycov``, whose own docstring reads "This is
a slow version of the filtering code used for testing purposes only. Use
cythonized version in the future!!" -- the cythonized version never arrived, so
that Python double loop is what `pairtools filterbycov` runs. It is 95s of the
110s this tool used to take on 5.6M pairs. It is still here as
``--backend python``, and is the reference the DuckDB backend is tested against.

Upstream buffers the whole file before computing coverage -- the algorithm needs
every pair at once -- so this does too, and the memory profile is no worse than
`pairtools filterbycov`.
"""

import numpy as np
import pyarrow as pa
from pairtools.lib import headerops, pairsam_format
from pairtools.lib.filterbycov import _filterbycov
from pairtools.lib.stats import PairCounter

from . import duckdb_utils
from .arrowio import PairsWriter, open_pairs
from .markasdup import SAM_COLUMNS, mark_sam_column

UTIL_NAME = "pairtools_parquet_filterbycov"

#: Pair type recorded in the statistics for a high-coverage pair.
HIGH_COVERAGE_PAIR_TYPE = "FF"

#: Backends for the coverage calculation.
SUPPORTED_BACKENDS = ("duckdb", "python")


def _encode(values, mapping):
    """Map labels to integers in first-seen order, as pairtools' fetchadd does."""
    out = np.empty(len(values), dtype=np.int64)
    for i, value in enumerate(values):
        if value not in mapping:
            mapping[value] = len(mapping)
        out[i] = mapping[value]
    return out


def _ends_table(table, c1, p1, c2, p2):
    """The 2N (chromosome, position) ends of N pairs, side 1 first.

    Same layout as the stacked array pairtools builds, so end `i` is side 1 of
    pair `i` and end `N + i` is side 2 -- which is what lets the two counts be
    recombined per pair by slicing.
    """
    n = table.num_rows
    chrom = pa.concat_arrays(
        [
            table.column(c1).combine_chunks().cast(pa.string()),
            table.column(c2).combine_chunks().cast(pa.string()),
        ]
    )
    pos = pa.concat_arrays(
        [
            table.column(p1).combine_chunks().cast(pa.int64()),
            table.column(p2).combine_chunks().cast(pa.int64()),
        ]
    )
    return pa.table(
        {
            "eid": pa.array(np.arange(2 * n, dtype=np.int64)),
            "chrom": chrom,
            "pos": pos,
        }
    )


NEIGHBOUR_SQL = """
WITH b AS (
  SELECT eid, chrom, pos, pos // {width} AS blk FROM ends
), p AS (
              SELECT eid, chrom, pos, blk     AS kb FROM b
    UNION ALL SELECT eid, chrom, pos, blk - 1       FROM b
    UNION ALL SELECT eid, chrom, pos, blk + 1       FROM b
)
SELECT a.eid AS eid, count(*) AS neighbours
FROM b a JOIN p ON a.chrom = p.chrom AND a.blk = p.kb
WHERE a.eid <> p.eid AND abs(a.pos - p.pos) <= {max_dist}
GROUP BY a.eid
"""


def coverage_duckdb(ends, max_dist, method, n_proc=4, tmpdir=None, memory=None):
    """Neighbours-within-`max_dist` per pair, as an equi-join over bucketed ends.

    A bucket is `max_dist` wide, so every end within `max_dist` of another is in
    that end's bucket or one either side; the probe side is emitted three times,
    once per neighbouring bucket, and each ordered pair is therefore generated
    exactly once.

    Counts exclude the end itself but include the pair's *other* end when the
    two are close enough, which is what the double loop does too.
    """
    if method not in ("sum", "max"):
        raise ValueError("Unknown method: {}".format(method))

    # `ends` holds both ends of every pair, side 1 first, so there are half as
    # many pairs as rows.
    n = ends.num_rows // 2
    if n == 0:
        return np.array([], dtype=np.int64)

    con = duckdb_utils.setup_duckdb_connection(
        temp_directory=tmpdir or None,
        memory_limit=memory or None,
        enable_progress_bar=False,
        enable_profiling="no_output",
        numb_threads=max(int(n_proc), 1),
    )
    try:
        con.register("ends", ends)
        found = con.execute(
            NEIGHBOUR_SQL.format(width=max(int(max_dist), 1), max_dist=int(max_dist))
        ).fetchnumpy()
    finally:
        con.close()

    counts = np.zeros(2 * n, dtype=np.int64)
    counts[found["eid"].astype(np.int64)] = found["neighbours"].astype(np.int64)

    first, second = counts[:n], counts[n:]
    # The +1 for the pair itself lands differently in the two methods: `sum`
    # adds it once to the total, `max` adds it to each side before taking the
    # larger. Not symmetric, but it is what `_filterbycov` does.
    if method == "sum":
        return first + second + 1
    return np.maximum(first + 1, second + 1)


def _add_pairs_to_stats(counter, table, columns, pair_types=None):
    """Feed rows to a PairCounter one at a time, as pairtools does here.

    `pair_types`, if given, replaces the column per row -- high-coverage pairs
    are counted as `FF`. Rows are fed in file order rather than grouped by
    outcome, because the counter records pair types in first-seen order and
    grouping them would reorder the output.
    """
    c1, p1, s1, c2, p2, s2, pt = (
        table.column(columns[name]).to_pylist()
        for name in ("c1", "p1", "s1", "c2", "p2", "s2", "pt")
    )
    for i in range(table.num_rows):
        counter.add_pair(
            c1[i], int(p1[i]), s1[i], c2[i], int(p2[i]), s2[i],
            pt[i] if pair_types is None else pair_types[i],
        )


def filterbycov_pairs(
    input_path,
    output,
    output_highcov=None,
    output_unmapped=None,
    output_stats=None,
    max_cov=8,
    max_dist=500,
    method="max",
    backend="duckdb",
    n_proc=4,
    mark_multi=False,
    unmapped_chrom=pairsam_format.UNMAPPED_CHROM,
    send_header_to="both",
    c1="chrom1",
    c2="chrom2",
    p1="pos1",
    p2="pos2",
    s1="strand1",
    s2="strand2",
    compress_program="auto",
    row_group_size=None,
    util_name=UTIL_NAME,
    **kwargs
):
    """Split `input_path` into low-coverage, high-coverage and unmapped pairs."""
    if backend not in SUPPORTED_BACKENDS:
        raise ValueError(
            "unknown backend {!r}; expected one of {}".format(
                backend, ", ".join(SUPPORTED_BACKENDS)
            )
        )

    header, reader = open_pairs(
        input_path,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )
    new_header = headerops.append_new_pg(header, ID=util_name, PN=util_name)
    column_names = headerops.extract_column_names(new_header)
    stat_columns = {
        "c1": c1, "c2": c2, "p1": p1, "p2": p2, "s1": s1, "s2": s2,
        "pt": "pair_type",
    }

    counter = PairCounter() if output_stats else None

    writer_kwargs = dict(
        compress_program=compress_program,
        row_group_size=row_group_size,
        nproc_out=kwargs.get("nproc_out", 8),
    )

    mapped_tables = []
    writers = {}
    try:
        if output_unmapped:
            writers["unmapped"] = PairsWriter(
                output_unmapped, new_header, **writer_kwargs
            )

        # Unmapped pairs stream straight out; the rest have to be held, because
        # coverage cannot be known until every pair has been seen.
        for batch in reader:
            table = pa.Table.from_batches([batch])
            is_unmapped = pa.compute.or_(
                pa.compute.equal(table.column(c1), unmapped_chrom),
                pa.compute.equal(table.column(c2), unmapped_chrom),
            )
            unmapped = table.filter(is_unmapped)
            if unmapped.num_rows:
                if counter is not None:
                    _add_pairs_to_stats(counter, unmapped, stat_columns)
                if "unmapped" in writers:
                    writers["unmapped"].write(unmapped)

            mapped = table.filter(pa.compute.invert(is_unmapped))
            if mapped.num_rows:
                mapped_tables.append(mapped)
    finally:
        for writer in writers.values():
            writer.close()

    mapped = (
        pa.concat_tables(mapped_tables)
        if mapped_tables
        else pa.table({name: pa.array([], type=pa.string()) for name in column_names})
    )

    coverage = np.array([], dtype=np.int64)
    if mapped.num_rows:
        if backend == "duckdb":
            coverage = coverage_duckdb(
                _ends_table(mapped, c1, p1, c2, p2),
                max_dist,
                method,
                n_proc=n_proc,
                tmpdir=kwargs.get("tmpdir"),
                memory=kwargs.get("memory"),
            )
        else:
            chrom_ids = {}
            coverage = np.asarray(
                _filterbycov(
                    _encode(mapped.column(c1).to_pylist(), chrom_ids),
                    mapped.column(p1).to_numpy(zero_copy_only=False),
                    _encode(mapped.column(c2).to_pylist(), chrom_ids),
                    mapped.column(p2).to_numpy(zero_copy_only=False),
                    max_dist,
                    method,
                )
            )

    is_high = coverage > max_cov
    low_rows = mapped.filter(pa.array(~is_high)) if mapped.num_rows else mapped
    high_rows = mapped.filter(pa.array(is_high)) if mapped.num_rows else mapped

    if counter is not None and mapped.num_rows:
        pair_types = [
            HIGH_COVERAGE_PAIR_TYPE if high else actual
            for high, actual in zip(is_high, mapped.column("pair_type").to_pylist())
        ]
        _add_pairs_to_stats(counter, mapped, stat_columns, pair_types=pair_types)

    with PairsWriter(
        output,
        new_header,
        write_header=send_header_to in ("both", "lowcov"),
        **writer_kwargs
    ) as writer:
        writer.write(low_rows)

    if output_highcov:
        if mark_multi and high_rows.num_rows:
            high_rows = _mark_as_duplicates(high_rows, column_names)
        with PairsWriter(
            output_highcov,
            new_header,
            write_header=send_header_to in ("both", "highcov"),
            **writer_kwargs
        ) as writer:
            writer.write(high_rows)

    if counter is not None:
        with open(output_stats, "w") as f:
            counter.save(f)


def _mark_as_duplicates(table, column_names):
    """Apply `--mark-multi`: tag the high-coverage pairs as duplicates."""
    columns = {name: table.column(name) for name in table.column_names}
    if "pair_type" in columns:
        columns["pair_type"] = pa.chunked_array(
            [pa.array(["DD"] * table.num_rows, type=pa.string())]
        )
    for name in SAM_COLUMNS:
        if name in columns:
            columns[name] = pa.chunked_array(
                [mark_sam_column(columns[name].combine_chunks())]
            )
    return pa.Table.from_arrays(
        [columns[name] for name in table.column_names], names=table.column_names
    )

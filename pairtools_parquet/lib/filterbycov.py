"""Removing pairs from high-coverage regions.

The coverage calculation is pairtools' ``_filterbycov``, called unchanged.
Upstream already buffers the whole file before calling it once -- the algorithm
needs every pair at the same time -- so this does the same, and the memory
profile is no worse than `pairtools filterbycov`.

Only the plumbing differs: pairs arrive as Arrow batches and leave through
:class:`~pairtools_parquet.lib.arrowio.PairsWriter`, so any input and output
format works.
"""

import numpy as np
import pyarrow as pa
from pairtools.lib import headerops, pairsam_format
from pairtools.lib.filterbycov import _filterbycov
from pairtools.lib.stats import PairCounter

from .arrowio import PairsWriter, open_pairs
from .markasdup import SAM_COLUMNS, mark_sam_column

UTIL_NAME = "pairtools_parquet_filterbycov"

#: Pair type recorded in the statistics for a high-coverage pair.
HIGH_COVERAGE_PAIR_TYPE = "FF"


def _encode(values, mapping):
    """Map labels to integers in first-seen order, as pairtools' fetchadd does."""
    out = np.empty(len(values), dtype=np.int64)
    for i, value in enumerate(values):
        if value not in mapping:
            mapping[value] = len(mapping)
        out[i] = mapping[value]
    return out


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
    mark_multi=False,
    unmapped_chrom=pairsam_format.UNMAPPED_CHROM,
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
        chrom_ids, strand_ids = {}, {}
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

    with PairsWriter(output, new_header, **writer_kwargs) as writer:
        writer.write(low_rows)

    if output_highcov:
        if mark_multi and high_rows.num_rows:
            high_rows = _mark_as_duplicates(high_rows, column_names)
        with PairsWriter(output_highcov, new_header, **writer_kwargs) as writer:
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

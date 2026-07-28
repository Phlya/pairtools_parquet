"""Flipping pairs onto the upper triangle.

Side 1 and side 2 are swapped where needed so that
``(order(chrom1), pos1) <= (order(chrom2), pos2)``, with the chromosome order
taken from a .chromsizes file. Unannotated chromosomes -- ones absent from that
file -- sort after annotated ones and compare to each other by name, matching
`pairtools flip`.

Unlike the other tools here this one is vectorized rather than row-wise: the
comparison is a few numpy operations over whole columns, so it is checked
against `pairtools flip` column by column in tests/test_smalltools.py.
"""

import numpy as np
import pyarrow as pa
from pairtools.lib import headerops, pairsam_format

from .._logging import get_logger
from .arrowio import PairsWriter, open_pairs

UTIL_NAME = "pairtools_parquet_flip"

logger = get_logger()


def chrom_order_map(chroms_path):
    """Map chromosome name to its rank, unmapped first, as `pairtools flip` does."""
    chromosomes = headerops.get_chrom_order(chroms_path)
    names = [pairsam_format.UNMAPPED_CHROM] + list(chromosomes)
    return {name: rank for rank, name in enumerate(names)}


def columns_to_flip(column_names):
    """Pairs of column indices to swap: every `<name>1` with its `<name>2`."""
    return [
        (col, col[:-1] + "2")
        for col in column_names
        if col.endswith("1") and (col[:-1] + "2") in column_names
    ]


def _reverse_pair_type(values):
    """Swap the two characters of each pair type, e.g. `UR` -> `RU`."""
    return pa.compute.binary_join_element_wise(
        pa.compute.utf8_slice_codeunits(values, 1, 2),
        pa.compute.utf8_slice_codeunits(values, 0, 1),
        "",
    )


def needs_flip(chrom1, chrom2, pos1, pos2, order):
    """Boolean mask of the rows whose two sides must be swapped."""
    ranks1 = np.array([order.get(c, -1) for c in chrom1], dtype=np.int64)
    ranks2 = np.array([order.get(c, -1) for c in chrom2], dtype=np.int64)
    annotated1 = ranks1 >= 0
    annotated2 = ranks2 >= 0

    # Both annotated: compare (rank, pos) lexicographically.
    both = annotated1 & annotated2
    correct = np.empty(len(ranks1), dtype=bool)
    correct[both] = (ranks1[both] < ranks2[both]) | (
        (ranks1[both] == ranks2[both]) & (pos1[both] <= pos2[both])
    )

    # Exactly one annotated: the annotated side goes first.
    correct[annotated1 & ~annotated2] = True
    correct[annotated2 & ~annotated1] = False

    # Neither annotated: fall back to comparing the names.
    neither = ~annotated1 & ~annotated2
    if neither.any():
        correct[neither] = np.array(
            [c1 < c2 for c1, c2 in zip(chrom1[neither], chrom2[neither])], dtype=bool
        )

    return ~correct


def flip_pairs(
    input_path,
    output,
    chroms_path,
    compress_program="auto",
    row_group_size=None,
    util_name=UTIL_NAME,
    **kwargs
):
    """Flip the pairs of `input_path` onto the upper triangle."""
    order = chrom_order_map(chroms_path)

    header, reader = open_pairs(
        input_path,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )
    new_header = headerops.append_new_pg(header, ID=util_name, PN=util_name)

    column_names = headerops.extract_column_names(new_header) or list(
        pairsam_format.COLUMNS
    )
    flip_pairs_of_columns = columns_to_flip(column_names)
    has_pair_type = "pair_type" in column_names
    warned = False

    with PairsWriter(
        output,
        new_header,
        schema=reader.schema,
        compress_program=compress_program,
        row_group_size=row_group_size,
        nproc_out=kwargs.get("nproc_out", 8),
    ) as writer:
        for batch in reader:
            table = pa.Table.from_batches([batch])
            chrom1 = np.asarray(table.column("chrom1").to_pylist(), dtype=object)
            chrom2 = np.asarray(table.column("chrom2").to_pylist(), dtype=object)
            pos1 = table.column("pos1").to_numpy(zero_copy_only=False)
            pos2 = table.column("pos2").to_numpy(zero_copy_only=False)

            if not warned:
                unannotated = [c for c in set(chrom1) | set(chrom2) if c not in order]
                if unannotated:
                    logger.warning("Unannotated chromosomes in the pairs file!")
                    warned = True

            mask = pa.array(needs_flip(chrom1, chrom2, pos1, pos2, order))

            columns = {name: table.column(name) for name in table.column_names}
            for col1, col2 in flip_pairs_of_columns:
                first, second = columns[col1], columns[col2]
                columns[col1] = pa.compute.if_else(mask, second, first)
                columns[col2] = pa.compute.if_else(mask, first, second)
            if has_pair_type:
                pair_type = table.column("pair_type")
                columns["pair_type"] = pa.compute.if_else(
                    mask, _reverse_pair_type(pair_type.combine_chunks()), pair_type
                )

            writer.write(
                pa.Table.from_arrays(
                    [columns[name] for name in table.column_names],
                    names=table.column_names,
                )
            )

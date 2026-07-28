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
import pyarrow.compute as pc
from pairtools.lib import headerops, pairsam_format

from .._logging import get_logger
from .arrowio import PairsWriter, open_pairs

UTIL_NAME = "pairtools_parquet_flip"

logger = get_logger()


def chrom_order(chroms_path):
    """The chromosome names in rank order, unmapped first, as `pairtools flip` has them."""
    chromosomes = headerops.get_chrom_order(chroms_path)
    return pa.array(
        [pairsam_format.UNMAPPED_CHROM] + list(chromosomes), type=pa.string()
    )


def chrom_ranks(column, order):
    """Rank of each chromosome; -1 for one the chromsizes file does not name.

    ``pc.index_in`` is the whole reason this is fast: looking the names up in a
    Python dict means materializing the column as Python strings first, which
    on 5.6M rows costs 7s per column against 0.4s here.
    """
    ranks = pc.fill_null(pc.index_in(column, value_set=order), -1)
    if isinstance(ranks, pa.ChunkedArray):
        ranks = ranks.combine_chunks()
    return ranks.to_numpy(zero_copy_only=False).astype(np.int64)


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


def needs_flip(ranks1, ranks2, pos1, pos2, column1=None, column2=None):
    """Boolean mask of the rows whose two sides must be swapped.

    `column1`/`column2` are only read for rows where *neither* side is in the
    chromsizes file, which is the one case that compares names rather than
    ranks. Passing the Arrow columns rather than Python strings keeps that cost
    proportional to the number of such rows, which is normally none.
    """
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

    # Neither annotated: compare (name, pos), as the both-annotated branch
    # compares (rank, pos). pairtools compares only the names here, which is
    # false for two sides of the *same* unannotated chromosome, so it swaps
    # them on every run and never looks at the positions -- flipping such a
    # file repeatedly oscillates instead of settling on the upper triangle.
    # See UPSTREAM.md; this is a deliberate divergence.
    neither = ~annotated1 & ~annotated2
    if neither.any():
        rows = pa.array(np.flatnonzero(neither))
        names1 = column1.take(rows).to_pylist()
        names2 = column2.take(rows).to_pylist()
        correct[neither] = np.array(
            [
                (c1, p1) <= (c2, p2)
                for c1, p1, c2, p2 in zip(
                    names1, pos1[neither], names2, pos2[neither]
                )
            ],
            dtype=bool,
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
    order = chrom_order(chroms_path)

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
            column1 = table.column("chrom1").combine_chunks()
            column2 = table.column("chrom2").combine_chunks()
            ranks1 = chrom_ranks(column1, order)
            ranks2 = chrom_ranks(column2, order)
            pos1 = table.column("pos1").to_numpy(zero_copy_only=False)
            pos2 = table.column("pos2").to_numpy(zero_copy_only=False)

            if not warned and ((ranks1 < 0).any() or (ranks2 < 0).any()):
                logger.warning("Unannotated chromosomes in the pairs file!")
                warned = True

            mask = pa.array(
                needs_flip(ranks1, ranks2, pos1, pos2, column1, column2)
            )

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

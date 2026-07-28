"""Tagging every pair in a file as a duplicate.

`pair_type` becomes `DD` for every row. When the file carries SAM records
(`sam1`/`sam2`), each record is tagged too, via pairtools' own
``mark_sam_as_dup`` -- that is bit-twiddling on SAM flags and tags, and is not
worth reimplementing.
"""

import pyarrow as pa
from pairtools.lib import headerops, pairsam_format
from pairtools.lib.dedup import mark_sam_as_dup

from .arrowio import PairsWriter, open_pairs

UTIL_NAME = "pairtools_parquet_markasdup"

SAM_COLUMNS = ("sam1", "sam2")


def mark_sam_column(values):
    """Tag every SAM record in one column as a duplicate."""
    marked = []
    for entry in values.to_pylist():
        if entry is None:
            marked.append(entry)
            continue
        marked.append(
            pairsam_format.INTER_SAM_SEP.join(
                mark_sam_as_dup(sam)
                for sam in entry.split(pairsam_format.INTER_SAM_SEP)
            )
        )
    return pa.array(marked, type=pa.string())


def markasdup_pairs(
    input_path,
    output,
    compress_program="auto",
    row_group_size=None,
    util_name=UTIL_NAME,
    **kwargs
):
    """Write `input_path` to `output` with every pair marked as a duplicate."""
    header, reader = open_pairs(
        input_path,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )
    new_header = headerops.append_new_pg(header, ID=util_name, PN=util_name)

    column_names = headerops.extract_column_names(new_header)
    sam_columns = [c for c in SAM_COLUMNS if c in column_names]

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
            columns = {name: table.column(name) for name in table.column_names}

            if "pair_type" in columns:
                columns["pair_type"] = pa.chunked_array(
                    [pa.array(["DD"] * table.num_rows, type=pa.string())]
                )
            for name in sam_columns:
                columns[name] = pa.chunked_array(
                    [mark_sam_column(columns[name].combine_chunks())]
                )

            writer.write(
                pa.Table.from_arrays(
                    [columns[name] for name in table.column_names],
                    names=table.column_names,
                )
            )

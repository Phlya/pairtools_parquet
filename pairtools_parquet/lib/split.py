"""Splitting a .pairsam back into .pairs and .sam.

The two outputs are different kinds of thing, so they are written differently:
the pairs go through :class:`PairsWriter` and can be .pairs, .pairs.gz or
.parquet, while the SAM records are plain text handed to pairtools' own
``auto_open``, which knows how to make a .bam or write to stdout.

Only the columns an output actually needs are read, which is where Parquet pays
off here: ``--output-pairs`` alone never touches the ``sam1``/``sam2`` columns,
and in a .pairsam those are most of the file.
"""

from pairtools.lib import fileio, headerops, pairsam_format

from .arrowio import PairsWriter, open_pairs

UTIL_NAME = "pairtools_parquet_split"

SAM_COLUMNS = ("sam1", "sam2")

#: What a .pairsam puts in a sam column when there are no records for that side.
NO_SAM_RECORDS = "."


def split_columns(columns):
    """Split a column list into (pairs columns, sam columns present).

    Raises if only one of the two sam columns is declared, as upstream does --
    a .pairsam with one sam column is not a format anything can read.
    """
    present = [name for name in SAM_COLUMNS if name in columns]
    if len(present) == 1:
        raise ValueError(
            "According to the #columns header field only one sam entry is present"
        )
    return [name for name in columns if name not in SAM_COLUMNS], present


def sam_records(value):
    """The tab-separated SAM lines packed into one sam1/sam2 field."""
    if value == NO_SAM_RECORDS:
        return []
    return [
        record.replace(pairsam_format.SAM_SEP, "\t")
        for record in value.split(pairsam_format.INTER_SAM_SEP)
    ]


def split_pairs(
    input_path,
    output_pairs="",
    output_sam="",
    compress_program="auto",
    row_group_size=None,
    util_name=UTIL_NAME,
    **kwargs
):
    """Write the pairs of `input_path` without their sam columns, and the SAM.

    Either output may be omitted, in which case that half is dropped. The SAM
    output may be ``-`` for stdout; the pairs output may not, since Parquet has
    no meaningful stdout form.
    """
    if not output_pairs and not output_sam:
        raise ValueError("At least one output (pairs and/or sam) must be specified!")
    if output_pairs == "-":
        raise ValueError(
            "The pairs output cannot be written to stdout; give it a path ending "
            "in .pairs, .pairs.gz or .parquet"
        )

    header = open_pairs(
        input_path,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )[0]
    columns = headerops.extract_column_names(header)
    pairs_columns, sam_columns = split_columns(columns)

    new_header = headerops.append_new_pg(list(header), ID=util_name, PN=util_name)
    new_header = headerops._update_header_entry(
        new_header, "columns", " ".join(pairs_columns)
    )

    # Read only what the requested outputs need, in the file's own order.
    wanted = set(pairs_columns if output_pairs else []) | set(
        sam_columns if output_sam else []
    )
    projection = [name for name in columns if name in wanted]

    _, reader = open_pairs(
        input_path,
        columns=projection,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )

    sam_stream = None
    pairs_writer = None
    try:
        if output_sam:
            sam_stream = fileio.auto_open(
                output_sam,
                mode="w",
                nproc=kwargs.get("nproc_out", 8),
                command=kwargs.get("cmd_out", None),
            )
            sam_stream.writelines(
                line[11:].strip() + "\n"
                for line in new_header
                if line.startswith("#samheader:")
            )

        if output_pairs:
            pairs_writer = PairsWriter(
                output_pairs,
                new_header,
                compress_program=compress_program,
                row_group_size=row_group_size,
                nproc_out=kwargs.get("nproc_out", 8),
            )

        for batch in reader:
            if pairs_writer is not None:
                pairs_writer.write(batch.select(pairs_columns))
            if sam_stream is not None and sam_columns:
                _write_sam(sam_stream, batch, sam_columns)
    finally:
        if pairs_writer is not None:
            pairs_writer.close()
        if sam_stream is not None and output_sam != "-":
            sam_stream.close()


def _write_sam(stream, batch, sam_columns):
    """Write the SAM records of one batch, side 1 before side 2 for each pair."""
    sides = [batch.column(name).to_pylist() for name in sam_columns]
    for row in zip(*sides):
        for value in row:
            for record in sam_records(value):
                stream.write(record)
                stream.write("\n")

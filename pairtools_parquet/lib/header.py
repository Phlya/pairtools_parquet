"""Manipulating the .pairs header.

pairtools' `header` subcommands change the header and copy the body through
untouched, which for text means shelling out to `cat`. Here the body goes
through the same Arrow reader and writer as every other tool, so these commands
also cross formats: a header can be transferred from a .pairs file onto a
.parquet, or the columns of a .parquet renamed and the result written back as
text.

The header itself is built by pairtools -- ``make_standard_pairsheader``,
``set_columns``, ``insert_samheader_pysam``, ``append_new_pg`` -- so a header
this package writes is the header pairtools would have written.

Two deliberate differences from upstream, both recorded in UPSTREAM.md:

* the output path is required, because Parquet cannot be written to a pipe and
  every other command in this package takes a required ``-o``;
* ``set-columns`` on a file with no ``#columns:`` line adds one. Upstream's
  ``headerops.set_columns`` only rewrites a line that already exists, so
  set-columns on a headerless file -- the case it exists for -- silently
  produces output with no header at all.
"""

import sys
import warnings

import pyarrow.parquet as pq
from pairtools.lib import fileio, headerops, pairsam_format

from . import arrowio
from .arrowio import PairsWriter, open_pairs

UTIL_NAME = "pairtools_parquet_header"


def columns_line_index(header):
    """Index of the ``#columns:`` line in `header`, or None."""
    for i, line in enumerate(header):
        if line.startswith("#columns:"):
            return i
    return None


def set_columns(header, columns):
    """``headerops.set_columns``, but adding the line when there is none."""
    if columns_line_index(header) is None:
        # Laid out exactly as headerops.set_columns lays it out, separator
        # after the colon included.
        return list(header) + [
            "#columns:" + headerops.SEP_COLS + headerops.SEP_COLS.join(columns)
        ]
    return headerops.set_columns(list(header), columns)


def check_body_columns(path, columns, nproc_in=3, cmd_in=None):
    """Raise unless the first data row has one field per name in `columns`.

    Mirrors pairtools' check, blind spot included: only the first row is looked
    at, so a file whose rows disagree with each other still passes.
    """
    found = None
    if arrowio.is_parquet(path):
        found = len(pq.ParquetFile(path).schema_arrow.names)
        matches = found == len(columns)
    else:
        instream = fileio.auto_open(path, mode="r", nproc=nproc_in, command=cmd_in)
        try:
            headerops.get_header(instream, ignore_warning=True)
            matches = headerops.validate_cols(instream, list(columns))
        finally:
            if instream is not sys.stdin:
                instream.close()

    if not matches:
        raise ValueError(
            "Number of columns mismatch:\n\t#columns: {}{}".format(
                headerops.SEP_COLS.join(columns),
                "" if found is None else "\n\tfile has {} columns".format(found),
            )
        )


def rewrite(input_path, output, new_header, **kwargs):
    """Stream the body of `input_path` into `output` under `new_header`.

    The columns are read positionally under the new header's names, so this both
    renames them and gives each one the type its new name implies.
    """
    columns = headerops.extract_column_names(new_header)
    if not columns:
        raise ValueError("the new header has no '#columns:' line")

    _, reader = open_pairs(
        input_path,
        column_names=columns,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )
    with PairsWriter(
        output,
        new_header,
        compress_program=kwargs.get("compress_program", "auto"),
        row_group_size=kwargs.get("row_group_size", None),
        nproc_out=kwargs.get("nproc_out", 8),
    ) as writer:
        writer.write_all(reader)


def generate_header(
    input_path,
    output,
    chroms_path=None,
    sam_path=None,
    columns="",
    extra_columns="",
    assembly="",
    no_flip=False,
    pairs=True,
    **kwargs
):
    """Put a standard header on a file, taking chromosome sizes from elsewhere.

    Parameters
    ----------
    chroms_path : str, optional
        A chrom.sizes file to take the chromosome sizes from.
    sam_path : str, optional
        A .sam/.bam file to inherit both the chromosome sizes and the @PG/@SQ
        records from. One of `chroms_path` and `sam_path` is required.
    """
    samheader = None
    if sam_path:
        # Imported here so pysam is only needed by the commands that read SAM.
        from pairtools.lib.parse_pysam import AlignmentFilePairtoolized

        input_sam = AlignmentFilePairtoolized(
            sam_path, "r", threads=kwargs.get("nproc_in", 1)
        )
        samheader = input_sam.header
        chromsizes = headerops.get_chromsizes_from_pysam_header(samheader)
    elif chroms_path:
        chromsizes = headerops.get_chromsizes_from_file(chroms_path)
    else:
        raise ValueError(
            "Either --chroms-path or --sam-path is required, to store the "
            "chromosome sizes in the header"
        )

    if columns:
        column_names = columns.split(",")
    else:
        column_names = list(
            pairsam_format.COLUMNS_PAIRS if pairs else pairsam_format.COLUMNS_PAIRSAM
        )
    if extra_columns:
        column_names += extra_columns.split(",")

    new_header = headerops.make_standard_pairsheader(
        assembly=assembly,
        chromsizes=chromsizes,
        columns=column_names,
        shape="whole matrix" if no_flip else "upper triangle",
    )
    if samheader is not None:
        new_header = headerops.insert_samheader_pysam(new_header, samheader)
    new_header = headerops.append_new_pg(new_header, ID=UTIL_NAME, PN=UTIL_NAME)

    check_body_columns(
        input_path,
        column_names,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )
    rewrite(input_path, output, new_header, **kwargs)


def transfer_header(input_path, output, reference_file, **kwargs):
    """Replace the header of `input_path` with the one from `reference_file`.

    The reference may itself be a .pairs or a .parquet file, so a header can be
    moved in either direction between the two formats.
    """
    reference_header = arrowio.read_header(
        reference_file,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )
    reference_columns = headerops.extract_column_names(reference_header)

    check_body_columns(
        input_path,
        reference_columns,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )
    new_header = headerops.append_new_pg(
        list(reference_header), ID=UTIL_NAME, PN=UTIL_NAME
    )
    rewrite(input_path, output, new_header, **kwargs)


def set_columns_header(input_path, output, columns, **kwargs):
    """Rename the columns of `input_path` to `columns`."""
    header = arrowio.read_header(
        input_path,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )
    column_names = columns.split(",")
    rewrite(input_path, output, set_columns(header, column_names), **kwargs)


def validate_columns_header(
    input_path, output, reference_file="", reference_columns="", **kwargs
):
    """Check a file's columns, and pass it through unchanged if they check out.

    With `reference_file` or `reference_columns`, the file's column names must
    match those. Either way the first data row must have one field per column.
    """
    header = arrowio.read_header(
        input_path,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )
    pairs_columns = headerops.extract_column_names(header)

    if reference_columns:
        reference_columns = reference_columns.split(",")

    if reference_file:
        if reference_columns:
            warnings.warn(
                "--reference-columns are ignored, as --reference-file is provided"
            )
        reference_columns = headerops.extract_column_names(
            arrowio.read_header(
                reference_file,
                nproc_in=kwargs.get("nproc_in", 3),
                cmd_in=kwargs.get("cmd_in", None),
            )
        )

    if reference_columns and pairs_columns != reference_columns:
        raise ValueError(
            "Pairs columns differ from reference columns:\n\t{}\n\t{}".format(
                pairs_columns, reference_columns
            )
        )

    check_body_columns(
        input_path,
        pairs_columns,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )
    new_header = headerops.append_new_pg(list(header), ID=UTIL_NAME, PN=UTIL_NAME)
    rewrite(input_path, output, new_header, **kwargs)

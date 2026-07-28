#!/usr/bin/env python
import functools

import click
from pairtools.lib import pairsam_format

from ..lib.header import (
    generate_header,
    set_columns_header,
    transfer_header,
    validate_columns_header,
)
from . import cli, common_io_options


@cli.group()
def header():
    """Manipulate the .pairs/.pairsam/.parquet header"""
    pass


def common_header_options(func):
    """The arguments and options every `header` subcommand takes.

    Unlike upstream, ``--output`` is required: Parquet has no meaningful stdout
    form, and every other command in this package writes to a named file.
    """

    @click.argument("pairs_path", type=str, required=True)
    @click.option(
        "-o",
        "--output",
        type=str,
        required=True,
        help="output file. The format follows the extension: .pairs, .pairs.gz"
        " or .parquet.",
    )
    @click.option(
        "--compress-program",
        type=str,
        default="auto",
        show_default=True,
        help="A binary to compress text output. Ignored for .parquet output.",
    )
    @common_io_options
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


@header.command()
@click.option(
    "--chroms-path",
    type=str,
    default=None,
    required=False,
    help="Chromosome order used to flip interchromosomal mates: "
    "path to a chromosomes file (e.g. UCSC chrom.sizes or similar) whose "
    "first column lists scaffold names. Any scaffolds not listed will be "
    "ordered lexicographically following the names provided.",
)
@click.option(
    "--sam-path",
    type=str,
    default=None,
    required=False,
    help="Input sam file to inherit the header."
    " Either --sam or --chroms-path should be provided to store the chromosome"
    " sizes in the header.",
)
@click.option(
    "--columns",
    type=click.STRING,
    default="",
    help="Report columns describing alignments "
    "Can take multiple values as a comma-separated list."
    "By default, assign standard .pairs columns: {}".format(
        ",".join(pairsam_format.COLUMNS)
    ),
)
@click.option(
    "--extra-columns",
    type=click.STRING,
    default="",
    help="Report extra columns describing alignments "
    "Can take multiple values as a comma-separated list.",
)
@click.option(
    "--assembly",
    type=str,
    default="",
    help="Name of genome assembly (e.g. hg19, mm10) to store in the pairs header.",
)
@click.option(
    "--no-flip",
    is_flag=True,
    help="If specified, assume that the pairs are not flipped in genomic order"
    " and instead preserve the order in which they were sequenced.",
)
@click.option(
    "--pairs/--pairsam",
    is_flag=True,
    default=True,
    help="If pairs, then the default columns will be set to: {}\nif pairsam,"
    " then to: {}".format(
        ",".join(pairsam_format.COLUMNS_PAIRS),
        ",".join(pairsam_format.COLUMNS_PAIRSAM),
    ),
)
@common_header_options
def generate(pairs_path, output, **kwargs):
    """Generate the header

    PAIRS_PATH : input .pairs/.pairsam/.parquet file, with or without a header
    of its own.
    """
    generate_header(pairs_path, output, **kwargs)


@header.command()
@click.option(
    "--reference-file", "-r", help="Header file for transfer", type=str, required=True
)
@common_header_options
def transfer(pairs_path, output, reference_file, **kwargs):
    """Transfer the header from one pairs file to another

    PAIRS_PATH : input .pairs/.pairsam/.parquet file.
    """
    transfer_header(pairs_path, output, reference_file, **kwargs)


@header.command(name="set-columns")
@click.option(
    "--columns",
    "-c",
    help="Comma-separated list of columns to be set, e.g.: {}".format(
        ",".join(pairsam_format.COLUMNS)
    ),
    type=str,
    required=True,
)
@common_header_options
def set_columns(pairs_path, output, columns, **kwargs):
    """Set the columns of the .pairs/.pairsam/.parquet file

    PAIRS_PATH : input .pairs/.pairsam/.parquet file.
    """
    set_columns_header(pairs_path, output, columns, **kwargs)


@header.command(name="validate-columns")
@click.option(
    "--reference-file",
    "-r",
    help="Header file for comparison (optional)",
    type=str,
    required=False,
    default="",
)
@click.option(
    "--reference-columns",
    "-c",
    help="Comma-separated list of columns for the check (optional), e.g.: {}".format(
        ",".join(pairsam_format.COLUMNS)
    ),
    type=str,
    required=False,
    default="",
)
@common_header_options
def validate_columns(pairs_path, output, reference_file, reference_columns, **kwargs):
    """Validate the columns of the file against a reference or within the file.

    If the checks pass, then the full file is written to the output. Otherwise
    an exception is raised.

    If reference_file is provided, check:
        1) columns are the same between the input and reference_file
        2) the number of columns in the body is the same as the number of columns

    If reference_columns are provided, check:
        1) the input's columns are the same as provided
        2) the number of columns in the body is the same as the number of columns

    If neither is given, check only the number of columns in the body.
    Checks only the first row of the body!

    PAIRS_PATH : input .pairs/.pairsam/.parquet file.
    """
    validate_columns_header(
        pairs_path, output, reference_file, reference_columns, **kwargs
    )


if __name__ == "__main__":
    header()

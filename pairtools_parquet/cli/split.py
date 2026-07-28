#!/usr/bin/env python
import click

from ..lib.split import split_pairs
from . import cli, common_io_options


@cli.command()
@click.argument("pairsam_path", type=str, required=True)
@click.option(
    "--output-pairs",
    type=str,
    default="",
    help="output pairs file. The format follows the extension: .pairs,"
    " .pairs.gz or .parquet."
    " If not specified, pairs are dropped.",
)
@click.option(
    "--output-sam",
    type=str,
    default="",
    help="output sam file."
    " If the path ends with .bam, the output is compressed into a bam file."
    " If -, sam entries are printed to stdout."
    " If not specified, sam entries are dropped.",
)
@click.option(
    "--compress-program",
    type=str,
    default="auto",
    show_default=True,
    help="A binary to compress text pairs output. Ignored for .parquet output.",
)
@common_io_options
def split(pairsam_path, output_pairs, output_sam, compress_program, **kwargs):
    """Split a .pairsam file into .pairs and .sam.

    Restore a .sam file from sam1 and sam2 fields of a .pairsam file. Create
    a .pairs file without sam1/sam2 fields.

    PAIRSAM_PATH : input .pairsam/.parquet file.
    """
    split_pairs(
        pairsam_path,
        output_pairs=output_pairs,
        output_sam=output_sam,
        compress_program=compress_program,
        **kwargs,
    )


if __name__ == "__main__":
    split()

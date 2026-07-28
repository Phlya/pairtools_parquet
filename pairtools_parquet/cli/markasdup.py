#!/usr/bin/env python
import click

from ..lib.markasdup import markasdup_pairs
from . import cli, common_io_options


@cli.command()
@click.argument("pairsam_path", type=str, required=False)
@click.option(
    "-o",
    "--output",
    type=str,
    default="",
    help="output file. The format follows the extension: .pairs, .pairs.gz or .parquet.",
)
@click.option(
    "--compress-program",
    type=str,
    default="auto",
    show_default=True,
    help="A binary to compress text output. Ignored for .parquet output.",
)
@common_io_options
def markasdup(pairsam_path, output, compress_program, **kwargs):
    """Tag all pairs in the input file as duplicates.

    Change the type of all pairs inside a .pairs/.pairsam/.parquet file to DD.
    If sam entries are present, change the pair type in the Yt SAM tag to
    'Yt:Z:DD'.

    PAIRSAM_PATH : input .pairs/.pairsam/.parquet file.
    """
    markasdup_pairs(
        pairsam_path, output, compress_program=compress_program, **kwargs
    )


if __name__ == "__main__":
    markasdup()

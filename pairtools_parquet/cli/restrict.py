#!/usr/bin/env python
import click

from ..lib.restrict import restrict_pairs
from . import cli, common_io_options


@cli.command()
@click.argument("pairs_path", type=str, required=False)
@click.option(
    "-f",
    "--frags",
    type=str,
    required=True,
    help="a tab-separated BED file with the positions of restriction fragments"
    " (chrom, start, end). Can be generated using cooler digest.",
)
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
def restrict(pairs_path, frags, output, compress_program, **kwargs):
    """Assign restriction fragments to pairs.

    Identify the restriction fragments that got ligated into a Hi-C molecule.

    Note: rfrags are 0-indexed.

    PAIRS_PATH : input .pairs/.pairsam/.parquet file.
    """
    restrict_pairs(
        pairs_path, output, frags, compress_program=compress_program, **kwargs
    )


if __name__ == "__main__":
    restrict()

#!/usr/bin/env python
import click

from ..lib.flip import flip_pairs
from . import cli, common_io_options


@cli.command()
@click.argument("pairs_path", type=str, required=True)
@click.option(
    "-c",
    "--chrom-subset",
    "chroms_path",
    type=str,
    required=True,
    help="Chromosome order used to flip interchromosomal mates: "
    "path to a chromosomes file (e.g. UCSC chrom.sizes or similar) whose "
    "first column lists scaffold names.",
)
@click.option(
    "-o",
    "--output",
    type=str,
    required=True,
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
def flip(pairs_path, chroms_path, output, compress_program, **kwargs):
    """Flip pairs to get an upper-triangular matrix.

    Change the order of side1 and side2 in pairs, such that
    (order(chrom1) < order(chrom2) or (order(chrom1) == order(chrom2)) and
    (pos1 <= pos2)). Equivalent to reflecting the lower triangle of a Hi-C
    matrix onto its upper triangle. The order of chromosomes must be provided
    via a .chromsizes file.

    PAIRS_PATH : input .pairs/.pairsam/.parquet file.
    """
    flip_pairs(
        pairs_path, output, chroms_path,
        compress_program=compress_program, **kwargs
    )


if __name__ == "__main__":
    flip()

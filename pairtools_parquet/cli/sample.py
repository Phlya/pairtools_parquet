#!/usr/bin/env python
import click

from ..lib.sample import sample_pairs
from . import cli, common_io_options


@cli.command()
@click.argument("fraction", type=float, required=True)
@click.argument("pairs_path", type=str, required=True)
@click.option(
    "-o",
    "--output",
    type=str,
    required=True,
    help="output file. The format follows the extension: .pairs, .pairs.gz or .parquet.",
)
@click.option(
    "-s",
    "--seed",
    type=int,
    default=None,
    help="the seed of the random number generator.",
)
@click.option(
    "--compress-program",
    type=str,
    default="auto",
    show_default=True,
    help="A binary to compress text output. Ignored for .parquet output.",
)
@common_io_options
def sample(fraction, pairs_path, output, seed, compress_program, **kwargs):
    """Select a random subset of pairs in a pairs file.

    FRACTION: the fraction of the randomly selected pairs subset

    PAIRS_PATH : input .pairs/.pairsam/.parquet file.
    """
    sample_pairs(
        pairs_path, output, fraction, seed=seed,
        compress_program=compress_program, **kwargs
    )


if __name__ == "__main__":
    sample()

#!/usr/bin/env python
import click

from ..lib.scaling import scaling_pairs
from . import cli, common_io_options


@cli.command()
@click.argument("input_path", type=str, nargs=-1, required=False)
@click.option(
    "-o", "--output", type=str, default="", help="output .tsv file with summary."
)
@click.option(
    "--view",
    "--regions",
    help="Path to a BED file which defines which regions (viewframe) of the"
    " chromosomes to use. By default, this is parsed from .pairs header. ",
    type=str,
    required=False,
    default=None,
)
@click.option(
    "--chunksize",
    type=int,
    default=100_000,
    show_default=True,
    required=False,
    help="Number of pairs in each chunk. Reduce for lower memory footprint.",
)
@click.option(
    "--dist-range",
    type=click.Tuple([int, int]),
    default=(1, 1_000_000_000),
    show_default=True,
    required=False,
    help="Distance range. ",
)
@click.option(
    "--n-dist-bins-decade",
    type=int,
    default=8,
    show_default=True,
    required=False,
    help="Number of bins to split the distance range in log10-space, specified"
    " per a factor of 10 difference.",
)
@common_io_options
def scaling(
    input_path, output, view, chunksize, dist_range, n_dist_bins_decade, **kwargs
):
    """Calculate pairs scalings.

    INPUT_PATH : a .pairs/.pairsam/.parquet file to calculate statistics.

    The files with paths ending with .gz/.lz4 are decompressed by bgzip/lz4c.

    Output is .tsv file with scaling stats (both cis scalings and trans levels).
    """
    if len(input_path) != 1:
        raise click.UsageError(
            "exactly one input file is expected, got {}".format(len(input_path))
        )

    scaling_pairs(
        input_path[0],
        output,
        view=view,
        chunksize=chunksize,
        dist_range=dist_range,
        n_dist_bins_decade=n_dist_bins_decade,
        **kwargs,
    )


if __name__ == "__main__":
    scaling()

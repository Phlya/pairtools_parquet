#!/usr/bin/env python
import click

from ..lib.stats import merge_stats, stats_pairs
from . import cli, common_io_options


@cli.command()
@click.argument("input_path", type=str, nargs=-1, required=True)
@click.option("-o", "--output", type=str, required=True, help="output stats tsv file.")
@click.option(
    "--merge",
    is_flag=True,
    help="If specified, merge multiple input stats files instead of calculating"
    " statistics of a .pairs file. Merging is performed via summation of"
    " all overlapping statistics. Non-overlapping statistics are appended to"
    " the end of the file.",
)
@click.option(
    "--n-dist-bins-decade",
    type=int,
    default=8,
    show_default=True,
    help="Number of bins to split the distance range in log10 space.",
)
@click.option(
    "--with-chromsizes/--no-chromsizes",
    default=True,
    show_default=True,
    help="If enabled, will store sizes of chromosomes from the header of the"
    " pairs file in the stats file.",
)
@click.option(
    "--yaml/--no-yaml",
    default=False,
    show_default=True,
    help="Output stats in yaml format instead of table.",
)
@click.option(
    "--bytile-dups/--no-bytile-dups",
    default=False,
    show_default=True,
    help="If enabled, will analyse by-tile duplication statistics to estimate"
    " library complexity more accurately. Requires parent_readID column.",
)
@click.option(
    "--output-bytile-stats",
    default="",
    type=str,
    help="If specified, will analyse by-tile duplication statistics and write"
    " them to this file.",
)
@click.option(
    "--filter",
    "filters",
    multiple=True,
    help="Filters with conditions to apply to the data, e.g."
    ' \'unique:(pair_type=="UU")\'. Multiple filters can be provided.',
)
@click.option(
    "--engine",
    type=str,
    default="pandas",
    show_default=True,
    help="Engine to use for filter evaluation.",
)
@click.option(
    "--startup-code",
    type=str,
    default="",
    help="An auxiliary code to execute before filtering.",
)
@click.option(
    "-t",
    "--type-cast",
    type=(str, str),
    default=(),
    multiple=True,
    help="Cast a given column to a given type for filter evaluation.",
)
@click.option(
    "--chunksize",
    type=int,
    default=100_000,
    show_default=True,
    help="Number of pairs read at a time.",
)
@common_io_options
def stats(
    input_path,
    output,
    merge,
    n_dist_bins_decade,
    with_chromsizes,
    yaml,
    bytile_dups,
    output_bytile_stats,
    filters,
    engine,
    startup_code,
    type_cast,
    chunksize,
    **kwargs,
):
    """Calculate pairs statistics.

    INPUT_PATH : by default, a .pairs/.pairsam/.parquet file to calculate
    statistics. If --merge is specified, then INPUT_PATH is interpreted as an
    arbitrary number of stats files to merge.
    """
    if merge:
        merge_stats(
            output, input_path, n_dist_bins_decade=n_dist_bins_decade,
            yaml=yaml, **kwargs
        )
        return

    if len(input_path) != 1:
        raise click.UsageError(
            "exactly one input file is expected without --merge, got {}".format(
                len(input_path)
            )
        )

    stats_pairs(
        input_path[0],
        output,
        n_dist_bins_decade=n_dist_bins_decade,
        with_chromsizes=with_chromsizes,
        yaml=yaml,
        bytile_dups=bytile_dups,
        output_bytile_stats=output_bytile_stats or None,
        filters=list(filters),
        startup_code=startup_code,
        type_cast=type_cast,
        engine=engine,
        chunksize=chunksize,
        **kwargs,
    )


if __name__ == "__main__":
    stats()

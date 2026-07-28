#!/usr/bin/env python
import click

from pairtools.lib import pairsam_format

from ..lib.filterbycov import filterbycov_pairs
from . import cli, common_io_options


@cli.command()
@click.argument("pairs_path", type=str, required=False)
@click.option(
    "-o", "--output", type=str, default="",
    help="output file for pairs from low coverage regions."
    " The format follows the extension: .pairs, .pairs.gz or .parquet.",
)
@click.option(
    "--output-highcov", type=str, default="",
    help="output file for pairs from high coverage regions."
    " If not provided, such pairs are dropped.",
)
@click.option(
    "--output-unmapped", type=str, default="",
    help="output file for unmapped pairs. If not provided, they are dropped.",
)
@click.option(
    "--output-stats", type=str, default="",
    help="output file for coverage filter statistics.",
)
@click.option(
    "--max-cov", type=int, default=8, show_default=True,
    help="The maximum allowed coverage per region.",
)
@click.option(
    "--max-dist", type=int, default=500, show_default=True,
    help="The resolution for calculating coverage. For each pair, the local "
    "coverage around each end is calculated as (1 + the number of neighbouring "
    "pairs within +/- max_dist bp).",
)
@click.option(
    "--method", type=click.Choice(["max", "sum"]), default="max", show_default=True,
    help='calculate the number of neighbouring pairs as either the "max" or the '
    '"sum" of the pairs within +/- max_dist bp of each end.',
)
@click.option(
    "--backend",
    type=click.Choice(["duckdb", "python"]),
    default="duckdb",
    show_default=True,
    help="How to compute coverage. duckdb counts neighbouring ends with a "
    "bucketed equi-join; python is pairtools' own double loop, kept as the "
    "reference duckdb is tested against.",
)
@click.option(
    "-p", "--n-proc", type=int, default=4, show_default=True,
    help="Number of cores to use. Applies to the duckdb backend.",
)
@click.option(
    "--mark-multi", is_flag=True,
    help="If specified, duplicate pairs are marked as FF in pair_type and as a "
    "duplicate in the sam entries.",
)
@click.option(
    "--unmapped-chrom", type=str, default=pairsam_format.UNMAPPED_CHROM,
    show_default=True, help="Placeholder for a chromosome on an unmapped side.",
)
@click.option(
    "--send-header-to",
    type=click.Choice(["lowcov", "highcov", "both", "none"]),
    default="both",
    show_default=True,
    help="Which of the outputs should receive header and comment lines. "
    "Applies to text outputs only: a .parquet keeps its header either way, "
    "since there it is metadata rather than leading lines and a file without "
    "it cannot be read back as pairs.",
)
@click.option("--c1", type=str, default="chrom1", show_default=True, help="Chrom 1 column.")
@click.option("--c2", type=str, default="chrom2", show_default=True, help="Chrom 2 column.")
@click.option("--p1", type=str, default="pos1", show_default=True, help="Position 1 column.")
@click.option("--p2", type=str, default="pos2", show_default=True, help="Position 2 column.")
@click.option("--s1", type=str, default="strand1", show_default=True, help="Strand 1 column.")
@click.option("--s2", type=str, default="strand2", show_default=True, help="Strand 2 column.")
@click.option(
    "--compress-program", type=str, default="auto", show_default=True,
    help="A binary to compress text output. Ignored for .parquet output.",
)
@common_io_options
def filterbycov(
    pairs_path, output, output_highcov, output_unmapped, output_stats,
    max_cov, max_dist, method, backend, n_proc, mark_multi, unmapped_chrom,
    send_header_to,
    c1, c2, p1, p2, s1, s2, compress_program, **kwargs,
):
    """Remove pairs from regions of high coverage.

    Find and remove pairs from regions with high coverage, which are likely to
    be spurious.

    PAIRS_PATH : input .pairs/.pairsam/.parquet file.
    """
    filterbycov_pairs(
        pairs_path, output,
        output_highcov=output_highcov or None,
        output_unmapped=output_unmapped or None,
        output_stats=output_stats or None,
        max_cov=max_cov, max_dist=max_dist, method=method,
        backend=backend, n_proc=n_proc,
        mark_multi=mark_multi, unmapped_chrom=unmapped_chrom,
        send_header_to=send_header_to,
        c1=c1, c2=c2, p1=p1, p2=p2, s1=s1, s2=s2,
        compress_program=compress_program, **kwargs,
    )


if __name__ == "__main__":
    filterbycov()

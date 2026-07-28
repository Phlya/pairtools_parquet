#!/usr/bin/env python
import click

from ..lib.phase import phase_pairs
from . import cli, common_io_options


@cli.command()
@click.argument("pairs_path", type=str, required=False)
@click.option(
    "-o", "--output", type=str, default="",
    help="output file. The format follows the extension: .pairs, .pairs.gz or .parquet.",
)
@click.option(
    "--phase-suffixes", nargs=2, type=str, required=True,
    help="Suffixes of chromosome names that correspond to two parental haplotypes.",
)
@click.option(
    "--clean-output", is_flag=True,
    help="Drop all columns besides the standard ones and phase1/2.",
)
@click.option(
    "--tag-mode", type=click.Choice(["XA", "XB"]), default="XB", show_default=True,
    help="Specifies the mode of bwa reporting. XB is bwa-mem with -M flag, XA "
    "is the default bwa-mem reporting.",
)
@click.option(
    "--report-scores/--no-report-scores", default=False, show_default=True,
    help="Report alignment scores of the alignments used for phasing.",
)
@click.option(
    "--compress-program", type=str, default="auto", show_default=True,
    help="A binary to compress text output. Ignored for .parquet output.",
)
@common_io_options
def phase(
    pairs_path, output, phase_suffixes, clean_output, tag_mode, report_scores,
    compress_program, **kwargs,
):
    """Phase pairs mapped to a diploid genome.

    PAIRS_PATH : input .pairs/.pairsam/.parquet file, parsed with the alignment
    tags kept (--add-columns XB,AS,XS or --add-columns XA,NM,AS,XS).
    """
    phase_pairs(
        pairs_path, output, list(phase_suffixes), tag_mode=tag_mode,
        clean_output=clean_output, report_scores=report_scores,
        compress_program=compress_program, **kwargs,
    )


if __name__ == "__main__":
    phase()

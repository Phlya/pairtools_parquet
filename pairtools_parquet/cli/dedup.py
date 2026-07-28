#!/usr/bin/env python
import click

from pairtools.lib import pairsam_format

from ..lib.dedup import dedup_pairs
from . import cli, common_io_options


@cli.command()
@click.argument("pairs_path", type=str, required=True)
@click.option(
    "-o",
    "--output",
    type=str,
    required=True,
    help="output file for pairs after duplicate removal."
    " The format follows the extension: .pairs, .pairs.gz or .parquet.",
)
@click.option(
    "--output-dups",
    type=str,
    default="",
    help="output file for duplicated pairs. "
    "If the same as --output, duplicates are marked and kept in the main "
    "output. If not provided, duplicates are dropped.",
)
@click.option(
    "--output-unmapped",
    type=str,
    default="",
    help="output file for unmapped pairs. "
    "If not provided, unmapped pairs are dropped.",
)
@click.option(
    "--output-stats",
    type=str,
    default="",
    help="output file for duplicate statistics. "
    "If not provided, statistics are not printed.",
)
@click.option(
    "--max-mismatch",
    type=int,
    default=3,
    show_default=True,
    help="Pairs with both sides mapped within this distance (bp) from each "
    "other are considered duplicates. [dedup option]",
)
@click.option(
    "--method",
    type=click.Choice(["max", "sum"]),
    default="max",
    show_default=True,
    help='define the mismatch as either the "max" or the "sum" of the '
    "mismatches of the genomic locations of the both sides of the two "
    "compared molecules. [dedup option]",
)
@click.option(
    "--backend",
    type=click.Choice(["scipy", "sklearn"]),
    default="scipy",
    show_default=True,
    help="What backend to use: scipy and sklearn are based on KD-trees. "
    "pairtools' cython backend is not available here, as it works "
    "line-by-line on text streams. [dedup option]",
)
@click.option(
    "--chunksize",
    type=int,
    default=10_000,
    show_default=True,
    help="Number of pairs in each chunk. Reduce for lower memory footprint. "
    "[dedup option]",
)
@click.option(
    "--carryover",
    type=int,
    default=100,
    show_default=True,
    help="Number of deduped pairs carried over from bin to bin to do dedup at "
    "the bin margins. [dedup option]",
)
@click.option(
    "-p",
    "--n-proc",
    type=int,
    default=1,
    show_default=True,
    help="Number of cores to use. Only applies with sklearn backend. "
    "[dedup option]",
)
@click.option(
    "--mark-dups/--no-mark-dups",
    default=True,
    show_default=True,
    help='If specified, duplicate pairs are marked as DD in "pair_type" and '
    "as a duplicate in the sam entries. [output format option]",
)
@click.option(
    "--keep-parent-id",
    is_flag=True,
    help="If specified, duplicate pairs are marked with the readID of the "
    "retained deduped read in the 'parent_readID' field. [output format option]",
)
@click.option(
    "--extra-col-pair",
    nargs=2,
    type=str,
    multiple=True,
    help="Extra columns that also must match for two pairs to be marked as "
    "duplicates. Can be either provided as 0-based column indices or as column "
    'names (requires the "#columns" header field). The option can be provided '
    "more than once. [dedup option]",
)
@click.option(
    "--unmapped-chrom",
    type=str,
    default=pairsam_format.UNMAPPED_CHROM,
    show_default=True,
    help="Placeholder for a chromosome on an unmapped side. [input format option]",
)
@click.option(
    "--c1",
    type=str,
    default=pairsam_format.COLUMNS_PAIRS[1],
    show_default=True,
    help="Chrom 1 column; default: chrom1. [input format option]",
)
@click.option(
    "--c2",
    type=str,
    default=pairsam_format.COLUMNS_PAIRS[3],
    show_default=True,
    help="Chrom 2 column; default: chrom2. [input format option]",
)
@click.option(
    "--p1",
    type=str,
    default=pairsam_format.COLUMNS_PAIRS[2],
    show_default=True,
    help="Position 1 column; default: pos1. [input format option]",
)
@click.option(
    "--p2",
    type=str,
    default=pairsam_format.COLUMNS_PAIRS[4],
    show_default=True,
    help="Position 2 column; default: pos2. [input format option]",
)
@click.option(
    "--s1",
    type=str,
    default=pairsam_format.COLUMNS_PAIRS[5],
    show_default=True,
    help="Strand 1 column; default: strand1. [input format option]",
)
@click.option(
    "--s2",
    type=str,
    default=pairsam_format.COLUMNS_PAIRS[6],
    show_default=True,
    help="Strand 2 column; default: strand2. [input format option]",
)
@click.option(
    "--yaml/--no-yaml",
    default=False,
    show_default=True,
    help="Output stats in yaml format instead of table. [output format option]",
)
@click.option(
    "--compress-program",
    type=str,
    default="auto",
    show_default=True,
    help="A binary to compress text output. Ignored for .parquet output.",
)
@common_io_options
def dedup(
    pairs_path,
    output,
    output_dups,
    output_unmapped,
    output_stats,
    max_mismatch,
    method,
    backend,
    chunksize,
    carryover,
    n_proc,
    mark_dups,
    keep_parent_id,
    extra_col_pair,
    unmapped_chrom,
    c1,
    c2,
    p1,
    p2,
    s1,
    s2,
    yaml,
    compress_program,
    **kwargs,
):
    """Find and remove PCR/optical duplicates.

    Find PCR/optical duplicates in an upper-triangular flipped sorted
    .pairs/.pairsam/.parquet file. Allow for a +/-N bp mismatch at each side of
    duplicated molecules.

    PAIRS_PATH : input triu-flipped sorted .pairs/.pairsam/.parquet file.
    """
    dedup_pairs(
        pairs_path,
        output,
        output_dups=output_dups or None,
        output_unmapped=output_unmapped or None,
        output_stats=output_stats or None,
        max_mismatch=max_mismatch,
        method=method,
        backend=backend,
        chunksize=chunksize,
        carryover=carryover,
        n_proc=n_proc,
        mark_dups=mark_dups,
        keep_parent_id=keep_parent_id,
        extra_col_pair=extra_col_pair,
        unmapped_chrom=unmapped_chrom,
        c1=c1,
        c2=c2,
        p1=p1,
        p2=p2,
        s1=s1,
        s2=s2,
        yaml=yaml,
        compress_program=compress_program,
        **kwargs,
    )


if __name__ == "__main__":
    dedup()

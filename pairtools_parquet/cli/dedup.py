#!/usr/bin/env python
import click

from pairtools.lib import pairsam_format

from ..lib.dedup import dedup_pairs
from . import cli, common_io_options


@cli.command()
@click.argument("pairs_path", type=str, required=False)
@click.option(
    "-o",
    "--output",
    type=str,
    default="",
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
    "--output-bytile-stats",
    type=str,
    default="",
    help="output file for by-tile duplicate statistics. The readID must "
    "contain tile information for this to work, and --keep-parent-id is "
    "forced on, since the analysis reads the parent's tile. "
    "[output stats filtering option]",
)
@click.option(
    "--filter",
    "filters",
    multiple=True,
    help="Filter stats with conditions to apply to the data (similar to "
    "`pairtools select` or `pairtools stats`). For non-YAML output only the "
    "first filter is reported. Example: --yaml --filter "
    "'unique:(pair_type==\"UU\")'. [output stats filtering option]",
)
@click.option(
    "--chrom-subset",
    type=str,
    default=None,
    help="A path to a chromosomes file (tab-separated, 1st column contains "
    "chromosome names) containing a chromosome subset of interest for the "
    "stats filter. Only pairs with both sides in the subset are counted; what "
    "is written out is unaffected. Note that `pairtools dedup` accepts this "
    "option and then ignores it. [output stats filtering option]",
)
@click.option(
    "--engine",
    type=str,
    default="pandas",
    show_default=True,
    help="Engine to use for stats filter evaluation. "
    "[output stats filtering option]",
)
@click.option(
    "--startup-code",
    type=str,
    default="",
    help="An auxiliary code to execute before stats filtering. "
    "[output stats filtering option]",
)
@click.option(
    "-t",
    "--type-cast",
    type=(str, str),
    default=(),
    multiple=True,
    help="Cast a given column to a given type for stats filter evaluation. "
    "[output stats filtering option]",
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
    type=click.Choice(["duckdb", "scipy", "sklearn"]),
    default="duckdb",
    show_default=True,
    help="What backend to use. duckdb finds duplicate candidates with a "
    "blocked equi-join and is several times faster; scipy and sklearn are "
    "pairtools' KD-tree implementations, kept as the reference duckdb is "
    "tested against. pairtools' cython backend is not available here, as it "
    "works line-by-line on text streams. [dedup option]",
)
@click.option(
    "--chunksize",
    type=int,
    default=None,
    help="Number of pairs in each chunk. Reduce for lower memory footprint. "
    "Ignored with --backend duckdb --max-mismatch 0, which decides the whole "
    "file in one pass and has no window to size. "
    "[default: 20000000 with --backend duckdb, 10000 otherwise] [dedup option]",
)
@click.option(
    "--tmpdir",
    type=str,
    default="",
    help="Custom temporary folder for duckdb to spill to when the key table "
    "does not fit in --memory. [duckdb backend option]",
)
@click.option(
    "--memory",
    type=str,
    default="",
    help="Memory limit for the duckdb backend, e.g. '8G'. Above it duckdb "
    "spills to --tmpdir rather than failing, and the answer is unchanged. "
    "[default: duckdb's own, 80% of RAM] [duckdb backend option]",
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
    default=4,
    show_default=True,
    help="Number of cores to use. Applies to the duckdb and sklearn backends. "
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
    "--send-header-to",
    type=click.Choice(["dups", "dedup", "both", "none"]),
    default="both",
    show_default=True,
    help="Which of the outputs should receive header and comment lines. "
    "Applies to text outputs only: a .parquet keeps its header either way, "
    "since there it is metadata rather than leading lines and a file without "
    "it cannot be read back as pairs. [input format option]",
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
    output_bytile_stats,
    filters,
    chrom_subset,
    engine,
    startup_code,
    type_cast,
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
    send_header_to,
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

    With --backend duckdb --max-mismatch 0 the input need not be sorted: exact
    equality is transitive, so the answer does not depend on the order rows
    arrive in. Every other combination requires sorted input, as pairtools does.

    PAIRS_PATH : input triu-flipped sorted .pairs/.pairsam/.parquet file.
    """
    dedup_pairs(
        pairs_path,
        output,
        output_dups=output_dups or None,
        output_unmapped=output_unmapped or None,
        output_stats=output_stats or None,
        output_bytile_stats=output_bytile_stats or None,
        filters=list(filters),
        chrom_subset=chrom_subset,
        engine=engine,
        startup_code=startup_code,
        type_cast=type_cast,
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
        send_header_to=send_header_to,
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

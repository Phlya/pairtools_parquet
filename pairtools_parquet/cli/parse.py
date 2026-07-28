#!/usr/bin/env python
import click
from pairtools.lib import pairsam_format

from ..lib.parse import parse_pairs
from . import cli, common_io_options


@cli.command()
@click.argument("sam_path", type=str, required=False)
@click.option(
    "-c",
    "--chroms-path",
    type=str,
    required=True,
    help="Chromosome order used to flip interchromosomal mates: "
    "path to a chromosomes file (e.g. UCSC chrom.sizes or similar) whose "
    "first column lists scaffold names. Any scaffolds not listed will be "
    "ordered lexicographically following the names provided.",
)
@click.option(
    "-o",
    "--output",
    type=str,
    required=True,
    help="output file. The format follows the extension: .pairs, .pairs.gz,"
    " .pairsam or .parquet.",
)
@click.option(
    "--assembly",
    type=str,
    help="Name of genome assembly (e.g. hg19, mm10) to store in the pairs header.",
)
@click.option(
    "--min-mapq",
    type=int,
    default=1,
    show_default=True,
    help="The minimal MAPQ score to consider a read as uniquely mapped",
)
@click.option(
    "--max-molecule-size",
    type=int,
    default=750,
    show_default=True,
    help="The maximal size of a Hi-C molecule; used to rescue single ligations"
    "(from molecules with three alignments) and to rescue complex ligations."
    "The default is based on oriented P(s) at short ranges of multiple Hi-C."
    "Not used with walks-policy all.",
)
@click.option(
    "--drop-readid",
    is_flag=True,
    help="If specified, do not add read ids to the output",
)
@click.option(
    "--drop-seq",
    is_flag=True,
    help="If specified, remove sequences and PHREDs from the sam fields",
)
@click.option(
    "--drop-sam", is_flag=True, help="If specified, do not add sams to the output"
)
@click.option(
    "--add-pair-index",
    is_flag=True,
    help="If specified, each pair will have pair index in the molecule",
)
@click.option(
    "--add-columns",
    type=click.STRING,
    default="",
    help="Report extra columns describing alignments "
    "Possible values (can take multiple values as a comma-separated "
    "list): a SAM tag (any pair of uppercase letters) or {}.".format(
        ", ".join(pairsam_format.EXTRA_COLUMNS)
    ),
)
@click.option(
    "--output-parsed-alignments",
    type=str,
    default="",
    help="output file for all parsed alignments, including walks."
    " Useful for debugging and analysis of walks."
    " If file exists, it will be open in the append mode."
    " If the path ends with .gz or .lz4, the output is bgzip-/lz4-compressed."
    " By default, not used.",
)
@click.option(
    "--output-stats",
    type=str,
    default="",
    help="output file for various statistics of pairs file. "
    " By default, statistics is not generated.",
)
@click.option(
    "--report-alignment-end",
    type=click.Choice(["5", "3"]),
    default="5",
    help="specifies whether the 5' or 3' end of the alignment is reported as"
    " the position of the Hi-C read.",
)
@click.option(
    "--max-inter-align-gap",
    type=int,
    default=20,
    show_default=True,
    help="read segments that are not covered by any alignment and"
    ' longer than the specified value are treated as "null" alignments.'
    " These null alignments convert otherwise linear alignments into walks,"
    " and affect how they get reported as a Hi-C pair (see --walks-policy).",
)
@click.option(
    "--walks-policy",
    type=click.Choice(["mask", "5any", "5unique", "3any", "3unique", "all"]),
    default="5unique",
    help="the policy for reporting unrescuable walks (reads containing more"
    " than one alignment on one or both sides, that can not be explained by a"
    " single ligation between two mappable DNA fragments)."
    ' "mask" - mask walks (chrom="!", pos=0, strand="-"); '
    ' "5any" - report the 5\'-most alignment on each side;'
    ' "5unique" - report the 5\'-most unique alignment on each side, if present;'
    ' "3any" - report the 3\'-most alignment on each side;'
    ' "3unique" - report the 3\'-most unique alignment on each side, if present;'
    ' "all" - report all available unique alignments on each side.',
    show_default=True,
)
@click.option(
    "--readid-transform",
    type=str,
    default=None,
    help="A Python expression to modify read IDs. Useful when read IDs differ "
    "between the two reads of a pair. Must be a valid Python expression that "
    "uses variables called readID and/or i (the 0-based index of the read pair "
    "in the bam file) and returns a new value, e.g. \"readID[:-2]+'_'+str(i)\". "
    "Make sure that transformed readIDs remain unique!",
    show_default=True,
)
@click.option(
    "--flip/--no-flip",
    is_flag=True,
    default=True,
    help="If specified, do not flip pairs in genomic order and instead preserve "
    "the order in which they were sequenced.",
)
@click.option(
    "--compress-program",
    type=str,
    default="auto",
    show_default=True,
    help="A binary to compress text output. Ignored for .parquet output.",
)
@common_io_options
def parse(
    sam_path,
    chroms_path,
    output,
    output_parsed_alignments,
    output_stats,
    compress_program,
    **kwargs,
):
    """Find ligation pairs in .sam data, make .pairs.

    SAM_PATH : an input .sam/.bam file with paired-end sequence alignments of
    Hi-C molecules. If the path ends with .bam, the input is decompressed from
    bam with samtools. By default, the input is read from stdin.
    """
    parse_pairs(
        sam_path,
        chroms_path,
        output,
        output_parsed_alignments=output_parsed_alignments,
        output_stats=output_stats,
        parse2=False,
        compress_program=compress_program,
        **kwargs,
    )


if __name__ == "__main__":
    parse()

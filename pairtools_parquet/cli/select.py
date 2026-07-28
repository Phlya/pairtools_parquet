import click

from pairtools.lib import pairsam_format

from ..lib.select import select_pairs
from . import cli, common_io_options


@cli.command()
@click.argument("condition", type=str)
@click.argument("pairs_path", type=str, required=False)
@click.option(
    "-o",
    "--output",
    type=str,
    default="",
    help="Output file for selected pairs."
    " The format follows the extension: .pairs, .pairs.gz or .parquet.",
)
@click.option(
    "--output-rest",
    type=str,
    default=None,
    help="Output file for pairs of other types."
    " By default, such pairs are dropped.",
)
@click.option(
    "--chrom-subset",
    type=str,
    default=None,
    help="A path to a chromosomes file (tab-separated, 1st column contains "
    "chromosome names) containing a chromosome subset of interest. "
    "If provided, additionally filter pairs with both sides originating from "
    "the provided subset of chromosomes. This operation modifies the "
    "#chromosomes: and #chromsize: header fields accordingly.",
)
@click.option(
    "--startup-code",
    type=str,
    default=None,
    help="An auxiliary code to execute before filtering. "
    "Use to define functions that can be evaluated in the CONDITION statement",
)
@click.option(
    "-t",
    "--type-cast",
    type=(str, str),
    default=(),
    multiple=True,
    help="Cast a given column to a given type. By default, only pos and mapq "
    "are cast to int, other columns are kept as str. Provide as "
    "-t <column_name> <type>, e.g. -t read_len1 int. Multiple entries are allowed.",
)
@click.option(
    "--remove-columns",
    "-r",
    type=str,
    default="",
    help="Comma-separated list of columns to be removed, e.g.: {}".format(
        ",".join(pairsam_format.COLUMNS)
    ),
)
@click.option(
    "--compress-program",
    type=str,
    default="auto",
    show_default=True,
    help="A binary to compress text output. "
    "Suggested alternatives: gzip, lzop, lz4c, snzip. "
    'If "auto", then use pigz or lz4c if available, and gzip otherwise. '
    "Ignored for .parquet output.",
)
@common_io_options
def select(
    condition,
    pairs_path,
    output,
    output_rest,
    chrom_subset,
    startup_code,
    type_cast,
    remove_columns,
    compress_program,
    **kwargs,
):
    """Select pairs according to some condition.

    CONDITION : A Python expression; if it returns True, select the read pair.
    Any column declared in the #columns line of the pairs header can be
    accessed by its name.

    PAIRS_PATH : input .pairs/.pairsam/.parquet file. If the path ends with .gz
    or .lz4, the input is decompressed by bgzip/lz4c.

    The following functions can be used in CONDITION besides the standard
    Python functions:

    - csv_match(x, csv) - True if variable x is contained in a list of
    comma-separated values, e.g. csv_match(chrom1, 'chr1,chr2')

    - wildcard_match(x, wildcard) - True if variable x matches a wildcard,
    e.g. wildcard_match(pair_type, 'C*')

    - regex_match(x, regex) - True if variable x matches a Python-flavor regex,
    e.g. regex_match(chrom1, 'chr\\d')

    - region_match(chrom, pos, region_chrom, region_start, region_end) - True
    if the position (chrom, pos) lies within the region.

    Examples:
    pairtools_parquet select '(pair_type=="UU") or (pair_type=="UR")' -o out.parquet in.parquet
    """
    select_pairs(
        pairs_path,
        output,
        condition,
        output_rest=output_rest,
        chrom_subset=chrom_subset,
        startup_code=startup_code,
        type_cast=type_cast,
        remove_columns=remove_columns,
        compress_program=compress_program,
        **kwargs,
    )


if __name__ == "__main__":
    select()

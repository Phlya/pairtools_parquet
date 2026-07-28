#!/usr/bin/env python
import click

from ..lib.merge import merge_pairs
from . import cli, common_io_options


@cli.command()
@click.argument("pairs_path", nargs=-1, type=str)
@click.option(
    "-o",
    "--output",
    type=str,
    default="",
    help="output file."
    " The format follows the extension: .pairs, .pairs.gz or .parquet.",
)
@click.option(
    "--tmpdir",
    type=str,
    default="",
    help="Custom temporary folder for merged intermediates.",
)
@click.option(
    "--memory",
    type=str,
    default="2G",
    show_default=True,
    help="The amount of memory used by default.",
)
@click.option(
    "--compress-program",
    type=str,
    default="auto",
    show_default=True,
    help="A binary to compress text output. "
    "Suggested alternatives: gzip, lzop, lz4c, snzip. "
    "Ignored for .parquet output.",
)
@click.option(
    "--nproc",
    type=int,
    default=8,
    show_default=True,
    help="Number of threads for merging.",
)
@click.option(
    "--max-nmerge",
    type=int,
    default=8,
    show_default=True,
    help="The maximal number of inputs merged at once. For more, merged "
    "intermediates are stored in temporary .parquet files. Pass 0 to merge "
    "every input in one pass.",
)
@click.option(
    "--keep-first-header/--no-keep-first-header",
    default=False,
    show_default=True,
    help="Keep the first header or merge the headers together. Default: merge headers.",
)
@click.option(
    "--concatenate/--no-concatenate",
    default=False,
    show_default=True,
    help="Simple concatenate instead of merging sorted files.",
)
@common_io_options
def merge(
    pairs_path,
    output,
    tmpdir,
    memory,
    compress_program,
    nproc,
    max_nmerge,
    keep_first_header,
    concatenate,
    **kwargs,
):
    """Merge .pairs/.pairsam/.parquet files.

    By default, assumes that the files are sorted and maintains the sorting.

    If present, the @SQ records of the SAM header must be identical; the
    sorting order of these lines is taken from the first file in the list.
    The ID fields of the @PG records of the SAM header are modified with a
    numeric suffix to produce unique records. The other unique SAM and non-SAM
    header lines are copied into the output header.

    PAIRS_PATH : upper-triangular flipped sorted .pairs/.pairsam/.parquet files
    to merge, or a group/groups of files specified by a wildcard. Input formats
    may be mixed.
    """
    merge_pairs(
        list(pairs_path),
        output,
        concatenate=concatenate,
        keep_first_header=keep_first_header,
        max_nmerge=max_nmerge,
        nproc=nproc,
        tmpdir=tmpdir,
        memory=memory,
        compress_program=compress_program,
        **kwargs,
    )


if __name__ == "__main__":
    merge()

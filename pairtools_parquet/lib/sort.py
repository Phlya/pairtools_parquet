"""Sorting pairs, in any input format, to any output format.

DuckDB does the sort -- it is parallel and spills to disk, so the input need not
fit in memory. The sort ordering matches ``pairtools sort``: lexicographic on
the chromosome and pair-type columns, numeric on the positions.
"""

from pairtools.lib import headerops

from . import arrowio
from .csv_parquet_converter import resolve_keys
from .transform import quote_identifier, scan_transform_write

UTIL_NAME = "pairtools_parquet_sort"

#: Default sort order, matching `pairtools sort` and `pairtools merge`.
DEFAULT_SORT_KEYS = ["chrom1", "chrom2", "pos1", "pos2", "pair_type"]

__all__ = ["UTIL_NAME", "DEFAULT_SORT_KEYS", "quote_identifier", "sort_pairs"]


def sort_pairs(
    input_path,
    output_path,
    sort_keys,
    nproc=8,
    tmpdir=None,
    memory=None,
    compress_program="auto",
    row_group_size=None,
    util_name=UTIL_NAME,
    **kwargs
):
    """Sort a .pairs/.pairsam/.parquet file by `sort_keys`.

    Parameters
    ----------
    input_path, output_path : str
        Paths; the format of each is taken from its extension.
    sort_keys : list of str
        Columns to order by, in priority order. Each may be a column name or a
        numeric column index, as ``pairtools sort`` accepts.
    nproc : int
        Threads for DuckDB.
    tmpdir, memory : str
        DuckDB spill directory and memory limit.
    compress_program : str
        Which compressor to use for compressed text output. Whether the output
        is compressed at all follows its extension.
    row_group_size : int, optional
        Rows per Parquet row group.
    """
    header = arrowio.read_header(
        input_path,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )

    column_names = headerops.extract_column_names(header)
    sort_keys = resolve_keys(sort_keys, column_names)
    missing = [key for key in sort_keys if key not in column_names]
    if missing:
        raise ValueError(
            "sort column(s) {} not present in {}, which has columns {}".format(
                ", ".join(missing), input_path, ", ".join(column_names)
            )
        )

    new_header = headerops.append_new_pg(header, ID=util_name, PN=util_name)
    new_header = headerops.mark_header_as_sorted(new_header)

    scan_transform_write(
        input_path,
        output_path,
        header,
        new_header,
        order_by=sort_keys,
        nproc=nproc,
        tmpdir=tmpdir,
        memory=memory,
        compress_program=compress_program,
        row_group_size=row_group_size,
        description="sort",
        **kwargs
    )


def convert_pairs(
    input_path,
    output_path,
    util_name="pairtools_parquet_convert",
    nproc=8,
    tmpdir=None,
    memory=None,
    compress_program="auto",
    row_group_size=None,
    **kwargs
):
    """Copy a file between formats, preserving row order."""
    header = arrowio.read_header(
        input_path,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )
    new_header = headerops.append_new_pg(header, ID=util_name, PN=util_name)

    scan_transform_write(
        input_path,
        output_path,
        header,
        new_header,
        order_by=None,
        nproc=nproc,
        tmpdir=tmpdir,
        memory=memory,
        compress_program=compress_program,
        row_group_size=row_group_size,
        description="conversion",
        **kwargs
    )

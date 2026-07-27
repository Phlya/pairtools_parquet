"""Sorting pairs, in any input format, to any output format.

DuckDB does the sort -- it is parallel and spills to disk, so the input need not
fit in memory. The sort ordering matches ``pairtools sort``: lexicographic on
the chromosome and pair-type columns, numeric on the positions.

Where DuckDB can open the input itself it does, because its readers are
parallel; a stream it cannot open by path (stdin) is scanned through Arrow
instead. Both routes produce the same output.
"""

import os

from pairtools.lib import headerops

from .._logging import get_logger
from . import arrowio, duckdb_utils, duckdbio
from .arrowio import PairsWriter, open_pairs
from .csv_parquet_converter import resolve_keys

UTIL_NAME = "pairtools_parquet_sort"

logger = get_logger()


def quote_identifier(name):
    """Quote a column name for SQL, so a column called e.g. `order` is usable."""
    return '"{}"'.format(str(name).replace('"', '""'))


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
        Compressor for text output.
    row_group_size : int, optional
        Rows per Parquet row group.
    """
    nproc_in = kwargs.get("nproc_in", 3)
    cmd_in = kwargs.get("cmd_in", None)

    native = duckdbio.can_scan_natively(input_path) and cmd_in is None
    reader = None
    if native:
        header = arrowio.read_header(input_path, nproc_in=nproc_in, cmd_in=cmd_in)
    else:
        header, reader = open_pairs(input_path, nproc_in=nproc_in, cmd_in=cmd_in)

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

    order_by = ", ".join(quote_identifier(key) for key in sort_keys)

    def run(use_chrom_enum):
        con = duckdb_utils.setup_duckdb_connection(
            temp_directory=tmpdir or None,
            memory_limit=memory or None,
            enable_progress_bar=False,
            enable_profiling="no_output",
            numb_threads=nproc,
        )
        try:
            chrom_type = None
            if native:
                if use_chrom_enum:
                    chrom_type = duckdbio.declare_chrom_enum(con, header)
                source = duckdbio.scan_sql(
                    input_path, header, nproc_in=nproc_in, chrom_type=chrom_type
                )
            else:
                # Scanned lazily: DuckDB pulls batches off the reader as it sorts.
                con.register("pairs_input", reader)
                source = "pairs_input"

            query = "SELECT * FROM {} ORDER BY {}".format(source, order_by)

            if arrowio.is_parquet(output_path):
                duckdbio.copy_to_parquet(
                    con, query, output_path, new_header, row_group_size
                )
            else:
                result = con.execute(query)
                with PairsWriter(
                    output_path,
                    new_header,
                    compress_program=compress_program,
                    nproc_out=kwargs.get("nproc_out", 8),
                ) as writer:
                    for batch in duckdbio.result_batches(result):
                        writer.write(batch)
        finally:
            con.close()

    try:
        run(use_chrom_enum=True)
    except Exception as error:
        if not duckdbio.is_enum_domain_error(error):
            raise
        # A chromosome the header never declared. Redo the sort comparing
        # strings, which is slower but has no fixed domain, rather than fail on
        # a file whose header is merely incomplete.
        logger.warning(
            "%s contains a chromosome absent from its header; re-sorting "
            "without the chromosome ENUM optimization",
            input_path,
        )
        if os.path.exists(output_path):
            os.remove(output_path)
        run(use_chrom_enum=False)

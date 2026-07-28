"""The scan -> optionally reorder -> write pipeline shared by the SQL tools.

`sort`, `csv-to-parquet` and `parquet-to-csv` are the same operation with and
without an ORDER BY, so they share one implementation. Where DuckDB can open
the input itself it does, because its readers are parallel; a stream it cannot
open by path is scanned through Arrow instead.
"""

from . import arrowio, duckdb_utils, duckdbio
from .arrowio import PairsWriter, open_pairs


def quote_identifier(name):
    """Quote a column name for SQL, so a column called e.g. `order` is usable."""
    return '"{}"'.format(str(name).replace('"', '""'))


def scan_transform_write(
    input_path,
    output_path,
    in_header,
    out_header,
    order_by=None,
    nproc=8,
    tmpdir=None,
    memory=None,
    compress_program="auto",
    row_group_size=None,
    description="operation",
    **kwargs
):
    """Copy `input_path` to `output_path`, optionally ordering by `order_by`.

    Parameters
    ----------
    in_header : list of str
        The input's header; the text scanner needs it to know how many lines to
        skip and what the columns are.
    out_header : list of str
        The header to write, with provenance already applied.
    order_by : list of str, optional
        Columns to ORDER BY. Omit to copy in input order.
    """
    nproc_in = kwargs.get("nproc_in", 3)
    cmd_in = kwargs.get("cmd_in", None)
    native = duckdbio.can_scan_natively(input_path) and cmd_in is None

    def run(use_chrom_enum):
        reader = None
        con = duckdb_utils.setup_duckdb_connection(
            temp_directory=tmpdir or None,
            memory_limit=memory or None,
            enable_progress_bar=False,
            enable_profiling="no_output",
            numb_threads=nproc,
        )
        try:
            if native:
                chrom_type = (
                    duckdbio.declare_chrom_enum(con, in_header)
                    if use_chrom_enum
                    else None
                )
                source = duckdbio.scan_sql(
                    input_path, in_header, nproc_in=nproc_in, chrom_type=chrom_type
                )
            else:
                _, reader = open_pairs(
                    input_path, nproc_in=nproc_in, cmd_in=cmd_in
                )
                # Scanned lazily: DuckDB pulls batches off the reader as it works.
                con.register("pairs_input", reader)
                source = "pairs_input"

            query = "SELECT * FROM {}".format(source)
            if order_by:
                query += " ORDER BY {}".format(
                    ", ".join(quote_identifier(key) for key in order_by)
                )

            if arrowio.is_parquet(output_path):
                duckdbio.copy_to_parquet(
                    con, query, output_path, out_header, row_group_size
                )
            else:
                result = con.execute(query)
                with PairsWriter(
                    output_path,
                    out_header,
                    compress_program=compress_program,
                    nproc_out=kwargs.get("nproc_out", 8),
                ) as writer:
                    for batch in duckdbio.result_batches(result):
                        writer.write(batch)
        finally:
            con.close()

    duckdbio.run_with_chrom_enum_fallback(run, output_path, description=description)

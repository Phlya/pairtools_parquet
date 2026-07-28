"""Bridging .pairs files into DuckDB and back out again.

An operation that is expressible entirely in SQL -- sorting, filtering, merging
-- should let DuckDB read and write the files itself: its readers and its
Parquet writer are both parallel, whereas feeding batches through a Python
generator serializes the scan and holds the GIL while doing it. Measured on 4M
pairs, letting DuckDB own both ends is about 3x faster than streaming the same
query through Arrow.

Tools that need Python-side processing use :mod:`pairtools_parquet.lib.arrowio`
instead; this module is for the ones that do not.
"""

import os

from pairtools.lib import headerops, pairsam_format

from .._logging import get_logger
from . import arrowio, headerio
from .schema import duckdb_types_from_columns

logger = get_logger()


def sql_string(value):
    """Render a Python string as a DuckDB string literal."""
    return "'" + str(value).replace("'", "''") + "'"


def quote_identifier(name):
    """Quote a column name for SQL, so a column called e.g. `order` is usable."""
    return '"{}"'.format(str(name).replace('"', '""'))


#: Columns whose domain the header declares, so they can be DuckDB ENUMs.
CHROM_COLUMNS = ("chrom1", "chrom2")

CHROM_ENUM = "PAIRS_CHROM"

#: A value no .pairs field will hold, used to disable DuckDB's NULL literal.
#: DuckDB rejects an empty nullstr, and a NUL byte breaks its SQL parser, so
#: this has to be a printable token instead.
NEVER_A_PAIRS_VALUE = "__pairtools_parquet_never_null__"


class EnumDomainError(Exception):
    """A value fell outside a declared ENUM domain."""


def declare_chrom_enum(con, header):
    """Declare the chromosome ENUM for `header`, returning its name or None.

    Comparing small integers rather than strings makes the sort roughly a third
    faster, and chrom1/chrom2 are the leading sort keys.

    The values are declared in lexicographic order, because DuckDB orders an
    ENUM by declaration index and ``pairtools sort`` orders chromosomes by byte
    value (it shells out to `sort` under LC_ALL=C). Declaring them in any other
    order silently changes the sort -- which is exactly the bug this replaces,
    where the pair_type ENUM was declared in itertools.product order and so
    sorted tied pairs as UU before DD instead of DD before UU.
    """
    # extract_chromsizes returns a pandas Series keyed by chromosome, so the
    # names are its keys -- iterating it would yield the sizes.
    chromsizes = headerops.extract_chromsizes(header)
    if chromsizes is None or len(chromsizes) == 0:
        return None

    values = sorted(set(chromsizes.keys()) | {pairsam_format.UNMAPPED_CHROM})
    con.execute("DROP TYPE IF EXISTS {}".format(CHROM_ENUM))
    con.execute(
        "CREATE TYPE {} AS ENUM ({})".format(
            CHROM_ENUM, ", ".join(sql_string(value) for value in values)
        )
    )
    return CHROM_ENUM


def is_enum_domain_error(error):
    """Whether a DuckDB error is a value falling outside an ENUM domain.

    A chromosome the header never declared is a conversion error, not a null,
    so the caller has to retry without the ENUM rather than lose rows.
    """
    message = str(error)
    return "Conversion Error" in message and "Could not convert string" in message


#: Column name :func:`scan_sql` exposes a file-order row number under.
ROW_NUMBER_COLUMN = "__rid"


def scan_sql(path, header, nproc_in=3, chrom_type=None, row_number=False):
    """Return the SQL FROM-clause that reads `path` natively in DuckDB.

    Raises ValueError for inputs DuckDB cannot open by path, e.g. stdin; the
    caller should fall back to scanning an Arrow reader for those.

    With `row_number`, the scan also yields :data:`ROW_NUMBER_COLUMN`, the
    row's position in the file. Only Parquet can do this natively, via
    ``file_row_number``; for text the caller has to count the rows itself with
    :func:`arrowio.with_row_ids`, so this raises rather than returning a scan
    whose ordinal would be wrong.
    """
    if arrowio.is_stdio(path):
        raise ValueError("DuckDB cannot read a .pairs stream from stdin by path")
    path = str(path)

    columns = headerops.extract_column_names(header)

    if arrowio.is_parquet(path):
        scan = "read_parquet({}{})".format(
            sql_string(path), ", file_row_number=true" if row_number else ""
        )
        if chrom_type is None and not row_number:
            return scan
        # Parquet carries its own types, so the ENUM has to be applied on top.
        projected = ", ".join(
            "{col}::{type} AS {col}".format(col=quote_identifier(col), type=chrom_type)
            if chrom_type is not None and col in CHROM_COLUMNS
            else quote_identifier(col)
            for col in columns
        )
        if row_number:
            projected += ", file_row_number AS {}".format(
                quote_identifier(ROW_NUMBER_COLUMN)
            )
        return "(SELECT {} FROM {})".format(projected, scan)

    if row_number:
        raise ValueError(
            "DuckDB's CSV scanner cannot number rows in file order; "
            "count them with arrowio.with_row_ids instead"
        )

    if not columns:
        raise ValueError(
            "{} has no '#columns:' header line, so its columns are unknown".format(
                path
            )
        )

    column_types = duckdb_types_from_columns(columns)
    if chrom_type is not None:
        for col in CHROM_COLUMNS:
            if col in column_types:
                column_types[col] = chrom_type

    # The header is skipped positionally: auto_detect would otherwise take the
    # last header line for column names.
    #
    # nullstr is set to a value no .pairs field can hold, because DuckDB
    # otherwise reads an empty field as NULL. .pairs has no null: an empty
    # field is the empty string, and tools depend on that -- `phase` reads an
    # empty XB tag as "no alternative alignment" and calls .split() on it.
    return (
        "read_csv({path}, delim='\\t', skip={skip}, columns={columns}, "
        "header=false, auto_detect=false, nullstr={nullstr})"
    ).format(
        path=sql_string(path),
        skip=len(header),
        columns=column_types,
        nullstr=sql_string(NEVER_A_PAIRS_VALUE),
    )


def can_scan_natively(path):
    """Whether :func:`scan_sql` can handle `path`."""
    path = str(path)
    return not arrowio.is_stdio(path) and (
        arrowio.is_parquet(path)
        or path.endswith((".pairs", ".pairs.gz", ".pairsam", ".pairsam.gz"))
    )


def kv_metadata_sql(header):
    """Render a .pairs header as a DuckDB ``KV_METADATA`` literal.

    Values are escaped as SQL string literals, so tabs, backslashes and quotes
    inside a ``@PG`` record survive the round trip verbatim.
    """
    metadata = headerio.header_to_metadata(header)
    return "{" + ", ".join(
        "{}: {}".format(sql_string(key.decode("utf-8")), sql_string(value.decode("utf-8")))
        for key, value in metadata.items()
    ) + "}"


def copy_to_parquet(con, query, output_path, header, row_group_size=None):
    """Write the result of `query` to a Parquet file, header included."""
    options = ["FORMAT PARQUET", "KV_METADATA {}".format(kv_metadata_sql(header))]
    if row_group_size:
        options.append("ROW_GROUP_SIZE {:d}".format(row_group_size))

    con.execute(
        "COPY ({query}) TO {path} ({options})".format(
            query=query, path=sql_string(output_path), options=", ".join(options)
        )
    )


def run_with_chrom_enum_fallback(run, output_path, description="operation"):
    """Call ``run(use_chrom_enum=True)``, retrying without the ENUM if needed.

    The chromosome ENUM makes comparisons integer-wide, but its domain comes
    from the header and DuckDB raises rather than nulling on a value outside
    it. A file whose header simply does not list every chromosome it contains
    is not corrupt, so it must still be processed -- just more slowly.
    """
    try:
        run(True)
    except Exception as error:
        if not is_enum_domain_error(error):
            raise
        logger.warning(
            "input contains a chromosome absent from its header; redoing the %s "
            "without the chromosome ENUM optimization",
            description,
        )
        # A partial output may already exist from the failed attempt.
        if output_path and os.path.exists(str(output_path)):
            os.remove(str(output_path))
        run(False)


def result_batches(result):
    """Iterate a DuckDB result as Arrow record batches, across duckdb versions."""
    if hasattr(result, "to_arrow_reader"):
        return result.to_arrow_reader()
    return result.fetch_record_batch()

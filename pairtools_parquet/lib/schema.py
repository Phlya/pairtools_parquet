"""Arrow schemas for .pairs data.

The .pairs format carries its column names in the ``#columns:`` header line and
its types nowhere at all -- every field is text on disk. This module is the one
place that decides what each column becomes in Arrow, so the text reader, the
Parquet writer and the DuckDB type map cannot drift apart.

Types follow ``pairtools.lib.pairsam_format``, extended with the columns the
tools append (restriction fragments from ``restrict``, phases from ``phase``,
duplicate bookkeeping from ``dedup``).
"""

import pyarrow as pa
from pairtools.lib import pairsam_format

#: Columns appended by `pairtools restrict`, which pairsam_format does not list.
DTYPES_RFRAG_COLUMNS = {
    "rfrag1": int,
    "rfrag_start1": int,
    "rfrag_end1": int,
    "rfrag2": int,
    "rfrag_start2": int,
    "rfrag_end2": int,
}

#: Columns appended by `pairtools phase` and `pairtools dedup`.
DTYPES_ANNOTATION_COLUMNS = {
    "phase1": str,
    "phase2": str,
    "parent_readID": str,
}

#: Every column whose type we know, merged in order of increasing specificity.
DTYPES = dict(pairsam_format.DTYPES_PAIRSAM)
DTYPES.update(pairsam_format.DTYPES_EXTRA_COLUMNS)
DTYPES.update(DTYPES_RFRAG_COLUMNS)
DTYPES.update(DTYPES_ANNOTATION_COLUMNS)

#: Columns drawn from a small vocabulary, worth dictionary-encoding.
LOW_CARDINALITY_COLUMNS = ("chrom1", "chrom2", "strand1", "strand2", "pair_type")

#: Positions are genomic coordinates: int32 covers every known assembly.
INT_TYPE = pa.int32()


def declared_type(column):
    """The Python type pairtools declares for a column, or None.

    ``pairtools parse --add-columns`` names its columns by appending the side to
    the base name -- ``mapq`` becomes ``mapq1`` and ``mapq2`` -- while
    ``DTYPES_EXTRA_COLUMNS`` is keyed by the base name alone. Looking up the
    full name only, as this did, silently typed every one of those columns as a
    string, which is why ``select 'mapq1>=30'`` used to fail on a comparison
    between str and int.
    """
    if column in DTYPES:
        return DTYPES[column]
    if column[-1:] in ("1", "2"):
        return DTYPES.get(column[:-1])
    return None


def arrow_type(column, dict_encode=False):
    """Return the Arrow type for a single .pairs column.

    Unknown columns fall back to string, which is what the .pairs format
    guarantees and what a custom column produced by a third-party tool will be.
    """
    if dict_encode and column in LOW_CARDINALITY_COLUMNS:
        index_type = pa.int16() if column in ("chrom1", "chrom2") else pa.int8()
        return pa.dictionary(index_type, pa.string())

    if declared_type(column) is int:
        return INT_TYPE

    return pa.string()


def schema_from_columns(columns, dict_encode=False):
    """Build an Arrow schema for a list of .pairs column names.

    Parameters
    ----------
    columns : list of str
        Column names, in file order -- normally from
        ``headerops.extract_column_names(header)``.
    dict_encode : bool
        Dictionary-encode the low-cardinality columns. Off by default: Parquet
        already applies RLE_DICTIONARY encoding to these columns on its own, so
        this buys compute rather than space, and nothing consumes that yet.

    Returns
    -------
    pyarrow.Schema
    """
    return pa.schema([pa.field(col, arrow_type(col, dict_encode)) for col in columns])


def duckdb_types_from_columns(columns):
    """Return the DuckDB ``columns=`` type map matching :func:`schema_from_columns`.

    Used when DuckDB, rather than pyarrow, reads a .pairs text file.
    """
    return {
        col: "INTEGER" if declared_type(col) is int else "VARCHAR" for col in columns
    }

"""Lossless round-tripping of a .pairs header through Parquet key-value metadata.

A .pairs header is an ordered list of ``#key: value`` lines. Parquet stores
arbitrary bytes-to-bytes key-value metadata in its schema. The mapping used here
keeps both:

``pairs_header``
    The header, verbatim and in order, as a JSON array. This is the source of
    truth on the way back out, so a header round-trips byte-for-byte no matter
    which keys it uses -- including the ``@PG`` provenance chain that every tool
    appends and any key a third-party tool invented.

The parsed convenience keys (``columns``, ``chromsize``, ``samheader``, ...)
are written alongside it, unchanged from the 0.2.0 layout, so that
``parquet_kv_metadata()`` queries in DuckDB and files read by older versions
keep working.

Files written before ``pairs_header`` existed have only the parsed keys; reading
those falls back to reconstructing the header from them, which is lossy in
exactly the ways it always was (fixed key order, unknown keys dropped).
"""

import json

import pyarrow.parquet as pq

from . import header_metadata, json_transform

#: Metadata key holding the verbatim header lines.
HEADER_KEY = b"pairs_header"


def header_to_metadata(header):
    """Convert a .pairs header into Parquet key-value metadata.

    Parameters
    ----------
    header : list of str
        Header lines, as returned by ``headerops.get_header``.

    Returns
    -------
    dict of bytes to bytes
        Suitable for ``pyarrow.Schema.with_metadata``.
    """
    metadata = {HEADER_KEY: json.dumps(list(header)).encode("utf-8")}

    # The parsed keys are best-effort: a malformed or unusual header must not
    # stop us writing the file, because `pairs_header` already preserves it.
    try:
        field_names = header_metadata.extract_field_names(header)
        parsed = json_transform.json_dict_to_json_str(
            json_transform.header_to_json_dict(header, field_names)
        )
    except (IndexError, KeyError, ValueError):
        return metadata

    for key, value in parsed.items():
        key = key.encode("utf-8")
        if key != HEADER_KEY:
            metadata[key] = value.encode("utf-8")

    return metadata


def metadata_to_header(metadata):
    """Reconstruct a .pairs header from Parquet key-value metadata.

    Parameters
    ----------
    metadata : dict of bytes to bytes, or None
        As found on ``pyarrow.Schema.metadata``.

    Returns
    -------
    list of str
    """
    if not metadata:
        return []

    if HEADER_KEY in metadata:
        return json.loads(metadata[HEADER_KEY].decode("utf-8"))

    # Written before pairs_header existed: rebuild from the parsed keys.
    decoded = {
        key.decode("utf-8"): json_transform.decode_and_parse_json(
            value.decode("utf-8")
        )
        for key, value in metadata.items()
    }
    if "columns" not in decoded:
        raise ValueError(
            "Parquet file carries neither a 'pairs_header' nor a 'columns' "
            "metadata key, so its .pairs header cannot be recovered"
        )
    return header_metadata.metadata_dict_to_header_list(decoded)


def read_header(path):
    """Read the .pairs header from a Parquet file without reading its rows."""
    return metadata_to_header(pq.read_schema(path).metadata)


def apply_to_schema(schema, header):
    """Return `schema` carrying `header` as key-value metadata."""
    return schema.with_metadata(header_to_metadata(header))

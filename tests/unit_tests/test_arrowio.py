# -*- coding: utf-8 -*-
import gzip
import os
import sys

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."))

from pairtools_parquet.lib import headerio  # noqa: E402
from pairtools_parquet.lib.arrowio import (  # noqa: E402
    PairsWriter,
    is_parquet,
    open_pairs,
    read_header,
)
from pairtools_parquet.lib.schema import schema_from_columns  # noqa: E402

DATADIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "..", "data"
)
MOCK_PAIRS = os.path.join(DATADIR, "mock.pairs")
MOCK_PARQUET = os.path.join(DATADIR, "mock.parquet")

HEADER = [
    "## pairs format v1.0.0",
    "#shape: upper triangle",
    "#chromsize: chr1 1000",
    "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type",
]


def write_pairs(path, header, rows):
    with open(path, "w") as f:
        f.write("".join(l + "\n" for l in header))
        f.write("".join("\t".join(str(v) for v in r) + "\n" for r in rows))
    return str(path)


def read_all(reader):
    return [row for batch in reader for row in batch.to_pylist()]


def test_is_parquet():
    assert is_parquet("a.parquet")
    assert is_parquet("a.PQ")
    assert not is_parquet("a.pairs")
    assert not is_parquet("a.pairs.gz")


def test_pairs_to_parquet_to_pairs_is_byte_identical(tmp_path):
    """The full file, header included, must survive a Parquet round trip."""
    parquet_path = tmp_path / "out.parquet"
    pairs_path = tmp_path / "out.pairs"

    header, reader = open_pairs(MOCK_PAIRS)
    with PairsWriter(parquet_path, header) as writer:
        writer.write_all(reader)

    header2, reader2 = open_pairs(str(parquet_path))
    assert header2 == header
    with PairsWriter(pairs_path, header2, compress_program="none") as writer:
        writer.write_all(reader2)

    with open(MOCK_PAIRS, "rb") as f:
        original = f.read()
    with open(pairs_path, "rb") as f:
        assert f.read() == original


def test_gzip_round_trip(tmp_path):
    gz_path = tmp_path / "out.pairs.gz"

    header, reader = open_pairs(MOCK_PAIRS)
    with PairsWriter(gz_path, header, compress_program="gzip") as writer:
        writer.write_all(reader)

    with open(MOCK_PAIRS, "rb") as f:
        original = f.read()
    with gzip.open(gz_path, "rb") as f:
        assert f.read() == original

    header2, reader2 = open_pairs(str(gz_path))
    assert header2 == header
    assert len(read_all(reader2)) == 9


def test_header_round_trips_unknown_keys(tmp_path):
    """A header key this package has never heard of must survive."""
    header = HEADER + [
        "#custom_key: something a third-party tool wrote",
        "#samheader: @PG\tID:x\tCL:tool -o 'quoted' --path C:\\weird\\path",
    ]
    # `#columns:` must stay last only by convention; keep the real order intact
    header = [l for l in header if not l.startswith("#columns")] + [
        l for l in header if l.startswith("#columns")
    ]

    path = tmp_path / "custom.parquet"
    with PairsWriter(path, header) as writer:
        pass

    assert headerio.read_header(str(path)) == header


def test_legacy_parquet_without_header_key(tmp_path):
    """Files written before `pairs_header` existed fall back to parsed keys."""
    header = read_header(MOCK_PARQUET)
    assert header[0] == "## pairs format v1.0.0"
    assert any(l.startswith("#columns:") for l in header)

    schema = pq.read_schema(MOCK_PARQUET)
    assert headerio.HEADER_KEY not in (schema.metadata or {})


def test_metadata_to_header_rejects_unrecoverable(tmp_path):
    with pytest.raises(ValueError, match="cannot be recovered"):
        headerio.metadata_to_header({b"unrelated": b'"x"'})


def test_empty_input_produces_valid_output(tmp_path):
    source = write_pairs(tmp_path / "empty.pairs", HEADER, [])
    parquet_path = tmp_path / "empty.parquet"

    header, reader = open_pairs(source)
    with PairsWriter(parquet_path, header) as writer:
        writer.write_all(reader)

    assert pq.ParquetFile(str(parquet_path)).metadata.num_rows == 0
    assert headerio.read_header(str(parquet_path)) == header


def test_column_projection(tmp_path):
    for path in (MOCK_PAIRS, MOCK_PARQUET):
        _, reader = open_pairs(path, columns=["chrom1", "pos1"])
        assert reader.schema.names == ["chrom1", "pos1"]
        assert len(read_all(reader)) == 9


def test_values_that_look_like_nulls_survive(tmp_path):
    """.pairs has no null literals: "NA" and "null" are ordinary readIDs."""
    rows = [
        ("NA", "chr1", 1, "chr1", 2, "+", "+", "UU"),
        ("null", "chr1", 3, "chr1", 4, "+", "+", "UU"),
        ("", "chr1", 5, "chr1", 6, "+", "+", "UU"),
        ("unmapped", "!", 0, "!", 0, "-", "-", "WW"),
    ]
    source = write_pairs(tmp_path / "odd.pairs", HEADER, rows)

    header, reader = open_pairs(source)
    read_back = read_all(reader)
    assert [r["readID"] for r in read_back] == [r[0] for r in rows]
    assert read_back[3]["chrom1"] == "!"
    assert read_back[3]["pos1"] == 0

    out = tmp_path / "odd.pairs.out"
    header, reader = open_pairs(source)
    with PairsWriter(out, header, compress_program="none") as writer:
        writer.write_all(reader)
    with open(source, "rb") as a, open(out, "rb") as b:
        assert a.read() == b.read()


def test_quote_in_value_reads_but_cannot_be_written_as_text(tmp_path):
    """Known limitation, asserted so it fails loudly rather than silently.

    pyarrow's CSV writer rejects a structural character when quoting is off,
    and .pairs is unquoted. Reading is fine and Parquet output is fine; only
    text output refuses. pairtools itself writes such a value without complaint.
    """
    rows = [('read"with"quotes', "chr1", 5, "chr1", 6, "+", "+", "UU")]
    source = write_pairs(tmp_path / "quoted.pairs", HEADER, rows)

    header, reader = open_pairs(source)
    assert read_all(reader)[0]["readID"] == 'read"with"quotes'

    # Parquet is unaffected
    header, reader = open_pairs(source)
    with PairsWriter(tmp_path / "quoted.parquet", header) as writer:
        writer.write_all(reader)
    assert (
        pq.read_table(str(tmp_path / "quoted.parquet"))["readID"][0].as_py()
        == 'read"with"quotes'
    )

    header, reader = open_pairs(source)
    with pytest.raises(pa.ArrowInvalid, match="structural characters"):
        with PairsWriter(
            tmp_path / "quoted.out.pairs", header, compress_program="none"
        ) as writer:
            writer.write_all(reader)


def test_empty_fields_stay_empty_strings_in_parquet(tmp_path):
    """.pairs has no NULL: an empty field is the empty string.

    DuckDB's read_csv reads an empty field as NULL by default, which turned
    e.g. an empty XB tag into None on conversion -- `phase` reads that tag as
    "no alternative alignment" and calls .split() on it.
    """
    import subprocess
    import sys

    columns = ["readID", "chrom1", "pos1", "chrom2", "pos2", "strand1",
               "strand2", "pair_type", "XB1"]
    header = ["## pairs format v1.0.0", "#chromsize: chr1 1000",
              "#columns: " + " ".join(columns)]
    rows = [
        ("r1", "chr1", 1, "chr1", 2, "+", "+", "UU", ""),
        ("r2", "chr1", 3, "chr1", 4, "+", "+", "UU", "chr1,+5,100M,0,90,60;"),
    ]
    source = write_pairs(tmp_path / "in.pairs", header, rows)
    out = tmp_path / "out.parquet"

    subprocess.run(
        [sys.executable, "-m", "pairtools_parquet", "csv-to-parquet",
         "-o", str(out), str(source)],
        check=True, capture_output=True,
    )

    values = pq.read_table(str(out)).column("XB1").to_pylist()
    assert values == ["", "chr1,+5,100M,0,90,60;"]
    assert None not in values


def test_header_keeps_trailing_whitespace(tmp_path):
    """Only the newline is stripped from a header line, not trailing spaces.

    `pairtools header generate` with no --assembly emits `#genome_assembly: `
    with a trailing space, and `headerops.get_header` keeps it, so stripping it
    on the way out would break byte parity for a file pairtools wrote.
    """
    header = ["## pairs format v1.0.0", "#genome_assembly: "] + HEADER[2:]
    source = write_pairs(tmp_path / "spaced.pairs", header, [])
    out = tmp_path / "out.pairs"

    read_back, reader = open_pairs(source)
    assert read_back == header
    with PairsWriter(out, read_back, compress_program="none") as writer:
        writer.write_all(reader)

    with open(source, "rb") as a, open(out, "rb") as b:
        assert a.read() == b.read()


def test_column_names_override(tmp_path):
    """`column_names` names the columns positionally, header or no header."""
    rows = [("r1", "chr1", 1, "chr1", 2, "+", "+", "UU")]
    names = ["readID", "chrom1", "pos5", "chrom2", "pos3", "strand1", "strand2",
             "pair_type"]

    headerless = write_pairs(tmp_path / "headerless.pairs", [], rows)
    _, reader = open_pairs(headerless, column_names=names)
    assert reader.schema.names == names
    # the type follows the new name, not the old one
    assert reader.schema.field("pos5").type == pa.int32()
    assert read_all(reader)[0]["pos5"] == 1

    # the same override applies to Parquet, whose columns are already named
    _, reader = open_pairs(MOCK_PARQUET, column_names=names)
    assert reader.schema.names == names
    assert len(read_all(reader)) == 9


def test_column_names_must_match_the_column_count(tmp_path):
    with pytest.raises(ValueError, match="8 columns, but 2 names"):
        open_pairs(MOCK_PARQUET, column_names=["a", "b"])


def test_writer_casts_mismatched_batches(tmp_path):
    """A caller may hand over batches whose types differ from the schema."""
    header = HEADER
    columns = ["readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2",
               "pair_type"]
    wide = pa.schema(
        [pa.field(c, pa.int64() if c.startswith("pos") else pa.string())
         for c in columns]
    )
    batch = pa.RecordBatch.from_pylist(
        [{"readID": "r1", "chrom1": "chr1", "pos1": 1, "chrom2": "chr1", "pos2": 2,
          "strand1": "+", "strand2": "+", "pair_type": "UU"}],
        schema=wide,
    )

    path = tmp_path / "cast.parquet"
    with PairsWriter(path, header, schema=schema_from_columns(columns)) as writer:
        writer.write(batch)

    table = pq.read_table(str(path))
    assert table.schema.field("pos1").type == pa.int32()
    assert table.num_rows == 1

# -*- coding: utf-8 -*-
"""header against pairtools header.

The four subcommands change the header and pass the body through, so parity
means the whole file matching -- header and body both -- apart from the @PG
record every tool stamps with its own name.
"""

import pytest

from conftest import (
    read_pairs_body,
    read_pairs_header,
    read_parquet_body,
    read_parquet_header,
    run_cli,
    run_pairtools,
)


@pytest.fixture
def headerless(tmp_path, mock_pairs_path):
    """The mock pairs with the header stripped off, as `generate` expects."""
    path = tmp_path / "headerless.pairs"
    with open(path, "w") as f:
        f.write("".join(l + "\n" for l in read_pairs_body(mock_pairs_path)))
    return path


def read(path):
    with open(path) as f:
        return f.read()


def without_pg(lines):
    """Header lines other than the @PG records, which name the tool that ran."""
    return [l for l in lines if "@PG" not in l]


def assert_matches(reference, ours):
    assert without_pg(read_pairs_header(reference)) == without_pg(
        read_pairs_header(ours)
    )
    assert read_pairs_body(reference) == read_pairs_body(ours)


def test_generate_matches_pairtools(tmp_path, headerless, mock_chromsizes_path):
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"

    run_pairtools(
        "header", "generate", "--chroms-path", mock_chromsizes_path,
        "-o", reference, headerless
    )
    run_cli(
        "header", "generate", "--chroms-path", mock_chromsizes_path,
        "-o", ours, headerless
    )

    # generate adds no @PG record when the input has no samheader, so these are
    # byte-identical rather than merely equivalent
    assert read(reference) == read(ours)


def test_generate_options_match_pairtools(tmp_path, headerless, mock_chromsizes_path):
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"
    options = ["--assembly", "hg38", "--no-flip"]

    run_pairtools(
        "header", "generate", "--chroms-path", mock_chromsizes_path,
        "-o", reference, *options, headerless
    )
    run_cli(
        "header", "generate", "--chroms-path", mock_chromsizes_path,
        "-o", ours, *options, headerless
    )

    assert read(reference) == read(ours)


def test_generate_needs_chromsizes(tmp_path, headerless):
    with pytest.raises(AssertionError, match="chroms-path"):
        run_cli("header", "generate", "-o", tmp_path / "out.pairs", headerless)


def test_generate_rejects_wrong_column_count(
    tmp_path, headerless, mock_chromsizes_path
):
    with pytest.raises(AssertionError, match="columns mismatch"):
        run_cli(
            "header", "generate", "--chroms-path", mock_chromsizes_path,
            "--columns", "a,b,c", "-o", tmp_path / "out.pairs", headerless
        )


def test_transfer_matches_pairtools(tmp_path, headerless, mock_pairs_path):
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"

    run_pairtools(
        "header", "transfer", "-r", mock_pairs_path, "-o", reference, headerless
    )
    run_cli("header", "transfer", "-r", mock_pairs_path, "-o", ours, headerless)

    assert_matches(reference, ours)


def test_set_columns_matches_pairtools(tmp_path, mock_pairs_path):
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"
    columns = "readID,chrom1,pos1,chrom2,pos2,strand1,strand2,pair_type"

    run_pairtools(
        "header", "set-columns", "-c", columns, "-o", reference, mock_pairs_path
    )
    run_cli("header", "set-columns", "-c", columns, "-o", ours, mock_pairs_path)

    assert read(reference) == read(ours)


def test_set_columns_renames(tmp_path, mock_pairs_path):
    """The names really change, and the body follows them."""
    ours = tmp_path / "ours.parquet"
    columns = "readID,chrom1,pos5,chrom2,pos3,strand1,strand2,pair_type"
    run_cli("header", "set-columns", "-c", columns, "-o", ours, mock_pairs_path)

    import pyarrow.parquet as pq

    assert pq.read_schema(str(ours)).names == columns.split(",")
    assert read_parquet_body(ours) == read_pairs_body(mock_pairs_path)


def test_set_columns_adds_a_missing_columns_line(tmp_path, headerless):
    """Upstream leaves a headerless file headerless; we add the line.

    `headerops.set_columns` only rewrites a `#columns:` line that already
    exists, so upstream's set-columns is a no-op on exactly the input it is
    meant for.
    """
    ours = tmp_path / "ours.pairs"
    run_cli(
        "header", "set-columns", "-c", "a,b,c,d,e,f,g,h", "-o", ours, headerless
    )

    assert read_pairs_header(ours) == ["#columns: a b c d e f g h"]
    assert read_pairs_body(ours) == read_pairs_body(headerless)


def test_validate_columns_matches_pairtools(tmp_path, mock_pairs_path):
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"

    run_pairtools("header", "validate-columns", "-o", reference, mock_pairs_path)
    run_cli("header", "validate-columns", "-o", ours, mock_pairs_path)

    assert_matches(reference, ours)


def test_validate_columns_rejects_a_mismatched_reference(tmp_path, mock_pairs_path):
    with pytest.raises(AssertionError, match="differ from reference"):
        run_cli(
            "header", "validate-columns", "-c", "a,b,c",
            "-o", tmp_path / "out.pairs", mock_pairs_path
        )


def test_headers_cross_formats(tmp_path, headerless, mock_chromsizes_path):
    """A header generated onto Parquet transfers back onto a text file."""
    generated = tmp_path / "generated.parquet"
    run_cli(
        "header", "generate", "--chroms-path", mock_chromsizes_path,
        "-o", generated, headerless
    )

    transferred = tmp_path / "transferred.pairs"
    run_cli("header", "transfer", "-r", generated, "-o", transferred, headerless)

    assert without_pg(read_parquet_header(generated)) == without_pg(
        read_pairs_header(transferred)
    )
    assert read_parquet_body(generated) == read_pairs_body(transferred)


def test_validate_columns_on_parquet(tmp_path, mock_parquet_path):
    ours = tmp_path / "ours.parquet"
    run_cli("header", "validate-columns", "-o", ours, mock_parquet_path)

    assert read_parquet_body(ours) == read_parquet_body(mock_parquet_path)

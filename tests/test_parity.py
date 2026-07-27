# -*- coding: utf-8 -*-
"""Parity of pairtools_parquet against pairtools itself.

The whole premise of this package is "same results, faster", so every tool is
checked against the upstream implementation rather than against a hand-written
expectation: run both, convert whatever came out back to text, and require the
bodies to be identical. The header may differ only by the @PG provenance record
each tool appends to record itself.
"""

import pytest

from conftest import (
    read_pairs_body,
    read_pairs_header,
    read_parquet_body,
    run_cli,
    run_pairtools,
    write_pairs,
)

SORT_KEYS = ["chrom1", "chrom2", "pos1", "pos2", "pair_type"]


def assert_headers_agree(reference_header, our_header):
    """Headers must match except for each tool's own @PG record."""

    def without_provenance(header):
        return [
            l
            for l in header
            if not (l.startswith("#samheader: @PG") and "pairtools" in l and "sort" in l)
        ]

    assert without_provenance(reference_header) == without_provenance(our_header)


def test_sort_matches_pairtools(tmp_path, mock_pairs_path):
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"

    run_pairtools("sort", "-o", reference, mock_pairs_path)
    run_cli("sort", "-o", ours, "--compress-program", "none", mock_pairs_path)

    assert read_pairs_body(reference) == read_pairs_body(ours)
    assert_headers_agree(read_pairs_header(reference), read_pairs_header(ours))


def test_sort_matches_pairtools_through_parquet(tmp_path, mock_pairs_path):
    """Routing through Parquet must not change the result."""
    reference = tmp_path / "ref.pairs"
    intermediate = tmp_path / "sorted.parquet"
    ours = tmp_path / "ours.pairs"

    run_pairtools("sort", "-o", reference, mock_pairs_path)
    run_cli("sort", "-o", intermediate, mock_pairs_path)
    run_cli("sort", "-o", ours, "--compress-program", "none", intermediate)

    assert read_pairs_body(reference) == read_pairs_body(ours)
    assert read_pairs_body(reference) == read_parquet_body(intermediate)


def test_sort_orders_tied_pairs_by_pair_type(tmp_path):
    """Pairs identical but for pair_type must order lexicographically.

    Regression test: chrom/strand/pair_type used to be sorted as DuckDB ENUMs
    declared in itertools.product order, so `UU` sorted before `DD` where
    pairtools -- which sorts bytewise -- puts `DD` first. pair_type is the last
    sort key, so this only shows up on pairs tied on every position.
    """
    header = [
        "## pairs format v1.0.0",
        "#chromsize: chr1 1000",
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type",
    ]
    rows = [
        ("r_uu", "chr1", 1, "chr1", 50, "+", "+", "UU"),
        ("r_dd", "chr1", 1, "chr1", 50, "+", "+", "DD"),
        ("r_ww", "chr1", 1, "chr1", 50, "+", "+", "WW"),
        ("r_nu", "chr1", 1, "chr1", 50, "+", "+", "NU"),
    ]
    source = write_pairs(tmp_path / "tied.pairs", header, rows)

    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"
    run_pairtools("sort", "-o", reference, source)
    run_cli("sort", "-o", ours, "--compress-program", "none", source)

    body = read_pairs_body(ours)
    assert [l.split("\t")[7] for l in body] == ["DD", "NU", "UU", "WW"]
    assert read_pairs_body(reference) == body


def test_sort_handles_chromosome_missing_from_header(tmp_path):
    """A chromosome the header never declares must not be fatal.

    The chromosome ENUM that makes the sort fast has a fixed domain taken from
    the header, and DuckDB raises rather than nulling on a value outside it, so
    this has to fall back to comparing strings.
    """
    header = [
        "## pairs format v1.0.0",
        "#chromsize: chr1 1000",
        "#chromsize: chr2 1000",
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type",
    ]
    rows = [
        ("r1", "chr2", 5, "chr2", 9, "+", "+", "UU"),
        ("r2", "chrUNDECLARED", 3, "chr1", 7, "+", "-", "UU"),
        ("r3", "chr1", 1, "chr1", 2, "-", "+", "UU"),
    ]
    source = write_pairs(tmp_path / "undeclared.pairs", header, rows)

    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"
    run_pairtools("sort", "-o", reference, source)
    run_cli("sort", "-o", ours, "--compress-program", "none", source)

    assert read_pairs_body(reference) == read_pairs_body(ours)


@pytest.mark.parametrize("compress_program", ["none", "gzip"])
def test_sort_to_compressed_output(tmp_path, mock_pairs_path, compress_program):
    suffix = ".pairs" if compress_program == "none" else ".pairs.gz"
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / ("ours" + suffix)

    run_pairtools("sort", "-o", reference, mock_pairs_path)
    run_cli(
        "sort", "-o", ours, "--compress-program", compress_program, mock_pairs_path
    )

    # reading back through our own reader exercises the decompression path
    from pairtools_parquet.lib.arrowio import open_pairs

    _, reader = open_pairs(str(ours))
    rows = sum(batch.num_rows for batch in reader)
    assert rows == len(read_pairs_body(reference))

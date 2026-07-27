# -*- coding: utf-8 -*-
from pairtools.lib import pairsam_format

from conftest import read_parquet_body, read_parquet_header, run_cli


def test_preserve(tmp_path, mock_parquet_path):
    output_path = tmp_path / "select_preserve_mock.parquet"
    run_cli("select", "True", mock_parquet_path, "-o", output_path)

    original_body = read_parquet_body(mock_parquet_path)
    output_body = read_parquet_body(output_path)

    # a "True" condition must be a no-op, not an empty result
    assert sorted(output_body) == sorted(original_body)


def test_equal(tmp_path, mock_parquet_path):
    output_path = tmp_path / "select_pairTypeEqual_mock.parquet"
    run_cli(
        "select",
        '(pair_type == "RU") or (pair_type == "UR") or (pair_type == "UU")',
        mock_parquet_path,
        "-o",
        output_path,
    )

    original_body = read_parquet_body(mock_parquet_path)
    output_body = read_parquet_body(output_path)

    assert all(l.split("\t")[7] in ["RU", "UR", "UU"] for l in output_body)
    assert all(
        l in output_body for l in original_body if l.split("\t")[7] in ["RU", "UR", "UU"]
    )


def test_csv(tmp_path, mock_parquet_path):
    output_path = tmp_path / "select_csvMatch_mock.parquet"
    run_cli(
        "select",
        'csv_match(pair_type, "RU,UR,UU")',
        mock_parquet_path,
        "-o",
        output_path,
    )

    original_body = read_parquet_body(mock_parquet_path)
    output_body = read_parquet_body(output_path)

    assert all(l.split("\t")[7] in ["RU", "UR", "UU"] for l in output_body)
    assert all(
        l in output_body for l in original_body if l.split("\t")[7] in ["RU", "UR", "UU"]
    )


def test_wildcard(tmp_path, mock_parquet_path):
    output_path = tmp_path / "select_wildcard_mock.parquet"
    run_cli(
        "select",
        'wildcard_match(pair_type, "*U")',
        mock_parquet_path,
        "-o",
        output_path,
    )

    original_body = read_parquet_body(mock_parquet_path)
    output_body = read_parquet_body(output_path)

    assert all(l.split("\t")[7] in ["NU", "MU", "RU", "UU"] for l in output_body)
    assert all(
        l in output_body
        for l in original_body
        if l.split("\t")[7] in ["NU", "MU", "RU", "UU"]
    )


def test_regex(tmp_path, mock_parquet_path):
    output_path = tmp_path / "select_regex_mock.parquet"
    run_cli(
        "select",
        'regex_match(pair_type, "[NM]U")',
        mock_parquet_path,
        "-o",
        output_path,
    )

    original_body = read_parquet_body(mock_parquet_path)
    output_body = read_parquet_body(output_path)

    assert all(l.split("\t")[7] in ["NU", "MU"] for l in output_body)
    assert all(
        l in output_body for l in original_body if l.split("\t")[7] in ["NU", "MU"]
    )


def test_chrom_subset(tmp_path, mock_parquet_path, mock_chromsizes_path):
    output_path = tmp_path / "select_chromSubset_mock.parquet"
    run_cli(
        "select",
        "True",
        "--chrom-subset",
        mock_chromsizes_path,
        mock_parquet_path,
        "-o",
        output_path,
    )

    output_header = read_parquet_header(output_path)

    chroms_from_chrom_field = [
        l.split()[1:] for l in output_header if l.startswith("#chromosomes:")
    ][0]

    assert set(chroms_from_chrom_field) == set(["chr1", "chr2"])

    chroms_from_chrom_sizes = [
        l.split()[1] for l in output_header if l.startswith("#chromsize:")
    ]

    assert set(chroms_from_chrom_sizes) == set(["chr1", "chr2"])

    # the body must be restricted to the same subset as the header
    for l in read_parquet_body(output_path):
        fields = l.split("\t")
        assert fields[1] in ["chr1", "chr2"]
        assert fields[3] in ["chr1", "chr2"]


def test_remove_columns(tmp_path, mock_parquet_path):
    """Test removal of columns from the file
    Example run:
    pairtools_parquet select True --remove-columns sam1,sam2 tests/data/mock.parquet
    """
    output_path = tmp_path / "select_remove_columns_mock.parquet"
    run_cli(
        "select",
        "True",
        "--remove-columns",
        "sam1,sam2",
        mock_parquet_path,
        "-o",
        output_path,
    )

    # check if the columns are removed properly:
    output_header = read_parquet_header(output_path)
    output_body = read_parquet_body(output_path)

    for l in output_header:
        if l.startswith("#columns:"):
            assert (
                l
                == "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type"
            )

    # check that the pairs got assigned properly
    for l in output_body:
        assert len(l.split(pairsam_format.PAIRSAM_SEP)) == 8


def test_region_match(tmp_path, mock_parquet_path):
    output_path = tmp_path / "select_regionMatch_mock.parquet"
    run_cli(
        "select",
        'region_match(chrom1, pos1, "chr1", 0, 50)',
        mock_parquet_path,
        "-o",
        output_path,
    )

    original_body = read_parquet_body(mock_parquet_path)
    output_body = read_parquet_body(output_path)

    # Verify all output rows have chrom1="chr1" and pos1 within range
    for l in output_body:
        fields = l.split("\t")
        chrom1, pos1 = fields[1], int(fields[2])
        assert chrom1 == "chr1"
        assert 0 <= pos1 <= 50

    # Verify all matching rows from input are in output
    for l in original_body:
        fields = l.split("\t")
        chrom1, pos1 = fields[1], int(fields[2])
        if chrom1 == "chr1" and 0 <= pos1 <= 50:
            assert l in output_body


def test_region_match_no_end(tmp_path, mock_parquet_path):
    # Every chr1 row in the mock data sits at pos1 == 1, so the start must be
    # low enough to actually select something -- otherwise both loops below
    # iterate over nothing and the test passes without checking anything.
    start = 1
    output_path = tmp_path / "select_region_match_no_end_mock.parquet"
    run_cli(
        "select",
        'region_match(chrom1, pos1, "chr1", {})'.format(start),
        mock_parquet_path,
        "-o",
        output_path,
    )

    original_body = read_parquet_body(mock_parquet_path)
    output_body = read_parquet_body(output_path)

    def matches(line):
        fields = line.split("\t")
        return fields[1] == "chr1" and int(fields[2]) >= start

    expected = [l for l in original_body if matches(l)]
    assert len(expected) > 0, "the region must select some rows to be a useful test"

    # the open-ended region has no upper bound, so it selects exactly the rows
    # on chr1 at or past the start
    assert sorted(output_body) == sorted(expected)

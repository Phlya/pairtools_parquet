# -*- coding: utf-8 -*-
from conftest import (
    read_pairs_body,
    read_pairs_header,
    read_parquet_body,
    read_parquet_header,
    run_cli,
)


def assert_sorted(body):
    """Assert chrom1 -> chrom2 -> pos1 -> pos2 ordering over tab-separated lines."""
    prev_pair = None
    for l in body:
        cur_pair = l.split("\t")[1:8]
        if prev_pair is not None:
            assert cur_pair[0] >= prev_pair[0]
            if cur_pair[0] == prev_pair[0]:
                assert cur_pair[2] >= prev_pair[2]
                if cur_pair[2] == prev_pair[2]:
                    assert int(cur_pair[1]) >= int(prev_pair[1])
                    if int(cur_pair[1]) == int(prev_pair[1]):
                        assert int(cur_pair[3]) >= int(prev_pair[3])

        prev_pair = cur_pair


def assert_header_only_gained_provenance(input_header, output_header):
    """Assert the output header only adds a @PG record, #sorted and #chromosomes."""
    for l in output_header:
        if not any([l in l2 for l2 in input_header]):
            assert (
                l.startswith("#samheader: @PG")
                or l.startswith("#sorted")
                or l.startswith("#chromosomes")
            )


def test_mock_pairs(tmp_path, mock_pairs_path):
    output_path = tmp_path / "sorted_mock.pairs"
    run_cli("sort", "-o", output_path, "--compress-program", "none", mock_pairs_path)

    assert_header_only_gained_provenance(
        read_pairs_header(mock_pairs_path), read_pairs_header(output_path)
    )

    input_body = read_pairs_body(mock_pairs_path)
    output_body = read_pairs_body(output_path)

    # check that all pairs entries survived sorting:
    assert len(input_body) == len(output_body)
    assert sorted(input_body) == sorted(output_body)

    assert_sorted(output_body)


def test_mock_pairs_parquet(tmp_path, mock_pairs_path):
    output_path = tmp_path / "sorted_mock.parquet"
    run_cli("sort", "-o", output_path, "--compress-program", "none", mock_pairs_path)

    assert_header_only_gained_provenance(
        read_pairs_header(mock_pairs_path), read_parquet_header(output_path)
    )

    input_body = read_pairs_body(mock_pairs_path)
    output_body = read_parquet_body(output_path)

    assert len(input_body) == len(output_body)
    assert sorted(input_body) == sorted(output_body)

    assert_sorted(output_body)


def test_pairs_and_parquet_outputs_agree(tmp_path, mock_pairs_path):
    """Sorting the same input to either format must produce the same rows."""
    pairs_output = tmp_path / "sorted_mock.pairs"
    parquet_output = tmp_path / "sorted_mock.parquet"

    run_cli("sort", "-o", pairs_output, "--compress-program", "none", mock_pairs_path)
    run_cli("sort", "-o", parquet_output, "--compress-program", "none", mock_pairs_path)

    assert read_pairs_body(pairs_output) == read_parquet_body(parquet_output)

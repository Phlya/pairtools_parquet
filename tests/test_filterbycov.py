# -*- coding: utf-8 -*-
"""filterbycov against pairtools filterbycov."""

import pytest

from conftest import read_pairs_body, read_parquet_body, run_cli, run_pairtools
from test_dedup import make_duplicated_pairs


@pytest.fixture
def sorted_pairs(tmp_path):
    """Pairs with clustered positions, so some regions exceed the coverage cap."""
    raw = make_duplicated_pairs(tmp_path / "raw.pairs")
    path = tmp_path / "sorted.pairs"
    run_pairtools("sort", "-o", path, raw)
    return path


def read(path):
    with open(path) as f:
        return f.read()


OPTIONS = [
    [],
    ["--max-cov", "2"],
    ["--max-cov", "20", "--max-dist", "5000"],
    ["--method", "sum"],
    ["--mark-multi"],
]


@pytest.mark.parametrize("options", OPTIONS, ids=lambda o: " ".join(o) or "defaults")
def test_filterbycov_matches_pairtools(tmp_path, sorted_pairs, options):
    ref = {p: tmp_path / "ref_{}.pairs".format(p) for p in ("low", "high", "unmapped")}
    ours = {p: tmp_path / "our_{}.pairs".format(p) for p in ref}

    run_pairtools(
        "filterbycov", "--output", ref["low"], "--output-highcov", ref["high"],
        "--output-unmapped", ref["unmapped"], *options, sorted_pairs
    )
    run_cli(
        "filterbycov", "-o", ours["low"], "--output-highcov", ours["high"],
        "--output-unmapped", ours["unmapped"], *options, sorted_pairs
    )

    for part in ref:
        assert read_pairs_body(ref[part]) == read_pairs_body(ours[part]), part
    # some pairs must actually be filtered, or this proves nothing
    assert len(read_pairs_body(ours["high"])) > 0


def test_filterbycov_stats_match_pairtools(tmp_path, sorted_pairs):
    """Counted in file order, not grouped by outcome.

    The counter records pair types in first-seen order, so feeding the
    low-coverage pairs before the high-coverage ones would reorder the output
    even though every count matched.
    """
    reference = tmp_path / "ref.stats"
    ours = tmp_path / "ours.stats"

    run_pairtools(
        "filterbycov", "--output", tmp_path / "r.pairs", "--output-stats", reference,
        sorted_pairs
    )
    run_cli(
        "filterbycov", "-o", tmp_path / "o.pairs", "--output-stats", ours, sorted_pairs
    )

    assert read(reference) == read(ours)


def test_filterbycov_through_parquet(tmp_path, sorted_pairs):
    as_parquet = tmp_path / "in.parquet"
    run_cli("csv-to-parquet", "-o", as_parquet, sorted_pairs)

    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.parquet"
    run_pairtools("filterbycov", "--output", reference, sorted_pairs)
    run_cli("filterbycov", "-o", ours, as_parquet)

    assert read_pairs_body(reference) == read_parquet_body(ours)


def test_outputs_partition_the_input(tmp_path, sorted_pairs):
    low = tmp_path / "low.pairs"
    high = tmp_path / "high.pairs"
    unmapped = tmp_path / "unmapped.pairs"
    run_cli(
        "filterbycov", "-o", low, "--output-highcov", high,
        "--output-unmapped", unmapped, sorted_pairs
    )

    total = (
        len(read_pairs_body(low))
        + len(read_pairs_body(high))
        + len(read_pairs_body(unmapped))
    )
    assert total == len(read_pairs_body(sorted_pairs))

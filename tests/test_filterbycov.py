# -*- coding: utf-8 -*-
"""filterbycov against pairtools filterbycov."""

import pytest

from conftest import read_pairs_body, read_parquet_body, run_cli, run_pairtools
from test_dedup import make_duplicated_pairs


@pytest.fixture(scope="session")
def sorted_pairs(tmp_path_factory):
    """Pairs with clustered positions, so some regions exceed the coverage cap.

    Built once for the session; every test only reads it.
    """
    directory = tmp_path_factory.mktemp("filterbycov_input")
    raw = make_duplicated_pairs(directory / "raw.pairs")
    path = directory / "sorted.pairs"
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


@pytest.mark.parametrize("options", OPTIONS, ids=lambda o: " ".join(o) or "defaults")
def test_duckdb_backend_matches_the_python_one(tmp_path, sorted_pairs, options):
    """The bucketed equi-join must agree with pairtools' double loop.

    The python backend calls `_filterbycov` unchanged, so this is the same
    comparison as the pairtools parity tests -- but it runs on the two backends
    of one implementation, which is what catches a divergence introduced on the
    duckdb side alone.
    """
    outputs = {}
    for backend in ("duckdb", "python"):
        low = tmp_path / "{}_low.pairs".format(backend)
        high = tmp_path / "{}_high.pairs".format(backend)
        run_cli(
            "filterbycov", "--backend", backend, "-o", low,
            "--output-highcov", high, *options, sorted_pairs
        )
        outputs[backend] = (read_pairs_body(low), read_pairs_body(high))

    assert outputs["duckdb"] == outputs["python"]
    assert outputs["duckdb"][0], "the fixture should leave some low-coverage pairs"


def test_coverage_counts_the_pairs_other_end(tmp_path):
    """A pair whose own two ends are close counts each against the other.

    An easy thing to get wrong when the neighbour search is a self-join: the
    obvious "exclude this pair's rows" is wrong, only "exclude this *end*" is.
    `_filterbycov` counts the other end because it works on a flat list of ends
    with no memory of which pair they came from.
    """
    import numpy as np
    import pyarrow as pa

    from pairtools_parquet.lib.filterbycov import _ends_table, coverage_duckdb

    lone = pa.table(
        {
            "chrom1": pa.array(["chr1"]),
            "pos1": pa.array([1000]),
            "chrom2": pa.array(["chr1"]),
            "pos2": pa.array([1100]),
        }
    )
    ends = _ends_table(lone, "chrom1", "pos1", "chrom2", "pos2")

    # 100bp apart: each end sees the other, so each side counts 2.
    assert coverage_duckdb(ends, 500, "max").tolist() == [2]
    assert coverage_duckdb(ends, 500, "sum").tolist() == [3]
    # 100bp apart is outside a 50bp window, so neither end sees anything.
    assert coverage_duckdb(ends, 50, "max").tolist() == [1]
    assert coverage_duckdb(ends, 50, "sum").tolist() == [1]

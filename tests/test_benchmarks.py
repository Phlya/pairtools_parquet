# -*- coding: utf-8 -*-
"""The benchmark harness itself.

A benchmark nobody runs rots, and this one has two ways to rot quietly: the
runner can stop working, and the generated dataset can drift into something
that does not exercise what it claims to -- pairs with no duplicates left in
them would still produce a table, just a meaningless one. So the harness gets
run end to end at a size where it takes seconds, and the fixture is checked for
the properties the benchmark depends on.

The timings themselves are not asserted on. A test that fails when a shared CI
runner is busy is worse than no test.
"""

import json
import os
import subprocess
import sys

import pytest

TESTDIR = os.path.dirname(os.path.realpath(__file__))
ROOT = os.path.dirname(TESTDIR)
BENCHMARKS = os.path.join(ROOT, "benchmarks")
sys.path.insert(0, BENCHMARKS)

import make_data  # noqa: E402

from conftest import read_pairs_body, run_pairtools  # noqa: E402


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    """A dataset small enough to build in a couple of seconds."""
    directory = tmp_path_factory.mktemp("benchmark-data")
    return make_data.build(str(directory), n_pairs=4000, dup_rate=0.2,
                           log=lambda *a: None)


def test_generated_pairs_still_contain_duplicates(dataset, tmp_path):
    """The whole point of generating the data rather than downloading it.

    `dedup` is the benchmark's headline tool, so a fixture whose pairs are
    already unique would make its column measure nothing.
    """
    stats = tmp_path / "dedup.stats"
    run_pairtools("dedup", "--output-stats", stats, "-o", os.devnull,
                  dataset.pairs)
    reported = dict(
        line.split("\t")[:2] for line in stats.read_text().splitlines()
        if "\t" in line
    )

    duplicates = int(reported["total_dups"])
    mapped = int(reported["total_mapped"])
    assert 0.1 < duplicates / mapped < 0.35, (
        "asked for 20% duplicates, pairtools found {:.1%}".format(
            duplicates / mapped)
    )


def test_generated_pairs_are_sorted_and_mixed(dataset):
    """`dedup` and `filterbycov` require sorted input; `stats` wants both
    cis and trans pairs and something unmapped to skip."""
    rows = [line.split("\t") for line in read_pairs_body(dataset.pairs)]
    keys = [(r[1], r[3], int(r[2]), int(r[4])) for r in rows]
    assert keys == sorted(keys)

    pair_types = {r[7] for r in rows}
    assert "UU" in pair_types and len(pair_types) > 1
    assert any(r[1] == r[3] for r in rows) and any(r[1] != r[3] for r in rows)


def test_the_dataset_is_reused_rather_than_rebuilt(dataset):
    """Regenerating is slow, so the manifest has to actually match."""
    assert dataset.is_current()
    assert make_data.Dataset(dataset.directory, dataset.n_pairs + 1,
                             dataset.dup_rate, dataset.seed).is_current() is False


def test_the_harness_runs_and_checks_its_outputs(dataset, tmp_path):
    """End to end, on the tools that cover each shape of benchmark.

    `parse` reads the BAM, `merge` takes two inputs, `stats` writes something
    that is not pairs, and `dedup` is the one with three outputs.
    """
    results = tmp_path / "results.json"
    proc = subprocess.run(
        [sys.executable, os.path.join(BENCHMARKS, "run.py"),
         "-d", dataset.directory, "-n", str(dataset.n_pairs),
         "--dup-rate", str(dataset.dup_rate), "--seed", str(dataset.seed),
         "-t", "parse,merge,stats,dedup", "--json", str(results)],
        capture_output=True,
    )
    assert proc.returncode == 0, proc.stdout.decode() + proc.stderr.decode()

    reported = {r["name"]: r for r in json.loads(results.read_text())}
    assert set(reported) == {"parse", "merge", "stats", "dedup"}
    for name, result in reported.items():
        assert result["formats_agree"] is True, (
            "{}: text and Parquet input gave different answers".format(name))
        assert result["reference"] > 0 and result["parquet"] > 0


def test_unknown_tool_is_rejected(dataset):
    proc = subprocess.run(
        [sys.executable, os.path.join(BENCHMARKS, "run.py"),
         "-d", dataset.directory, "-t", "nosuchtool"],
        capture_output=True,
    )
    assert proc.returncode != 0
    assert "unknown tools: nosuchtool" in proc.stderr.decode()

# -*- coding: utf-8 -*-
"""stats against pairtools stats."""

import pytest

from conftest import run_cli, run_pairtools
from test_dedup import make_duplicated_pairs


@pytest.fixture
def sorted_pairs(tmp_path):
    raw = make_duplicated_pairs(tmp_path / "raw.pairs")
    path = tmp_path / "sorted.pairs"
    run_pairtools("sort", "-o", path, raw)
    return path


@pytest.fixture
def deduped_with_parents(tmp_path, sorted_pairs):
    """Duplicates carrying parent_readID, which by-tile stats need."""
    dups = tmp_path / "dups.pairs"
    run_pairtools(
        "dedup", "--output", tmp_path / "nd.pairs", "--output-dups", dups,
        "--keep-parent-id", sorted_pairs
    )
    return dups


def read(path):
    with open(path) as f:
        return f.read()


STATS_OPTIONS = [
    [],
    ["--no-chromsizes"],
    ["--n-dist-bins-decade", "4"],
    ["--yaml"],
    ["--filter", 'unique:(pair_type=="UU")'],
    ["--yaml", "--filter", 'unique:(pair_type=="UU")', "--filter", "cis:(chrom1==chrom2)"],
]


@pytest.mark.parametrize(
    "options", STATS_OPTIONS, ids=lambda o: " ".join(o) or "defaults"
)
def test_stats_matches_pairtools(tmp_path, sorted_pairs, options):
    reference = tmp_path / "ref.stats"
    ours = tmp_path / "ours.stats"

    run_pairtools("stats", "-o", reference, *options, sorted_pairs)
    run_cli("stats", "-o", ours, *options, sorted_pairs)

    assert read(reference) == read(ours)


def test_stats_from_parquet(tmp_path, sorted_pairs):
    as_parquet = tmp_path / "in.parquet"
    run_cli("csv-to-parquet", "-o", as_parquet, sorted_pairs)

    reference = tmp_path / "ref.stats"
    ours = tmp_path / "ours.stats"
    run_pairtools("stats", "-o", reference, sorted_pairs)
    run_cli("stats", "-o", ours, as_parquet)

    assert read(reference) == read(ours)


def test_stats_are_independent_of_chunksize(tmp_path, sorted_pairs):
    """Counts are additive, so how the input is chunked must not matter.

    `pairtools stats` has no --chunksize (it hard-codes 100_000), so this is
    checked against ourselves rather than upstream.
    """
    small = tmp_path / "small.stats"
    large = tmp_path / "large.stats"
    run_cli("stats", "--chunksize", "137", "-o", small, sorted_pairs)
    run_cli("stats", "--chunksize", "1000000", "-o", large, sorted_pairs)

    assert read(small) == read(large)


def test_stats_merge_matches_pairtools(tmp_path, sorted_pairs):
    """Compared as sets of lines: upstream's line order is not reproducible.

    `pairtools stats --merge` emits its keys in an order that varies between
    processes, so two runs of pairtools itself disagree. We delegate to its
    `do_merge`, so we inherit that exactly -- which means the values can be
    compared but the ordering cannot.
    """
    part = tmp_path / "part.stats"
    run_pairtools("stats", "-o", part, sorted_pairs)

    reference = tmp_path / "ref.stats"
    ours = tmp_path / "ours.stats"
    run_pairtools("stats", "--merge", "-o", reference, part, part)
    run_cli("stats", "--merge", "-o", ours, part, part)

    assert sorted(read(reference).splitlines()) == sorted(read(ours).splitlines())


def test_bytile_dups_go_into_the_main_output(tmp_path, deduped_with_parents):
    """Without a separate path, the by-tile table is prepended to the stats."""
    reference = tmp_path / "ref.stats"
    ours = tmp_path / "ours.stats"

    run_pairtools("stats", "--bytile-dups", "-o", reference, deduped_with_parents)
    run_cli("stats", "--bytile-dups", "-o", ours, deduped_with_parents)

    assert read(reference) == read(ours)
    assert "dup_count" in read(ours)


def test_bytile_dups_to_separate_file(tmp_path, deduped_with_parents):
    ref_stats = tmp_path / "ref.stats"
    ref_tile = tmp_path / "ref.tile.tsv"
    our_stats = tmp_path / "ours.stats"
    our_tile = tmp_path / "ours.tile.tsv"

    run_pairtools(
        "stats", "--output-bytile-stats", ref_tile, "-o", ref_stats,
        deduped_with_parents
    )
    run_cli(
        "stats", "--output-bytile-stats", our_tile, "-o", our_stats,
        deduped_with_parents
    )

    assert read(ref_stats) == read(our_stats)
    assert read(ref_tile) == read(our_tile)


def test_dedup_stats_equal_standalone_stats(tmp_path, sorted_pairs):
    """`dedup --output-stats` and `stats` on the same pairs must agree.

    They accumulate through the same counter, so a discrepancy would mean one
    of the two pipelines is feeding it different rows.
    """
    from_dedup = tmp_path / "dedup.stats"
    run_cli(
        "dedup", "-o", tmp_path / "nd.pairs", "--output-stats", from_dedup,
        sorted_pairs
    )

    reference = tmp_path / "ref.stats"
    run_pairtools(
        "dedup", "--output", tmp_path / "r_nd.pairs", "--output-stats", reference,
        sorted_pairs
    )
    assert read(reference) == read(from_dedup)

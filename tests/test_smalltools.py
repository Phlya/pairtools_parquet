# -*- coding: utf-8 -*-
"""flip, markasdup and sample against pairtools."""

import random

import pytest

from conftest import (
    read_pairs_body,
    read_parquet_body,
    run_cli,
    run_pairtools,
    write_pairs,
)

HEADER = (
    ["## pairs format v1.0.0"]
    + ["#chromsize: chr{} 100000".format(i) for i in (1, 2, 3)]
    + ["#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type"]
)

PAIR_TYPES = ["UU", "UR", "RU", "NU", "UN", "MU", "UM", "WW", "DD"]


@pytest.fixture
def chromsizes(tmp_path):
    """Chromosome order that is deliberately not alphabetical.

    Flipping compares chromosomes by their rank in this file, not by name, so
    an alphabetical file would not tell the two apart.
    """
    path = tmp_path / "order.chromsizes"
    path.write_text("chr3\t100000\nchr1\t100000\nchr2\t100000\n")
    return path


@pytest.fixture
def mixed_pairs(tmp_path):
    """Pairs covering every branch of the flip decision.

    Includes chromosomes absent from the chromsizes file, the unmapped `!`
    sentinel, equal positions, and same-chromosome pairs -- each of which is a
    separate case in `pairtools flip`.
    """
    rng = random.Random(11)
    chroms = ["chr1", "chr2", "chr3", "chrUNKNOWN", "chrOTHER", "!"]
    rows = []
    for i in range(1500):
        c1, c2 = rng.choice(chroms), rng.choice(chroms)
        p1, p2 = rng.randrange(0, 100000), rng.randrange(0, 100000)
        if rng.random() < 0.15:
            p2 = p1
        if rng.random() < 0.15:
            c2 = c1
        rows.append(
            ("r{}".format(i), c1, p1, c2, p2, rng.choice("+-"), rng.choice("+-"),
             rng.choice(PAIR_TYPES))
        )
    return write_pairs(tmp_path / "mixed.pairs", HEADER, rows)


def test_flip_matches_pairtools(tmp_path, mixed_pairs, chromsizes):
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"

    run_pairtools("flip", "-c", chromsizes, "-o", reference, mixed_pairs)
    run_cli("flip", "-c", chromsizes, "-o", ours, mixed_pairs)

    assert read_pairs_body(reference) == read_pairs_body(ours)

    # the test is only meaningful if a good share of rows actually moved
    flipped = sum(
        1
        for before, after in zip(read_pairs_body(mixed_pairs), read_pairs_body(ours))
        if before != after
    )
    assert flipped > len(read_pairs_body(mixed_pairs)) // 4


def test_flip_through_parquet(tmp_path, mixed_pairs, chromsizes):
    as_parquet = tmp_path / "in.parquet"
    run_cli("csv-to-parquet", "-o", as_parquet, mixed_pairs)

    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.parquet"
    run_pairtools("flip", "-c", chromsizes, "-o", reference, mixed_pairs)
    run_cli("flip", "-c", chromsizes, "-o", ours, as_parquet)

    assert read_pairs_body(reference) == read_parquet_body(ours)


def test_flip_twice_matches_pairtools_twice(tmp_path, mixed_pairs, chromsizes):
    """Flipping is not idempotent upstream, and we match that too.

    For a pair whose two sides are on the same chromosome *absent from the
    chromsizes file*, `pairtools flip` decides by `chrom1 < chrom2`, which is
    false when the names are equal -- so it swaps the sides on every run,
    ignoring the positions, and oscillates. For annotated chromosomes it is
    idempotent, as intended. See UPSTREAM.md.
    """
    ours = [tmp_path / "ours1.pairs", tmp_path / "ours2.pairs"]
    theirs = [tmp_path / "ref1.pairs", tmp_path / "ref2.pairs"]

    run_cli("flip", "-c", chromsizes, "-o", ours[0], mixed_pairs)
    run_cli("flip", "-c", chromsizes, "-o", ours[1], ours[0])
    run_pairtools("flip", "-c", chromsizes, "-o", theirs[0], mixed_pairs)
    run_pairtools("flip", "-c", chromsizes, "-o", theirs[1], theirs[0])

    assert read_pairs_body(theirs[1]) == read_pairs_body(ours[1])


def test_flip_is_idempotent_for_annotated_chromosomes(tmp_path, chromsizes):
    """Where the chromosomes are known, flipping twice is a no-op."""
    rows = [
        ("r1", "chr1", 900, "chr1", 100, "-", "+", "RU"),
        ("r2", "chr2", 500, "chr1", 100, "+", "+", "UU"),
        ("r3", "chr3", 100, "chr1", 900, "+", "-", "UR"),
    ]
    source = write_pairs(tmp_path / "annotated.pairs", HEADER, rows)

    once = tmp_path / "once.pairs"
    twice = tmp_path / "twice.pairs"
    run_cli("flip", "-c", chromsizes, "-o", once, source)
    run_cli("flip", "-c", chromsizes, "-o", twice, once)

    assert read_pairs_body(once) == read_pairs_body(twice)


def test_markasdup_matches_pairtools(tmp_path, mock_pairs_path):
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"

    run_pairtools("markasdup", "-o", reference, mock_pairs_path)
    run_cli("markasdup", "-o", ours, mock_pairs_path)

    body = read_pairs_body(ours)
    assert read_pairs_body(reference) == body
    assert all(l.split("\t")[7] == "DD" for l in body)


def test_markasdup_through_parquet(tmp_path, mock_pairs_path):
    as_parquet = tmp_path / "in.parquet"
    run_cli("csv-to-parquet", "-o", as_parquet, mock_pairs_path)

    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.parquet"
    run_pairtools("markasdup", "-o", reference, mock_pairs_path)
    run_cli("markasdup", "-o", ours, as_parquet)

    assert read_pairs_body(reference) == read_parquet_body(ours)


@pytest.mark.parametrize("seed", [1, 42, 12345])
@pytest.mark.parametrize("fraction", [0.1, 0.5, 0.9])
def test_sample_matches_pairtools(tmp_path, mixed_pairs, seed, fraction):
    """The same seed must select the same pairs as pairtools.

    pairtools draws once per row from Python's `random` in file order, so
    matching it means reproducing that sequence of draws, not just the
    expected count.
    """
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"

    run_pairtools("sample", str(fraction), "-s", str(seed), "-o", reference, mixed_pairs)
    run_cli("sample", str(fraction), mixed_pairs, "-s", str(seed), "-o", ours)

    assert read_pairs_body(reference) == read_pairs_body(ours)
    assert len(read_pairs_body(ours)) > 0


def test_sample_is_independent_of_batching(tmp_path, mixed_pairs):
    """Reading Parquet in different batch sizes must not change the selection."""
    as_parquet = tmp_path / "in.parquet"
    run_cli("csv-to-parquet", "-o", as_parquet, mixed_pairs)

    from_text = tmp_path / "text.pairs"
    from_parquet = tmp_path / "parquet.pairs"
    run_cli("sample", "0.5", mixed_pairs, "-s", "7", "-o", from_text)
    run_cli("sample", "0.5", as_parquet, "-s", "7", "-o", from_parquet)

    assert read_pairs_body(from_text) == read_pairs_body(from_parquet)

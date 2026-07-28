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


ANNOTATED = {"chr1", "chr2", "chr3"}


def same_unannotated_chromosome(line):
    """Rows where we deliberately differ from pairtools (see UPSTREAM.md)."""
    fields = line.split("\t")
    return fields[1] == fields[3] and fields[1] not in ANNOTATED


def test_flip_matches_pairtools(tmp_path, mixed_pairs, chromsizes):
    """Identical to pairtools except for the case whose bug we fixed."""
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"

    run_pairtools("flip", "-c", chromsizes, "-o", reference, mixed_pairs)
    run_cli("flip", "-c", chromsizes, "-o", ours, mixed_pairs)

    theirs_body = read_pairs_body(reference)
    ours_body = read_pairs_body(ours)
    source_body = read_pairs_body(mixed_pairs)

    differing = [
        (src, a, b)
        for src, a, b in zip(source_body, theirs_body, ours_body)
        if a != b
    ]
    assert all(same_unannotated_chromosome(src) for src, _, _ in differing), (
        "diverged from pairtools outside the same-unannotated-chromosome case"
    )
    # and that case must actually be present, or this proves nothing
    assert any(same_unannotated_chromosome(l) for l in source_body)

    # the test is only meaningful if a good share of rows actually moved
    flipped = sum(1 for before, after in zip(source_body, ours_body) if before != after)
    assert flipped > len(source_body) // 4


def test_flip_orders_same_unannotated_chromosome_by_position(tmp_path, chromsizes):
    """Our fix: an unannotated chromosome is ordered by position, like any other.

    pairtools compares only the chromosome names here, which are equal, so it
    swaps the sides unconditionally. See UPSTREAM.md.
    """
    rows = [
        ("lower", "chrUNKNOWN", 900, "chrUNKNOWN", 100, "-", "+", "RU"),
        ("upper", "chrUNKNOWN", 100, "chrUNKNOWN", 900, "+", "-", "UR"),
    ]
    source = write_pairs(tmp_path / "unannotated.pairs", HEADER, rows)

    ours = tmp_path / "ours.pairs"
    run_cli("flip", "-c", chromsizes, "-o", ours, source)

    body = read_pairs_body(ours)
    for line in body:
        fields = line.split("\t")
        assert int(fields[2]) <= int(fields[4]), line
    # the already-ordered pair is left alone
    assert body[1].startswith("upper\tchrUNKNOWN\t100")


def test_flip_through_parquet(tmp_path, mixed_pairs, chromsizes):
    """Flipping must not depend on the format it reads and writes."""
    as_parquet = tmp_path / "in.parquet"
    run_cli("csv-to-parquet", "-o", as_parquet, mixed_pairs)

    from_text = tmp_path / "text.pairs"
    from_parquet = tmp_path / "parquet.parquet"
    run_cli("flip", "-c", chromsizes, "-o", from_text, mixed_pairs)
    run_cli("flip", "-c", chromsizes, "-o", from_parquet, as_parquet)

    assert read_pairs_body(from_text) == read_parquet_body(from_parquet)


def test_flip_is_idempotent(tmp_path, mixed_pairs, chromsizes):
    """Flipping is a projection onto the upper triangle, so it must settle.

    This is the property `pairtools flip` loses for pairs on the same
    unannotated chromosome, which oscillate there forever.
    """
    once = tmp_path / "once.pairs"
    twice = tmp_path / "twice.pairs"
    run_cli("flip", "-c", chromsizes, "-o", once, mixed_pairs)
    run_cli("flip", "-c", chromsizes, "-o", twice, once)

    assert read_pairs_body(once) == read_pairs_body(twice)


def test_pairtools_flip_oscillates(tmp_path, chromsizes):
    """Pins the upstream behaviour our fix departs from.

    If a future pairtools fixes this, this test fails and tells us the
    divergence -- and the note in UPSTREAM.md -- can be retired.
    """
    rows = [("r1", "chrUNKNOWN", 100, "chrUNKNOWN", 900, "+", "-", "UR")]
    source = write_pairs(tmp_path / "unannotated.pairs", HEADER, rows)

    once = tmp_path / "ref1.pairs"
    twice = tmp_path / "ref2.pairs"
    run_pairtools("flip", "-c", chromsizes, "-o", once, source)
    run_pairtools("flip", "-c", chromsizes, "-o", twice, once)

    assert read_pairs_body(once) != read_pairs_body(twice), (
        "pairtools flip no longer oscillates; our divergence can be removed"
    )


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

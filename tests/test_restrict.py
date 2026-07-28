# -*- coding: utf-8 -*-
"""restrict against pairtools restrict."""

import random

import numpy as np
import pytest

from conftest import (
    read_pairs_body,
    read_pairs_header,
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


@pytest.fixture
def frags(tmp_path):
    rng = random.Random(5)
    path = tmp_path / "frags.bed"
    with open(path, "w") as f:
        for chrom in ("chr1", "chr2", "chr3"):
            pos = 0
            while pos < 100000:
                nxt = pos + rng.randint(500, 4000)
                f.write("{}\t{}\t{}\n".format(chrom, pos, min(nxt, 100000)))
                pos = nxt
    return path


@pytest.fixture
def annotated_pairs(tmp_path):
    """Pairs only on chromosomes the fragment file covers.

    Includes fragment boundaries, position 0, and positions past the last
    fragment, which are the edges of the searchsorted lookup.
    """
    rng = random.Random(21)
    rows = []
    for i in range(1500):
        c1, c2 = rng.choice(["chr1", "chr2", "chr3"]), rng.choice(["chr1", "chr2", "chr3"])
        p1 = rng.choice([0, rng.randrange(0, 100000), 99999, 100000, 150000])
        p2 = rng.choice([0, rng.randrange(0, 100000), 99999, 100000, 150000])
        rows.append(("r{}".format(i), c1, p1, c2, p2, rng.choice("+-"), rng.choice("+-"), "UU"))
    for i in range(30):
        rows.append(("u{}".format(i), "!", 0, "chr1", rng.randrange(0, 100000), "-", "+", "NU"))
    return write_pairs(tmp_path / "in.pairs", HEADER, rows)


def test_restrict_matches_pairtools(tmp_path, annotated_pairs, frags):
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"

    run_pairtools("restrict", "-f", frags, "-o", reference, annotated_pairs)
    run_cli("restrict", "-f", frags, "-o", ours, annotated_pairs)

    assert read_pairs_body(reference) == read_pairs_body(ours)
    assert [l for l in read_pairs_header(reference) if l.startswith("#columns:")] == [
        l for l in read_pairs_header(ours) if l.startswith("#columns:")
    ]


def test_restrict_through_parquet(tmp_path, annotated_pairs, frags):
    as_parquet = tmp_path / "in.parquet"
    run_cli("csv-to-parquet", "-o", as_parquet, annotated_pairs)

    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.parquet"
    run_pairtools("restrict", "-f", frags, "-o", reference, annotated_pairs)
    run_cli("restrict", "-f", frags, "-o", ours, as_parquet)

    assert read_pairs_body(reference) == read_parquet_body(ours)


def test_unmapped_sides_are_unannotated(tmp_path, annotated_pairs, frags):
    ours = tmp_path / "ours.pairs"
    run_cli("restrict", "-f", frags, "-o", ours, annotated_pairs)

    unmapped = [l for l in read_pairs_body(ours) if l.split("\t")[1] == "!"]
    assert unmapped
    for line in unmapped:
        rfrag1, start1, end1 = line.split("\t")[8:11]
        assert (rfrag1, start1, end1) == ("-1", "0", "0")


def test_chromosome_without_fragments_is_not_fatal(tmp_path, frags):
    """`pairtools restrict` crashes here; we do what its code intends.

    `find_rfrag` catches ValueError around a dict lookup, which raises
    KeyError, so the intended "warn and return empty" path is unreachable and
    the command dies. See UPSTREAM.md.
    """
    rows = [
        ("r1", "chr1", 100, "chr1", 900, "+", "-", "UU"),
        ("r2", "chrABSENT", 100, "chr1", 900, "+", "-", "UU"),
    ]
    source = write_pairs(tmp_path / "absent.pairs", HEADER, rows)

    ours = tmp_path / "ours.pairs"
    run_cli("restrict", "-f", frags, "-o", ours, source)

    body = read_pairs_body(ours)
    assert len(body) == 2
    # the row on the missing chromosome gets the unannotated sentinels
    assert body[1].split("\t")[8:11] == ["-1", "0", "0"]
    # the row on a known chromosome is still annotated
    assert body[0].split("\t")[8] != "-1"


def test_pairtools_restrict_crashes_on_missing_chromosome(tmp_path, frags):
    """Pins the upstream behaviour our fix departs from."""
    rows = [("r1", "chrABSENT", 100, "chr1", 900, "+", "-", "UU")]
    source = write_pairs(tmp_path / "absent.pairs", HEADER, rows)

    with pytest.raises(AssertionError, match="KeyError"):
        run_pairtools("restrict", "-f", frags, "-o", tmp_path / "r.pairs", source)


def test_fragment_loader_matches_genfromtxt(frags):
    """Our pyarrow loader must produce exactly what upstream's genfromtxt does."""
    from pairtools_parquet.lib.restrict import load_rfrags

    ours = load_rfrags(str(frags))

    upstream = np.genfromtxt(
        str(frags), delimiter="\t", comments="#", dtype=None, encoding="ascii",
        names=["chrom", "start", "end"],
    )
    upstream.sort(order=["chrom", "start", "end"])
    borders = np.r_[
        0, 1 + np.where(upstream["chrom"][:-1] != upstream["chrom"][1:])[0],
        upstream.shape[0],
    ]
    expected = {
        upstream["chrom"][i]: np.concatenate([[0], upstream["end"][i:j] + 1])
        for i, j in zip(borders[:-1], borders[1:])
    }

    assert set(ours) == set(expected)
    for chrom in expected:
        assert np.array_equal(ours[chrom], expected[chrom]), chrom

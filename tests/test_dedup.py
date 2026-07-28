# -*- coding: utf-8 -*-
"""dedup against pairtools dedup, across its options.

Duplicate detection depends on which pairs share a chunk, so these run over a
few thousand pairs with planted near-duplicates rather than the 9-row mock.
"""

import random

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
    ["## pairs format v1.0.0", "#shape: upper triangle"]
    + ["#chromsize: chr{} 1000000".format(i) for i in range(1, 4)]
    + ["#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type"]
)


def make_duplicated_pairs(path, n_unique=800, seed=7):
    """Write pairs with planted near-duplicates and some unmapped reads."""
    rng = random.Random(seed)
    rows = []
    n = 0
    for _ in range(n_unique):
        c1 = "chr{}".format(rng.randint(1, 3))
        c2 = "chr{}".format(rng.randint(1, 3))
        if c1 > c2:
            c1, c2 = c2, c1
        p1 = rng.randrange(1, 1000000)
        p2 = rng.randrange(1, 1000000)
        s1, s2 = rng.choice("+-"), rng.choice("+-")
        rows.append(("r{}".format(n), c1, p1, c2, p2, s1, s2, "UU"))
        n += 1
        # near-duplicates, within the default --max-mismatch of 3bp
        for _ in range(rng.choice([0, 0, 1, 2])):
            rows.append(
                (
                    "r{}".format(n),
                    c1,
                    p1 + rng.randint(0, 3),
                    c2,
                    p2 + rng.randint(0, 3),
                    s1,
                    s2,
                    "UU",
                )
            )
            n += 1
    for i in range(30):
        rows.append(("u{}".format(i), "!", 0, "chr1", rng.randrange(1, 1000), "-", "+", "NU"))
    return write_pairs(path, HEADER, rows)


@pytest.fixture
def sorted_duplicated_pairs(tmp_path):
    raw = make_duplicated_pairs(tmp_path / "raw.pairs")
    path = tmp_path / "sorted.pairs"
    run_pairtools("sort", "-o", path, raw)
    return path


def has_sklearn():
    try:
        import sklearn  # noqa: F401
    except ImportError:
        return False
    return True


DEDUP_OPTIONS = [
    [],
    ["--max-mismatch", "0"],
    ["--max-mismatch", "10", "--method", "sum"],
    ["--chunksize", "137"],
    ["--chunksize", "100000"],
    ["--carryover", "5"],
    ["--no-mark-dups"],
    pytest.param(
        ["--backend", "sklearn"],
        marks=pytest.mark.skipif(
            not has_sklearn(), reason="the sklearn backend needs scikit-learn"
        ),
    ),
]


@pytest.mark.parametrize("options", DEDUP_OPTIONS, ids=lambda o: " ".join(o) or "defaults")
def test_dedup_matches_pairtools(tmp_path, sorted_duplicated_pairs, options):
    ref = {p: tmp_path / "ref_{}.pairs".format(p) for p in ("nodups", "dups", "unmapped")}
    ours = {p: tmp_path / "our_{}.pairs".format(p) for p in ref}

    run_pairtools(
        "dedup", "--output", ref["nodups"], "--output-dups", ref["dups"],
        "--output-unmapped", ref["unmapped"], *options, sorted_duplicated_pairs
    )
    run_cli(
        "dedup", "-o", ours["nodups"], "--output-dups", ours["dups"],
        "--output-unmapped", ours["unmapped"], *options, sorted_duplicated_pairs
    )

    for part in ref:
        assert read_pairs_body(ref[part]) == read_pairs_body(ours[part]), part
    # the planted duplicates must actually be found, or this proves nothing
    assert len(read_pairs_body(ours["dups"])) > 0


def test_dedup_stats_match_pairtools(tmp_path, sorted_duplicated_pairs):
    ref_stats = tmp_path / "ref.stats"
    our_stats = tmp_path / "our.stats"

    run_pairtools(
        "dedup", "--output", tmp_path / "r.pairs", "--output-stats", ref_stats,
        sorted_duplicated_pairs
    )
    run_cli(
        "dedup", "-o", tmp_path / "o.pairs", "--output-stats", our_stats,
        sorted_duplicated_pairs
    )

    with open(ref_stats) as a, open(our_stats) as b:
        assert a.read() == b.read()


def test_dedup_through_parquet(tmp_path, sorted_duplicated_pairs):
    """Deduplicating a Parquet file must give the same pairs as the text file."""
    as_parquet = tmp_path / "in.parquet"
    run_cli("csv-to-parquet", "-o", as_parquet, sorted_duplicated_pairs)

    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.parquet"
    run_pairtools("dedup", "--output", reference, sorted_duplicated_pairs)
    run_cli("dedup", "-o", ours, as_parquet)

    assert read_pairs_body(reference) == read_parquet_body(ours)


def test_dedup_dups_into_main_output(tmp_path, sorted_duplicated_pairs):
    """--output-dups equal to --output keeps duplicates in the main file."""
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"

    run_pairtools(
        "dedup", "--output", reference, "--output-dups", reference,
        sorted_duplicated_pairs
    )
    run_cli("dedup", "-o", ours, "--output-dups", ours, sorted_duplicated_pairs)

    assert read_pairs_body(reference) == read_pairs_body(ours)


def test_dedup_keep_parent_id(tmp_path, sorted_duplicated_pairs):
    """parent_readID is added to the duplicates, but not to the non-duplicates."""
    ref = {p: tmp_path / "ref_{}.pairs".format(p) for p in ("nodups", "dups", "unmapped")}
    ours = {p: tmp_path / "our_{}.pairs".format(p) for p in ref}

    run_pairtools(
        "dedup", "--output", ref["nodups"], "--output-dups", ref["dups"],
        "--output-unmapped", ref["unmapped"], "--keep-parent-id",
        sorted_duplicated_pairs
    )
    run_cli(
        "dedup", "-o", ours["nodups"], "--output-dups", ours["dups"],
        "--output-unmapped", ours["unmapped"], "--keep-parent-id",
        sorted_duplicated_pairs
    )

    for part in ref:
        assert read_pairs_body(ref[part]) == read_pairs_body(ours[part]), part

    def columns_of(path):
        line = [l for l in read_pairs_header(path) if l.startswith("#columns:")][0]
        return line[len("#columns: ") :].split()

    assert "parent_readID" in columns_of(ours["dups"])
    assert "parent_readID" not in columns_of(ours["nodups"])


def test_unmapped_header_matches_its_own_rows(tmp_path, sorted_duplicated_pairs):
    """A deliberate divergence from pairtools, asserted so it stays deliberate.

    Upstream writes the unmapped stream before dropping parent_readID but with
    the header from before it was added, declaring 8 columns for 9-field rows.
    Parquet cannot represent that at all -- its schema is real -- so our
    unmapped header counts the columns its rows actually have.
    """
    ours_unmapped = tmp_path / "our_unmapped.pairs"
    run_cli(
        "dedup", "-o", tmp_path / "o.pairs", "--output-unmapped", ours_unmapped,
        "--keep-parent-id", sorted_duplicated_pairs
    )

    line = [l for l in read_pairs_header(ours_unmapped) if l.startswith("#columns:")][0]
    n_columns = len(line[len("#columns: ") :].split())

    # not read_pairs_body: parent_readID is empty for unmapped pairs, so the
    # rows end in a tab that .strip() would eat along with the field it delimits
    with open(ours_unmapped) as f:
        body = [l.rstrip("\n") for l in f if not l.startswith("#") and l.strip()]

    assert body, "expected some unmapped pairs"
    assert all(len(l.split("\t")) == n_columns for l in body)


def test_cython_backend_reports_clearly(tmp_path, sorted_duplicated_pairs):
    """pairtools' cython backend is line-based, so it is not offered here."""
    with pytest.raises(AssertionError, match="backend"):
        run_cli(
            "dedup", "-o", tmp_path / "o.pairs", "--backend", "cython",
            sorted_duplicated_pairs
        )

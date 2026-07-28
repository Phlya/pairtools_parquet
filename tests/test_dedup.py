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


@pytest.fixture(scope="session")
def sorted_duplicated_pairs(tmp_path_factory):
    """Built once for the session: every test only reads it.

    Sorting it per test meant a `pairtools sort` subprocess for each of the
    thirty-odd tests here, which dominated the suite's runtime.
    """
    directory = tmp_path_factory.mktemp("dedup_input")
    raw = make_duplicated_pairs(directory / "raw.pairs")
    path = directory / "sorted.pairs"
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


# --------------------------------------------------------------------------
# The duckdb backend, against the pandas one it replaces
# --------------------------------------------------------------------------

BACKEND_OPTIONS = [
    [],
    ["--max-mismatch", "0"],
    ["--max-mismatch", "10"],
    ["--method", "sum"],
    ["--method", "sum", "--max-mismatch", "8"],
    ["--no-mark-dups"],
    ["--keep-parent-id"],
]


@pytest.mark.parametrize(
    "options", BACKEND_OPTIONS, ids=lambda o: " ".join(o) or "defaults"
)
def test_duckdb_backend_matches_scipy(tmp_path, sorted_duplicated_pairs, options):
    """The fast path has to agree with the implementation it replaces."""
    out = {}
    for backend in ("scipy", "duckdb"):
        out[backend] = {
            part: tmp_path / "{}_{}.pairs".format(backend, part)
            for part in ("nodups", "dups", "unmapped")
        }
        run_cli(
            "dedup", "--backend", backend,
            "-o", out[backend]["nodups"],
            "--output-dups", out[backend]["dups"],
            "--output-unmapped", out[backend]["unmapped"],
            *options, sorted_duplicated_pairs
        )

    for part in out["scipy"]:
        assert read_pairs_body(out["scipy"][part]) == read_pairs_body(
            out["duckdb"][part]
        ), part
    assert len(read_pairs_body(out["duckdb"]["dups"])) > 0


# Small windows mean one query per window, so these stay modest; the
# single-row windows are exercised on the tiny fixtures below instead.
@pytest.mark.parametrize("chunksize", [11, 101, 100000])
def test_duckdb_windows_do_not_change_the_answer(
    tmp_path, sorted_duplicated_pairs, chunksize
):
    """--chunksize is a memory knob for this backend, not a semantic one.

    Each window is deduplicated with a lookback of --carryover rows, and a
    cluster whose smallest row falls in the lookback inherits that row's own
    keeper -- so a family split across a window boundary still resolves to one
    keeper however small the window is.
    """
    whole = tmp_path / "whole.pairs"
    windowed = tmp_path / "windowed.pairs"

    run_cli("dedup", "-o", whole, sorted_duplicated_pairs)
    run_cli(
        "dedup", "-o", windowed, "--chunksize", str(chunksize),
        sorted_duplicated_pairs
    )

    assert read_pairs_body(whole) == read_pairs_body(windowed)


def test_duckdb_backend_rejects_mismatched_extra_col_pair(
    tmp_path, sorted_duplicated_pairs
):
    """An extra column pair naming two different columns is not a join key."""
    with pytest.raises(AssertionError, match="backend scipy"):
        run_cli(
            "dedup", "-o", tmp_path / "o.pairs",
            "--extra-col-pair", "strand1", "strand2",
            sorted_duplicated_pairs
        )


def test_duckdb_backend_keeps_chains_across_a_window_boundary(tmp_path):
    """A chain of near-duplicates must resolve to one keeper, boundary or not.

    `pairtools dedup` loses this case: its carryover holds only the previous
    chunk's *non*-duplicates, so a chain that reaches the boundary through a
    duplicate is cut and the next read is kept as unique. See UPSTREAM.md.
    """
    header = (
        ["## pairs format v1.0.0", "#shape: upper triangle", "#chromsize: chr1 10000"]
        + ["#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type"]
    )
    # each row is 3bp from the previous one and 6bp from the one before that,
    # so they are one cluster only by chaining
    rows = [
        ("r{}".format(i), "chr1", 1000 + 3 * i, "chr1", 5000, "+", "+", "UU")
        for i in range(4)
    ]
    path = write_pairs(tmp_path / "chain.pairs", header, rows)

    for chunksize in (1, 2, 3, 4, 100):
        out = tmp_path / "out_{}.pairs".format(chunksize)
        run_cli("dedup", "-o", out, "--chunksize", str(chunksize), path)
        body = read_pairs_body(out)
        assert len(body) == 1, "chunksize {} split the chain".format(chunksize)
        assert body[0].split("\t")[0] == "r0"


def test_duckdb_window_joins_two_families_through_a_boundary(tmp_path):
    """A window can be the first to see that two earlier families are one.

    Rows 0 and 1 start separate; row 2 (in the next window) is within reach of
    both, so all three are one family and only row 0 survives. Resolving the
    component to the keeper of its own smallest member -- row 1 here -- would
    keep row 0 *and* row 1, because row 1's family reaches back to row 0 only
    through the row that joined them.
    """
    header = (
        ["## pairs format v1.0.0", "#shape: upper triangle", "#chromsize: chr1 10000"]
        + ["#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type"]
    )
    # pos2 apart by 4 (> --max-mismatch 3) for r0/r1, but r2 sits between them
    rows = [
        ("r0", "chr1", 1000, "chr1", 5000, "+", "+", "UU"),
        ("r1", "chr1", 1000, "chr1", 5004, "+", "+", "UU"),
        ("r2", "chr1", 1000, "chr1", 5002, "+", "+", "UU"),
    ]
    path = write_pairs(tmp_path / "join.pairs", header, rows)

    for chunksize in (1, 2, 3, 100):
        out = tmp_path / "out_{}.pairs".format(chunksize)
        run_cli("dedup", "-o", out, "--chunksize", str(chunksize), path)
        body = read_pairs_body(out)
        assert [l.split("\t")[0] for l in body] == ["r0"], (
            "chunksize {} kept {}".format(chunksize, [l.split("\t")[0] for l in body])
        )


@pytest.mark.parametrize("chunksize", [11, 100000])
def test_duckdb_parent_ids_do_not_depend_on_the_window(
    tmp_path, sorted_duplicated_pairs, chunksize
):
    """A keeper can move earlier as a window revises its lookback.

    The recorded parent has to end up the same as a single-window run, or
    --keep-parent-id would name a different member of the family depending on
    a memory setting.
    """
    whole = tmp_path / "whole.dups"
    windowed = tmp_path / "windowed.dups"

    run_cli(
        "dedup", "-o", tmp_path / "w.pairs", "--output-dups", whole,
        "--keep-parent-id", sorted_duplicated_pairs
    )
    run_cli(
        "dedup", "-o", tmp_path / "n.pairs", "--output-dups", windowed,
        "--keep-parent-id", "--chunksize", str(chunksize), sorted_duplicated_pairs
    )

    assert read_pairs_body(whole) == read_pairs_body(windowed)

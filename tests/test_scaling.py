# -*- coding: utf-8 -*-
"""scaling against pairtools scaling.

The chromosome sizes come from the file's own `#chromsize:` header lines, as
they do in pairtools master. pairtools 1.1.x parses them and then drops them,
leaving every region's end at -1 and making `n_bp2` -- the area P(s) is
normalised by -- meaningless, so a bare run means different things on the two
versions. The parity checks here therefore compare against `pairtools scaling
--view <the header's chromsizes>`, which both versions answer identically and
which is what our output equals exactly.
"""

import pandas as pd
import pytest

from conftest import run_cli, run_pairtools
from test_dedup import make_duplicated_pairs


@pytest.fixture
def pairs(tmp_path):
    return make_duplicated_pairs(tmp_path / "pairs.pairs")


@pytest.fixture
def view(tmp_path):
    """A viewframe restricting the calculation to part of each chromosome.

    `pairtools scaling --view` reads the file with `pd.read_table`, so it wants
    a named header row rather than a bare BED.
    """
    path = tmp_path / "view.bed"
    with open(path, "w") as f:
        f.write("chrom\tstart\tend\n")
        for i in range(1, 4):
            f.write("chr{}\t0\t500000\n".format(i))
    return path


def read(path):
    with open(path) as f:
        return f.read()


def header_view(pairs_path, destination):
    """A viewframe spelling out what the file's header already declares."""
    from pairtools.lib import headerops

    from pairtools_parquet.lib import arrowio

    chromsizes = headerops.extract_chromsizes(arrowio.read_header(str(pairs_path)))
    with open(destination, "w") as f:
        f.write("chrom\tstart\tend\tname\n")
        for chrom, size in chromsizes.items():
            f.write("{}\t0\t{}\t{}\n".format(chrom, size, chrom))
    return destination


def counts_only(text):
    """The table without the columns the chromosome sizes decide.

    `n_pairs` is what the binning counts; the region bounds and `n_bp2` are
    what it is normalised by. The counts must match pairtools whatever the
    normalisation, so they are compared separately.
    """
    rows = [line.split("\t") for line in text.splitlines()]
    keep = [i for i, name in enumerate(rows[0])
            if name not in ("end1", "end2", "n_bp2")]
    return ["\t".join(row[i] for i in keep) for row in rows]


OPTIONS = [
    [],
    ["--dist-range", "10", "100000"],
    ["--n-dist-bins-decade", "3"],
    ["--chunksize", "100"],
]


@pytest.mark.parametrize("options", OPTIONS, ids=lambda o: " ".join(o) or "defaults")
def test_scaling_matches_pairtools(tmp_path, pairs, options):
    reference = tmp_path / "ref.tsv"
    ours = tmp_path / "ours.tsv"

    view = header_view(pairs, tmp_path / "header_view.tsv")
    run_pairtools("scaling", "--view", view, "-o", reference, *options, pairs)
    run_cli("scaling", "-o", ours, *options, pairs)

    assert read(reference) == read(ours)


def test_scaling_with_view_matches_pairtools(tmp_path, pairs, view):
    reference = tmp_path / "ref.tsv"
    ours = tmp_path / "ours.tsv"

    run_pairtools("scaling", "--view", view, "-o", reference, pairs)
    run_cli("scaling", "--view", view, "-o", ours, pairs)

    assert read(reference) == read(ours)


def test_scaling_on_mock_matches_pairtools(tmp_path, mock_pairs_path):
    reference = tmp_path / "ref.tsv"
    ours = tmp_path / "ours.tsv"

    view = header_view(mock_pairs_path, tmp_path / "header_view.tsv")
    run_pairtools("scaling", "--view", view, "-o", reference, mock_pairs_path)
    run_cli("scaling", "-o", ours, mock_pairs_path)

    assert read(reference) == read(ours)


def test_scaling_through_parquet(tmp_path, pairs):
    as_parquet = tmp_path / "in.parquet"
    run_cli("csv-to-parquet", "-o", as_parquet, pairs)

    reference = tmp_path / "ref.tsv"
    ours = tmp_path / "ours.tsv"
    view = header_view(pairs, tmp_path / "header_view.tsv")
    run_pairtools("scaling", "--view", view, "-o", reference, pairs)
    run_cli("scaling", "-o", ours, as_parquet)

    assert read(reference) == read(ours)


def test_chunked_loop_matches_compute_scaling(pairs):
    """The restated chunk loop must agree with pairtools' own.

    `compute_scaling` cannot be handed an iterator of DataFrames, so the loop
    around its per-chunk kernel lives here instead. This pins that loop against
    the upstream one, run on the whole file at once.
    """
    from pairtools.lib.pairsio import read_pairs
    from pairtools.lib.scaling import compute_scaling

    from pairtools_parquet.lib.scaling import (
        accumulate_scaling,
        dist_bin_edges,
        scaling_chunks,
    )

    whole_df, _, _ = read_pairs(str(pairs))
    ref_cis, ref_trans = compute_scaling(whole_df)

    cis, trans = accumulate_scaling(
        scaling_chunks(str(pairs), chunksize=100), dist_bin_edges((1, int(1e9)), 8)
    )

    for reference, ours in ((ref_cis, cis), (ref_trans, trans)):
        assert list(reference.columns) == list(ours.columns)
        pd.testing.assert_frame_equal(
            reference.reset_index(drop=True),
            ours.reset_index(drop=True),
            check_dtype=False,
        )


def test_scaling_of_a_file_with_no_pairs(tmp_path, pairs):
    """A header with no data rows still has a scaling: a table of zeros."""
    from conftest import read_pairs_header, write_pairs

    empty = write_pairs(tmp_path / "empty.pairs", read_pairs_header(pairs), [])
    reference = tmp_path / "ref.tsv"
    ours = tmp_path / "ours.tsv"

    view = header_view(empty, tmp_path / "header_view.tsv")
    run_pairtools("scaling", "--view", view, "-o", reference, empty)
    run_cli("scaling", "-o", ours, empty)

    assert read(reference) == read(ours)


def test_scaling_writes_to_stdout(pairs):
    """An empty --output means stdout, as it does upstream."""
    assert run_cli("scaling", pairs).startswith("chrom1\t")


KEYS = ["chrom1", "chrom2", "strand1", "strand2", "min_dist", "max_dist"]


def test_chromsizes_come_from_the_header(tmp_path, pairs):
    """The header's chromosome sizes reach the binning.

    pairtools master takes them from `read_pairs`' third return value; 1.1.x
    takes them and drops them, leaving every region end at the -1 sentinel, so
    the area each distance bin covers is computed from a region one base long
    and `n_bp2` -- what P(s) is divided by -- is not the genome's.

    So this asserts our answer outright, and then, against whichever pairtools
    is installed, either full equality or exactly that difference.
    """
    reference = tmp_path / "ref.tsv"
    ours = tmp_path / "ours.tsv"
    run_pairtools("scaling", "-o", reference, pairs)
    run_cli("scaling", "-o", ours, pairs)

    theirs = pd.read_csv(reference, sep="\t")
    mine = pd.read_csv(ours, sep="\t")

    assert (mine["end1"] > 0).all()
    assert mine["n_bp2"].sum() > 0
    assert "!" not in set(mine["chrom1"]) | set(mine["chrom2"])

    if not (set(theirs["end1"]) == {-1}):  # pairtools has the fix
        assert read(reference) == read(ours)
        return

    # A -1 end leaves a region one base long, so 1.1.x's areas are a rounding
    # error next to the genome's.
    assert mine["n_bp2"].sum() > 1000 * theirs["n_bp2"].sum()

    # The counting itself is untouched: every bin the two tables share holds
    # the same number of pairs.
    both = mine.merge(theirs, on=KEYS, suffixes=("_mine", "_theirs"))
    assert len(both) == len(mine)
    assert list(both["n_pairs_mine"]) == list(both["n_pairs_theirs"])

    # Knowing the chromosome sizes drops rows, and none that a P(s) curve
    # wants: bins reaching past the end of a chromosome, which are empty, and
    # the unmapped `!` "region" 1.1.x invents because with no sizes to go on it
    # takes its regions from the data.
    dropped = theirs.merge(mine[KEYS], on=KEYS, how="left", indicator=True)
    dropped = dropped[dropped["_merge"] == "left_only"]
    assert len(dropped) > 0
    counted = dropped[dropped["n_pairs"] > 0]
    assert len(counted) > 0
    assert ((counted["chrom1"] == "!") | (counted["chrom2"] == "!")).all()


def test_a_header_without_chromsizes_falls_back_to_sentinel_regions(tmp_path, pairs):
    """There is nothing better to give it, so the -1 sentinel stands.

    `pairtools scaling` cannot read such a file at all, on any version --
    `read_pairs` calls `extract_chromsizes` unconditionally, which zips the
    parsed lines and indexes the result without checking there were any -- so
    there is no reference run to compare against.
    """
    from conftest import read_pairs_body, read_pairs_header, write_pairs

    header = [l for l in read_pairs_header(pairs) if not l.startswith("#chromsize")]
    rows = [line.split("\t") for line in read_pairs_body(pairs)]
    stripped = write_pairs(tmp_path / "no_chromsizes.pairs", header, rows)

    with pytest.raises(Exception):
        run_pairtools("scaling", "-o", tmp_path / "ref.tsv", stripped)

    ours = tmp_path / "ours.tsv"
    run_cli("scaling", "-o", ours, stripped)
    mine = pd.read_csv(ours, sep="\t")

    assert set(mine["end1"]) == {-1}
    assert mine["n_pairs"].sum() > 0

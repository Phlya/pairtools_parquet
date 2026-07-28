# -*- coding: utf-8 -*-
"""phase against pairtools phase.

Phasing needs pairs carrying the aligner's tag columns, which `pairtools parse
--add-columns` produces. Since parse is not ported yet, the fixtures here
synthesize those columns directly, covering the branches phase_side_XB and
phase_side_XA take: no alternative alignment, one, two, an alternative on the
homologous chromosome, and an unmapped side.
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

CHROMS = ["chr1_pat", "chr1_mat", "chr2_pat", "chr2_mat", "!"]
CHROMSIZES = [
    "#chromsize: chr{}_{} 100000".format(i, s) for i in (1, 2) for s in ("pat", "mat")
]


def header_for(columns):
    return ["## pairs format v1.0.0"] + CHROMSIZES + ["#columns: " + " ".join(columns)]


@pytest.fixture
def xb_pairs(tmp_path):
    """Pairs with bwa-mem -M style XB tags."""
    rng = random.Random(31)
    columns = [
        "readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2",
        "pair_type", "XB1", "XB2", "AS1", "AS2", "XS1", "XS2",
    ]

    def alignment(chrom, score):
        other = (
            chrom.replace("_pat", "_mat")
            if chrom.endswith("_pat")
            else chrom.replace("_mat", "_pat")
        )
        # chrom,pos,CIGAR,NM,AS,mapq
        return "{},+{},100M,{},{},{}".format(
            rng.choice([other, chrom, "chr2_pat"]),
            rng.randrange(1, 100000), rng.randint(0, 3),
            rng.randint(score - 10, score), rng.choice([0, 30, 60]),
        )

    def xb(chrom, score):
        if chrom == "!":
            return ""
        n = rng.choice([0, 1, 1, 2])
        return ";".join(alignment(chrom, score) for _ in range(n)) + (";" if n else "")

    rows = []
    for i in range(500):
        c1, c2 = rng.choice(CHROMS), rng.choice(CHROMS)
        p1 = 0 if c1 == "!" else rng.randrange(1, 100000)
        p2 = 0 if c2 == "!" else rng.randrange(1, 100000)
        s1 = "-" if c1 == "!" else rng.choice("+-")
        s2 = "-" if c2 == "!" else rng.choice("+-")
        as1, as2 = rng.randint(80, 100), rng.randint(80, 100)
        # XS >= AS is what sends the decision into the XB-inspecting branch
        xs1 = rng.choice([0, as1, as1 - rng.randint(1, 20)])
        xs2 = rng.choice([0, as2, as2 - rng.randint(1, 20)])
        rows.append(
            ("r{}".format(i), c1, p1, c2, p2, s1, s2, "UU",
             xb(c1, as1), xb(c2, as2), as1, as2, xs1, xs2)
        )
    return write_pairs(tmp_path / "xb.pairs", header_for(columns), rows)


@pytest.fixture
def xa_pairs(tmp_path):
    """Pairs with default bwa-mem XA tags."""
    rng = random.Random(41)
    columns = [
        "readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2",
        "pair_type", "XA1", "XA2", "NM1", "NM2", "AS1", "AS2", "XS1", "XS2",
    ]

    def xa(chrom):
        if chrom == "!":
            return ""
        other = (
            chrom.replace("_pat", "_mat")
            if chrom.endswith("_pat")
            else chrom.replace("_mat", "_pat")
        )
        n = rng.choice([0, 1, 1, 2])
        return ";".join(
            "{},+{},100M,{}".format(
                rng.choice([other, chrom]), rng.randrange(1, 100000), rng.randint(0, 3)
            )
            for _ in range(n)
        ) + (";" if n else "")

    rows = []
    for i in range(400):
        c1, c2 = rng.choice(CHROMS), rng.choice(CHROMS)
        p1 = 0 if c1 == "!" else rng.randrange(1, 100000)
        p2 = 0 if c2 == "!" else rng.randrange(1, 100000)
        s1 = "-" if c1 == "!" else rng.choice("+-")
        s2 = "-" if c2 == "!" else rng.choice("+-")
        as1, as2 = rng.randint(80, 100), rng.randint(80, 100)
        rows.append(
            ("r{}".format(i), c1, p1, c2, p2, s1, s2, "UU", xa(c1), xa(c2),
             rng.randint(0, 3), rng.randint(0, 3), as1, as2,
             rng.choice([0, as1, as1 - 5]), rng.choice([0, as2, as2 - 5]))
        )
    return write_pairs(tmp_path / "xa.pairs", header_for(columns), rows)


def columns_line(path):
    return [l for l in read_pairs_header(path) if l.startswith("#columns:")][0]


@pytest.mark.parametrize(
    "options",
    [[], ["--clean-output"], ["--report-scores"], ["--clean-output", "--report-scores"]],
    ids=lambda o: " ".join(o) or "defaults",
)
def test_phase_xb_matches_pairtools(tmp_path, xb_pairs, options):
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"

    run_pairtools(
        "phase", "--phase-suffixes", "_pat", "_mat", "-o", reference, *options, xb_pairs
    )
    run_cli(
        "phase", xb_pairs, "--phase-suffixes", "_pat", "_mat", "-o", ours, *options
    )

    assert read_pairs_body(reference) == read_pairs_body(ours)
    assert columns_line(reference) == columns_line(ours)


def test_phase_xa_matches_pairtools(tmp_path, xa_pairs):
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"

    run_pairtools(
        "phase", "--phase-suffixes", "_pat", "_mat", "--tag-mode", "XA",
        "-o", reference, xa_pairs
    )
    run_cli(
        "phase", xa_pairs, "--phase-suffixes", "_pat", "_mat", "--tag-mode", "XA",
        "-o", ours
    )

    assert read_pairs_body(reference) == read_pairs_body(ours)
    assert columns_line(reference) == columns_line(ours)


def test_all_phase_outcomes_are_exercised(tmp_path, xb_pairs):
    """The fixture must reach every branch, or parity proves little."""
    ours = tmp_path / "ours.pairs"
    run_cli("phase", xb_pairs, "--phase-suffixes", "_pat", "_mat", "-o", ours)

    phases = {l.split("\t")[14] for l in read_pairs_body(ours)}
    assert phases == {"!", ".", "0", "1"}


def test_phase_through_parquet(tmp_path, xb_pairs):
    as_parquet = tmp_path / "in.parquet"
    run_cli("csv-to-parquet", "-o", as_parquet, xb_pairs)

    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.parquet"
    run_pairtools("phase", "--phase-suffixes", "_pat", "_mat", "-o", reference, xb_pairs)
    run_cli("phase", as_parquet, "--phase-suffixes", "_pat", "_mat", "-o", ours)

    assert read_pairs_body(reference) == read_parquet_body(ours)


def test_missing_tag_columns_are_reported(tmp_path, xb_pairs):
    """XA mode on an XB-tagged file must say which columns are missing."""
    with pytest.raises(AssertionError, match="XA1"):
        run_cli(
            "phase", xb_pairs, "--phase-suffixes", "_pat", "_mat",
            "--tag-mode", "XA", "-o", tmp_path / "out.pairs"
        )

# -*- coding: utf-8 -*-
"""split against pairtools split.

A .pairsam packs whole SAM records into the sam1/sam2 columns, with fields
joined by \\x19 and records by \\x19NEXT_SAM\\x19. There is no .pairsam fixture
in the repo and `parse` is not ported yet, so the fixture here builds those
columns the way `pairtools parse` does: one record per side for a simple pair,
several for a walk, and `.` for a side with no alignment.
"""

import random

import pytest
from pairtools.lib import pairsam_format

from conftest import (
    read_pairs_body,
    read_pairs_header,
    read_parquet_body,
    run_cli,
    run_pairtools,
    write_pairs,
)

CHROMS = ["chr1", "chr2", "chr3"]

HEADER = (
    ["## pairs format v1.0.0", "#shape: upper triangle"]
    + ["#chromsize: {} 1000000".format(c) for c in CHROMS]
    + ["#samheader: @HD\tVN:1.6\tSO:unknown"]
    + ["#samheader: @SQ\tSN:{}\tLN:1000000".format(c) for c in CHROMS]
    + ["#samheader: @PG\tID:bwa\tPN:bwa\tVN:0.7.17\tCL:bwa mem -SP genome reads.fq"]
    + [
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type "
        "sam1 sam2"
    ]
)


def sam_record(read_id, flag, chrom, pos, strand):
    """One SAM line, packed as a .pairsam field packs it."""
    fields = [
        read_id, str(flag), chrom, str(pos), "60", "100M", "*", "0", "0",
        "A" * 100, "I" * 100, "NM:i:0", "AS:i:100",
    ]
    return pairsam_format.SAM_SEP.join(fields)


def sam_field(records):
    if not records:
        return "."
    return pairsam_format.INTER_SAM_SEP.join(records)


@pytest.fixture
def pairsam(tmp_path):
    rng = random.Random(11)
    rows = []
    for i in range(300):
        read_id = "r{}".format(i)
        c1, c2 = rng.choice(CHROMS), rng.choice(CHROMS)
        p1, p2 = rng.randrange(1, 1000000), rng.randrange(1, 1000000)
        s1, s2 = rng.choice("+-"), rng.choice("+-")

        kind = rng.choice(["simple", "simple", "walk", "unmapped"])
        if kind == "unmapped":
            c1, p1, s1 = "!", 0, "-"
            rows.append(
                (read_id, c1, p1, c2, p2, s1, s2, "NU",
                 sam_field([]), sam_field([sam_record(read_id, 16, c2, p2, s2)]))
            )
            continue

        # a walk puts several supplementary alignments on one side
        n1 = 2 if kind == "walk" else 1
        rows.append(
            (read_id, c1, p1, c2, p2, s1, s2, "UU",
             sam_field([
                 sam_record(read_id, 65 if j else 1, c1, p1 + j, s1)
                 for j in range(n1)
             ]),
             sam_field([sam_record(read_id, 129, c2, p2, s2)]))
        )
    return write_pairs(tmp_path / "in.pairsam", HEADER, rows)


def read(path):
    with open(path) as f:
        return f.read()


def without_pg(lines):
    """Drop the @PG records, which name the tool that produced the file."""
    return [l for l in lines if "@PG\tID:pairtools" not in l]


def sam_lines(path):
    """The SAM output, minus the @PG records our tools and pairtools' differ on."""
    with open(path) as f:
        return without_pg([l.rstrip("\n") for l in f])


def test_split_pairs_matches_pairtools(tmp_path, pairsam):
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"

    run_pairtools("split", "--output-pairs", reference, pairsam)
    run_cli("split", "--output-pairs", ours, pairsam)

    assert read_pairs_body(reference) == read_pairs_body(ours)
    assert without_pg(read_pairs_header(reference)) == without_pg(
        read_pairs_header(ours)
    )


def test_split_drops_the_sam_columns(tmp_path, pairsam):
    ours = tmp_path / "ours.pairs"
    run_cli("split", "--output-pairs", ours, pairsam)

    columns = [l for l in read_pairs_header(ours) if l.startswith("#columns:")][0]
    assert columns.split()[1:] == [
        "readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2",
        "pair_type",
    ]
    assert all(len(l.split("\t")) == 8 for l in read_pairs_body(ours))


def test_split_sam_matches_pairtools(tmp_path, pairsam):
    reference = tmp_path / "ref.sam"
    ours = tmp_path / "ours.sam"

    run_pairtools("split", "--output-sam", reference, pairsam)
    run_cli("split", "--output-sam", ours, pairsam)

    assert sam_lines(reference) == sam_lines(ours)


def test_split_both_outputs_match_pairtools(tmp_path, pairsam):
    ref_pairs, ref_sam = tmp_path / "ref.pairs", tmp_path / "ref.sam"
    our_pairs, our_sam = tmp_path / "ours.pairs", tmp_path / "ours.sam"

    run_pairtools(
        "split", "--output-pairs", ref_pairs, "--output-sam", ref_sam, pairsam
    )
    run_cli("split", "--output-pairs", our_pairs, "--output-sam", our_sam, pairsam)

    assert read_pairs_body(ref_pairs) == read_pairs_body(our_pairs)
    assert sam_lines(ref_sam) == sam_lines(our_sam)


def test_split_from_parquet(tmp_path, pairsam):
    """The sam columns survive a Parquet round trip and still split correctly."""
    as_parquet = tmp_path / "in.parquet"
    run_cli("csv-to-parquet", "-o", as_parquet, pairsam)

    ref_pairs, ref_sam = tmp_path / "ref.pairs", tmp_path / "ref.sam"
    our_pairs, our_sam = tmp_path / "ours.parquet", tmp_path / "ours.sam"

    run_pairtools(
        "split", "--output-pairs", ref_pairs, "--output-sam", ref_sam, pairsam
    )
    run_cli(
        "split", "--output-pairs", our_pairs, "--output-sam", our_sam, as_parquet
    )

    assert read_pairs_body(ref_pairs) == read_parquet_body(our_pairs)
    assert sam_lines(ref_sam) == sam_lines(our_sam)


def test_split_sam_to_stdout(pairsam):
    out = run_cli("split", "--output-sam", "-", pairsam)
    assert out.startswith("@HD\tVN:1.6")
    assert "\x19" not in out


def test_split_needs_an_output(tmp_path, pairsam):
    with pytest.raises(AssertionError, match="At least one output"):
        run_cli("split", pairsam)


def test_split_rejects_pairs_on_stdout(tmp_path, pairsam):
    with pytest.raises(AssertionError, match="cannot be written to stdout"):
        run_cli("split", "--output-pairs", "-", pairsam)


def test_split_rejects_a_lone_sam_column(tmp_path):
    header = [
        "## pairs format v1.0.0",
        "#chromsize: chr1 1000",
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type sam1",
    ]
    path = write_pairs(
        tmp_path / "lone.pairsam", header,
        [("r1", "chr1", 1, "chr1", 2, "+", "+", "UU", ".")],
    )
    with pytest.raises(AssertionError, match="only one sam entry"):
        run_cli("split", "--output-pairs", tmp_path / "out.pairs", path)

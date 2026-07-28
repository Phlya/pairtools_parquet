# -*- coding: utf-8 -*-
"""parse and parse2 against pairtools parse and parse2.

The parsing itself is pairtools' `streaming_classify`, called unchanged; what
is ported is where its output goes. So the fixture's job is to reach the
branches that decide what a row looks like: simple pairs, a walk with a
supplementary alignment, an unmapped read pair, and a multi-mapper below the
default MAPQ threshold.
"""

import random
import subprocess
import sys

import pytest

from conftest import read_pairs_body, read_pairs_header, run_cli

CHROMS = [("chr1", 200000), ("chr2", 150000), ("chr3", 100000)]


def sam_line(read_id, flag, chrom, pos, mapq, cigar, seq, tags=()):
    return "\t".join(
        [read_id, str(flag), chrom, str(pos), str(mapq), cigar, "*", "0", "0",
         seq, "I" * len(seq)]
        + list(tags)
    )


@pytest.fixture
def sam_path(tmp_path):
    rng = random.Random(5)
    lines = ["@HD\tVN:1.6\tSO:unknown"]
    lines += ["@SQ\tSN:{}\tLN:{}".format(c, n) for c, n in CHROMS]
    lines += ["@PG\tID:bwa\tPN:bwa\tVN:0.7.17\tCL:bwa mem -SP genome r1.fq r2.fq"]

    for i in range(400):
        rid = "read{:04d}".format(i)
        (c1, l1), (c2, l2) = rng.choice(CHROMS), rng.choice(CHROMS)
        p1, p2 = rng.randrange(1, l1 - 200), rng.randrange(1, l2 - 200)
        kind = rng.choice(["pair", "pair", "pair", "walk", "unmapped", "lowmapq"])

        if kind == "unmapped":
            lines.append(sam_line(rid, 77, "*", 0, 0, "*", "A" * 80))
            lines.append(sam_line(rid, 141, "*", 0, 0, "*", "C" * 80))
            continue

        tags = ["NM:i:{}".format(rng.randint(0, 3)),
                "AS:i:{}".format(rng.randint(60, 80)),
                "XS:i:{}".format(rng.choice([0, 20, 70]))]
        lines.append(
            sam_line(rid, rng.choice([65, 81]), c1, p1,
                     0 if kind == "lowmapq" else rng.choice([30, 60]),
                     "80M", "A" * 80, tags)
        )
        if kind == "walk":
            # a supplementary alignment on R1 is what makes this a walk
            lines.append(
                sam_line(rid, 2113, c2, p2, 60, "40S40M", "G" * 80,
                         ["NM:i:0", "AS:i:40", "XS:i:0",
                          "SA:Z:{},{},+,80M,60,0;".format(c1, p1)])
            )
        lines.append(
            sam_line(rid, rng.choice([129, 145]), c2, p2, rng.choice([30, 60]),
                     "80M", "T" * 80, tags)
        )

    path = tmp_path / "reads.sam"
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def chromsizes_path(tmp_path):
    path = tmp_path / "genome.chrom.sizes"
    path.write_text("".join("{}\t{}\n".format(c, n) for c, n in CHROMS))
    return path


def reference(*args):
    """Run pairtools, skipping the test if it crashes rather than fails.

    pairtools 1.1.2's compiled `AlignedSegmentPairtoolized` segfaults against
    pysam 0.24 on any path that touches a read's sequence -- `--drop-seq` and
    `--add-columns seq`. That is an ABI mismatch in the environment, not
    something this package can be right or wrong about, so those option
    combinations have no reference to compare against here.
    """
    cmd = ["pairtools"] + [str(a) for a in args]
    proc = subprocess.run(cmd, capture_output=True)
    if proc.returncode < 0:
        pytest.skip(
            "pairtools crashed with signal {} on `{}`; pysam ABI mismatch".format(
                -proc.returncode, " ".join(cmd)
            )
        )
    if proc.returncode != 0:
        raise AssertionError(
            "command failed: {}\n{}".format(" ".join(cmd), proc.stderr.decode())
        )
    return proc.stdout.decode()


def without_pg(lines):
    return [l for l in lines if "@PG\tID:pairtools" not in l]


def assert_matches(ref, ours):
    assert read_pairs_body(ref) == read_pairs_body(ours)
    assert without_pg(read_pairs_header(ref)) == without_pg(read_pairs_header(ours))


PARSE_OPTIONS = [
    [],
    ["--drop-sam"],
    ["--drop-sam", "--add-columns", "mapq,pos5,pos3,cigar,read_len,matched_bp,"
     "algn_ref_span,algn_read_span,dist_to_5,dist_to_3,mismatches"],
    ["--drop-sam", "--add-columns", "NM,AS,XS"],
    ["--drop-sam", "--add-pair-index"],
    ["--drop-sam", "--walks-policy", "all"],
    ["--drop-sam", "--walks-policy", "mask"],
    ["--drop-sam", "--no-flip"],
    ["--drop-sam", "--report-alignment-end", "3"],
    ["--drop-sam", "--drop-readid"],
    ["--drop-sam", "--min-mapq", "30"],
    ["--drop-sam", "--max-molecule-size", "200"],
]


@pytest.mark.parametrize(
    "options", PARSE_OPTIONS, ids=lambda o: " ".join(o) or "defaults"
)
def test_parse_matches_pairtools(tmp_path, sam_path, chromsizes_path, options):
    ref = tmp_path / "ref.pairsam"
    ours = tmp_path / "ours.pairsam"

    reference("parse", "-c", chromsizes_path, "-o", ref, *options, sam_path)
    run_cli("parse", "-c", chromsizes_path, "-o", ours, *options, sam_path)

    assert_matches(ref, ours)
    # `--walks-policy all` reports more than one pair for a walk, so only the
    # count being non-trivial is common to every case
    assert len(read_pairs_body(ours)) >= 400


def test_parse_reports_one_pair_per_read_by_default(
    tmp_path, sam_path, chromsizes_path
):
    ours = tmp_path / "ours.pairsam"
    run_cli("parse", "-c", chromsizes_path, "-o", ours, sam_path)

    body = read_pairs_body(ours)
    assert len(body) == 400
    # the fixture must reach the mapped, unmapped and multi-mapping branches
    assert {l.split("\t")[7] for l in body} >= {"UU", "NN", "MU"}


PARSE2_OPTIONS = [
    [],
    ["--drop-sam"],
    ["--drop-sam", "--add-pair-index"],
    ["--drop-sam", "--expand"],
    ["--drop-sam", "--expand", "--max-expansion-depth", "2"],
    ["--drop-sam", "--report-position", "junction"],
    ["--drop-sam", "--report-orientation", "walk"],
    ["--drop-sam", "--single-end"],
    ["--drop-sam", "--flip"],
    ["--drop-sam", "--add-columns", "NM,AS"],
]


@pytest.mark.parametrize(
    "options", PARSE2_OPTIONS, ids=lambda o: " ".join(o) or "defaults"
)
def test_parse2_matches_pairtools(tmp_path, sam_path, chromsizes_path, options):
    ref = tmp_path / "ref.pairsam"
    ours = tmp_path / "ours.pairsam"

    reference("parse2", "-c", chromsizes_path, "-o", ref, *options, sam_path)
    run_cli("parse2", "-c", chromsizes_path, "-o", ours, *options, sam_path)

    assert_matches(ref, ours)


def test_parse_side_outputs_match_pairtools(tmp_path, sam_path, chromsizes_path):
    """--output-stats and --output-parsed-alignments are plain text either way."""
    ref_stats, our_stats = tmp_path / "ref.stats", tmp_path / "ours.stats"
    ref_algn, our_algn = tmp_path / "ref.algn", tmp_path / "ours.algn"

    reference(
        "parse", "-c", chromsizes_path, "--drop-sam",
        "--output-stats", ref_stats, "--output-parsed-alignments", ref_algn,
        "-o", tmp_path / "ref.pairs", sam_path
    )
    run_cli(
        "parse", "-c", chromsizes_path, "--drop-sam",
        "--output-stats", our_stats, "--output-parsed-alignments", our_algn,
        "-o", tmp_path / "ours.pairs", sam_path
    )

    assert ref_stats.read_text() == our_stats.read_text()
    assert ref_algn.read_text() == our_algn.read_text()


def test_parse_to_parquet_equals_parse_then_convert(
    tmp_path, sam_path, chromsizes_path
):
    """The whole point: parsing straight to Parquet skips a text round trip.

    The result has to be what the round trip would have produced, sam columns
    and their \\x19 separators included.
    """
    as_text = tmp_path / "ours.pairsam"
    direct = tmp_path / "direct.parquet"
    run_cli("parse", "-c", chromsizes_path, "-o", as_text, sam_path)
    run_cli("parse", "-c", chromsizes_path, "-o", direct, sam_path)

    back = tmp_path / "back.pairsam"
    run_cli("parquet-to-csv", "-o", back, direct)

    assert read_pairs_body(as_text) == read_pairs_body(back)
    assert without_pg(read_pairs_header(as_text)) == without_pg(
        read_pairs_header(back)
    )


def test_parse_rejects_stdout(tmp_path, sam_path, chromsizes_path):
    with pytest.raises(AssertionError, match="output path is required"):
        run_cli("parse", "-c", chromsizes_path, "-o", "", sam_path)


def test_parse_rejects_an_unknown_extra_column(tmp_path, sam_path, chromsizes_path):
    with pytest.raises(AssertionError, match="not a valid extra column"):
        run_cli(
            "parse", "-c", chromsizes_path, "--add-columns", "nonsense",
            "-o", tmp_path / "out.pairs", sam_path
        )


def test_extra_columns_get_their_declared_types(tmp_path, sam_path, chromsizes_path):
    """`mapq1` is an int because pairtools declares `mapq` as one.

    DTYPES_EXTRA_COLUMNS is keyed by the base name, so a lookup on the full
    column name misses -- which used to leave every --add-columns column a
    string, and `select 'mapq1>=30'` comparing str to int.
    """
    import pyarrow.parquet as pq

    out = tmp_path / "out.parquet"
    run_cli(
        "parse", "-c", chromsizes_path, "--drop-sam", "--add-columns",
        "mapq,cigar,read_len", "-o", out, sam_path
    )

    schema = pq.read_schema(str(out))
    types = dict(zip(schema.names, (str(t) for t in schema.types)))
    assert types["mapq1"] == "int32"
    assert types["read_len2"] == "int32"
    assert types["cigar1"] == "string"

    selected = tmp_path / "selected.pairs"
    run_cli("select", "mapq1>=30", "-o", selected, out)

    ref = tmp_path / "ref.pairs"
    ref_input = tmp_path / "ref_in.pairs"
    reference(
        "parse", "-c", chromsizes_path, "--drop-sam", "--add-columns",
        "mapq,cigar,read_len", "-o", ref_input, sam_path
    )
    reference("select", "mapq1>=30", "-o", ref, ref_input)

    assert len(read_pairs_body(selected)) == len(read_pairs_body(ref))
    assert len(read_pairs_body(selected)) > 0

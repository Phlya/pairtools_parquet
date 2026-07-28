# -*- coding: utf-8 -*-
"""Parity of pairtools_parquet against pairtools itself.

The whole premise of this package is "same results, faster", so every tool is
checked against the upstream implementation rather than against a hand-written
expectation: run both, convert whatever came out back to text, and require the
bodies to be identical. The header may differ only by the @PG provenance record
each tool appends to record itself.
"""

import shutil

import pytest

from conftest import (
    read_pairs_body,
    read_pairs_header,
    read_parquet_body,
    run_cli,
    run_pairtools,
    write_pairs,
)

SORT_KEYS = ["chrom1", "chrom2", "pos1", "pos2", "pair_type"]


def assert_headers_agree(reference_header, our_header):
    """Headers must match except for each tool's own @PG record."""

    def without_provenance(header):
        return [
            l
            for l in header
            if not (l.startswith("#samheader: @PG") and "pairtools" in l and "sort" in l)
        ]

    assert without_provenance(reference_header) == without_provenance(our_header)


def test_sort_matches_pairtools(tmp_path, mock_pairs_path):
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"

    run_pairtools("sort", "-o", reference, mock_pairs_path)
    run_cli("sort", "-o", ours, "--compress-program", "none", mock_pairs_path)

    assert read_pairs_body(reference) == read_pairs_body(ours)
    assert_headers_agree(read_pairs_header(reference), read_pairs_header(ours))


def test_sort_matches_pairtools_through_parquet(tmp_path, mock_pairs_path):
    """Routing through Parquet must not change the result."""
    reference = tmp_path / "ref.pairs"
    intermediate = tmp_path / "sorted.parquet"
    ours = tmp_path / "ours.pairs"

    run_pairtools("sort", "-o", reference, mock_pairs_path)
    run_cli("sort", "-o", intermediate, mock_pairs_path)
    run_cli("sort", "-o", ours, "--compress-program", "none", intermediate)

    assert read_pairs_body(reference) == read_pairs_body(ours)
    assert read_pairs_body(reference) == read_parquet_body(intermediate)


def test_sort_orders_tied_pairs_by_pair_type(tmp_path):
    """Pairs identical but for pair_type must order lexicographically.

    Regression test: chrom/strand/pair_type used to be sorted as DuckDB ENUMs
    declared in itertools.product order, so `UU` sorted before `DD` where
    pairtools -- which sorts bytewise -- puts `DD` first. pair_type is the last
    sort key, so this only shows up on pairs tied on every position.
    """
    header = [
        "## pairs format v1.0.0",
        "#chromsize: chr1 1000",
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type",
    ]
    rows = [
        ("r_uu", "chr1", 1, "chr1", 50, "+", "+", "UU"),
        ("r_dd", "chr1", 1, "chr1", 50, "+", "+", "DD"),
        ("r_ww", "chr1", 1, "chr1", 50, "+", "+", "WW"),
        ("r_nu", "chr1", 1, "chr1", 50, "+", "+", "NU"),
    ]
    source = write_pairs(tmp_path / "tied.pairs", header, rows)

    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"
    run_pairtools("sort", "-o", reference, source)
    run_cli("sort", "-o", ours, "--compress-program", "none", source)

    body = read_pairs_body(ours)
    assert [l.split("\t")[7] for l in body] == ["DD", "NU", "UU", "WW"]
    assert read_pairs_body(reference) == body


def test_sort_handles_chromosome_missing_from_header(tmp_path):
    """A chromosome the header never declares must not be fatal.

    The chromosome ENUM that makes the sort fast has a fixed domain taken from
    the header, and DuckDB raises rather than nulling on a value outside it, so
    this has to fall back to comparing strings.
    """
    header = [
        "## pairs format v1.0.0",
        "#chromsize: chr1 1000",
        "#chromsize: chr2 1000",
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type",
    ]
    rows = [
        ("r1", "chr2", 5, "chr2", 9, "+", "+", "UU"),
        ("r2", "chrUNDECLARED", 3, "chr1", 7, "+", "-", "UU"),
        ("r3", "chr1", 1, "chr1", 2, "-", "+", "UU"),
    ]
    source = write_pairs(tmp_path / "undeclared.pairs", header, rows)

    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"
    run_pairtools("sort", "-o", reference, source)
    run_cli("sort", "-o", ours, "--compress-program", "none", source)

    assert read_pairs_body(reference) == read_pairs_body(ours)


SELECT_CONDITIONS = [
    'pair_type == "UU"',
    '(pair_type == "RU") or (pair_type == "UR") or (pair_type == "UU")',
    'csv_match(pair_type, "RU,UR,UU")',
    'wildcard_match(pair_type, "*U")',
    # '?' is fnmatch's single-character wildcard. The SQL translation this
    # replaces passed it through to LIKE, where it is a literal, and so
    # silently selected nothing.
    'wildcard_match(chrom1, "chr?")',
    'regex_match(chrom1, "chr[0-9]+")',
    "chrom1 == chrom2",
    "chrom1 == chrom2 and abs(pos1 - pos2) < 10",
    'not (pair_type == "UU")',
    'pair_type in ("UU", "DD")',
    # a Python method call: no SQL translation of the condition language can
    # express this, which is why the condition is evaluated by pairtools
    'readID.endswith("9")',
    "pos1 < pos2 < 100",
]


@pytest.mark.parametrize("condition", SELECT_CONDITIONS)
def test_select_matches_pairtools(tmp_path, mock_pairs_path, condition):
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"

    run_pairtools("select", condition, "-o", reference, mock_pairs_path)
    run_cli(
        "select", condition, mock_pairs_path, "-o", ours, "--compress-program", "none"
    )

    assert read_pairs_body(reference) == read_pairs_body(ours)


@pytest.mark.parametrize("condition", SELECT_CONDITIONS)
def test_select_matches_pairtools_via_parquet(tmp_path, mock_pairs_path, condition):
    """Selecting from Parquet must give the same rows as selecting from text."""
    reference = tmp_path / "ref.pairs"
    as_parquet = tmp_path / "in.parquet"
    ours = tmp_path / "ours.parquet"

    run_pairtools("select", condition, "-o", reference, mock_pairs_path)
    run_cli("csv-to-parquet", "-o", as_parquet, mock_pairs_path)
    run_cli("select", condition, as_parquet, "-o", ours)

    assert read_pairs_body(reference) == read_parquet_body(ours)


def test_select_string_literal_containing_operator_words(tmp_path):
    """A string literal in CONDITION must not be rewritten.

    The SQL translation this replaces uppercased ' and ' / ' or ' everywhere in
    the condition, including inside quoted values, turning a readID of
    "a and b" into a comparison against 'a AND b'.
    """
    header = [
        "## pairs format v1.0.0",
        "#chromsize: chr1 1000",
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type",
    ]
    rows = [
        ("a and b", "chr1", 1, "chr1", 2, "+", "+", "UU"),
        ("a AND b", "chr1", 3, "chr1", 4, "+", "+", "UU"),
        ("plain", "chr1", 5, "chr1", 6, "+", "+", "UU"),
    ]
    source = write_pairs(tmp_path / "literal.pairs", header, rows)
    condition = 'readID == "a and b"'

    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"
    run_pairtools("select", condition, "-o", reference, source)
    run_cli(
        "select", condition, source, "-o", ours, "--compress-program", "none"
    )

    body = read_pairs_body(ours)
    assert [l.split("\t")[0] for l in body] == ["a and b"]
    assert read_pairs_body(reference) == body


def test_select_output_rest_partitions_the_input(tmp_path, mock_pairs_path):
    selected = tmp_path / "sel.pairs"
    rest = tmp_path / "rest.pairs"
    run_cli(
        "select",
        'pair_type == "UU"',
        mock_pairs_path,
        "-o",
        selected,
        "--output-rest",
        rest,
        "--compress-program",
        "none",
    )

    body = read_pairs_body(mock_pairs_path)
    assert sorted(read_pairs_body(selected) + read_pairs_body(rest)) == sorted(body)
    assert all(l.split("\t")[7] == "UU" for l in read_pairs_body(selected))
    assert all(l.split("\t")[7] != "UU" for l in read_pairs_body(rest))


def test_select_startup_code(tmp_path, mock_pairs_path):
    """--startup-code defines helpers usable in CONDITION, as in pairtools."""
    # kept to one line: both tools record the full command line in a @PG
    # header record, and an embedded newline would split that header line
    startup = "def is_cis(c1, c2): return c1 == c2"
    condition = "is_cis(chrom1, chrom2)"

    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"
    run_pairtools(
        "select", condition, "--startup-code", startup, "-o", reference, mock_pairs_path
    )
    run_cli(
        "select",
        condition,
        mock_pairs_path,
        "-o",
        ours,
        "--startup-code",
        startup,
        "--compress-program",
        "none",
    )

    assert read_pairs_body(reference) == read_pairs_body(ours)
    assert len(read_pairs_body(ours)) > 0


def test_select_matching_nothing_still_feeds_sort(tmp_path, mock_pairs_path):
    """An empty selection is a valid .pairs file and must stay usable."""
    empty = tmp_path / "empty.pairs"
    run_cli(
        "select",
        'pair_type == "NOPE"',
        mock_pairs_path,
        "-o",
        empty,
        "--compress-program",
        "none",
    )
    assert read_pairs_body(empty) == []

    sorted_empty = tmp_path / "empty.sorted.parquet"
    run_cli("sort", "-o", sorted_empty, empty)
    assert read_parquet_body(sorted_empty) == []


def split_into_sorted_inputs(tmp_path, mock_pairs_path):
    """Split the mock file into two separately-sorted .pairs files."""
    header = read_pairs_header(mock_pairs_path)
    body = read_pairs_body(mock_pairs_path)
    halves = []
    for i, rows in enumerate((body[: len(body) // 2], body[len(body) // 2 :])):
        raw = tmp_path / "part{}.pairs".format(i)
        with open(raw, "w") as f:
            f.write("".join(l + "\n" for l in header + rows))
        part = tmp_path / "part{}.sorted.pairs".format(i)
        run_pairtools("sort", "-o", part, raw)
        halves.append(part)
    return halves


def test_merge_matches_pairtools(tmp_path, mock_pairs_path):
    parts = split_into_sorted_inputs(tmp_path, mock_pairs_path)
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"

    run_pairtools("merge", "-o", reference, *parts)
    run_cli("merge", "-o", ours, *parts)

    assert read_pairs_body(reference) == read_pairs_body(ours)


def test_merge_header_matches_pairtools(tmp_path, mock_pairs_path):
    """Merging must not re-mark the header as sorted.

    `headerops.mark_header_as_sorted` rewrites `#chromosomes: a b c` as
    `#chromosomes: : a b c`; `pairtools merge` never calls it, so neither may
    we, or merged headers drift from upstream's.
    """
    parts = split_into_sorted_inputs(tmp_path, mock_pairs_path)
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"

    run_pairtools("merge", "-o", reference, *parts)
    run_cli("merge", "-o", ours, *parts)

    def normalize(header):
        # the @PG records differ by tool name and command line
        return [l for l in header if not l.startswith("#samheader: @PG")]

    assert normalize(read_pairs_header(reference)) == normalize(
        read_pairs_header(ours)
    )
    assert "#chromosomes: : chr1 chr2 chr3" not in read_pairs_header(ours)


def test_merge_concatenate_matches_pairtools(tmp_path, mock_pairs_path):
    parts = split_into_sorted_inputs(tmp_path, mock_pairs_path)
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.pairs"

    run_pairtools("merge", "--concatenate", "-o", reference, *parts)
    run_cli("merge", "--concatenate", "-o", ours, *parts)

    assert read_pairs_body(reference) == read_pairs_body(ours)


def test_merge_mixed_input_formats(tmp_path, mock_pairs_path):
    """Inputs may be a mix of .pairs and .parquet."""
    parts = split_into_sorted_inputs(tmp_path, mock_pairs_path)
    as_parquet = tmp_path / "part0.parquet"
    run_cli("csv-to-parquet", "-o", as_parquet, parts[0])

    reference = tmp_path / "ref.pairs"
    ours = tmp_path / "ours.parquet"
    run_pairtools("merge", "-o", reference, *parts)
    run_cli("merge", "-o", ours, as_parquet, parts[1])

    assert read_pairs_body(reference) == read_parquet_body(ours)


@pytest.mark.parametrize("max_nmerge", [2, 3, 7])
def test_merge_staging_does_not_change_the_result(
    tmp_path, mock_pairs_path, max_nmerge
):
    """--max-nmerge bounds how many inputs are open at once, nothing else.

    Beyond it the merge runs in rounds through temporary files, so the header
    is decided once from the real inputs and handed to every round -- otherwise
    each round would append its own @PG record and a staged merge would not
    look like a single-pass one.
    """
    parts = split_into_sorted_inputs(tmp_path, mock_pairs_path)
    many = []
    for i in range(8):
        copy = tmp_path / "copy{}.pairs".format(i)
        shutil.copy(parts[i % len(parts)], copy)
        many.append(copy)

    one_pass = tmp_path / "one.pairs"
    staged = tmp_path / "staged.pairs"
    run_cli("merge", "--max-nmerge", "0", "-o", one_pass, *many)
    run_cli("merge", "--max-nmerge", str(max_nmerge), "-o", staged, *many)

    assert read_pairs_body(one_pass) == read_pairs_body(staged)
    # The headers match except for the recorded command line, which honestly
    # differs -- it is the invocation, and the invocations differ.
    assert _without_command_lines(one_pass) == _without_command_lines(staged)


def _without_command_lines(path):
    return [line.split("\tCL:")[0] for line in read_pairs_header(path)]


def test_merge_single_input_leaves_header_alone(tmp_path, mock_pairs_path):
    """As in pairtools, a lone input is passed through without a new @PG."""
    parts = split_into_sorted_inputs(tmp_path, mock_pairs_path)
    ours = tmp_path / "ours.pairs"
    run_cli("merge", "-o", ours, parts[0])

    assert read_pairs_header(ours) == read_pairs_header(parts[0])
    assert read_pairs_body(ours) == read_pairs_body(parts[0])


@pytest.mark.parametrize(
    "command", [("sort",), ("merge",), ("parquet-to-csv",), ("csv-to-parquet",)]
)
def test_text_output_is_compressed_only_when_the_extension_says_so(
    tmp_path, mock_pairs_path, command
):
    """`-o out.pairs` must be text, whatever --compress-program defaults to.

    Compression follows the extension, as pairtools' auto_open does; the flag
    only chooses which compressor.
    """
    source = mock_pairs_path
    if command[0] == "parquet-to-csv":
        source = tmp_path / "in.parquet"
        run_cli("csv-to-parquet", "-o", source, mock_pairs_path)

    plain = tmp_path / "out.pairs"
    run_cli(*command, "-o", plain, source)
    with open(plain, "rb") as f:
        assert f.read(2) != b"\x1f\x8b", "plain .pairs output was gzip-compressed"
    assert read_pairs_header(plain)[0] == "## pairs format v1.0.0"

    compressed = tmp_path / "out.pairs.gz"
    run_cli(*command, "-o", compressed, source)
    with open(compressed, "rb") as f:
        assert f.read(2) == b"\x1f\x8b", ".gz output was not compressed"


@pytest.mark.parametrize("compress_program", ["none", "gzip"])
def test_sort_to_compressed_output(tmp_path, mock_pairs_path, compress_program):
    suffix = ".pairs" if compress_program == "none" else ".pairs.gz"
    reference = tmp_path / "ref.pairs"
    ours = tmp_path / ("ours" + suffix)

    run_pairtools("sort", "-o", reference, mock_pairs_path)
    run_cli(
        "sort", "-o", ours, "--compress-program", compress_program, mock_pairs_path
    )

    # reading back through our own reader exercises the decompression path
    from pairtools_parquet.lib.arrowio import open_pairs

    _, reader = open_pairs(str(ours))
    rows = sum(batch.num_rows for batch in reader)
    assert rows == len(read_pairs_body(reference))

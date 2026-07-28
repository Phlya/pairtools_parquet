# -*- coding: utf-8 -*-
"""Arrow IPC as a third format, alongside .pairs and .parquet.

Parquet writes its footer last, so it cannot be read until it is complete --
which is why a `.parquet` cannot go through a pipe. The Arrow IPC *stream*
format has no footer: each batch is self-describing and can be read as it
arrives. That makes it the one binary format the tools can pipe between them,
and it carries the .pairs header in its schema metadata exactly as Parquet
does, so nothing is lost on the way.

`-.arrow` means the standard stream carrying Arrow, the way `-` means it
carrying text. On input nothing has to be declared: a stream beginning with
Arrow's continuation token is Arrow, and .pairs text always begins with '#'.
"""

import os
import subprocess
import sys

import pytest

from conftest import read_pairs_body, read_pairs_header, run_cli

TESTDIR = os.path.dirname(os.path.realpath(__file__))
DATADIR = os.path.join(TESTDIR, "data")


def body(text):
    return [l for l in text.splitlines() if l and not l.startswith("#")]


def pipeline(steps, cwd):
    """Run `steps` joined by shell pipes and return the finished process."""
    joined = " | ".join(
        "{} -m pairtools_parquet {}".format(sys.executable, step) for step in steps
    )
    return subprocess.run(joined, shell=True, capture_output=True, cwd=str(cwd))


def test_arrow_file_round_trips(tmp_path, mock_pairs_path):
    """.pairs -> .arrow -> .pairs must change nothing but the @PG record."""
    as_arrow = tmp_path / "mid.arrow"
    back = tmp_path / "back.pairs"
    run_cli("markasdup", "-o", as_arrow, mock_pairs_path)
    run_cli("markasdup", "-o", back, as_arrow)

    reference = tmp_path / "ref.pairs"
    run_cli("markasdup", "-o", reference, mock_pairs_path)
    run_cli("markasdup", "-o", reference, reference)

    assert read_pairs_body(back) == read_pairs_body(reference)


def test_arrow_carries_the_header(tmp_path, mock_pairs_path):
    """The header rides in the schema metadata, as it does in Parquet.

    Compared as a set, and against Parquet: both binary formats put the header
    through the same metadata round trip, which groups the lines by key and so
    moves the `#samheader:` records after the `#chromsize:` ones. That is the
    existing behaviour of Parquet output, not something Arrow introduces, so
    what is checked here is that Arrow loses nothing Parquet keeps.
    """
    from_arrow = tmp_path / "from_arrow.pairs"
    from_parquet = tmp_path / "from_parquet.pairs"
    for extension, destination in [("arrow", from_arrow), ("parquet", from_parquet)]:
        binary = tmp_path / ("mid." + extension)
        run_cli("markasdup", "-o", binary, mock_pairs_path)
        run_cli("markasdup", "-o", destination, binary)

    def without_provenance(path):
        return [l for l in read_pairs_header(path) if "@PG" not in l]

    assert set(without_provenance(from_arrow)) == set(
        without_provenance(mock_pairs_path)
    )
    assert without_provenance(from_arrow) == without_provenance(from_parquet)


def test_a_tool_that_reads_the_header_first_can_read_an_arrow_file(
    tmp_path, mock_pairs_path
):
    """`sort` resolves its keys from the header before reading any rows."""
    as_arrow = tmp_path / "mid.arrow"
    run_cli("markasdup", "-o", as_arrow, mock_pairs_path)

    ours = tmp_path / "sorted.pairs"
    reference = tmp_path / "ref.pairs"
    run_cli("sort", "-o", ours, as_arrow)
    run_cli("markasdup", "-o", reference, mock_pairs_path)
    run_cli("sort", "-o", reference, reference)

    assert read_pairs_body(ours) == read_pairs_body(reference)


def test_arrow_through_a_pipe_matches_text_through_a_pipe(
    tmp_path, mock_pairs_path
):
    """The format on the wire must not change the answer."""
    steps_text = [
        "select 'pair_type==\"UU\"' -o - {}".format(mock_pairs_path),
        "sort -o - -",
        "markasdup -o - -",
    ]
    steps_arrow = [
        "select 'pair_type==\"UU\"' -o -.arrow {}".format(mock_pairs_path),
        "sort -o -.arrow -",
        "markasdup -o - -",
    ]

    as_text = pipeline(steps_text, tmp_path)
    as_arrow = pipeline(steps_arrow, tmp_path)

    assert as_text.returncode == 0, as_text.stderr.decode()
    assert as_arrow.returncode == 0, as_arrow.stderr.decode()
    assert body(as_arrow.stdout.decode()) == body(as_text.stdout.decode())
    assert body(as_text.stdout.decode()), "the pipeline produced no pairs"


def test_arrow_on_stdout_is_not_a_file(tmp_path, mock_pairs_path):
    proc = subprocess.run(
        [sys.executable, "-m", "pairtools_parquet", "markasdup",
         "-o", "-.arrow", str(mock_pairs_path)],
        capture_output=True, cwd=str(tmp_path),
    )

    assert proc.returncode == 0, proc.stderr.decode()
    assert proc.stdout.startswith(b"\xff\xff\xff\xff"), "not an Arrow stream"
    assert not (tmp_path / "-.arrow").exists()
    assert not (tmp_path / "-").exists()


def test_an_arrow_stream_is_recognised_without_being_declared(
    tmp_path, mock_pairs_path
):
    """Input needs no flag: Arrow starts with its continuation token, .pairs
    with '#'."""
    as_arrow = tmp_path / "mid.arrow"
    run_cli("markasdup", "-o", as_arrow, mock_pairs_path)

    out = tmp_path / "out.pairs"
    with open(as_arrow, "rb") as f:
        proc = subprocess.run(
            [sys.executable, "-m", "pairtools_parquet", "sort", "-o", str(out), "-"],
            stdin=f, capture_output=True,
        )
    assert proc.returncode == 0, proc.stderr.decode()

    reference = tmp_path / "ref.pairs"
    run_cli("sort", "-o", reference, mock_pairs_path)
    assert len(read_pairs_body(out)) == len(read_pairs_body(reference))


@pytest.mark.parametrize("extension", ["arrow", "ipc"])
def test_both_extensions_work(tmp_path, mock_pairs_path, extension):
    path = tmp_path / ("mid." + extension)
    back = tmp_path / "back.pairs"
    run_cli("markasdup", "-o", path, mock_pairs_path)
    run_cli("markasdup", "-o", back, path)

    assert len(read_pairs_body(back)) == len(read_pairs_body(mock_pairs_path))

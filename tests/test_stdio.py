# -*- coding: utf-8 -*-
"""Reading stdin and writing stdout, so the tools can be composed with pipes.

`-` means stdin as an input and stdout as an output, as it does in pairtools.
Both directions were broken in ways that produced no error: `-o -` wrote a file
literally named `-` and left stdout empty, and any tool that read the header
before opening the input for real -- `sort` resolving its keys, `dedup` its
columns -- found the header already consumed, because a pipe cannot be reopened.

Text only. Parquet keeps its footer at the end of the file, so it cannot be
written to a stream at all; `-o -` therefore always means .pairs.
"""

import os
import subprocess
import sys

import pytest

TESTDIR = os.path.dirname(os.path.realpath(__file__))
DATADIR = os.path.join(TESTDIR, "data")


def run(args, stdin=None, cwd=None, check=True):
    """Run a subcommand, optionally feeding it stdin, and return its stdout."""
    proc = subprocess.run(
        [sys.executable, "-m", "pairtools_parquet"] + [str(a) for a in args],
        input=stdin, capture_output=True, cwd=str(cwd) if cwd else None,
    )
    if check and proc.returncode != 0:
        raise AssertionError(
            "failed: {}\n{}".format(" ".join(str(a) for a in args),
                                    proc.stderr.decode())
        )
    return proc


def body(text):
    return [l for l in text.splitlines() if l and not l.startswith("#")]


@pytest.fixture
def mock_bytes(mock_pairs_path):
    with open(mock_pairs_path, "rb") as f:
        return f.read()


# Each entry is a tool and the arguments that make it do something, chosen so
# the output is pairs rather than statistics.
PIPEABLE = [
    ("select", ['pair_type=="UU"']),
    ("sort", []),
    ("markasdup", []),
    ("flip", ["-c", os.path.join(DATADIR, "mock.chrom.sizes")]),
    ("sample", ["--seed", "1", "0.5"]),
    ("merge", []),
]


@pytest.mark.parametrize("tool,args", PIPEABLE, ids=[t for t, _ in PIPEABLE])
def test_output_to_stdout_matches_output_to_a_file(
    tmp_path, mock_pairs_path, tool, args
):
    """`-o -` must write the pairs to stdout, and not to a file called `-`.

    Bodies are compared, not whole files: every tool records its own command
    line in an `@PG` header record, and `-o -` is not `-o out.pairs`.
    """
    to_file = tmp_path / "out.pairs"
    run([tool] + args + ["-o", to_file, mock_pairs_path])
    with open(to_file) as f:
        expected = f.read()

    proc = run([tool] + args + ["-o", "-", mock_pairs_path], cwd=tmp_path)
    written = proc.stdout.decode()

    assert body(written) == body(expected)
    assert written.startswith("## pairs format"), "header missing from stdout"
    assert not (tmp_path / "-").exists(), "wrote a file named '-' instead of stdout"


@pytest.mark.parametrize("tool,args", PIPEABLE, ids=[t for t, _ in PIPEABLE])
def test_input_from_stdin_matches_input_from_a_file(
    tmp_path, mock_pairs_path, mock_bytes, tool, args
):
    """A tool must give the same answer whether it read a file or a pipe."""
    to_file = tmp_path / "out.pairs"
    run([tool] + args + ["-o", to_file, mock_pairs_path])
    with open(to_file) as f:
        from_file = f.read()

    from_stdin = tmp_path / "piped.pairs"
    run([tool] + args + ["-o", from_stdin, "-"], stdin=mock_bytes)
    with open(from_stdin) as f:
        assert body(f.read()) == body(from_file)


def test_a_pipeline_of_three_tools(tmp_path, mock_pairs_path, mock_bytes):
    """The point of all of it: composing without touching the disk."""
    expected = tmp_path / "staged.pairs"
    step1 = tmp_path / "s1.pairs"
    step2 = tmp_path / "s2.pairs"
    run(["select", 'pair_type=="UU"', "-o", step1, mock_pairs_path])
    run(["sort", "-o", step2, step1])
    run(["markasdup", "-o", expected, step2])

    piped = subprocess.run(
        "{py} -m pairtools_parquet select 'pair_type==\"UU\"' -o - {inp} "
        "| {py} -m pairtools_parquet sort -o - - "
        "| {py} -m pairtools_parquet markasdup -o - -".format(
            py=sys.executable, inp=mock_pairs_path
        ),
        shell=True, capture_output=True, cwd=str(tmp_path),
    )
    assert piped.returncode == 0, piped.stderr.decode()

    with open(expected) as f:
        assert body(piped.stdout.decode()) == body(f.read())


def test_dedup_refuses_stdin_rather_than_emitting_nothing(tmp_path, mock_bytes):
    """dedup needs two passes, which a pipe cannot give it.

    It used to read the spent stream and write an empty file, reporting
    success -- the worst way to be wrong.
    """
    proc = run(["dedup", "-o", tmp_path / "out.pairs", "-"],
               stdin=mock_bytes, check=False)

    assert proc.returncode != 0
    assert "cannot read from stdin" in proc.stderr.decode()


def test_merge_accepts_a_mix_of_stdin_and_files(tmp_path, mock_pairs_path,
                                                mock_bytes):
    out = tmp_path / "merged.pairs"
    run(["merge", "-o", out, mock_pairs_path, "-"], stdin=mock_bytes)

    with open(mock_pairs_path) as f:
        one = len(body(f.read()))
    with open(out) as f:
        assert len(body(f.read())) == 2 * one


def test_stdout_output_is_not_compressed_by_a_compress_program(
    tmp_path, mock_pairs_path
):
    """A stream has no extension, so it is plain text whatever the flag says.

    `pairtools` decides compression by the output's extension too, so `-o -`
    is uncompressed there as well; pipe through a compressor if you want one.
    """
    proc = run(["sort", "--compress-program", "lz4", "-o", "-", mock_pairs_path],
               cwd=tmp_path)

    assert proc.stdout.startswith(b"## pairs format")

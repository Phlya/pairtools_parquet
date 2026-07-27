# -*- coding: utf-8 -*-
import os
import subprocess
import sys

import pyarrow.parquet as pq
import pytest

from pairtools_parquet.lib.duckdb_utils import duckdb_kv_metadata_to_header

TESTDIR = os.path.dirname(os.path.realpath(__file__))
DATADIR = os.path.join(TESTDIR, "data")


@pytest.fixture
def mock_pairs_path():
    return os.path.join(DATADIR, "mock.pairs")


@pytest.fixture
def mock_parquet_path():
    return os.path.join(DATADIR, "mock.parquet")


@pytest.fixture
def mock_chromsizes_path():
    return os.path.join(DATADIR, "mock.chrom.sizes")


def run_cli(*args):
    """Run a pairtools_parquet subcommand in a subprocess.

    Raises with the captured stdout/stderr on a non-zero exit, so a failing
    command reports what actually went wrong instead of a bare CalledProcessError.
    """
    cmd = [sys.executable, "-m", "pairtools_parquet"] + [str(a) for a in args]
    proc = subprocess.run(cmd, capture_output=True)
    if proc.returncode != 0:
        raise AssertionError(
            "command failed: {}\n--- stdout ---\n{}\n--- stderr ---\n{}".format(
                " ".join(cmd), proc.stdout.decode(), proc.stderr.decode()
            )
        )
    return proc.stdout.decode()


def run_pairtools(*args):
    """Run a pairtools subcommand, for comparing our output against upstream."""
    cmd = ["pairtools"] + [str(a) for a in args]
    proc = subprocess.run(cmd, capture_output=True)
    if proc.returncode != 0:
        raise AssertionError(
            "command failed: {}\n--- stdout ---\n{}\n--- stderr ---\n{}".format(
                " ".join(cmd), proc.stdout.decode(), proc.stderr.decode()
            )
        )
    return proc.stdout.decode()


def write_pairs(path, header, rows):
    """Write a small .pairs file from header lines and row tuples."""
    with open(path, "w") as f:
        f.write("".join(line + "\n" for line in header))
        f.write("".join("\t".join(str(v) for v in row) + "\n" for row in rows))
    return path


def read_pairs_header(path):
    with open(path, "r") as f:
        return [l.strip() for l in f if l.startswith("#")]


def read_pairs_body(path):
    with open(path, "r") as f:
        return [l.strip() for l in f if not l.startswith("#") and l.strip()]


def read_parquet_header(path):
    return [l.strip() for l in duckdb_kv_metadata_to_header(path)]


def read_parquet_body(path):
    """Read a .parquet file as tab-separated lines, matching the .pairs layout.

    Lets assertions index fields the same way for both output formats.
    """
    table = pq.read_table(path)
    columns = table.column_names
    lines = [
        "\t".join("" if row[col] is None else str(row[col]) for col in columns)
        for row in table.to_pylist()
    ]
    return [l.strip() for l in lines if not l.startswith("#") and l.strip()]

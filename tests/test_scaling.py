# -*- coding: utf-8 -*-
"""scaling against pairtools scaling."""

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

    run_pairtools("scaling", "-o", reference, *options, pairs)
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

    run_pairtools("scaling", "-o", reference, mock_pairs_path)
    run_cli("scaling", "-o", ours, mock_pairs_path)

    assert read(reference) == read(ours)


def test_scaling_through_parquet(tmp_path, pairs):
    as_parquet = tmp_path / "in.parquet"
    run_cli("csv-to-parquet", "-o", as_parquet, pairs)

    reference = tmp_path / "ref.tsv"
    ours = tmp_path / "ours.tsv"
    run_pairtools("scaling", "-o", reference, pairs)
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


def test_scaling_writes_to_stdout(pairs):
    """An empty --output means stdout, as it does upstream."""
    assert run_cli("scaling", pairs).startswith("chrom1\t")

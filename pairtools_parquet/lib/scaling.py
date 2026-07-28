"""Contact-frequency-vs-distance curves.

Every number here is computed by pairtools: ``bins_pairs_by_distance`` bins one
chunk, ``contact_areas_same_reg`` turns bin edges into areas, and ``geomspace``
lays out the bins. What this module adds is the loop that feeds Arrow batches
into them.

That loop exists upstream as ``compute_scaling``, and it already iterates over
chunks -- but its entry point accepts only a DataFrame or a path, and raises
``ValueError`` for anything else, so an iterator of DataFrames cannot be handed
to it. UPSTREAM.md records the one-line widening that would let this call into
pairtools instead of restating the loop.

Only six columns take part in the binning, so Parquet input reads six columns
instead of the whole file -- the projection is where the format pays off here.
"""

import numpy as np
import pandas as pd
from pairtools.lib import fileio
from pairtools.lib.scaling import (
    bins_pairs_by_distance,
    contact_areas_same_reg,
    geomspace,
)

from .arrowio import open_pairs
from .chunking import rechunk

UTIL_NAME = "pairtools_parquet_scaling"

#: The only columns the binning looks at.
SCALING_COLUMNS = ["chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2"]

#: Pairs per chunk, matching `pairtools scaling --chunksize`.
DEFAULT_CHUNKSIZE = 100_000


def dist_bin_edges(dist_range, n_dist_bins_decade):
    """Log-spaced distance bin edges, laid out as ``compute_scaling`` lays them."""
    return geomspace(
        dist_range[0],
        dist_range[1],
        int(np.round(np.log10(dist_range[1] / dist_range[0]) * n_dist_bins_decade)),
    )


def accumulate_scaling(
    chunks,
    dist_bins,
    regions=None,
    chromsizes=None,
    ignore_trans=False,
    keep_unassigned=False,
):
    """Bin every chunk and sum the results, returning (cis scalings, trans levels).

    This is the body of ``pairtools.lib.scaling.compute_scaling`` past its input
    dispatch, restated over an iterator of DataFrames.
    """
    sc, trans_counts = None, None

    for chunk in chunks:
        sc_chunk, trans_chunk = bins_pairs_by_distance(
            chunk,
            dist_bins,
            regions=regions,
            chromsizes=chromsizes,
            ignore_trans=ignore_trans,
            keep_unassigned=keep_unassigned,
        )

        sc = sc_chunk if sc is None else sc.add(sc_chunk, fill_value=0)
        trans_counts = (
            trans_chunk
            if trans_counts is None
            else trans_counts.add(trans_chunk, fill_value=0)
        )

    if sc is None:
        raise ValueError("no pairs to compute scalings from")

    sc.reset_index(inplace=True)
    sc["n_bp2"] = contact_areas_same_reg(
        sc["min_dist"], sc["max_dist"], sc["end1"] - sc["start1"]
    )

    if not ignore_trans:
        trans_counts.reset_index(inplace=True)
        trans_counts["n_bp2"] = (trans_counts["end1"] - trans_counts["start1"]) * (
            trans_counts["end2"] - trans_counts["start2"]
        )

    return sc, trans_counts


def scaling_chunks(input_path, chunksize=DEFAULT_CHUNKSIZE, nproc_in=3, cmd_in=None):
    """Yield the pairs of `input_path` as DataFrames of exactly `chunksize` rows.

    The chunk boundaries have to match pairtools' because per-chunk results are
    combined with ``DataFrame.add``, which promotes the integer pair counts to
    floats -- so a file read in one chunk and the same file read in two do not
    format identically, even though every count agrees.
    """
    _, reader = open_pairs(
        input_path,
        columns=SCALING_COLUMNS,
        batch_size=chunksize,
        nproc_in=nproc_in,
        cmd_in=cmd_in,
    )
    schema = reader.schema
    return _at_least_one_chunk(
        rechunk((batch.to_pandas() for batch in reader), chunksize), schema
    )


def _at_least_one_chunk(chunks, schema):
    """Yield `chunks`, or a single empty frame if there were none.

    A file with a header and no data rows still has a scaling: `pd.read_csv`
    hands pairtools one empty chunk for it, and pairtools answers with a table
    of zeros. An Arrow reader yields no batches at all for that file, so the
    empty frame has to be put back or the loop would have nothing to bin.
    """
    empty = True
    for chunk in chunks:
        empty = False
        yield chunk
    if empty:
        yield schema.empty_table().to_pandas()


def scaling_pairs(
    input_path,
    output,
    view=None,
    chunksize=DEFAULT_CHUNKSIZE,
    dist_range=(int(1e0), int(1e9)),
    n_dist_bins_decade=8,
    **kwargs
):
    """Write the scaling curves of `input_path` to `output` as a .tsv table.

    Parameters
    ----------
    input_path : str
        A .pairs/.pairsam/.parquet file.
    output : str
        Where to write the table; empty writes to stdout.
    view : str, optional
        Path to a table of regions to restrict the calculation to, read with
        ``pd.read_table`` and so needing a named header row. Without it each
        chromosome present in the data is its own region -- pairtools describes
        this as taking the regions from the header, but it never reads the
        header's chromsizes.
    """
    regions = pd.read_table(view) if view is not None else None

    cis_scalings, trans_levels = accumulate_scaling(
        scaling_chunks(
            input_path,
            chunksize=chunksize,
            nproc_in=kwargs.get("nproc_in", 3),
            cmd_in=kwargs.get("cmd_in", None),
        ),
        dist_bin_edges(dist_range, n_dist_bins_decade),
        regions=regions,
    )
    summary_stats = pd.concat([cis_scalings, trans_levels])

    outstream = fileio.auto_open(
        output,
        mode="w",
        nproc=kwargs.get("nproc_out", 8),
        command=kwargs.get("cmd_out", None),
    )
    try:
        summary_stats.to_csv(outstream, sep="\t", index=False)
    finally:
        if output:
            outstream.close()

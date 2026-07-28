"""Randomly subsampling pairs.

`pairtools sample` seeds Python's `random` and draws once per row, in order, so
a given seed selects a specific set of pairs. Reproducing that selection means
reproducing that sequence of draws exactly -- drawing a vector with numpy would
be faster but would select different pairs for the same seed, which would make
`--seed` mean something different here than in pairtools.
"""

import random

import numpy as np
import pyarrow as pa
from pairtools.lib import headerops

from .arrowio import PairsWriter, open_pairs

UTIL_NAME = "pairtools_parquet_sample"


def selection_mask(n_rows, fraction, rng):
    """Draw `n_rows` decisions, one per row, in order."""
    return np.fromiter(
        (rng.random() <= fraction for _ in range(n_rows)), dtype=bool, count=n_rows
    )


def sample_pairs(
    input_path,
    output,
    fraction,
    seed=None,
    compress_program="auto",
    row_group_size=None,
    util_name=UTIL_NAME,
    **kwargs
):
    """Write a random `fraction` of the pairs of `input_path` to `output`."""
    header, reader = open_pairs(
        input_path,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )
    new_header = headerops.append_new_pg(header, ID=util_name, PN=util_name)

    # A private Random rather than the module-level one, so that sampling does
    # not depend on, or disturb, anything else that uses `random`.
    rng = random.Random(seed)

    with PairsWriter(
        output,
        new_header,
        schema=reader.schema,
        compress_program=compress_program,
        row_group_size=row_group_size,
        nproc_out=kwargs.get("nproc_out", 8),
    ) as writer:
        for batch in reader:
            mask = selection_mask(batch.num_rows, fraction, rng)
            writer.write(batch.filter(pa.array(mask)))

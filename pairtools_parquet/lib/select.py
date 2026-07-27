"""Filtering pairs by a `pairtools select` CONDITION.

The condition language is pairtools', evaluated by pairtools: this calls
``pairtools.lib.select.evaluate_df`` per batch rather than reimplementing the
expression semantics. That matters because the language is arbitrary Python
with a handful of helpers (``csv_match``, ``wildcard_match``, ``regex_match``,
``region_match``) in scope, and any translation of it into SQL is an
approximation -- the one this replaces silently mishandled ``not``, chained
comparisons, ``--startup-code`` and anything else it had no regex for.

Evaluation is row-wise, as it is in pairtools. What this saves over
``pairtools select`` is the text parsing and re-serialization on either side,
not the evaluation itself; vectorizing the common comparison-only conditions is
left to the point where the other tools get vectorized, so that it can be
checked against this implementation.
"""

import warnings

import numpy as np
import pyarrow as pa
import pairtools.lib.select as pairtools_select
from pairtools.lib import headerops, pairsam_format
from pairtools.lib.select import evaluate_df

from .arrowio import PairsWriter, open_pairs

UTIL_NAME = "pairtools_parquet_select"


def _backport_region_match():
    """Make `region_match` available on pairtools releases that lack it.

    CONDITION is evaluated inside pairtools' own module namespace, so the
    helper functions have to live there. `region_match` was added upstream in
    open2c/pairtools#278 and is not in 1.1.2, the oldest release we support,
    while this package has offered it since 0.2.0 -- so on an older pairtools
    we install upstream's implementation rather than drop the feature.

    Delete this once the declared pairtools minimum contains it.
    """
    if hasattr(pairtools_select, "region_match"):
        return

    def region_match(chrom, pos, region_chrom, region_start=-1, region_end=-1):
        if region_end == -1:
            region_end = np.inf
        return chrom == region_chrom and region_start <= pos <= region_end

    pairtools_select.region_match = region_match


_backport_region_match()


def read_chrom_subset(path):
    """Read chromosome names from the first column of a chromsizes file."""
    with open(path, "r") as f:
        return [line.strip().split("\t")[0] for line in f if line.strip()]


def header_update(header, util_name=UTIL_NAME, remove_columns="", chrom_subset=None):
    """Apply a tool's provenance, column removal and chromosome subset to a header."""
    new_header = headerops.append_new_pg(header, ID=util_name, PN=util_name)

    if remove_columns:
        input_columns = headerops.extract_column_names(header)
        removed = remove_columns.split(",")
        for col in removed:
            if col in pairsam_format.COLUMNS_PAIRS:
                warnings.warn(
                    "Removing required {} column for .pairs format. "
                    "Output is not .pairs anymore".format(col)
                )
            elif col in pairsam_format.COLUMNS_PAIRSAM:
                warnings.warn(
                    "Removing required {} column for .pairsam format. "
                    "Output is not .pairsam anymore".format(col)
                )
        updated_columns = [x for x in input_columns if x not in removed]

        if len(updated_columns) == len(input_columns):
            warnings.warn(
                "Some column(s) {} not in the file, the operation has no "
                "effect".format(",".join(removed))
            )
        elif not updated_columns:
            raise ValueError("--remove-columns removed every column")
        else:
            new_header = headerops.set_columns(new_header, updated_columns)

    if chrom_subset is not None:
        new_header = headerops.subset_chroms_in_pairsheader(
            new_header, read_chrom_subset(chrom_subset)
        )

    return new_header


def evaluate_batch(batch, condition, type_cast=(), startup_code=None):
    """Return a boolean mask of the rows of `batch` satisfying `condition`."""
    df = batch.to_pandas()
    # engine="python" is the one that matches `pairtools select`: its pandas
    # engine goes through DataFrame.eval, which cannot call the helper
    # functions and rejects `and`/`or`, so almost every real condition would
    # fail or, worse, mean something different.
    mask = evaluate_df(df, condition, type_cast, startup_code, engine="python")
    return np.asarray(mask, dtype=bool)


def select_pairs(
    input_path,
    output,
    condition,
    output_rest=None,
    chrom_subset=None,
    startup_code=None,
    type_cast=(),
    remove_columns="",
    compress_program="auto",
    row_group_size=None,
    util_name=UTIL_NAME,
    **kwargs
):
    """Write the pairs of `input_path` satisfying `condition` to `output`.

    Parameters
    ----------
    input_path, output : str
        Paths; the format of each is taken from its extension.
    condition : str
        A pairtools CONDITION expression.
    output_rest : str, optional
        Where to write the pairs that did not satisfy the condition.
    chrom_subset : str, optional
        Path to a chromsizes file; additionally requires both sides of a pair
        to be on those chromosomes, and subsets the header accordingly.
    startup_code : str, optional
        Executed before filtering, to define helpers used in `condition`.
    type_cast : tuple of (str, str)
        Extra column-to-type casts for the evaluation.
    remove_columns : str
        Comma-separated columns to drop from the output.
    """
    header, reader = open_pairs(
        input_path,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )

    new_header = header_update(header, util_name, remove_columns, chrom_subset)
    keep_columns = headerops.extract_column_names(new_header)
    chroms = set(read_chrom_subset(chrom_subset)) if chrom_subset else None

    writer_kwargs = dict(
        compress_program=compress_program,
        row_group_size=row_group_size,
        nproc_out=kwargs.get("nproc_out", 8),
    )
    selected_schema = pa.schema([reader.schema.field(c) for c in keep_columns])

    with PairsWriter(
        output, new_header, schema=selected_schema, **writer_kwargs
    ) as selected_writer:
        rest_writer = (
            PairsWriter(output_rest, new_header, schema=selected_schema,
                        **writer_kwargs)
            if output_rest
            else None
        )
        try:
            for batch in reader:
                mask = evaluate_batch(batch, condition, type_cast, startup_code)
                if chroms is not None:
                    mask &= np.isin(batch.column("chrom1").to_numpy(zero_copy_only=False), list(chroms))
                    mask &= np.isin(batch.column("chrom2").to_numpy(zero_copy_only=False), list(chroms))

                projected = batch.select(keep_columns)
                selected_writer.write(projected.filter(pa.array(mask)))
                if rest_writer is not None:
                    rest_writer.write(projected.filter(pa.array(~mask)))
        finally:
            if rest_writer is not None:
                rest_writer.close()

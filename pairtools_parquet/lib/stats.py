"""Summary statistics over pairs.

`pairtools.lib.stats.PairCounter` already accumulates from DataFrames, so this
is a thin adapter: Arrow batches in, the same counter, the same output file.
Nothing about the statistics themselves is reimplemented.

`--merge` operates on stats files rather than on pairs, so it is handed
straight to pairtools' `do_merge`.
"""

from pairtools.lib import headerops
from pairtools.lib.stats import PairCounter, do_merge

from .._logging import get_logger
from .arrowio import open_pairs

UTIL_NAME = "pairtools_parquet_stats"

logger = get_logger()

#: Rows per batch. pairtools reads stats chunks of this size.
DEFAULT_CHUNKSIZE = 100_000


def stats_pairs(
    input_path,
    output,
    n_dist_bins_decade=8,
    with_chromsizes=True,
    yaml=False,
    bytile_dups=False,
    output_bytile_stats=None,
    filters=None,
    startup_code="",
    type_cast=(),
    engine="pandas",
    chunksize=DEFAULT_CHUNKSIZE,
    **kwargs
):
    """Accumulate statistics over `input_path` and write them to `output`.

    Parameters
    ----------
    input_path : str
        A .pairs/.pairsam/.parquet file.
    output : str
        Where to write the statistics.
    filters : list of str, optional
        ``name:condition`` filters, as `pairtools stats --filter` takes.
    """
    header, reader = open_pairs(
        input_path,
        batch_size=chunksize,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )
    columns = headerops.extract_column_names(header)

    if bytile_dups and "parent_readID" not in columns:
        logger.warning(
            "No 'parent_readID' column in the file, not generating duplicate stats."
        )
        bytile_dups = False
    if output_bytile_stats:
        bytile_dups = True

    # Non-YAML output can only show one filter, so the first one names it.
    first_filter_name = "no_filter"
    filter_map = None
    if filters:
        first_filter_name = filters[0].split(":", 1)[0]
        if len(filters) > 1 and not yaml:
            logger.warning(
                "Output the first filter only in non-YAML output: %s",
                first_filter_name,
            )
        filter_map = dict(f.split(":", 1) for f in filters)

    counter = PairCounter(
        n_dist_bins_decade=n_dist_bins_decade,
        bytile_dups=bytile_dups,
        filters=filter_map,
        startup_code=startup_code,
        type_cast=type_cast,
        engine=engine,
    )

    for batch in reader:
        counter.add_pairs_from_dataframe(batch.to_pandas())

    if with_chromsizes:
        counter.add_chromsizes(headerops.extract_chromsizes(header))

    with open(output, "w") as outstream:
        # With --bytile-dups and no separate path, pairtools writes the by-tile
        # table into the main output, ahead of the statistics themselves.
        if bytile_dups:
            if output_bytile_stats:
                with open(output_bytile_stats, "w") as f:
                    counter.save_bytile_dups(f)
            else:
                counter.save_bytile_dups(outstream)

        counter.save(
            outstream, yaml=yaml, filter=None if yaml else first_filter_name
        )


def merge_stats(output, stats_paths, **kwargs):
    """Merge several stats files into one, via pairtools' own implementation."""
    do_merge(output, list(stats_paths), **kwargs)

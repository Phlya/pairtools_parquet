"""Summary statistics over pairs.

`pairtools.lib.stats.PairCounter` already accumulates from DataFrames, so this
is a thin adapter: Arrow batches in, the same counter, the same output file.
Nothing about the statistics themselves is reimplemented.

`--merge` operates on stats files rather than on pairs, so it is handed
straight to pairtools' `do_merge`.
"""

import pyarrow as pa
import pyarrow.compute as pc
from pairtools.lib import headerops
from pairtools.lib.stats import PairCounter, do_merge

from .._logging import get_logger
from .arrowio import open_pairs
from .chunking import rechunk
from .select import read_chrom_subset

UTIL_NAME = "pairtools_parquet_stats"

logger = get_logger()

#: Rows per batch. pairtools reads stats chunks of this size.
DEFAULT_CHUNKSIZE = 100_000


def parse_filters(filters, yaml):
    """Turn ``name:condition`` strings into what `PairCounter` wants.

    Shared with `dedup`, which takes the same `--filter` option and builds the
    same counter from it.

    Returns ``(filter map or None, the name to report in non-YAML output)``.
    """
    if not filters:
        return None, "no_filter"

    # Non-YAML output can only show one filter, so the first one names it.
    first_filter_name = filters[0].split(":", 1)[0]
    if len(filters) > 1 and not yaml:
        logger.warning(
            "Output the first filter only in non-YAML output: %s", first_filter_name
        )
    return dict(f.split(":", 1) for f in filters), first_filter_name


def restrict_to_chroms(df, chroms, c1="chrom1", c2="chrom2"):
    """The rows of `df` with both sides on a chromosome in `chroms`.

    For counting only -- `dedup --chrom-subset` narrows what the statistics
    describe, not what is written out.
    """
    if chroms is None:
        return df
    return df[df[c1].isin(chroms) & df[c2].isin(chroms)]


def chrom_subset_mask(table, chroms, c1="chrom1", c2="chrom2"):
    """Rows with both sides on a chromosome in `chroms`."""
    value_set = pa.array(sorted(chroms), type=pa.string())
    return pc.and_(
        pc.is_in(table.column(c1).cast(pa.string()), value_set=value_set),
        pc.is_in(table.column(c2).cast(pa.string()), value_set=value_set),
    )


def stats_pairs(
    input_path,
    output,
    n_dist_bins_decade=8,
    with_chromsizes=True,
    yaml=False,
    bytile_dups=False,
    output_bytile_stats=None,
    filters=None,
    chrom_subset=None,
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
    chrom_subset : str, optional
        A chromsizes file naming the chromosomes of interest. Only pairs with
        both sides on one of them are counted, and only those chromosomes are
        reported. `pairtools stats` declares this option but never reads it,
        so it silently counts everything; see UPSTREAM.md.
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

    chroms = None
    if chrom_subset:
        chroms = read_chrom_subset(chrom_subset)
        header = headerops.subset_chroms_in_pairsheader(header, chroms)

    filter_map, first_filter_name = parse_filters(filters, yaml)

    counter = PairCounter(
        n_dist_bins_decade=n_dist_bins_decade,
        bytile_dups=bytile_dups,
        filters=filter_map,
        startup_code=startup_code,
        type_cast=type_cast,
        engine=engine,
    )

    def frames():
        for batch in reader:
            if chroms is not None:
                batch = batch.filter(chrom_subset_mask(batch, chroms))
            yield batch.to_pandas()

    # `add_pairs_from_dataframe` groups each chunk by chromosome pair, and
    # pandas sorts group keys, so the order chromosome pairs appear in
    # `chrom_freq` follows the chunk boundaries. Text batches are cut by bytes
    # rather than rows, which would put those lines in a different order than
    # `pairtools stats --chunksize` puts them for the same file.
    for frame in rechunk(frames(), chunksize):
        counter.add_pairs_from_dataframe(frame)

    # `extract_chromsizes` zips the parsed `#chromsize:` lines and indexes the
    # result without checking there were any, so it raises IndexError rather
    # than returning nothing on a header that declares none -- which is what a
    # --chrom-subset naming nothing in the file leaves behind.
    if with_chromsizes and headerops.extract_fields(header, "chromsize"):
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

"""Duplicate detection.

Two backends live here, and they compute the same thing.

``scipy``/``sklearn`` are pairtools': every chunk goes through
``pairtools.lib.dedup._dedup_chunk`` unchanged. Only the loop around it is
reimplemented -- the ~30 lines that chunk the input, carry the tail of each
chunk into the next, and split the marked chunk into outputs. That loop exists
upstream as ``_dedup_stream``, but it hard-codes ``pd.read_table`` for input
and ``.to_csv()`` for output, so it cannot be handed Arrow batches.

``duckdb`` is the default and is 6-7x faster. It rests on one observation:
``--max-mismatch`` is 3bp by default and dedup input is sorted, so finding
duplicates is not a nearest-neighbour search over a plane -- it is a lookup in
a 3bp window. Bucketing on ``(chrom1, chrom2, strand1, strand2, pos1 // r)``
puts every pair within `r` of each other in the same bucket or an adjacent one,
which turns the search into an equi-join whose buckets hold exactly the
duplicate families. Profiling the pandas path showed the KD-trees were only 12%
of its runtime and pandas bookkeeping the rest, so this replaces both.

Clustering stays in Python, because it is not the expensive part: the edge list
is small and ``scipy.sparse.csgraph.connected_components`` -- the same function
pairtools uses -- labels 2M edges in under 0.1s.

The pandas backend is kept as the reference the DuckDB one is tested against.
"""

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
from pairtools.lib import headerops
from pairtools.lib.dedup import _dedup_chunk
from pairtools.lib.stats import PairCounter
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from . import arrowio, duckdb_utils
from .arrowio import PairsWriter, open_pairs
from .chunking import rechunk
from .duckdbio import quote_identifier

UTIL_NAME = "pairtools_parquet_dedup"

#: Backends that work on DataFrames, and so can be driven from Arrow batches.
SUPPORTED_BACKENDS = ("duckdb", "scipy", "sklearn")

#: Rows per window, per backend. The pandas backend's window is a semantic
#: knob -- which rows share a chunk decides what gets compared -- so it stays at
#: pairtools' default. The DuckDB backend's is only a memory knob, since its
#: lookback joins clusters back together across the boundary, so it can be
#: large enough that most files go through in one pass.
DEFAULT_CHUNKSIZE = {"duckdb": 20_000_000, "scipy": 10_000, "sklearn": 10_000}


def mark_duplicate_chunks(
    frames,
    colnames,
    carryover,
    max_mismatch,
    method,
    mark_dups,
    keep_parent_id,
    extra_col_pairs,
    backend,
    n_proc,
    unmapped_chrom,
    c1="chrom1",
    c2="chrom2",
    p1="pos1",
    p2="pos2",
    s1="strand1",
    s2="strand2",
):
    """Yield each chunk with a `duplicate` column added.

    Mirrors ``pairtools.lib.dedup._dedup_stream``, differing only in taking
    DataFrames rather than reading them from a text stream.
    """
    df_prev_nodups = pd.DataFrame([])

    for df in frames:
        df = df.copy()
        df["carryover"] = False
        input_chunk = pd.concat(
            [df_prev_nodups, df], axis=0, ignore_index=True
        ).reset_index(drop=True)

        df_marked = _dedup_chunk(
            input_chunk,
            r=max_mismatch,
            method=method,
            keep_parent_id=keep_parent_id,
            extra_col_pairs=list(extra_col_pairs),
            backend=backend,
            n_proc=n_proc,
            c1=c1,
            c2=c2,
            p1=p1,
            p2=p2,
            s1=s1,
            s2=s2,
            unmapped_chrom=unmapped_chrom,
        )

        df_marked = (
            df_marked[~df_marked["carryover"]]
            .drop(columns=["carryover"])
            .reset_index(drop=True)
        )

        mask_duplicated = df_marked["duplicate"]
        if mark_dups:
            df_marked.loc[mask_duplicated, "pair_type"] = "DD"

        yield df_marked

        # The tail of the non-duplicates is prepended to the next chunk, so a
        # duplicate pair straddling the boundary is still seen as a pair.
        df_nodups = df_marked.loc[~mask_duplicated, colnames]
        df_prev_nodups = df_nodups.tail(carryover).reset_index(drop=True)
        df_prev_nodups["carryover"] = True


def _write(writer, df, columns):
    if writer is None or not len(df):
        return
    writer.write(pa.Table.from_pandas(df[columns], preserve_index=False))


# --------------------------------------------------------------------------
# DuckDB backend
# --------------------------------------------------------------------------

#: Distance between two pairs, per --method.
METRICS = {
    "max": "greatest(abs(a.{p1} - p.{p1}), abs(a.{p2} - p.{p2}))",
    "sum": "abs(a.{p1} - p.{p1}) + abs(a.{p2} - p.{p2})",
}


def key_columns(extra_col_pairs, c1, c2, s1, s2):
    """The columns that must match exactly for two pairs to be duplicates.

    `extra_col_pairs` whose two names differ compare one column of a read
    against a *different* column of the other, which is not an equality and so
    cannot be a join key; those fall back to the pandas backend.
    """
    columns = [c1, c2, s1, s2]
    for left, right in extra_col_pairs:
        if left != right:
            raise ValueError(
                "the duckdb backend cannot match {!r} against {!r}: it groups "
                "candidates by equality, so an extra column pair has to name "
                "the same column on both sides. Use --backend scipy for "
                "this.".format(left, right)
            )
        columns.append(left)
    return columns


def edge_sql(keys, p1, p2, r, method, lo, hi):
    """SQL returning every pair of rows within `r` of each other, as (u, v).

    Rows are bucketed by position so the search is an equi-join rather than a
    scan of the whole (chrom, chrom, strand, strand) group, which for a cis
    chromosome would be millions of rows. A bucket is `r` bases wide, so two
    rows within `r` are either in the same bucket or in adjacent ones; the
    probe side is emitted three times, once per neighbouring bucket.

    Emitting all three -- rather than relying on the input being sorted and
    only looking forwards -- means every within-`r` pair is generated in both
    directions, and ``a.rid < p.rid`` keeps the one whose first row comes first
    in the file. That is pairtools' keeper, and it holds whatever order the
    input is in.

    Rows sharing an exact position are collapsed first, because they have
    identical neighbourhoods: whatever one of them is within `r` of, so are the
    others. The smallest row id of each position stands for the rest, joined to
    them by a star, and only those representatives take part in the pairwise
    search.

    That matters where families are deep. A family of ten otherwise contributes
    45 edges where nine span it, and the count grows quadratically with
    duplication: on a 90%-duplicate library the edge list drops from 47M to 16M
    and the clustering stage from 13.5s to 4.9s. It costs about 1.4s on 5.6M
    rows where there is nothing to collapse, which is the trade -- a fixed cost
    against a quadratic one, in both time and the memory the edge arrays take.

    The representative is found with a window rather than a GROUP BY joined
    back to the rows: the join form makes DuckDB build a plan that did not
    finish at all on a low-duplication file, while the window is one pass.
    """
    if method not in METRICS:
        raise ValueError('Unknown method, only "sum" or "max" allowed')

    columns = ", ".join(quote_identifier(name) for name in keys + [p1, p2])
    point = ", ".join(quote_identifier(name) for name in keys + [p1, p2])
    on = " AND ".join(
        "a.{k} = p.{k}".format(k=quote_identifier(name)) for name in keys
    )
    metric = METRICS[method].format(
        p1=quote_identifier(p1), p2=quote_identifier(p2)
    )
    # A zero radius means exact equality, but a bucket still needs a width.
    width = max(int(r), 1)

    return """
WITH w AS (
  SELECT rid, {columns}, min(rid) OVER (PARTITION BY {point}) AS rep
  FROM keys WHERE rid >= {lo} AND rid < {hi}
), stars AS (
  SELECT rep AS u, rid AS v FROM w WHERE rid <> rep
), b AS (
  SELECT rid, {columns}, {p1} // {width} AS blk FROM w WHERE rid = rep
), p AS (
              SELECT rid, {columns}, blk     AS kb FROM b
    UNION ALL SELECT rid, {columns}, blk - 1       FROM b
    UNION ALL SELECT rid, {columns}, blk + 1       FROM b
), near AS (
  SELECT a.rid AS u, p.rid AS v
  FROM b a JOIN p ON {on} AND a.blk = p.kb
  WHERE a.rid < p.rid AND {metric} <= {r}
)
SELECT u, v FROM stars UNION ALL SELECT u, v FROM near
""".format(
        columns=columns,
        point=point,
        p1=quote_identifier(p1),
        width=width,
        lo=int(lo),
        hi=int(hi),
        on=on,
        metric=metric,
        r=int(r),
    )


def exact_duplicates(con, keys, p1, p2):
    """Find duplicates when `--max-mismatch` is 0, which needs no graph at all.

    Exact equality is transitive, so each group of identical pairs is already a
    cluster: there is nothing to connect, no edge list to build or move, and no
    window to bound it -- one pass over the key columns answers the whole file
    at once. The first row of each group is the keeper, as everywhere else.

    Only the duplicates come back, so at low duplication almost nothing is
    transferred.

    Returns ``(duplicate row ids, their keepers' row ids)``.
    """
    point = ", ".join(quote_identifier(name) for name in keys + [p1, p2])
    found = con.execute(
        """
        SELECT rid, rep FROM (
          SELECT rid, min(rid) OVER (PARTITION BY {point}) AS rep FROM keys
        ) WHERE rid <> rep
        """.format(point=point)
    ).fetchnumpy()
    return found["rid"].astype(np.int64), found["rep"].astype(np.int64)


def cluster_window(u, v, lo, hi):
    """Label connected components over rows [lo, hi).

    The graph is built over the raw row-id range rather than over a compacted
    node list, which skips a ``np.unique`` -- 2.3s on 4M ids under numpy 2.4,
    against 0.07s for the component labelling itself.

    Returns ``(labels, keepers)``: the component of each row in the range, and
    the smallest row id in each component, which is the read pairtools keeps.
    """
    size = int(hi - lo)
    graph = coo_matrix(
        (np.ones(len(u), dtype=np.int8), (u - lo, v - lo)), shape=(size, size)
    )
    n_components, labels = connected_components(graph, directed=False)

    local = np.arange(size, dtype=np.int64)
    # Assigning in reverse leaves the smallest member of each component as the
    # last write, and NumPy keeps the last value for a repeated index.
    keepers = np.empty(n_components, dtype=np.int64)
    keepers[labels[::-1]] = local[::-1]
    return labels, keepers + lo


def duplicate_windows(con, keys, p1, p2, r, method, n_rows, window, carryover):
    """Yield ``(lo, end, keeper of each row in [lo, end))``, window by window.

    Each window is deduplicated together with a lookback of `carryover` rows
    from the one before, so a duplicate family straddling the boundary is still
    one family. The lookback rows are *re-decided*, not just consulted: a
    window is the first to see a family joined by one of its rows, and the
    answer for a row is only final once `carryover` rows have gone past it.
    Since a component can only grow as later rows arrive, a keeper can only
    move earlier, so revisions never take a duplicate back.

    This is where pairtools loses duplicates: its carryover holds only the
    previous chunk's non-duplicates, so a chain reaching the boundary through a
    duplicate is cut, and its decisions are never revisited.
    """
    window = max(int(window), 1)
    keeper, keeper_lo = None, 0

    for start in range(0, max(int(n_rows), 1), window):
        end = min(start + window, n_rows)
        # The lookback is exactly the range the previous window kept keepers
        # for, so every row in it can be resolved.
        lo = keeper_lo if keeper is not None else start

        edges = con.execute(edge_sql(keys, p1, p2, r, method, lo, end)).fetchnumpy()
        labels, keepers = cluster_window(
            edges["u"].astype(np.int64), edges["v"].astype(np.int64), lo, end
        )

        if keeper is not None:
            # A component reaching into the lookback inherits the keeper those
            # rows already resolved to, and the smallest of them wins. Taking
            # the keeper of the component's own smallest row instead would be
            # wrong where this window is the first to join two families that
            # were separate before: the join can run through a lookback row
            # whose family reaches further back than the smallest one's does.
            # There are at most `carryover` rows to fold in.
            np.minimum.at(keepers, labels[: len(keeper)], keeper)

        parents = keepers[labels]
        yield lo, end, parents

        take = min(int(carryover), int(end - lo))
        if take > 0:
            # Copied, so the whole window's array is not held alive by a slice.
            keeper, keeper_lo = parents[-take:].copy(), end - take
        else:
            keeper, keeper_lo = None, 0


def _batches_with_row_id(reader, columns, counter):
    """Re-emit `reader`'s batches carrying a file-order row id.

    DuckDB's own ``row_number() OVER ()`` is not usable for this: its CSV
    scanner is parallel, so the order rows reach the window operator is not the
    order they are in the file. The Arrow reader is sequential, so counting
    here is exact for every input format.
    """
    for batch in reader:
        rid = pa.array(
            np.arange(counter[0], counter[0] + batch.num_rows, dtype=np.int64)
        )
        counter[0] += batch.num_rows
        yield pa.RecordBatch.from_arrays(
            [rid] + [batch.column(name) for name in columns],
            names=["rid"] + list(columns),
        )


def load_key_table(con, input_path, columns, c1, c2, unmapped_chrom, **kwargs):
    """Build the DuckDB table the neighbour search runs over, and count rows.

    Only the columns that decide whether two pairs are duplicates are loaded.
    For Parquet that is a projected scan of a handful of columns; the rest of
    the file is never touched until the marked rows are written out.
    """
    mapped = "{} <> '{}' AND {} <> '{}'".format(
        quote_identifier(c1), unmapped_chrom, quote_identifier(c2), unmapped_chrom
    )
    projection = ", ".join(quote_identifier(name) for name in columns)

    if arrowio.is_parquet(input_path) and not kwargs.get("cmd_in"):
        import pyarrow.parquet as pq

        n_rows = pq.ParquetFile(str(input_path)).metadata.num_rows
        con.execute(
            """
            CREATE OR REPLACE TEMP TABLE keys AS
            SELECT file_row_number AS rid, {projection}
            FROM read_parquet('{path}', file_row_number=true)
            WHERE {mapped}
            """.format(
                projection=projection,
                path=str(input_path).replace("'", "''"),
                mapped=mapped,
            )
        )
        return n_rows

    counter = [0]
    _, reader = open_pairs(
        input_path,
        columns=list(columns),
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )
    indexed = pa.RecordBatchReader.from_batches(
        pa.schema([pa.field("rid", pa.int64())] + list(reader.schema)),
        _batches_with_row_id(reader, list(columns), counter),
    )
    con.register("pairs_keys", indexed)
    con.execute(
        "CREATE OR REPLACE TEMP TABLE keys AS "
        "SELECT * FROM pairs_keys WHERE {}".format(mapped)
    )
    con.unregister("pairs_keys")
    return counter[0]


def parent_read_ids(con, dup_rid, parent_rid, input_path, **kwargs):
    """The readID of each duplicate's keeper, ordered by the duplicate's row id.

    Resolved with one join rather than a row-by-row lookup during the output
    pass, so ``--keep-parent-id`` stays vectorized.
    """
    con.register(
        "dup_ids",
        pa.table({"rid": pa.array(dup_rid), "parent_rid": pa.array(parent_rid)}),
    )
    try:
        if arrowio.is_parquet(input_path) and not kwargs.get("cmd_in"):
            source = (
                "(SELECT file_row_number AS rid, readID "
                "FROM read_parquet('{}', file_row_number=true))".format(
                    str(input_path).replace("'", "''")
                )
            )
        else:
            source = "read_ids"
            counter = [0]
            _, reader = open_pairs(
                input_path,
                columns=["readID"],
                nproc_in=kwargs.get("nproc_in", 3),
                cmd_in=kwargs.get("cmd_in", None),
            )
            con.register(
                "read_ids",
                pa.RecordBatchReader.from_batches(
                    pa.schema([pa.field("rid", pa.int64())] + list(reader.schema)),
                    _batches_with_row_id(reader, ["readID"], counter),
                ),
            )

        return (
            con.execute(
                "SELECT s.readID FROM dup_ids d JOIN {} s ON s.rid = d.parent_rid "
                "ORDER BY d.rid".format(source)
            )
            .fetch_arrow_table()
            .column("readID")
            .combine_chunks()
        )
    finally:
        con.unregister("dup_ids")


def _scatter(values, positions, length):
    """A length-`length` string array with `values` at `positions`, "" elsewhere."""
    out = np.full(length, "", dtype=object)
    if len(positions):
        out[positions] = values
    return pa.array(out, type=pa.string())


def dedup_pairs(
    input_path,
    output,
    output_dups=None,
    output_unmapped=None,
    output_stats=None,
    max_mismatch=3,
    method="max",
    backend="duckdb",
    chunksize=None,
    carryover=100,
    n_proc=1,
    mark_dups=True,
    keep_parent_id=False,
    extra_col_pair=(),
    unmapped_chrom="!",
    c1="chrom1",
    c2="chrom2",
    p1="pos1",
    p2="pos2",
    s1="strand1",
    s2="strand2",
    yaml=False,
    compress_program="auto",
    row_group_size=None,
    util_name=UTIL_NAME,
    **kwargs
):
    """Remove PCR duplicates from a sorted .pairs/.parquet file.

    Parameters
    ----------
    input_path : str
        Input path; must be sorted, as `pairtools dedup` requires.
    output : str
        Where the non-duplicate pairs go.
    output_dups, output_unmapped : str, optional
        Where duplicates and unmapped pairs go. Passing the same path as
        `output` for `output_dups` keeps duplicates in the main output, as in
        pairtools.
    output_stats : str, optional
        Where to write the duplication statistics.
    """
    if backend not in SUPPORTED_BACKENDS:
        raise NotImplementedError(
            "the {!r} backend is not available in pairtools_parquet: it works "
            "line-by-line on text streams rather than on dataframes. Use "
            "--backend duckdb (the default), scipy or sklearn.".format(backend)
        )
    if chunksize is None:
        chunksize = DEFAULT_CHUNKSIZE[backend]

    header = arrowio.read_header(
        input_path,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )
    columns = headerops.extract_column_names(header)

    out_header = headerops.append_new_pg(header, ID=util_name, PN=util_name)
    dups_header = out_header
    dups_columns = columns
    if keep_parent_id:
        dups_columns = columns + ["parent_readID"]
        # append_columns rewrites the #columns line of the list it is given and
        # returns that same list, so it has to be handed a copy.
        dups_header = headerops.append_columns(list(out_header), ["parent_readID"])

    stats = PairCounter() if output_stats else None
    dups_to_main_output = bool(output_dups) and str(output_dups) == str(output)

    writers = {}
    try:
        # pairtools drops parent_readID only just before writing the
        # non-duplicates, so the dups and unmapped streams still carry it.
        writer_kwargs = dict(
            compress_program=compress_program,
            row_group_size=row_group_size,
            nproc_out=kwargs.get("nproc_out", 8),
        )
        writers["out"] = PairsWriter(
            output,
            dups_header if dups_to_main_output else out_header,
            **writer_kwargs
        )
        if output_dups and not dups_to_main_output:
            writers["dups"] = PairsWriter(output_dups, dups_header, **writer_kwargs)
        if output_unmapped:
            writers["unmapped"] = PairsWriter(
                output_unmapped, dups_header, **writer_kwargs
            )

        marking = dict(
            columns=columns,
            dups_columns=dups_columns,
            dups_to_main_output=dups_to_main_output,
            keep_parent_id=keep_parent_id,
            mark_dups=mark_dups,
            max_mismatch=max_mismatch,
            method=method,
            chunksize=chunksize,
            carryover=carryover,
            extra_col_pair=extra_col_pair,
            unmapped_chrom=unmapped_chrom,
            c1=c1, c2=c2, p1=p1, p2=p2, s1=s1, s2=s2,
        )
        run = _dedup_duckdb if backend == "duckdb" else _dedup_pandas
        run(
            input_path,
            writers,
            stats,
            backend=backend,
            n_proc=n_proc,
            marking=marking,
            **kwargs
        )
    finally:
        for writer in writers.values():
            writer.close()

    if stats is not None:
        with open(output_stats, "w") as f:
            stats.save(f, yaml=yaml)


def _dedup_pandas(input_path, writers, stats, backend, n_proc, marking, **kwargs):
    """pairtools' own chunked KD-tree implementation, over Arrow batches."""
    columns = marking["columns"]
    dups_columns = marking["dups_columns"]
    unmapped_chrom = marking["unmapped_chrom"]
    c1, c2 = marking["c1"], marking["c2"]

    _, reader = open_pairs(
        input_path,
        batch_size=marking["chunksize"],
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )
    frames = rechunk((batch.to_pandas() for batch in reader), marking["chunksize"])
    marked = mark_duplicate_chunks(
        frames,
        colnames=columns,
        carryover=marking["carryover"],
        max_mismatch=marking["max_mismatch"],
        method=marking["method"],
        mark_dups=marking["mark_dups"],
        keep_parent_id=marking["keep_parent_id"],
        extra_col_pairs=marking["extra_col_pair"],
        backend=backend,
        n_proc=n_proc,
        unmapped_chrom=unmapped_chrom,
        c1=c1,
        c2=c2,
        p1=marking["p1"],
        p2=marking["p2"],
        s1=marking["s1"],
        s2=marking["s2"],
    )

    for df in marked:
        if stats is not None:
            stats.add_pairs_from_dataframe(df, unmapped_chrom=unmapped_chrom)

        mask_mapped = np.logical_and(
            (df[c1] != unmapped_chrom), (df[c2] != unmapped_chrom)
        )
        mask_duplicates = df["duplicate"].to_numpy(dtype=bool)
        df = df.drop(columns=["duplicate"])

        _write(writers.get("unmapped"), df.loc[~mask_mapped], dups_columns)

        if marking["dups_to_main_output"]:
            _write(writers["out"], df.loc[mask_mapped], dups_columns)
        else:
            _write(writers.get("dups"), df.loc[mask_duplicates], dups_columns)
            _write(writers["out"], df.loc[mask_mapped & (~mask_duplicates)], columns)


def _dedup_duckdb(input_path, writers, stats, backend, n_proc, marking, **kwargs):
    """Find duplicates with a blocked equi-join, then stream the rows out.

    Two passes over the input. The first loads only the columns that decide
    whether two pairs are duplicates and produces a boolean per row; the second
    reads the file in order and splits it into the outputs. Reading twice is
    cheaper than sorting the join's output back into file order, and for
    Parquet the first pass touches a handful of columns.
    """
    keys = key_columns(
        marking["extra_col_pair"], marking["c1"], marking["c2"],
        marking["s1"], marking["s2"],
    )
    p1, p2 = marking["p1"], marking["p2"]
    keep_parent_id = marking["keep_parent_id"]

    con = duckdb_utils.setup_duckdb_connection(
        temp_directory=kwargs.get("tmpdir") or None,
        memory_limit=kwargs.get("memory") or None,
        enable_progress_bar=False,
        enable_profiling="no_output",
        numb_threads=max(int(n_proc), 1),
    )
    try:
        n_rows = load_key_table(
            con,
            input_path,
            keys + [p1, p2],
            marking["c1"],
            marking["c2"],
            marking["unmapped_chrom"],
            **kwargs
        )

        is_dup = np.zeros(n_rows, dtype=bool)
        parent_of = None

        if int(marking["max_mismatch"]) == 0:
            dup_rid, parent_rid = exact_duplicates(con, keys, p1, p2)
            is_dup[dup_rid] = True
            if keep_parent_id:
                parent_of = np.arange(n_rows, dtype=np.int64)
                parent_of[dup_rid] = parent_rid
        else:
            # A duplicate never becomes unique again as later windows revise
            # their lookback, so the flags can be OR-ed in; only the keeper can
            # move, and only when --keep-parent-id makes it visible.
            if keep_parent_id:
                parent_of = np.zeros(n_rows, dtype=np.int64)
            for lo, end, keepers in duplicate_windows(
                con,
                keys,
                p1,
                p2,
                marking["max_mismatch"],
                marking["method"],
                n_rows,
                marking["chunksize"],
                marking["carryover"],
            ):
                is_dup[lo:end] |= keepers != np.arange(lo, end, dtype=np.int64)
                if parent_of is not None:
                    parent_of[lo:end] = keepers

        parents = None
        if keep_parent_id:
            dup_rid = np.flatnonzero(is_dup)
            parents = parent_read_ids(
                con, dup_rid, parent_of[dup_rid], input_path, **kwargs
            )
    finally:
        con.close()

    _write_marked(input_path, writers, stats, is_dup, parents, marking, **kwargs)


def _write_marked(input_path, writers, stats, is_dup, parents, marking, **kwargs):
    """Second pass: read the input in order and split it into the outputs."""
    columns = marking["columns"]
    dups_columns = marking["dups_columns"]
    unmapped_chrom = marking["unmapped_chrom"]
    c1, c2 = marking["c1"], marking["c2"]

    _, reader = open_pairs(
        input_path,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )

    offset, consumed = 0, 0
    for batch in reader:
        n = batch.num_rows
        duplicate = is_dup[offset : offset + n]
        offset += n

        table = pa.Table.from_batches([batch])
        mask = pa.array(duplicate)

        if parents is not None:
            positions = np.flatnonzero(duplicate)
            values = parents.slice(consumed, len(positions)).to_pylist()
            consumed += len(positions)
            table = table.append_column(
                "parent_readID", _scatter(values, positions, n)
            )

        if marking["mark_dups"] and "pair_type" in table.column_names:
            index = table.column_names.index("pair_type")
            table = table.set_column(
                index,
                "pair_type",
                pc.if_else(mask, "DD", table.column("pair_type")),
            )

        if stats is not None:
            stats.add_pairs_from_dataframe(
                table.to_pandas(), unmapped_chrom=unmapped_chrom
            )

        mapped = pc.and_(
            pc.not_equal(table.column(c1), unmapped_chrom),
            pc.not_equal(table.column(c2), unmapped_chrom),
        )

        if writers.get("unmapped") is not None:
            writers["unmapped"].write(
                table.filter(pc.invert(mapped)).select(dups_columns)
            )

        if marking["dups_to_main_output"]:
            writers["out"].write(table.filter(mapped).select(dups_columns))
        else:
            if writers.get("dups") is not None:
                writers["dups"].write(table.filter(mask).select(dups_columns))
            writers["out"].write(
                table.filter(pc.and_(mapped, pc.invert(mask))).select(columns)
            )

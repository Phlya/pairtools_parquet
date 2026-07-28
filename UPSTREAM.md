# Notes for upstream (open2c/pairtools)

This package aims to merge into pairtools eventually, so it calls into pairtools
rather than copying it. Where that was not possible, or where pairtools looks
wrong, it is recorded here as something to raise or patch upstream.

## Bugs found while porting

### `headerops.mark_header_as_sorted` corrupts `#chromosomes:`

`mark_header_as_sorted` rewrites

    #chromosomes: chr1 chr2 chr3

as

    #chromosomes: : chr1 chr2 chr3

so every file produced by `pairtools sort` carries the doubled colon, and it is
doubled again on each re-sort. `pairtools merge` does not call the function, so
its output is unaffected — which is why merging sorted files silently "repairs"
the line.

We reproduce the behaviour in `sort` (matching upstream byte-for-byte matters
more than the line being right) and avoid it in `merge`, matching upstream
there too. Both are pinned by tests in `tests/test_parity.py`.

### `dedup --keep-parent-id` writes an unmapped file whose header undercounts its columns

`lib/dedup.py:streaming_dedup` writes the unmapped pairs before dropping the
`parent_readID` column:

```python
if outstream_unmapped:
    df_chunk.loc[~mask_mapped, :].to_csv(outstream_unmapped, ...)   # 9 fields
...
if keep_parent_id:
    df_chunk = df_chunk.drop(columns=["parent_readID"])
```

but the header written to that stream is the un-extended one. The result
declares 8 columns and contains 9 fields per row, so anything that parses by
`#columns` misreads it.

**We deliberately diverge here**: the unmapped output gets the 9-column header
that matches its rows. Bug-for-bug parity is not available in Parquet anyway —
a Parquet file has a real schema, so there is no way to write 9 columns while
declaring 8. The bodies are identical to pairtools; only the header differs,
and only under `--keep-parent-id --output-unmapped`.

### `flip` oscillates for pairs on the same unannotated chromosome

`pairtools flip` is meant to project pairs onto the upper triangle, so running
it twice should be the same as running it once. It is not, for a pair whose two
sides are on the same chromosome that is *absent from the chromsizes file*.
That case falls into the both-unannotated branch:

```python
elif not is_annotated1 and not is_annotated2:
    has_correct_order = cols[chrom1_col] < cols[chrom2_col]
```

which is `False` when the two names are equal, so the sides are swapped —
every run, regardless of position:

```
in    r1  chrUNKNOWN 100  chrUNKNOWN 900  + - UR
x1    r1  chrUNKNOWN 900  chrUNKNOWN 100  - + RU
x2    r1  chrUNKNOWN 100  chrUNKNOWN 900  + - UR
x3    r1  chrUNKNOWN 900  chrUNKNOWN 100  - + RU
```

The positions are never compared, so neither ordering is "upper triangular".
Comparing `(chrom, pos)` rather than `chrom` alone in that branch would fix it,
matching what the both-annotated branch does.

**We deliberately diverge here**: `lib/flip.py` compares `(chrom, pos)` in that
branch, so flipping settles after one pass for every input. The two agree on
every other case — the only rows that differ are those whose two sides are on
the *same* unannotated chromosome, where upstream's answer is not stable
anyway. `tests/test_smalltools.py` pins that the divergence is confined to
exactly those rows, that our flip is idempotent, and — via
`test_pairtools_flip_oscillates` — that upstream still has the bug, so the
test fails and tells us to drop the divergence once it is fixed.

### `restrict` crashes on a chromosome with no annotated fragments

`lib/restrict.py:find_rfrag` means to warn and return the unannotated
sentinels for a chromosome missing from the fragment file:

```python
try:
    rsites_chrom = rfrags[chrom]
except ValueError as e:
    warnings.warn(f"Chomosome {chrom} does not have annotated restriction fragments, return empty.")
    return (UNANNOTATED_RFRAG, UNMAPPED_POS, UNMAPPED_POS)
```

but `rfrags` is a dict, so the lookup raises `KeyError`, which that clause does
not catch. The recovery path is unreachable and `pairtools restrict` dies with
a traceback — reproducible with any pairs file mentioning a scaffold the
fragment BED omits. Catching `KeyError` would restore the intended behaviour.

**We deliberately diverge here**: we do what the code intends — warn once per
chromosome and emit `-1, 0, 0`. Matching a crash is not a useful form of
parity. `tests/test_restrict.py` pins both our behaviour and, via
`test_pairtools_restrict_crashes_on_missing_chromosome`, that upstream still
crashes, so we learn when it is fixed.

### `dedup` silently misses duplicates in input that is not sorted

`pairtools dedup` compares pairs only within a `--chunksize` chunk plus a
`--carryover` tail, which is sound for sorted input, where duplicates are
adjacent. It does not check that its input *is* sorted, and on input that is
not, most duplicate pairs never share a chunk. On a 200k-row shuffled file:

```
true unique reads (same file, sorted)   194189
pairtools dedup on the shuffled file    197802     (3613 duplicates missed)
```

It is also 18x slower doing it, because every chunk's KD-tree is built over
scattered coordinates.

This is not something we work around for `--max-mismatch > 0` — a windowed
algorithm has the same constraint. But **exact** deduplication has no such
constraint: equality is transitive, so a group of identical pairs is a cluster
regardless of where its members sit in the file. `--max-mismatch 0` takes a
single grouped pass here and gives the same answer on a shuffled file as on a
sorted one, pinned by
`test_exact_dedup_does_not_depend_on_input_order`. Upstream could do the same
by special-casing `r == 0` to a hash group-by, and would then also be able to
deduplicate before sorting rather than after.

### `dedup`'s carryover cuts duplicate chains at every chunk boundary

`_dedup_stream` carries the tail of each chunk into the next so a duplicate
family straddling the boundary is still found — but it carries only the
**non**-duplicates:

```python
df_nodups = df_marked.loc[~mask_duplicated, colnames]
df_prev_nodups = df_nodups.tail(carryover).reset_index(drop=True)
```

Duplicates are found by connected components, so a family can be a *chain*: A
and B within `r`, B and C within `r`, A and C not. If the chunk boundary falls
between B and C, B is a duplicate and so is dropped from the carryover — C
never sees it, is not within `r` of A, and is kept as a unique read.

Four pairs 3bp apart with `--max-mismatch 3`, which are one family:

```
              chunksize 1  2  3  4  100
pairtools keeps          2  2  2  1    1
```

It is rare with the defaults (0 rows in 4M with `--method max --max-mismatch
3`) but grows with the radius: `--method sum --max-mismatch 5` on the same 4M
pairs leaves 2 reads that should have been dropped, both at row indices that
are exact multiples of the 10000-row chunk size. Carrying over the last
`carryover` rows regardless of duplicate status, and keeping their resolved
parents, would fix it.

**We deliberately diverge here**: the `duckdb` backend's lookback keeps every
row and substitutes each lookback row's already-final keeper, so a chain
resolves to one read at any window size. `tests/test_dedup.py` pins that our
answer does not depend on `--chunksize` and that the chain above collapses to a
single read. The `scipy`/`sklearn` backends still call upstream's
`_dedup_chunk` through upstream's carryover rule, and so still reproduce it.

### `stats --merge` output ordering is not reproducible

`pairtools stats --merge` writes its keys in an order that varies between
processes — running it twice on the same inputs can put `pair_types/UU` before
or after `pair_types/NU`. The values are always the same; only the line order
moves, which is enough to break byte-comparison of pipeline outputs and to make
diffs noisy. It looks like a set iteration somewhere in `PairCounter.__add__`
picking up string hash randomization.

We delegate to `do_merge`, so we inherit the behaviour rather than papering
over it; `tests/test_stats.py` compares merged stats as sets of lines and says
why.

### `header set-columns` does nothing to a file with no `#columns:` line

`headerops.set_columns` only rewrites a line that is already there:

```python
for i in range(len(header)):
    if header[i].startswith("#columns:"):
        header[i] = "#columns:" + SEP_COLS + SEP_COLS.join(columns)
return header
```

so on a headerless file — the input the command exists for — it returns the
empty header unchanged, and `pairtools header set-columns` writes out a file
with no header at all, exit code 0. Appending the line when it is absent would
fix it.

**We deliberately diverge here**: `lib/header.py:set_columns` adds the line,
laid out exactly as `headerops.set_columns` lays it out. Pinned by
`test_set_columns_adds_a_missing_columns_line`.

### `headerops.append_columns` mutates its argument

It rewrites the `#columns` line of the list it is passed and returns that same
list, unlike the other `headerops` functions, which return new headers. Callers
that keep a reference to the original — e.g. to write a second output with the
unmodified columns — silently get the modified one.

### `select.evaluate_df` casts columns only when they are already the target type

`pairtools/lib/select.py`:

```python
if not str(df.dtypes[col]) != TYPES[col]:
    df[col] = df[col].astype(TYPES[col])
```

The double negative means the cast happens when the column *already* has the
target type, and is skipped when it does not — the opposite of the intent. It
is invisible to `pairtools select`, which uses `evaluate_stream`, and harmless
for us because Arrow already gives the columns their proper types.

### `evaluate_df(engine="pandas")` cannot evaluate the documented condition language

The pandas engine goes through `DataFrame.eval`, which cannot call
`csv_match`/`wildcard_match`/`regex_match`/`region_match` and rejects `and`/`or`.
Since that is the default, `evaluate_df` fails on most real conditions. We pass
`engine="python"`, which matches `evaluate_stream`'s semantics.

## Version constraints we work around

### `parse --drop-seq` and `--add-columns seq` segfault against pysam 0.24

pairtools declares `pysam>=0.15.0` with no upper bound, but ships
`AlignedSegmentPairtoolized`, a compiled Cython subclass of pysam's
`AlignedSegment`. Against pysam 0.24.0 any path that touches a read's sequence
crashes the interpreter:

```
Fatal Python error: Segmentation fault
  File "pairtools/lib/parse.py", line 361 in parse_pysam_entry     # algn["seq"] = sam.seq
  File "pairtools/lib/parse.py", line 1502 in write_pairsam        # sam.query_sequence = ""
```

`pairtools parse` crashes on its own, identically to ours — we call the same
code — so this is an ABI mismatch to pin, not a bug either package can fix in
Python. An upper bound on the pysam requirement, or a rebuild against the
installed pysam, would resolve it.

`tests/test_parse.py:reference` skips a parity case when the pairtools run dies
of a signal, so those option combinations are reported as skipped rather than
silently passing on two identical crashes.

### `region_match` is not in any release

Added in open2c/pairtools#278, after 1.1.2. `lib/select.py` installs upstream's
implementation into `pairtools.lib.select` when it is missing, so that a
condition this package has supported since 0.2.0 keeps working on the released
pairtools. Delete `_backport_region_match` once the declared pairtools minimum
contains it.

## Changes that would let us drop local code

### A column-list return from `write_pairsam`, and a function for `parse`'s columns

`lib/parse.py:write_pairsam` builds a list of fields and then joins and writes
it in one statement, so the only way to get the fields rather than the line is
to split the line again. We pass `streaming_classify` a sink that does exactly
that (via the CSV reader, so the conversion matches the one every .pairs file
goes through), rather than copy the ~40 lines that know how to pack SAM
records and drop sequences. Splitting the function in two — `pairsam_columns()`
returning the list, `write_pairsam()` joining and writing it — would remove the
round trip.

Likewise, the column bookkeeping at `cli/parse.py:parse_py` (and its verbatim
copy in `cli/parse2.py`) is inline in the command body, so a caller that writes
the pairs somewhere else has to restate it. As a function it would be reusable
by both commands and by us.

### A DataFrame-iterator entry point for `dedup`

`lib/dedup.py:_dedup_stream` hard-codes `pd.read_table(in_stream, chunksize=...)`
as its input and `.to_csv(outstream)` as its output. Taking an iterable of
DataFrames instead — with the current signature kept as a thin wrapper that
builds one from a stream — would let any caller feed it Arrow batches without
reimplementing the carryover logic that surrounds `_dedup_chunk`.

### `compute_scaling` should accept the chunk iterator it already loops over

`lib/scaling.py:compute_scaling` iterates over chunks internally:

```python
for pairs_chunk in [pairs_df] if isinstance(pairs_df, pd.DataFrame) else pairs_df:
```

but its entry point accepts only a DataFrame or a path/file-like object and
raises `ValueError` for anything else, so an iterable of DataFrames — which is
what the loop wants — cannot be passed in. Widening that dispatch by one branch
would let `lib/scaling.py` here call `compute_scaling` instead of restating the
loop around `bins_pairs_by_distance`.

`tests/test_scaling.py:test_chunked_loop_matches_compute_scaling` pins our loop
against `compute_scaling` run on a whole file, so the two cannot drift while
the local copy exists.

## Rough edges we match rather than fix

### Chunked `scaling` reports pair counts as floats

`compute_scaling` combines per-chunk results with `DataFrame.add(fill_value=0)`,
which promotes the integer `n_pairs` counts to `float64`. A file smaller than
`--chunksize` is read in one chunk and never goes through `add`, so it reports
`5`; the same file read in two chunks reports `5.0`. The counts are right either
way, but the output of a command is not supposed to depend on a memory knob.

We match it, which means matching pairtools' chunk boundaries exactly rather
than Arrow's — see `lib/chunking.py`. Casting `n_pairs` back to `int64` after
the loop would fix it upstream.

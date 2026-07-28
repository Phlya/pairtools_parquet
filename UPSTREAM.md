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

### `region_match` is not in any release

Added in open2c/pairtools#278, after 1.1.2. `lib/select.py` installs upstream's
implementation into `pairtools.lib.select` when it is missing, so that a
condition this package has supported since 0.2.0 keeps working on the released
pairtools. Delete `_backport_region_match` once the declared pairtools minimum
contains it.

## Changes that would let us drop local code

### A DataFrame-iterator entry point for `dedup`

`lib/dedup.py:_dedup_stream` hard-codes `pd.read_table(in_stream, chunksize=...)`
as its input and `.to_csv(outstream)` as its output. Taking an iterable of
DataFrames instead — with the current signature kept as a thin wrapper that
builds one from a stream — would let any caller feed it Arrow batches without
reimplementing the carryover logic that surrounds `_dedup_chunk`.

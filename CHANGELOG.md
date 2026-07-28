# Changelog


---

## [Unreleased]
### Added
- `dedup --backend duckdb`, now the default, is **6.4x faster than `pairtools
  dedup`** (62.2s -> 9.7s on 5.6M pairs, 4 threads) and produces byte-identical
  output. `--max-mismatch` is 3bp by default and dedup input is sorted, so
  finding duplicates is not a nearest-neighbour search over a plane but a
  lookup in a 3bp window: bucketing on
  `(chrom1, chrom2, strand1, strand2, pos1 // r)` turns it into an equi-join
  whose buckets hold exactly the duplicate families. Profiling the old path
  showed its KD-trees were only 12% of the runtime and pandas bookkeeping the
  rest, so this replaces both. `scipy` and `sklearn` still work and are the
  reference the new backend is tested against.
  With this backend `--chunksize` is a memory knob rather than a semantic one:
  its lookback resolves clusters across a window boundary, so the answer does
  not depend on it. The default rises to 20,000,000 rows accordingly, and
  `-p/--n-proc` now sets DuckDB's thread count (default 4).
- `parse` and `parse2`, ported from their pairtools counterparts. The parser is
  pairtools' own `streaming_classify`, called unchanged; only its output is
  redirected, from formatted text lines into Arrow batches. `parse -o x.parquet`
  produces exactly what `parse -o x.pairs` followed by `csv-to-parquet` would,
  without the text round trip. `--output-stats` and
  `--output-parsed-alignments` are unchanged.
- `split`, ported from `pairtools split`. The pairs go to any supported format
  and the SAM records to text, .bam or stdout. With `--output-pairs` alone and
  Parquet input, the `sam1`/`sam2` columns are never read; from text input
  `pairtools split` is faster, since it splits lines rather than typing fields.
- `header`, ported from `pairtools header`, with all four subcommands —
  `generate`, `transfer`, `set-columns` and `validate-columns`. Headers cross
  formats: the header of a `.parquet` file can be transferred onto a `.pairs`
  file and back. Unlike upstream, `--output` is required, since Parquet has no
  meaningful stdout form.
- `scaling`, ported from `pairtools scaling`. Only six columns take part in the
  binning, so Parquet input reads six columns instead of the whole file. The
  binning itself is pairtools' own `bins_pairs_by_distance`.
- `phase`, ported from `pairtools phase`, in both XB and XA tag modes, with
  `--clean-output` and `--report-scores`.
- `filterbycov`, ported from `pairtools filterbycov`. The coverage calculation
  is pairtools' own `_filterbycov`, called unchanged.
- `restrict`, ported from `pairtools restrict`. About 2.8x faster than
  upstream on a genome-scale fragment file (8.1s -> 2.9s for 540k pairs),
  mostly by not loading the fragment BED with `np.genfromtxt`.
- `flip`, `markasdup` and `sample`, ported from their pairtools counterparts.
  `sample` reproduces pairtools' per-row draw sequence, so a given `--seed`
  selects the same pairs as `pairtools sample` does.
- `stats`, ported from `pairtools stats`, including `--filter`, `--yaml`,
  `--bytile-dups` and `--merge`. Statistics come from pairtools' own
  `PairCounter`, so they are identical.
- `dedup`, ported from `pairtools dedup`. The algorithm is pairtools' own —
  every chunk goes through its `_dedup_chunk` unchanged — so output and
  statistics are identical. The `scipy` and `sklearn` backends are supported;
  pairtools' `cython` backend is not, as it works line-by-line on text streams.
- `merge`, ported from `pairtools merge`. Inputs may be a mix of `.pairs`,
  `.pairs.gz` and `.parquet`; the output is whichever the extension names.
- `lib/arrowio.py`, an Arrow record-batch I/O layer. Every tool reads and writes
  through it, so `.pairs`, `.pairs.gz` and `.parquet` are interchangeable
  wherever a path is accepted. A `.pairs` file round-trips through Parquet
  byte-for-byte, header included.
- `select` now reads and writes any of those formats, not just Parquet, and
  gained `--startup-code` and `--compress-program`.
- `tests/test_parity.py`, which runs each tool against real `pairtools` and
  requires identical output.

### Changed
- The header is now stored verbatim in Parquet under a `pairs_header` key, so
  arbitrary and unknown header lines survive a round trip. The 0.2.0 keys are
  still written, and files that have only those are still readable.
- `select` evaluates its CONDITION with pairtools itself instead of translating
  it to SQL, so the condition language is whatever pairtools supports.

### Fixed
- `sort` no longer reorders pairs tied on chrom1/chrom2/pos1/pos2. Chromosome,
  strand and pair type were cast to DuckDB ENUMs, which order by declaration
  index, and the pair-type ENUM was declared in `itertools.product` order — so
  `UU` sorted before `DD` where pairtools puts `DD` first.
- `sort` no longer aborts on a chromosome absent from the file's header.
- `select` conditions that the SQL translation could not express now work:
  Python method calls and chained comparisons raised, and `wildcard_match`
  with a `?` silently matched nothing.
- `flip` now settles after one pass for every input. For a pair whose two
  sides are on the same chromosome *absent from the chromsizes file*,
  `pairtools flip` compares only the chromosome names — equal, so it swaps the
  sides on every run and never consults the positions. This is a deliberate
  divergence from upstream, confined to exactly those rows; see UPSTREAM.md.
- `restrict` no longer dies on a chromosome missing from the fragment file.
  Upstream catches `ValueError` around a dict lookup that raises `KeyError`,
  so its intended warn-and-continue path is unreachable. Also a deliberate
  divergence; see UPSTREAM.md.
- A `.pairs` file with a header and no data rows can be read, so an empty
  `select` result can be fed to another tool.
- Empty fields no longer become NULL when converting text to Parquet. DuckDB's
  `read_csv` reads an empty field as NULL by default, but `.pairs` has no null:
  an empty field is the empty string, and `phase` reads an empty `XB` tag as
  "no alternative alignment".
- Columns added by `parse --add-columns` now get the type pairtools declares
  for them. `DTYPES_EXTRA_COLUMNS` is keyed by the base name (`mapq`) while the
  columns are named per side (`mapq1`, `mapq2`), so the lookup missed and every
  one of them was typed as a string — which made `select 'mapq1>=30'` fail
  comparing a str to an int, where `pairtools select` works.
- `select` is no longer 22x **slower** than `pairtools select` (201s -> 3.9s on
  5.6M pairs, against pairtools' 8.5s). `evaluate_df`, which 1.0.0 adopted for
  its exact condition semantics, evaluates the condition once per row through
  `DataFrame.iterrows`. Conditions are now rewritten to run over whole columns:
  `and`/`or`/`not` and chained comparisons become their array equivalents on
  the parse tree — not on the text, where `&` binding tighter than `==` would
  change the meaning — and `csv_match`/`wildcard_match`/`regex_match`/
  `region_match` get column-wise versions that evaluate once per distinct
  value. Anything the rewrite cannot express still goes to `evaluate_df`, which
  stays the definition of what a condition means.
- `dedup` no longer loses duplicates at chunk boundaries. `pairtools dedup`
  carries only the previous chunk's *non*-duplicates into the next, so a chain
  of near-duplicates that crosses the boundary through a duplicate is cut and
  the next read is kept as unique. Applies to `--backend duckdb`; the pandas
  backends call upstream's code and still reproduce it. See UPSTREAM.md.
- Header lines no longer lose their trailing whitespace when written as text.
  `pairtools header generate` with no `--assembly` emits `#genome_assembly: `
  with a trailing space, and `headerops.get_header` keeps it, so stripping it
  broke byte parity for a file pairtools wrote. Parquet was unaffected, which
  stores the header verbatim.
- `header set-columns` adds a `#columns:` line to a file that has none.
  Upstream's `headerops.set_columns` only rewrites a line already present, so
  its set-columns is a no-op on exactly the headerless input it exists for.
  Another deliberate divergence; see UPSTREAM.md.
- Text output is now compressed only when the output extension says so, as
  `pairtools` does. `-o out.pairs` previously produced a gzip-compressed file
  named `.pairs` unless `--compress-program none` was passed explicitly.

### Changed (breaking)
- **Renamed the package from `pairs_to_parquet` to `pairtools_parquet`**, along
  with the CLI entry point and the `@PG` provenance records it writes. This is a
  breaking change for imports and for the command name.

### Fixed
- `pairtools_parquet.lib` was missing an `__init__.py` and so was dropped from
  non-editable installs.
- Declared the `pandas` and `numpy` runtime dependencies; dropped the unused
  `fire` dependency.
- Tests no longer write into tracked files under `tests/data/`.

---

## [0.2.0] - 2025-11-26
### Added
- Includes CLI tool for selecting `.parquet` files only.
- Added testing.

---

## [0.1.0] - 2025-10-21
### Added
- First prototype version for internal testing of **pairs-to-parquet**.
- Includes CLI tool for sorting `.pairs` and `.parquet` files.
- Added documentation and optional dependencies for testing and docs.


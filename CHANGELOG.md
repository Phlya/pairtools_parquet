# Changelog


---

## [Unreleased]
### Added
- Arrow IPC as a third format, written and read wherever `.pairs` and
  `.parquet` are: `-o mid.arrow`, and `-o -.arrow` for the standard stream.
  Parquet writes its footer last, so it cannot be read until it is complete and
  therefore cannot go through a pipe; the Arrow IPC *stream* format has no
  footer, so it can. The .pairs header rides in the schema metadata exactly as
  it does in Parquet, and the two agree line for line.
  Input needs nothing declared. An Arrow stream opens with the continuation
  token `ff ff ff ff` and .pairs text with `#`, so the first four bytes settle
  it; `-` accepts either. `-.arrow` on output keeps the format in the path,
  where every other format already lives, rather than in a flag that would have
  had to be threaded through all eighteen tools.
  This turns out to be the fastest way to compose commands, which is not what
  the previous entry predicted. `select | sort | markasdup` over 5.6M pairs,
  same source and same final output:

      Arrow pipes     12.7-13.0s
      Parquet files   14.5-14.8s   116 MB of intermediates
      text pipes      14.9-16.2s
      Arrow files     15.9-16.4s   457 MB
      text files      ~40s cold    401 MB

  A pipe still cannot be seeked, so it still gives up column projection and
  parallel scanning -- but it also skips the text encoding and the disk, and
  the stages run at the same time rather than one after another, which more
  than pays for the loss. Arrow *files* are the worst of both: uncompressed on
  disk and no pipelining. Parquet remains the right choice for anything worth
  keeping.
- `benchmarks/`, so every speedup in the README can be reproduced rather than
  taken on trust. `python benchmarks/run.py` generates a dataset, times each
  tool three ways — `pairtools` text-to-text, ours text-to-text, ours
  Parquet-to-Parquet — and compares the outputs, so a speedup that came from
  doing less work fails the run instead of printing a good number. It takes
  `--pairs`/`--bam` for your own files, `-t` for a subset of tools, `-r` for
  repeats, and reports the interpreter startup cost so a table dominated by
  imports can be recognised as one.
  The dataset is generated rather than downloaded because published `.pairs`
  files are contact lists, which is to say already deduplicated — `dedup` on
  one measures how fast each implementation finds nothing — and the raw BAMs
  behind them are not generally published. `make_data.py` builds a bwa-mem-like
  BAM instead, with PCR duplicates introduced at the read level where PCR
  introduces them, a `--dup-rate` you set, and Illumina read IDs carrying tile
  coordinates; the pairs are `parse`d out of that BAM, so `parse` has real
  input and the duplicates in the pairs are the library's own.
  `tests/test_benchmarks.py` runs the harness end to end at a size that takes
  seconds, and checks the fixture still has the properties the benchmark
  depends on.
- `tests/test_cli_parity.py`, which requires every pairtools command and option
  to exist here under the same name. Nothing checked this before, so a renamed
  option passed the whole suite — the rest of the tests call options by
  whichever spelling this package chose. It found five gaps, all now closed
  (below). Options we deliberately do not accept are listed in the test with
  the reason, so an omission is a decision on the record.
- `filterbycov --backend duckdb`, now the default, is **7.4x faster than
  `pairtools filterbycov`** (136.9s -> 18.4s on 5.6M pairs) with byte-identical
  output. Coverage is a neighbour count — for each end of each pair, how many
  other ends lie within `--max-dist` on the same chromosome — so it gets the
  same treatment as `dedup`: bucket the ends by position and the search is an
  equi-join. pairtools computes it with `_filterbycov`, whose own docstring
  reads "This is a slow version of the filtering code used for testing purposes
  only. Use cythonized version in the future!!"; the cythonized version never
  arrived, and that double loop was 95s of the 110s this tool took. It is still
  available as `--backend python` and is the reference the DuckDB backend is
  tested against.
- `dedup --filter`, `--chrom-subset`, `--output-bytile-stats`, `--startup-code`,
  `-t/--type-cast` and `--engine`, matching `pairtools dedup`. The filters feed
  pairtools' own `PairCounter`, so `--yaml --filter 'unique:(pair_type=="UU")'`
  produces the same statistics file upstream does.
  `--output-bytile-stats` forces `--keep-parent-id`, as upstream does.
- `dedup --tmpdir` and `--memory`, which the DuckDB backend already honoured but
  which only `sort` exposed. They default to DuckDB's own limits rather than to
  `sort`'s 2G.
- `dedup` and `filterbycov` accept `--send-header-to`. It applies to text
  outputs; a `.parquet` keeps its header either way, since there the header is
  key-value metadata rather than leading lines and a file without it cannot be
  read back as pairs.
- `merge --max-nmerge`, matching upstream. Beyond it the merge runs in rounds
  through temporary Parquet files, so a fan-in of thousands does not need every
  input open at once. The header is decided once from the real inputs and handed
  to every round, so a staged merge writes what a single-pass one would.
- `stats --chrom-subset`, which upstream declares and never reads (see
  UPSTREAM.md).
- `dedup --max-mismatch 0` — exact duplicate detection — takes a path with no
  graph in it at all. Exact equality is transitive, so each group of identical
  pairs is already a cluster: no edge list, no connected components, no window.
  **50.8s -> 7.0s against `pairtools dedup`** on 5.6M pairs, of which 2.1s is
  the detection itself and the rest is reading and writing the file.
  It is also **independent of the input's order**, so it needs no sort. On a
  shuffled file it finds exactly the duplicates it finds on a sorted one, where
  `pairtools dedup` compares only within a chunk and misses most of them —
  197,802 reads reported unique against a true 194,189 on a 200k-row shuffled
  file, and 18x slower doing it.
  This holds however large the file is, because the exact path never chunks the
  detection: one aggregate answers the whole file, so there is no window for a
  duplicate family to fall across. `--chunksize` is ignored on this path — a
  5.6M-row shuffled file gives byte-identical output at 1,000, 250,000 and
  20,000,000. What bounds the file size is memory, and DuckDB spills rather than
  changing the answer: the same file capped at `--memory 200M` produces the same
  bytes. Those two knobs, `--memory` and `--tmpdir`, are now on `dedup` as they
  already were on `sort`.
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
  With a non-zero `--max-mismatch`, rows sharing an exact position are
  collapsed before the pairwise search, since they have identical
  neighbourhoods. That keeps the edge list from growing quadratically with
  duplication: on a 90%-duplicate library it is 16M edges instead of 47M, and
  the clustering stage 4.9s instead of 13.5s. It costs about 1.4s where there
  is nothing to collapse.
  With this backend `--chunksize` is a memory knob rather than a semantic one:
  its lookback resolves clusters across a window boundary, so the answer does
  not depend on it. The default rises to 20,000,000 rows accordingly, and
  `-p/--n-proc` now sets DuckDB's thread count (default 4).
  This is the one place the output is not byte-identical to pairtools', and
  the benchmark harness now shows it on real-sized input: `pairtools dedup`
  carries only *non*-duplicate rows into the next chunk, so a chain A~B~C split
  over a boundary loses the B~C link when B was marked a duplicate of A, and C
  is reported unique. Our lookback holds every row and re-decides them, so the
  chain survives. On a 1M-pair library at a 15% duplicate rate it is one row.
  `--backend scipy` and `--backend sklearn` reproduce pairtools exactly. See
  UPSTREAM.md.
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
- `scaling` normalises P(s) by the genome its input declares. Without a
  `--view`, the chromosome sizes now come from the `#chromsize:` header lines.
  `pairtools scaling` extracts them and then throws them away — `pairs_df, _, _
  = pairsio.read_pairs(...)`, with the sizes going into the second `_` — so
  every region's `end` stays at the `-1` sentinel, the area each distance bin
  covers is computed from a region one base long, and `n_bp2` (which is what
  the pair counts are divided by to get P(s)) is meaningless: 48 in total for a
  file whose real answer is 300000. The unmapped `!` chromosome also becomes a
  region of its own, so unmapped pairs land in a scaling curve.
  The `n_pairs` counts are unchanged; what changes is everything they are
  normalised by. The output is now exactly what `pairtools scaling --view <a
  viewframe spelling out the header's chromsizes>` produces, which is what the
  parity tests compare against. A `--view` still takes precedence, as upstream,
  and a header that declares no sizes still gets upstream's answer. See
  UPSTREAM.md.
- `-o` defaults to stdout and the input path to stdin, as they do in pairtools,
  so a pipeline needs no `-` at all:
  `select '...' in.pairs | sort | markasdup -o out.pairs`. Both were required
  arguments here, which meant every composed command had to spell out `-o -`
  and a trailing `-`. pairtools reads `not path or path == "-"`, and an empty
  path, no path and `-` now all mean the standard stream here too.
- `dedup` reads stdin, by spooling it. It needs two passes over the pairs --
  one to find the duplicates, one to write them out with the answer applied --
  which a pipe cannot give it, so the stream is written to a temporary Parquet
  file and deduplicated from there. Parquet because the second pass is then a
  native DuckDB scan over a compressed, seekable file, which is what the
  algorithm wants anyway; the temporary file lands in `--tmpdir` when one is
  given. The cost is one extra write and read, which is what writing a Parquet
  intermediate by hand would have cost: `markasdup | dedup` over 5.6M pairs
  takes 10.5-12.0s through a pipe against 11.0-11.4s through a Parquet file.
  It used to read the spent stream and write an empty file, reporting success.
- `-` works in both directions, so the tools can be composed with pipes:
  `select ... -o - in.pairs | sort -o - - | markasdup -o - -`. Text only —
  Parquet's footer is written at the end of the file, so it cannot go to a
  stream, and `-o -` always means `.pairs`.
  Neither half worked before, and neither failed loudly. `-o -` opened a file
  literally named `-` in the working directory and left stdout empty, exiting
  0, so a pipeline produced an empty result and a stray file rather than an
  error. And any tool that read the header to decide what to do before opening
  the input for real — `sort` resolving its sort keys, `dedup` its columns —
  reopened a pipe that had already given up its header, and stopped with "has
  no '#columns:' header line". stdin is now read once and the stream kept
  positioned after the header, and stdout is written directly, wrapped so that
  pyarrow closing its sink does not close the process's stdout.
- `dedup` refuses `-` as input rather than silently producing nothing. It reads
  the pairs twice — once to find the duplicates, once to write them out with
  the answer applied — which a pipe cannot serve. It used to read the spent
  stream and write an empty file, reporting success.
  Everything else pipes: `select`, `sort`, `markasdup`, `flip`, `sample`,
  `merge` (including a mix of files and stdin), `restrict`, `stats`, `scaling`
  and `filterbycov` all give the same answer from a pipe as from a file, which
  `tests/test_stdio.py` now checks by comparing the two.
  Worth knowing when composing: piping *text* is not the fast path. A pipe
  cannot be seeked, so it gives up both column projection and parallel
  scanning, and text costs an encode and a decode on top -- `select | sort |
  dedup` over 5.6M pairs takes 24.6s through Parquet files against 34.6s
  through text files. Piping Arrow is a different matter; see the Arrow IPC
  entry above.
- `merge` is now a merge rather than a sort, and gives the same answer whatever
  format its inputs are in. Rows tied on all five sort keys came out in a
  different order from Parquet input than from text — and in a different order
  between runs of the same command. Unmapped pairs all carry `! 0 ! 0`, so this
  was not a corner case: 364,003 rows of the 5.6M benchmark, every one of them
  a tie. No pair was ever lost, gained or altered, but the files did not
  compare equal.
  Both paths ran one `UNION ALL ... ORDER BY chrom1, chrom2, pos1, pos2,
  pair_type`. Nothing followed those five keys, and DuckDB's sort is not
  stable, so tied rows emerged in whatever order the parallel scans handed
  them over. Ordering by `(input index, row number within input)` after the
  keys reproduces what `sort --merge` does: `pairtools merge` passes GNU sort
  neither `-s` nor a last-resort comparison, and `--merge` does not re-sort
  tied rows, so upstream's order is the first file's ties followed by the
  second's, each in file order. Verified against the reference output rather
  than assumed.
  Getting the row number needs different machinery per format. Parquet numbers
  its own rows, so it keeps the native scan (`read_parquet(file_row_number=
  true)`) and its timing is unchanged. DuckDB's CSV scanner cannot: it is
  parallel, so `row_number() OVER ()` counts in arrival order, not file order —
  the same trap `dedup` documents. Text inputs are therefore read through a
  sequential Arrow reader, which costs them DuckDB's parallel CSV parsing:
  merging two 5.6M-row text files takes about 7.3s against `pairtools merge`'s
  5.9s, best of four runs each. Parquet in and Parquet out, the path this
  package exists for, is unaffected at 3.3s.
  Two cheaper-looking alternatives were measured and rejected, so they need not
  be tried again. Reading the CSV with `parallel=false` and numbering rows with
  `row_number() OVER ()` keeps upstream's semantics and looked much faster in
  isolation, but end to end it was slower than the Arrow reader on every path.
  Dropping the row numbers entirely and tie-breaking on every remaining column
  is deterministic by construction -- rows still tied are byte-identical, so
  their order cannot be observed -- and is the fastest option for text input,
  but it is slower for Parquet, more so for `.pairsam` where the SAM blobs
  bloat the sort payload, and it abandons byte-identity with `pairtools merge`.
  Nothing in the .pairs specification or pairtools' documentation constrains the
  order of tied rows, so that last point is a choice rather than a requirement;
  it can be revisited if text input ever matters more than parity.
  The row-numbering helper is now `arrowio.with_row_ids`, shared with `dedup`
  rather than duplicated.
  Found by `benchmarks/run.py`, which compares outputs as well as timing them
  and fails the run when our own two paths disagree. It reproduces only at
  scale, and only sometimes — the regression test asserts the stable order
  itself, which holds regardless of how the scan happened to be scheduled.
- `stats` re-cuts its input to `--chunksize` rows before counting, so text and
  Parquet input now produce identical output. Text batches are sized in bytes
  and Parquet batches in rows, and `PairCounter.add_pairs_from_dataframe`
  groups each chunk by chromosome pair with pandas — which sorts group keys —
  so which chunk a pair landed in decided the order the `chrom_freq` lines were
  written in. The counts were always right; the files did not `diff` clean
  against `pairtools stats`, and now do. `dedup` and `scaling` already re-cut
  their input for the same reason. Found by the new benchmark harness, which
  compares outputs as well as timing them.
  A consequence worth knowing: `--chunksize` now genuinely reaches text input,
  and since the ordering follows chunk boundaries, a non-default `--chunksize`
  reorders the `chrom_freq` lines (values unchanged). That is inherited from
  pairtools, which hard-codes 100,000 and exposes no way to vary it.
- `flip`'s chromosome-order file is `-c/--chroms-path`, as in pairtools. It was
  named `--chrom-subset` here, which is not only a different name but the name
  of an unrelated pairtools option — a filter, in `select`, `stats` and `dedup`.
  Every test used `-c`, so nothing caught it.
- `flip` is **2.7x faster** (15.9s -> 6.0s on 5.6M pairs) and no longer slower
  than `pairtools flip`. Ranking a chromosome column meant materializing it as
  Python strings to look the names up in a dict: 7.1s per column, twice per
  batch. `pc.index_in` against the chromosome order does it in Arrow in 0.4s.
  Only rows whose *both* sides are missing from the chromsizes file are still
  compared by name, and only those rows are materialized.
- `restrict` is **3.3x faster** (25.3s -> 7.6s), 5.9x against `pairtools
  restrict`. Same cause, plus `np.unique` over a column of Python strings — 7.4s
  of the 23s, and the same numpy pathology that shaped `dedup`'s clustering.
  Dictionary-encoding the chromosome column gives the distinct values for free
  and turns per-chromosome row selection into an integer comparison.
- `select --chrom-subset` tests membership with `pc.is_in` rather than rebuilding
  a Python list and calling `np.isin` on decoded strings once per batch.
- `sample` draws with `iter(rng.random, None)` instead of a generator
  expression. The draws have to stay one-per-row to keep `--seed` meaning what
  it means in pairtools, but the per-row Python frame does not: 1.6s of frame
  overhead against 0.38s of actual random number generation.
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
- Project URLs now point at `Phlya/pairtools_parquet` and
  `pairtools-parquet.readthedocs.io`, following the repository rename. They
  pointed at `ayaksvals/pairs_to_parquet`, the upstream this was forked from,
  which is neither where this work lives nor reachable under the new name. The
  `Changelog` URL was also broken independently of the rename — it pointed at
  `/issues/blob/master/CHANGES.md`, a path with an extra `issues/` segment and
  the wrong filename.

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


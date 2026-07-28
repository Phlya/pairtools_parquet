# Changelog


---

## [Unreleased]
### Added
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


# Benchmarks

Where the numbers in the top-level `README.md` come from, and how to get your
own:

```sh
pip install -e '.[test]' pairtools
python benchmarks/run.py
```

That generates a dataset under `benchmarks/data/` (about 50s and 200 MB for the
default million read pairs, reused on later runs), then times each tool three
ways and checks that all three produced the same pairs:

| column | what it runs |
|---|---|
| `pairtools` | `pairtools <tool>`, text in and text out — the baseline |
| ours (text) | `pairtools_parquet <tool>`, the same text in and out |
| ours (parquet) | `pairtools_parquet <tool>`, Parquet in and Parquet out |

The arguments are identical across all three columns — this package mirrors
`pairtools` option for option, so anything that needed different arguments for
the two engines would be measuring different work.

Outputs are compared as well as timed, so a speedup that came from skipping
work fails the run instead of printing a good number (`--no-check` turns that
off). The `match` column answers two questions at once:

- our text and Parquet columns must always agree with each other — the format a
  file is in must not change the answer, and if it does the run fails whatever
  else it found;
- agreeing with `pairtools` is also required, except where a benchmark declares
  a divergence and says why. Those print as a row count with a `*`, and the
  reason is printed under the table. Today there is one: `dedup` at the default
  `--max-mismatch 3` keeps near-duplicate chains that `pairtools dedup` cuts at
  its chunk boundaries — a handful of rows in a million, explained in
  `UPSTREAM.md`.

## Options worth knowing

```sh
python benchmarks/run.py -n 5000000        # a bigger dataset (~1 GB, ~4 min)
python benchmarks/run.py -t dedup,restrict # just these tools
python benchmarks/run.py -r 3              # best of three runs each
python benchmarks/run.py -p 8              # --nproc passed to both engines
python benchmarks/run.py --dup-rate 0.4    # a badly over-amplified library
python benchmarks/run.py --json out.json   # machine-readable results
```

Changing `-n`, `--dup-rate` or `--seed` regenerates the dataset; a manifest
records what the cached one was built from, so re-running the same settings
does not rebuild it. `--regenerate` forces a rebuild.

`filterbycov` is the slowest row by a wide margin, because `pairtools
filterbycov` is a per-pair Python loop over neighbours — seven minutes against
eight seconds on the default dataset. `-t` past it if you are iterating.

Every measurement is a fresh subprocess, so Python startup and the pyarrow /
duckdb / pandas imports are inside the number, exactly as they are when you run
the command yourself. That is around 1.5s on either engine, which is why small
datasets make every tool score 1.0x — the harness measures the startup cost and
warns when the table is dominated by it. A million pairs is comfortably past
that point; below ~200k pairs the numbers mean nothing.

## Your own data

```sh
python benchmarks/run.py \
    --pairs mine.pairs -c mine.chrom.sizes -f mine.frags.bed
```

`--pairs` has to be text (`.pairs` or `.pairs.gz`) — the `pairtools` column
needs a file pairtools can read. The Parquet copy is converted from it once,
before timing starts; pass `--parquet mine.parquet` instead if you already have
one and want that exact file measured. Add `--bam mine.bam` to include `parse`,
which is skipped without it; `--bam` on its own runs `parse` and nothing else.

Two things to watch for in your own file: `dedup` only means something on
input that still contains duplicates, and `dedup`/`filterbycov` both require
sorted input.

## Why the default data is synthetic

The awkward part of benchmarking this toolset is that `parse` needs an aligned
BAM and `dedup` needs pairs that have not been deduplicated, and public data
gives you neither:

- Published `.pairs` files — 4DN, ENCODE — are *contact lists*, which is to say
  the output of a pipeline that already ran `dedup`. Running `dedup` on one
  measures how fast each implementation finds nothing, and it also skews every
  other tool, since the duplicate-heavy positions are exactly the ones that
  make `filterbycov` and `restrict` work hard.
- The raw BAMs behind them are not generally published, and the FASTQs that
  are run to tens or hundreds of GB and would need `bwa` and an indexed genome
  before this harness had anything to time.

So the default dataset is generated instead, by `make_data.py`, from one seed:
a bwa-mem-like BAM of `-n` read pairs over a 590 Mb mock genome, two thirds
cis with log-uniform cis distances, a realistic scattering of unmapped and
multi-mapping reads, Illumina-style read IDs with tile coordinates (so `dedup
--output-bytile-stats` has something to group by), and PCR duplicates
introduced where PCR introduces them — at the read level before alignment,
with a few bp of wobble so `--max-mismatch` has something to do. The `.pairs`
and `.parquet` files are then `parse`d and `sort`ed out of that BAM, so the
duplicates in them are the library's duplicates rather than something written
in afterwards.

This is a benchmark fixture, not a simulation: the P(s) is approximate and
there are no restriction sites, so do not draw biology from it. What it does
give you is the one thing a public file cannot — a duplicate rate you set with
`--dup-rate`, which is the single number `dedup`'s cost depends on most.

If you want a real-data check anyway, the harness runs on anything you point
`--pairs` at, including a 4DN contact list. Expect `dedup` to report a
duplicate rate near zero on one.

# pairtools_parquet

[![Documentation Status](https://readthedocs.org/projects/pairtools-parquet/badge/?version=latest)](https://pairtools-parquet.readthedocs.io/en/latest/)
[![Join the chat on Slack](https://img.shields.io/badge/chat-slack-%233F0F3F?logo=slack)](https://bit.ly/2UaOpAe)


## Transform 3D contacts (.pairs) to parquet & process them

`pairtools_parquet` is a .parquet extention of .pairs file format.

The main purpose of this extension is to leverage the row groups and metadata features of the Parquet format in order to:

- speed up data selection, filtering and sorting
- address a minor limitation of the .pairs format, where metadata cannot be easily parsed by generic CSV readers
- reduce storage space required for the data 
- improve I/O performance



## Data formats

There are 2 main file formats, which are used by our converter & processor: 

1. `.pairs`: 
`pairtools` produce and operate on tab-separated files compliant with the [.pairs](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md) format defined by the [4D Nucleome Consortium](https://www.4dnucleome.org/). All pairtools properly manage file headers and keep track of the data processing history.

2. `.parquet`: 
a columnar or hybrid file format, which is highly optimized for big data processing. It supports features like predicate pushdown and column projection -> better query performance by minimizing data read from disk.
In our workflow, .parquet files are processed using [duckdb](https://duckdb.org/docs/stable/data/parquet/overview) - an open-source column-oriented Relational Database Management System. 
For more information, see [parquet metadata](https://duckdb.org/docs/stable/data/parquet/metadata) by duckdb

## Operations: 
- convert .pairs -> .parquet
- convert .parquet -> .pairs
- sort pairs in lexycographic order

## Installation

Requirements:
- Python 3.x

Currently there is only 1 option for installing `pairtools_parquet`:

And it is the same, when you want to modify `pairtools_parquet`: build `pairtools_parquet` from source via pip's "editable" mode:

```sh
$ git clone https://github.com/Phlya/pairtools_parquet
$ cd pairtools_parquet
$ pip install -e .
```

## Tools

- `csv-to-parquet`: transform standard .pairs files into the optimized .parquet format for faster querying and reduced storage. Header of .pairs becomes key-value metadata in a new parquet file

- `parquet-to-csv`: export Parquet data back into .pairs format for compatibility with existing pairtools pipelines.

- `sort`: sort .pairs or .parquet files(the lexicographic order for chromosomes, the numeric order for the positions, the lexicographic order for pair types)

- `select`: filter pairs by a `pairtools select` condition. 1.6x faster than `pairtools select`. Conditions are rewritten to evaluate over whole columns; anything the rewrite cannot express falls back to pairtools' own evaluator, so the full condition language works, including `--startup-code`.

- `merge`: merge sorted files, keeping them sorted. Inputs may be a mix of formats. Rows tied on every sort key keep their input order, as `sort --merge` gives them, so the result does not depend on which format the inputs were in.

- `dedup`: find and remove PCR/optical duplicates, with statistics. `--backend scipy` restores pairtools' KD-tree implementation and reproduces its output exactly; the default backend differs from it only by keeping near-duplicate chains that `pairtools dedup` cuts at its chunk boundaries — 1 row in a million-pair library, 16 in a 5.6M one, since chains get likelier as coverage rises (UPSTREAM.md).

  - **`--max-mismatch 0` (exact): 8.8x faster** — 60.7s → 6.9s on 5.6M pairs. Exact equality is transitive, so each group of identical pairs is already a cluster and no graph is needed; the detection itself is 2.1s and the rest is I/O. This path does not care what order the file is in, so it needs no sort — and unlike `pairtools dedup`, which only compares within a chunk, it gives the same answer on a shuffled file as on a sorted one. It never chunks the detection either, so file size does not reintroduce a boundary: `--chunksize` is ignored, and where the key table outgrows `--memory` DuckDB spills to `--tmpdir` rather than comparing fewer rows.
  - **Default `--max-mismatch 3`: 6.4x faster** — 57.7s → 9.0s. A 3bp tolerance makes this a lookup in a 3bp window rather than a nearest-neighbour search, so bucketing on `(chrom1, chrom2, strand1, strand2, pos1 // r)` turns it into an equi-join. Rows at identical positions are collapsed first, which keeps the edge list linear rather than quadratic in the duplication rate.

- `flip`: reflect pairs onto the upper triangle, given a chromosome order. 1.5x faster than `pairtools flip`.

- `markasdup`: tag every pair as a duplicate.

- `sample`: take a seeded random subset of pairs. A given `--seed` selects the same pairs `pairtools sample` does, which means keeping its one-draw-per-row sequence, so this is the one tool where text input is not faster than pairtools — it types every field where `pairtools sample` only splits lines.

- `filterbycov`: remove pairs from regions of unusually high coverage. **50x faster** (473.5s → 9.4s on 1M pairs; the ratio grows with coverage, and at 5.6M the reference run takes hours). Coverage is a count of neighbouring pair ends, so it becomes an equi-join over position buckets, exactly as in `dedup`. `--backend python` restores pairtools' own double loop, which its docstring describes as "a slow version of the filtering code used for testing purposes only".

- `phase`: assign pairs to parental haplotypes in a diploid genome.

- `restrict`: assign pairs to restriction fragments. **5.4x faster** (38.1s → 7.1s), by not loading the fragment BED with `np.genfromtxt` and by ranking chromosomes as dictionary codes rather than as Python strings.

- `stats`: summary statistics, and merging of stats files.

- `scaling`: contact frequency as a function of genomic distance. Reads only the six columns the binning looks at, so Parquet input skips the rest of the file.

- `header`: `generate`, `transfer`, `set-columns` and `validate-columns`. Headers move between formats, so one generated onto a `.parquet` can be transferred onto a `.pairs` and back.

- `parse` / `parse2`: extract pairs from `.sam`/`.bam`, writing Parquet directly. The parser is pairtools' own, and pysam dominates the runtime, so parsing itself is no faster — what goes is the pipeline's first text round trip (on 400k read pairs: 5.4s straight to Parquet against 5.5s + 1.2s for `pairtools parse` then `csv-to-parquet`).

- `split`: separate a `.pairsam` back into pairs and SAM. From Parquet, `--output-pairs` alone reads only the columns it writes, skipping `sam1`/`sam2` entirely. From text there is no such saving — every column has to be parsed either way, and `pairtools split` is faster there, since it splits lines rather than typing fields.

Every tool takes `.pairs`, `.pairs.gz`, `.parquet` or `.arrow` as input and writes whichever of those the output path's extension names, so the formats are interchangeable wherever a path is accepted.

## Composing commands

Omit `-o` and it writes to stdout; omit the input path and it reads stdin, as
in pairtools. So the tools pipe with nothing extra to say:

```sh
pairtools_parquet select 'pair_type=="UU"' in.pairs \
  | pairtools_parquet sort \
  | pairtools_parquet markasdup -o out.pairs
```

`-` says the same thing explicitly, and is what you need when a command takes
several inputs — `merge a.pairs -` reads the second from the pipe.

`-o -` writes .pairs. To keep the data binary on the wire, use `-.arrow`:
the Arrow IPC *stream* format has no footer, so unlike Parquet it can be read
as it arrives. Input needs no flag either way — a stream starting with Arrow's
continuation token is Arrow, and .pairs text always starts with `#`.

```sh
pairtools_parquet select 'pair_type=="UU"' -o -.arrow in.parquet \
  | pairtools_parquet sort -o -.arrow - \
  | pairtools_parquet markasdup -o out.pairs -
```

`dedup` reads the pairs twice, once to find the duplicates and once to write
them out, which a pipe cannot serve directly — so it spools the stream to a
temporary Parquet file and works from that. It costs about what writing a
Parquet intermediate yourself would cost, because that is what it is; the
difference is not having to manage it.

**Arrow through a pipe is the fastest way to compose**, because it skips both
the text encoding and the disk, and because the stages then run at the same
time instead of one after another. `select | sort | markasdup` over 5.6M pairs,
same source and same final output, differing only in what goes between:

| between steps | time | intermediate size |
|---|---|---|
| **Arrow pipes** | **12.7–13.0s** | — |
| Parquet files | 14.5–14.8s | 116 MB |
| Arrow files | 15.9–16.4s | 457 MB |
| text pipes | 14.9–16.2s | — |
| text files | ~40s cold | 401 MB |

Arrow *files* are the worst of both: uncompressed on disk and no pipelining.
Use Parquet when a step's output is worth keeping, and `-.arrow` when it is
not.


## Benchmarks

Every number here is reproducible:

```sh
python benchmarks/run.py -n 5600000 -r 2
```

which generates its own dataset, times each tool against `pairtools`, and
checks that both produced the same pairs — see [benchmarks/README.md](benchmarks/README.md)
for running it on your own files, and for why the default data is generated
rather than downloaded.

5.6M pairs over a 590 Mb genome with a 15% duplicate rate, 4 threads, best of
two runs. "ours" is Parquet in and Parquet out. Each measurement is a fresh
process, so Python startup — about 1.5s on either side — is inside every
number, as it is when you run the command yourself.

| Tool | `pairtools` | ours | | output |
|---|---|---|---|---|
| `filterbycov` † | 473.5s | **9.4s** | 50x | identical |
| `dedup --max-mismatch 0` | 60.7s | **6.9s** | 8.8x | identical |
| `dedup` (default 3bp) | 57.7s | **9.0s** | 6.4x | 16 rows differ ‡ |
| `restrict` | 38.1s | **7.1s** | 5.4x | identical |
| `sort` | 9.1s | **3.0s** | 3.1x | identical |
| `stats` | 14.3s | **5.1s** | 2.8x | identical |
| `markasdup` | 13.5s | **5.1s** | 2.7x | identical |
| `scaling` | 19.8s | **9.0s** | 2.2x | identical |
| `select` | 10.2s | **6.4s** | 1.6x | identical |
| `merge` | 5.9s | **3.3s** | 1.8x | identical ¶ |
| `flip` | 9.6s | **6.3s** | 1.5x | identical |
| `parse` | 97.0s | **104.7s** | 0.9x | identical |
| `sample` | 2.3s | **3.2s** | 0.7x | identical |

† `filterbycov` is measured on 1M pairs, not 5.6M. `pairtools filterbycov` is a
per-pair Python loop over neighbouring ends, so its cost grows with the square
of coverage: at 5.6M on this genome the reference run alone would take hours.
Its ratio is the one that depends most on your library — the denser the data,
the further ahead the equi-join gets.

‡ Not identical, and deliberately so.
`pairtools dedup` carries only *non*-duplicate rows into the next chunk, so a
near-duplicate chain A~B~C split across a boundary loses the B~C link when B
was marked a duplicate of A, and C is reported unique. Our lookback holds every
row and re-decides them, so the chain survives. It scales with density — 1 row
per million pairs at 1M, 16 at 5.6M. `--backend scipy` reproduces pairtools
exactly. See UPSTREAM.md.

¶ `merge` is the one tool whose *text* input is slower than pairtools —
about 7.3s against 5.9s, best of four runs each. Rows tied on all five sort
keys — every unmapped pair shares `! 0 ! 0` — need an explicit tie-break to
come out in upstream's order, since DuckDB's sort is not stable, and the row
numbers that tie-break needs are free in Parquet but not in DuckDB's parallel
CSV scanner. Text inputs are read through a sequential Arrow reader to get
them. Parquet is unaffected.

`sample` and `parse` are the two that are not faster, both structurally.
Reproducing `pairtools sample`'s `--seed` means keeping its one-draw-per-row
sequence, and the pairs then go through typed columns where `pairtools sample`
only decides whether to copy a line. `parse` runs pairtools' own parser with
pysam dominating; what a Parquet pipeline saves there is the text round trip
that would otherwise follow, not the parse itself.

These ratios are not constants — they depend on the library, `filterbycov` most
of all. Run the harness on your own data if the number matters to you.

## Why to use `.parquet` extention for sorting (and many more future processing tools)?
If we use the same 2.4 GB file, 35 GB of memory, 4 threads:



| Tool, input & output formats          | Memory (2.4GB)| real time | user time | sys time  |                              Comments                                 |
|---------------------------------------|---------------|-----------|-----------|-----------|-----------------------------------------------------------------------|
| Pairtools sort                        |     2.3GB     | 10min 10s | 20min 23s | 3min 14s  |                                                                       |
| pairtools_parquet csv-parquet sort     |     2.5GB     | 2min  24s | 6min  56s | 0min 47s  | 5x times faster in real time, 3x times faster in user time            |
| pairtools_parquet parquet-parquet sort |     2.6GB     | 2min  33s | 6min  18s | 2m   6s   | also a major speed up                                                 |
| pairtools_parquet csv-csv sort         |     2.2GB     | 4mim  39s | 15min 55s | 0min 53s  | worse, that first 2; better, than pairtools sort; better compression  |
| pairtools_parquet parquet-csv sort     |     2.6GB     | 5min  15s | 14min 10s | 1m   56s  | worse, that first 2, but still better, than pairtools sort            |

So pairtools_parquet with any input and output format will outperform pairtools sort on csv. Here csv-parquet and parquet-parquet show the best results. 
Spoiler alert: on bigger files, like 10GB compressed, the difference feels even more dramatic. pairtools sort ~25 min, pairtools_parquet sort csv-parquet ~12 minutes. 

Working directly with Parquet files (parquet → parquet sort) delivers performance close to the best case, confirming that the Parquet format maintains efficiency across repeated operations.

As a result, switching from .pairs (CSV) to .parquet for sorting (and we will show in the future other data processing) yields 3–4× faster runtimes, better I/O performance, and improved scalability for large datasets.

So welcome to the world of Parquet, it is been waiting for you! 
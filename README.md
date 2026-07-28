# pairtools_parquet

[![Documentation Status](https://readthedocs.org/projects/pairs_to_parquet/badge/?version=latest)](https://pairs-to-parquet.readthedocs.io/en/latest/)
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
$ git clone https://github.com/ayaksvals/pairs_to_parquet
$ cd pairs_to_parquet
$ pip install -e .
```

## Tools

- `csv-to-parquet`: transform standard .pairs files into the optimized .parquet format for faster querying and reduced storage. Header of .pairs becomes key-value metadata in a new parquet file

- `parquet-to-csv`: export Parquet data back into .pairs format for compatibility with existing pairtools pipelines.

- `sort`: sort .pairs or .parquet files(the lexicographic order for chromosomes, the numeric order for the positions, the lexicographic order for pair types)

- `select`: filter pairs by a `pairtools select` condition. 1.4x faster than `pairtools select`. Conditions are rewritten to evaluate over whole columns; anything the rewrite cannot express falls back to pairtools' own evaluator, so the full condition language works, including `--startup-code`.

- `merge`: merge sorted files, keeping them sorted. Inputs may be a mix of formats.

- `dedup`: find and remove PCR/optical duplicates, with statistics. Byte-identical output, `--backend scipy` restores pairtools' KD-tree implementation.

  - **`--max-mismatch 0` (exact): 7.3x faster** — 50.8s → 7.0s on 5.6M pairs. Exact equality is transitive, so each group of identical pairs is already a cluster and no graph is needed; the detection itself is 2.1s and the rest is I/O. This path does not care what order the file is in, so it needs no sort — and unlike `pairtools dedup`, which only compares within a chunk, it gives the same answer on a shuffled file as on a sorted one. It never chunks the detection either, so file size does not reintroduce a boundary: `--chunksize` is ignored, and where the key table outgrows `--memory` DuckDB spills to `--tmpdir` rather than comparing fewer rows.
  - **Default `--max-mismatch 3`: 6.4x faster** — 62.2s → 9.7s. A 3bp tolerance makes this a lookup in a 3bp window rather than a nearest-neighbour search, so bucketing on `(chrom1, chrom2, strand1, strand2, pos1 // r)` turns it into an equi-join. Rows at identical positions are collapsed first, which keeps the edge list linear rather than quadratic in the duplication rate.

- `flip`: reflect pairs onto the upper triangle, given a chromosome order. 1.5x faster than `pairtools flip`.

- `markasdup`: tag every pair as a duplicate.

- `sample`: take a seeded random subset of pairs. A given `--seed` selects the same pairs `pairtools sample` does, which means keeping its one-draw-per-row sequence, so this is the one tool where text input is not faster than pairtools — it types every field where `pairtools sample` only splits lines.

- `filterbycov`: remove pairs from regions of unusually high coverage. **7.4x faster** (136.9s → 18.4s on 5.6M pairs). Coverage is a count of neighbouring pair ends, so it becomes an equi-join over position buckets, exactly as in `dedup`. `--backend python` restores pairtools' own double loop, which its docstring describes as "a slow version of the filtering code used for testing purposes only".

- `phase`: assign pairs to parental haplotypes in a diploid genome.

- `restrict`: assign pairs to restriction fragments. **5.9x faster** (44.2s → 7.6s), by not loading the fragment BED with `np.genfromtxt` and by ranking chromosomes as dictionary codes rather than as Python strings.

- `stats`: summary statistics, and merging of stats files.

- `scaling`: contact frequency as a function of genomic distance. Reads only the six columns the binning looks at, so Parquet input skips the rest of the file.

- `header`: `generate`, `transfer`, `set-columns` and `validate-columns`. Headers move between formats, so one generated onto a `.parquet` can be transferred onto a `.pairs` and back.

- `parse` / `parse2`: extract pairs from `.sam`/`.bam`, writing Parquet directly. The parser is pairtools' own, and pysam dominates the runtime, so parsing itself is no faster — what goes is the pipeline's first text round trip (on 400k read pairs: 5.4s straight to Parquet against 5.5s + 1.2s for `pairtools parse` then `csv-to-parquet`).

- `split`: separate a `.pairsam` back into pairs and SAM. From Parquet, `--output-pairs` alone reads only the columns it writes, skipping `sam1`/`sam2` entirely. From text there is no such saving — every column has to be parsed either way, and `pairtools split` is faster there, since it splits lines rather than typing fields.

Every tool takes `.pairs`, `.pairs.gz` or `.parquet` as input and writes whichever of those the output path's extension names, so the formats are interchangeable wherever a path is accepted.


## Benchmarks

Every tool against `pairtools` itself, on 5.6M pairs, 4 threads. "ours" is
Parquet in and Parquet out; output is byte-identical to pairtools in every row.

| Tool | `pairtools` | ours | |
|---|---|---|---|
| `filterbycov` | 136.9s | **18.4s** | 7.4x |
| `restrict` | 44.2s | **7.6s** | 5.9x |
| `dedup --max-mismatch 0` | 50.8s | **7.0s** | 7.3x |
| `dedup` (default 3bp) | 62.2s | **9.7s** | 6.4x |
| `stats` | 10.1s | **3.7s** | 2.7x |
| `sort` | 6.2s | **2.6s** | 2.4x |
| `scaling` | 15.5s | **8.3s** | 1.9x |
| `flip` | 8.8s | **6.0s** | 1.5x |
| `select` | 9.0s | **6.4s** | 1.4x |
| `markasdup` | 6.6s | **5.2s** | 1.3x |
| `sample` | **2.3s** | 2.8s | 0.8x |

`sample` is the one tool that is slower, and structurally so: reproducing
`pairtools sample`'s `--seed` means keeping its one-draw-per-row sequence, and
the pairs then go through typed columns where `pairtools sample` only has to
decide whether to copy a line. `split` from text input is slower for the same
reason. Both are noted in their entries above.

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
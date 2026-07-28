# -*- coding: utf-8 -*-
"""Generate a synthetic Hi-C dataset to benchmark against.

Why synthetic rather than a public download: the benchmark needs an aligned
BAM (so `parse` has something to parse) whose pairs are *not* deduplicated (so
`dedup` has something to remove). Public `.pairs` files -- 4DN, ENCODE -- are
published as finished contact lists, which is to say post-dedup, so running
`dedup` on one measures the cost of finding nothing. The raw BAMs those were
built from are tens to hundreds of GB, which is not a thing to download to run
a benchmark. See benchmarks/README.md for pointing the harness at real data.

What is generated, from one seed:

  bench.bam           N read pairs, bwa-mem-like: two records per read, no
                      proper-pair flag, mates adjacent, Illumina read IDs
  bench.pairs         `parse` of that BAM, sorted -- so the pairs really are
                      the BAM's pairs, and the duplicates in them are the
                      library's duplicates
  bench.parquet       the same, converted
  bench.chrom.sizes   the genome
  bench.frags.bed     restriction fragments, for `restrict`

Duplicates are introduced where PCR puts them -- at the read level, before
alignment -- so a duplicate is a second read pair at the same coordinates
under a different read ID, with a few bp of wobble to give `--max-mismatch`
something to do. Family sizes are geometric with mean ``1 / (1 - dup_rate)``,
so most are pairs and a few are deep.
"""

import argparse
import json
import os
import random
import subprocess
import sys

# A read ID that pairtools can extract a tile from, so `dedup
# --output-bytile-stats` has something to group by. The fields are
# instrument:run:flowcell:lane:tile:x:y, and pairtools reads flowcell:lane:tile.
READ_ID = "A00111:222:HGKMNDSXX:{lane}:{tile}:{x}:{y}"

READ_LENGTH = 80
LANES = 4
TILES = 24

# Enough chromosomes that cis/trans and the chromosome-order comparisons in
# `flip` and `sort` are exercised, sized so a 1kb fragment map is not enormous.
CHROMS = [("chr{}".format(i), 30_000_000) for i in range(1, 20)] + [
    ("chrX", 20_000_000),
    ("chrM", 16_000),
]

FRAGMENT_LENGTH = 1000


def read_id(rng):
    return READ_ID.format(
        lane=rng.randrange(1, LANES + 1),
        tile=rng.randrange(1101, 1101 + TILES),
        x=rng.randrange(1000, 30000),
        y=rng.randrange(1000, 30000),
    )


def sam_record(rid, flag, chrom, pos, mapq, seq, tags):
    return "\t".join(
        [rid, str(flag), chrom, str(pos), str(mapq), "{}M".format(len(seq)),
         "*", "0", "0", seq, "I" * len(seq)] + list(tags)
    )


def ligation(rng):
    """One Hi-C molecule: two genomic positions and their strands.

    Two thirds cis, matching the cis fraction of a passable Hi-C library, and
    the cis distances are drawn log-uniformly so P(s) looks roughly like P(s)
    -- `scaling` bins by log distance, so a uniform draw would pile every pair
    into the last few bins.
    """
    chrom1, length1 = rng.choice(CHROMS)
    pos1 = rng.randrange(1, length1 - READ_LENGTH)
    if rng.random() < 0.66:
        chrom2, length2 = chrom1, length1
        distance = int(10 ** rng.uniform(2, 7))
        pos2 = pos1 + distance * rng.choice([1, -1])
        if not 1 <= pos2 < length2 - READ_LENGTH:
            pos2 = rng.randrange(1, length2 - READ_LENGTH)
    else:
        chrom2, length2 = rng.choice(CHROMS)
        pos2 = rng.randrange(1, length2 - READ_LENGTH)
    return chrom1, pos1, rng.choice("+-"), chrom2, pos2, rng.choice("+-")


def sequence_pool(rng, size=512):
    """A pool of read sequences to draw from.

    Drawing 80 bases per read costs more than everything else in the generator
    put together, and nothing downstream looks at the bases -- pairtools reads
    the length, the flags and the CIGAR. So they are generated once and reused.
    """
    return ["".join(rng.choice("ACGT") for _ in range(READ_LENGTH))
            for _ in range(size)]


def records_for(rng, rid, molecule, kind, sequences):
    """The two SAM records of one read pair, as bwa mem -SP would emit them."""
    chrom1, pos1, strand1, chrom2, pos2, strand2 = molecule
    seq1 = rng.choice(sequences)
    seq2 = rng.choice(sequences)

    if kind == "unmapped":
        # 0x4d = paired, both unmapped, first in pair; 0x8d likewise for second
        return [sam_record(rid, 77, "*", 0, 0, seq1, []),
                sam_record(rid, 141, "*", 0, 0, seq2, [])]

    tags = ["NM:i:{}".format(rng.randint(0, 3)), "AS:i:{}".format(rng.randint(60, 80))]
    mapq1 = mapq2 = 60
    if kind == "multimapper":
        # a low MAPQ on one side is what makes the pair MU/NU rather than UU
        mapq1 = rng.choice([0, 1])
        tags = tags + ["XS:i:75"]

    flag1 = 65 if strand1 == "+" else 81
    flag2 = 129 if strand2 == "+" else 145
    return [sam_record(rid, flag1, chrom1, pos1, mapq1, seq1, tags),
            sam_record(rid, flag2, chrom2, pos2, mapq2, seq2, tags)]


def write_sam(path, n_pairs, dup_rate, seed):
    """Write `n_pairs` read pairs, of which `dup_rate` are PCR duplicates."""
    rng = random.Random(seed)
    sequences = sequence_pool(rng)
    written = 0
    with open(path, "w") as f:
        f.write("@HD\tVN:1.6\tSO:unknown\n")
        for chrom, length in CHROMS:
            f.write("@SQ\tSN:{}\tLN:{}\n".format(chrom, length))
        f.write("@PG\tID:bwa\tPN:bwa\tVN:0.7.17\tCL:bwa mem -SP genome r1.fq r2.fq\n")

        molecule = None
        lines = []
        while written < n_pairs:
            if molecule is None or rng.random() >= dup_rate:
                molecule = ligation(rng)
                kind = rng.choices(
                    ["mapped", "multimapper", "unmapped"], [0.90, 0.07, 0.03]
                )[0]
            else:
                # A PCR duplicate: same molecule under a new read ID, a few bp
                # off where the aligner put the last one. The wobble is applied
                # to the previous copy rather than to the original, so a deep
                # family walks: members are near-duplicates of their neighbours
                # without all being within --max-mismatch of the first. Those
                # chains are the interesting case for dedup, and they are what
                # a real family of duplicates with soft-clipped ends looks
                # like once the aligner has had its say.
                chrom1, pos1, strand1, chrom2, pos2, strand2 = molecule
                molecule = (chrom1, pos1 + rng.randint(0, 3), strand1,
                            chrom2, pos2 + rng.randint(0, 3), strand2)

            lines.extend(records_for(rng, read_id(rng), molecule, kind, sequences))
            written += 1
            if len(lines) >= 20000:
                f.write("\n".join(lines) + "\n")
                lines = []
        if lines:
            f.write("\n".join(lines) + "\n")
    return path


def write_chromsizes(path):
    with open(path, "w") as f:
        f.write("".join("{}\t{}\n".format(c, n) for c, n in CHROMS))
    return path


def write_fragments(path):
    """A BED of evenly spaced fragments, standing in for a digest."""
    with open(path, "w") as f:
        for chrom, length in CHROMS:
            for start in range(0, length, FRAGMENT_LENGTH):
                f.write("{}\t{}\t{}\n".format(
                    chrom, start, min(start + FRAGMENT_LENGTH, length)))
    return path


def run(*args):
    proc = subprocess.run([str(a) for a in args], capture_output=True)
    if proc.returncode != 0:
        raise SystemExit("failed: {}\n{}".format(
            " ".join(str(a) for a in args), proc.stderr.decode()))


def sam_to_bam(sam_path, bam_path):
    import pysam

    pysam.view("-b", "-o", str(bam_path), str(sam_path), catch_stdout=False)
    return bam_path


class Dataset(object):
    """Where the generated files live, and what settings produced them."""

    def __init__(self, directory, n_pairs, dup_rate, seed):
        self.directory = directory
        self.n_pairs = n_pairs
        self.dup_rate = dup_rate
        self.seed = seed

    def path(self, name):
        return os.path.join(self.directory, name)

    bam = property(lambda self: self.path("bench.bam"))
    pairs = property(lambda self: self.path("bench.pairs"))
    parquet = property(lambda self: self.path("bench.parquet"))
    chromsizes = property(lambda self: self.path("bench.chrom.sizes"))
    fragments = property(lambda self: self.path("bench.frags.bed"))
    manifest = property(lambda self: self.path("manifest.json"))

    @property
    def settings(self):
        return {"n_pairs": self.n_pairs, "dup_rate": self.dup_rate,
                "seed": self.seed, "read_length": READ_LENGTH,
                "chroms": len(CHROMS), "format": 2}

    def is_current(self):
        """True when the files on disk were built from these settings."""
        try:
            with open(self.manifest) as f:
                if json.load(f) != self.settings:
                    return False
        except (IOError, OSError, ValueError):
            return False
        return all(os.path.exists(p) for p in
                   [self.bam, self.pairs, self.parquet,
                    self.chromsizes, self.fragments])


def build(directory, n_pairs=1_000_000, dup_rate=0.15, seed=17, force=False,
          nproc=4, log=print):
    """Generate the dataset, or reuse what is already there.

    Regenerating is not cheap, so a manifest records what the files were built
    from and generation is skipped when it matches.
    """
    dataset = Dataset(directory, n_pairs, dup_rate, seed)
    os.makedirs(directory, exist_ok=True)
    if dataset.is_current() and not force:
        log("reusing dataset in {}".format(directory))
        return dataset

    log("generating {:,} read pairs ({:.0%} duplicates) in {}".format(
        n_pairs, dup_rate, directory))
    write_chromsizes(dataset.chromsizes)
    write_fragments(dataset.fragments)

    sam = dataset.path("bench.sam")
    write_sam(sam, n_pairs, dup_rate, seed)
    log("  -> bam")
    sam_to_bam(sam, dataset.bam)
    os.remove(sam)

    # Parse with our own tool rather than pairtools': it is the same parser,
    # and this way generating the fixture does not require pairtools to be
    # installed. Sorting is what `dedup` and `filterbycov` require of input.
    log("  -> pairs")
    unsorted_pairs = dataset.path("unsorted.parquet")
    run(sys.executable, "-m", "pairtools_parquet", "parse",
        "-c", dataset.chromsizes, "--drop-sam", "--add-columns", "mapq",
        "--assembly", "synthetic", "-o", unsorted_pairs, dataset.bam)
    run(sys.executable, "-m", "pairtools_parquet", "sort",
        "--nproc", nproc, "-o", dataset.parquet, unsorted_pairs)
    os.remove(unsorted_pairs)
    run(sys.executable, "-m", "pairtools_parquet", "parquet-to-csv",
        "-o", dataset.pairs, dataset.parquet)

    with open(dataset.manifest, "w") as f:
        json.dump(dataset.settings, f, indent=2, sort_keys=True)
    log("done")
    return dataset


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("directory", help="where to write the dataset")
    parser.add_argument("-n", "--n-pairs", type=int, default=1_000_000)
    parser.add_argument("--dup-rate", type=float, default=0.15,
                        help="share of read pairs that are PCR duplicates")
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--force", action="store_true",
                        help="regenerate even if the dataset is current")
    args = parser.parse_args(argv)
    build(args.directory, n_pairs=args.n_pairs, dup_rate=args.dup_rate,
          seed=args.seed, force=args.force)


if __name__ == "__main__":
    main()

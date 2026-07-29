# -*- coding: utf-8 -*-
"""Time every tool against `pairtools`, and check the answers still match.

    python benchmarks/run.py                    # generated data, every tool
    python benchmarks/run.py -t dedup,sort -r 3
    python benchmarks/run.py --pairs mine.pairs -c mine.chrom.sizes \\
                             -f mine.frags.bed --bam mine.bam

Three columns are timed for each tool, all running the same arguments through
the same option names:

    pairtools        text in, text out -- the baseline
    ours (text)      the same text in and out, through this package
    ours (parquet)   Parquet in, Parquet out -- what a Parquet pipeline pays

Unless `--no-check` is given the outputs are compared as well as timed, so a
speedup that came from doing less work fails the run rather than printing a
good number. Every timing measures a fresh subprocess, so interpreter and
import time are in the total, as they are for a real invocation.
"""

import argparse
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import make_data  # noqa: E402

OURS = [sys.executable, "-m", "pairtools_parquet"]
PAIRTOOLS = ["pairtools"]


class Context(object):
    """The inputs a benchmark can ask for, and the knobs it should respect."""

    def __init__(self, pairs, parquet, bam, chromsizes, fragments, nproc, workdir):
        self.pairs = pairs
        self.parquet = parquet
        self.bam = bam
        self.chromsizes = chromsizes
        self.fragments = fragments
        self.nproc = nproc
        self.workdir = workdir
        # Filled in only when a benchmark needs several inputs; see `merge`.
        self.halves_pairs = None
        self.halves_parquet = None


class Benchmark(object):
    """One tool, and how to invoke it on either engine.

    `args` is shared between the two binaries on purpose: `pairtools_parquet`
    mirrors `pairtools` option for option, so a benchmark that needed
    different arguments for the two would be measuring different work.
    """

    def __init__(self, name, tool=None, args=lambda ctx: [], source="pairs",
                 output="pairs", note="", differs_from_pairtools=None):
        self.name = name
        self.tool = tool or name
        self.args = args
        self.source = source
        self.output = output
        self.note = note
        # Set where our answer is known to differ from pairtools' and we mean
        # it -- the string says why. Everything else disagreeing is a bug, and
        # fails the run.
        self.differs_from_pairtools = differs_from_pairtools

    def inputs(self, ctx, engine):
        """The positional input paths, in the format this engine reads."""
        if self.source == "bam":
            return [ctx.bam]
        if self.source == "halves":
            return ctx.halves_parquet if engine == "parquet" else ctx.halves_pairs
        return [ctx.parquet if engine == "parquet" else ctx.pairs]

    def output_path(self, ctx, engine):
        # Tools whose output is not pairs -- stats, scaling -- write the same
        # format whichever engine produced it, so there is nothing to switch on.
        if self.output == "pairs":
            extension = "parquet" if engine == "parquet" else "pairs"
        else:
            extension = self.output
        return os.path.join(
            ctx.workdir, "{}.{}.{}".format(self.name, engine, extension))

    def command(self, ctx, engine):
        binary = PAIRTOOLS if engine == "reference" else OURS
        return ([str(a) for a in binary + [self.tool] + list(self.args(ctx))]
                + ["-o", self.output_path(ctx, engine)]
                + [str(p) for p in self.inputs(ctx, engine)])


BENCHMARKS = [
    Benchmark(
        "parse", source="bam",
        args=lambda ctx: ["-c", ctx.chromsizes, "--drop-sam",
                          "--add-columns", "mapq", "--assembly", "synthetic"],
        note="the parser is pairtools' own; what differs is where the rows go",
    ),
    Benchmark("sort", args=lambda ctx: ["--nproc", ctx.nproc]),
    Benchmark("merge", source="halves", args=lambda ctx: ["--nproc", ctx.nproc]),
    Benchmark("select", args=lambda ctx: ['(pair_type=="UU") and (chrom1==chrom2)']),
    Benchmark(
        "dedup", args=lambda ctx: ["--nproc-in", ctx.nproc, "--nproc-out", ctx.nproc],
        note="default --max-mismatch 3",
        differs_from_pairtools=(
            "a few rows, where duplicates chain across a window boundary. "
            "`pairtools dedup` carries only non-duplicate rows into the next "
            "chunk, so a chain A~B~C split over a boundary loses the B~C link "
            "when B was marked a duplicate of A. Our lookback holds every row "
            "and re-decides them, so the chain survives -- the answer is the "
            "one a single chunk would give. See UPSTREAM.md."
        ),
    ),
    Benchmark(
        "dedup-exact", tool="dedup",
        args=lambda ctx: ["--max-mismatch", 0, "--nproc-in", ctx.nproc,
                          "--nproc-out", ctx.nproc],
    ),
    Benchmark("filterbycov", args=lambda ctx: []),
    Benchmark("restrict", args=lambda ctx: ["-f", ctx.fragments]),
    Benchmark("flip", args=lambda ctx: ["-c", ctx.chromsizes]),
    Benchmark("markasdup"),
    Benchmark("sample", args=lambda ctx: ["--seed", 1, 0.1]),
    Benchmark("stats", output="stats"),
    Benchmark(
        "scaling", output="tsv",
        differs_from_pairtools=(
            "the region bounds and the `n_bp2` areas, on every row. "
            "`pairtools scaling` extracts the header's chromosome sizes and "
            "then discards them, so every region end stays at the -1 sentinel "
            "and P(s) is normalised by an area computed from a negative "
            "length. We use the sizes, which gives exactly what `pairtools "
            "scaling --view <those sizes>` gives. The `n_pairs` counts are "
            "identical. See UPSTREAM.md."
        ),
    ),
]

BY_NAME = {b.name: b for b in BENCHMARKS}


def run(command):
    proc = subprocess.run(command, capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError("failed: {}\n{}".format(
            " ".join(command), proc.stderr.decode()[-4000:]))


def time_command(command, repeat):
    """Best of `repeat` runs.

    The minimum, not the mean: every source of noise on a shared machine adds
    time, so the fastest run is the closest to the cost of the work itself.
    """
    best = float("inf")
    for _ in range(repeat):
        started = time.perf_counter()
        run(command)
        best = min(best, time.perf_counter() - started)
    return best


def startup_cost():
    """How long each binary takes to do nothing.

    Both are Python programs importing pandas, pyarrow and duckdb, so a run on
    a small file measures imports rather than the work. Reported so a table
    where every tool scores 1.0x can be recognised for what it is.
    """
    costs = {}
    for label, binary in [("pairtools", PAIRTOOLS), ("ours", OURS)]:
        costs[label] = time_command(binary + ["--help"], 3)
    return costs


def as_text(path, workdir, label):
    """The result as text, so results in different formats can be compared.

    Parquet is converted back, so what gets compared is the pairs themselves
    rather than the encoding.
    """
    if not path.endswith(".parquet"):
        return path
    converted = os.path.join(workdir, "converted.{}.pairs".format(label))
    run(OURS + ["parquet-to-csv", "-o", converted, path])
    return converted


def body_digest(path):
    """A hash of a result's rows.

    Headers are excluded: every tool appends its own @PG record naming itself,
    so those are expected to differ.
    """
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for line in f:
            if not line.startswith(b"#"):
                digest.update(line)
    return digest.hexdigest()


def count_differing_rows(left, right):
    """How many body rows the two files disagree on, or None if unknown.

    Only ever called once a hash has already shown they differ, so the cost is
    paid only where there is something to explain. `diff` does the work
    because the alternative -- holding a hash of every row of both files in
    memory -- does not scale to the file sizes this harness is pointed at.
    """
    bash = shutil.which("bash")
    if bash is None or shutil.which("diff") is None:
        return None
    script = "diff <(grep -v '^#' {}) <(grep -v '^#' {}) | grep -c '^[<>]'".format(
        shlex.quote(left), shlex.quote(right))
    try:
        proc = subprocess.run([bash, "-c", script], capture_output=True)
        return int(proc.stdout.decode().strip())
    except (OSError, ValueError):
        return None


def split_in_half(path, workdir, name):
    """Two sorted halves of a pairs file, for benchmarking `merge`."""
    first = os.path.join(workdir, "half1." + name)
    second = os.path.join(workdir, "half2." + name)
    run(OURS + ["select", "pos1 % 2 == 0", "-o", first, "--output-rest", second,
                path])
    return [first, second]


def compare(benchmark, ctx):
    """Check the three outputs against each other.

    Two separate questions, because they mean different things. Our own text
    and Parquet paths disagreeing is always a bug -- the format a file is in
    must not change the answer. Disagreeing with `pairtools` is a bug too,
    unless the benchmark declares the divergence and says why.
    """
    outputs = {engine: benchmark.output_path(ctx, engine)
               for engine in ["reference", "text", "parquet"]}
    texts = {engine: as_text(path, ctx.workdir, engine)
             for engine, path in outputs.items()}
    digests = {engine: body_digest(path) for engine, path in texts.items()}

    result = {"formats_agree": digests["text"] == digests["parquet"],
              "matches_pairtools": digests["reference"] == digests["parquet"]}
    if not result["matches_pairtools"]:
        result["differing_rows"] = count_differing_rows(
            texts["reference"], texts["parquet"])
        result["expected"] = benchmark.differs_from_pairtools is not None

    # The outputs themselves are the caller's to clean up; only the text
    # conversions made here are ours.
    for engine, path in texts.items():
        if path != outputs[engine]:
            os.remove(path)
    return result


def measure(benchmark, ctx, repeat, check):
    """Time the three variants of one tool and compare what they produced."""
    timings = {}
    for engine in ["reference", "text", "parquet"]:
        try:
            timings[engine] = time_command(benchmark.command(ctx, engine), repeat)
        except RuntimeError as error:
            return {"name": benchmark.name, "error": str(error).split("\n")[0],
                    "note": benchmark.note}

    result = {"name": benchmark.name, "note": benchmark.note,
              "reference": timings["reference"], "text": timings["text"],
              "parquet": timings["parquet"],
              "speedup": timings["reference"] / timings["parquet"]}
    if check:
        result.update(compare(benchmark, ctx))
    for engine in ["reference", "text", "parquet"]:
        path = benchmark.output_path(ctx, engine)
        if os.path.exists(path):
            os.remove(path)
    return result


def verdict(result):
    """The 'match' cell: what the comparison found, in a few characters."""
    if result.get("formats_agree") is False:
        return "**FORMATS DISAGREE**"
    if result.get("matches_pairtools"):
        return "yes"
    rows = result.get("differing_rows")
    if rows is None:
        counted = "differs"
    else:
        counted = "{} row{}".format(rows, "" if rows == 1 else "s")
    return counted + (" *" if result.get("expected") else " **?**")


def is_failure(result):
    """Whether this result should fail the run."""
    if "error" in result:
        return True
    if result.get("formats_agree") is False:
        return True
    return result.get("matches_pairtools") is False and not result.get("expected")


def format_table(results, check):
    lines = ["| tool | `pairtools` | ours (text) | ours (parquet) | speedup |"
             + (" match |" if check else ""),
             "|---|---|---|---|---|" + ("---|" if check else "")]
    for r in results:
        if "error" in r:
            lines.append("| `{}` | {} |".format(r["name"], r["error"]))
            continue
        row = "| `{}` | {:.1f}s | {:.1f}s | **{:.1f}s** | {:.1f}x |".format(
            r["name"], r["reference"], r["text"], r["parquet"], r["speedup"])
        if check:
            row += " {} |".format(verdict(r))
        lines.append(row)
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    parser.add_argument("--pairs", help="benchmark this .pairs/.pairs.gz file "
                                        "instead of generating one")
    parser.add_argument("--parquet", help="the same pairs as Parquet; "
                                          "converted from --pairs if omitted")
    parser.add_argument("--bam", help="benchmark `parse` on this .bam/.sam")
    parser.add_argument("-c", "--chroms-path",
                        help="chromosome sizes, required with --pairs/--bam")
    parser.add_argument("-f", "--frags", help="restriction fragment BED, "
                                              "required with --pairs/--bam")
    parser.add_argument("-d", "--data-dir", default="benchmarks/data",
                        help="where the generated dataset is cached "
                             "[default: %(default)s]")
    parser.add_argument("-n", "--n-pairs", type=int, default=1_000_000,
                        help="size of the generated dataset [default: %(default)s]")
    parser.add_argument("--dup-rate", type=float, default=0.15,
                        help="duplicate share of the generated dataset "
                             "[default: %(default)s]")
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--regenerate", action="store_true")
    parser.add_argument("-t", "--tools", default="",
                        help="comma-separated subset of {}".format(
                            ",".join(BY_NAME)))
    parser.add_argument("-r", "--repeat", type=int, default=1,
                        help="runs per measurement, best taken [default: %(default)s]")
    parser.add_argument("-p", "--nproc", type=int, default=4)
    parser.add_argument("--no-check", dest="check", action="store_false",
                        help="time only, do not compare the outputs")
    parser.add_argument("--json", help="also write the results here")
    args = parser.parse_args(argv)

    if args.pairs or args.bam:
        if not (args.chroms_path and args.frags):
            parser.error("--chroms-path and --frags are required with "
                         "--pairs/--bam")
        if args.parquet and not args.pairs:
            parser.error("--parquet needs --pairs alongside it: the `pairtools` "
                         "column has to read the same pairs as text")
        if args.pairs and args.pairs.endswith(".parquet"):
            parser.error("--pairs has to be text: the `pairtools` column needs "
                         "a file pairtools can read. Pass the Parquet copy as "
                         "--parquet, or convert with `pairtools_parquet "
                         "parquet-to-csv`.")
        pairs, parquet, bam = args.pairs, args.parquet, args.bam
        chromsizes, fragments = args.chroms_path, args.frags
    else:
        dataset = make_data.build(args.data_dir, n_pairs=args.n_pairs,
                                  dup_rate=args.dup_rate, seed=args.seed,
                                  force=args.regenerate, nproc=args.nproc)
        pairs, parquet, bam = dataset.pairs, dataset.parquet, dataset.bam
        chromsizes, fragments = dataset.chromsizes, dataset.fragments

    names = [n.strip() for n in args.tools.split(",") if n.strip()] or list(BY_NAME)
    unknown = [n for n in names if n not in BY_NAME]
    if unknown:
        parser.error("unknown tools: {}".format(", ".join(unknown)))
    selected = [BY_NAME[n] for n in names]
    if not bam:
        selected = [b for b in selected if b.source != "bam"]
    if not pairs:
        selected = [b for b in selected if b.source == "bam"]
    if not selected:
        parser.error("nothing left to run: of {}, {} needs a --bam and the "
                     "rest need --pairs".format(", ".join(names), "parse"))

    workdir = tempfile.mkdtemp(prefix="pairs-bench-")
    try:
        if pairs and not parquet:
            parquet = os.path.join(workdir, "input.parquet")
            print("converting {} to Parquet ...".format(pairs), file=sys.stderr,
                  flush=True)
            run(OURS + ["csv-to-parquet", "-o", parquet, pairs])

        ctx = Context(pairs, parquet, bam, chromsizes, fragments, args.nproc,
                      workdir)
        if any(b.source == "halves" for b in selected):
            ctx.halves_pairs = split_in_half(pairs, workdir, "pairs")
            ctx.halves_parquet = split_in_half(parquet, workdir, "parquet")

        results = []
        for benchmark in selected:
            print("running {} ...".format(benchmark.name), file=sys.stderr,
                  flush=True)
            results.append(measure(benchmark, ctx, args.repeat, args.check))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    startup = startup_cost()
    results.sort(key=lambda r: -r.get("speedup", 0))
    print()
    print(format_table(results, args.check))
    print()
    print("interpreter startup: {:.1f}s for pairtools, {:.1f}s for ours".format(
        startup["pairtools"], startup["ours"]))
    floor = 5 * max(startup.values())
    dominated = [r["name"] for r in results if r.get("reference", floor) < floor]
    if dominated:
        print("WARNING: {} ran in less than 5x startup, so their rows measure "
              "imports more than work. Re-run those with a larger -n.".format(
                  ", ".join(dominated)))
    for r in results:
        if r.get("note"):
            print("`{}`: {}".format(r["name"], r["note"]))
    for benchmark in selected:
        if benchmark.differs_from_pairtools and not reported_match(
                results, benchmark.name):
            print("`{}` *: {}".format(benchmark.name,
                                      benchmark.differs_from_pairtools))
    if args.json:
        with open(args.json, "w") as f:
            json.dump(results, f, indent=2)

    failed = [r["name"] for r in results if is_failure(r)]
    if failed:
        print("\nFAILED: {}".format(", ".join(failed)))
    return 1 if failed else 0


def reported_match(results, name):
    for r in results:
        if r["name"] == name:
            return r.get("matches_pairtools", True)
    return True


if __name__ == "__main__":
    sys.exit(main())

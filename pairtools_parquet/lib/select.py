"""Filtering pairs by a `pairtools select` CONDITION.

The condition language is pairtools', evaluated by pairtools: this calls
``pairtools.lib.select.evaluate_df`` per batch rather than reimplementing the
expression semantics. That matters because the language is arbitrary Python
with a handful of helpers (``csv_match``, ``wildcard_match``, ``regex_match``,
``region_match``) in scope, and any translation of it into SQL is an
approximation -- the one this replaces silently mishandled ``not``, chained
comparisons, ``--startup-code`` and anything else it had no regex for.

``evaluate_df`` evaluates the condition once per row, through
``DataFrame.iterrows``, which costs about 36us a row -- 200s for 5.6M pairs,
where `pairtools select` itself takes 9s. So a condition is first rewritten to
run over whole columns at once: `and`/`or`/`not` and chained comparisons become
their array equivalents on the parse tree, and the helper functions get
column-wise versions that evaluate the predicate once per *distinct* value,
since a chromosome column has a handful of them and millions of rows.

Anything the rewrite cannot express falls back to ``evaluate_df``, which
remains the definition of what a condition means. The fast path never guesses:
it returns nothing unless it produced a boolean array of the right length, and
`tests/test_select.py` checks the two agree condition by condition.
"""

import ast
import warnings

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pairtools.lib.select as pairtools_select
from pairtools.lib import headerops, pairsam_format
from pairtools.lib.select import evaluate_df

from .arrowio import PairsWriter, open_pairs

UTIL_NAME = "pairtools_parquet_select"


def _backport_region_match():
    """Make `region_match` available on pairtools releases that lack it.

    CONDITION is evaluated inside pairtools' own module namespace, so the
    helper functions have to live there. `region_match` was added upstream in
    open2c/pairtools#278 and is not in 1.1.2, the oldest release we support,
    while this package has offered it since 0.2.0 -- so on an older pairtools
    we install upstream's implementation rather than drop the feature.

    Delete this once the declared pairtools minimum contains it.
    """
    if hasattr(pairtools_select, "region_match"):
        return

    def region_match(chrom, pos, region_chrom, region_start=-1, region_end=-1):
        if region_end == -1:
            region_end = np.inf
        return chrom == region_chrom and region_start <= pos <= region_end

    pairtools_select.region_match = region_match


_backport_region_match()


def read_chrom_subset(path):
    """Read chromosome names from the first column of a chromsizes file."""
    with open(path, "r") as f:
        return [line.strip().split("\t")[0] for line in f if line.strip()]


def header_update(header, util_name=UTIL_NAME, remove_columns="", chrom_subset=None):
    """Apply a tool's provenance, column removal and chromosome subset to a header."""
    new_header = headerops.append_new_pg(header, ID=util_name, PN=util_name)

    if remove_columns:
        input_columns = headerops.extract_column_names(header)
        removed = remove_columns.split(",")
        for col in removed:
            if col in pairsam_format.COLUMNS_PAIRS:
                warnings.warn(
                    "Removing required {} column for .pairs format. "
                    "Output is not .pairs anymore".format(col)
                )
            elif col in pairsam_format.COLUMNS_PAIRSAM:
                warnings.warn(
                    "Removing required {} column for .pairsam format. "
                    "Output is not .pairsam anymore".format(col)
                )
        updated_columns = [x for x in input_columns if x not in removed]

        if len(updated_columns) == len(input_columns):
            warnings.warn(
                "Some column(s) {} not in the file, the operation has no "
                "effect".format(",".join(removed))
            )
        elif not updated_columns:
            raise ValueError("--remove-columns removed every column")
        else:
            new_header = headerops.set_columns(new_header, updated_columns)

    if chrom_subset is not None:
        new_header = headerops.subset_chroms_in_pairsheader(
            new_header, read_chrom_subset(chrom_subset)
        )

    return new_header


_BOOL_OPS = {ast.And: ast.BitAnd, ast.Or: ast.BitOr}


class _ArrayRewriter(ast.NodeTransformer):
    """Rewrite a CONDITION so it evaluates over whole columns at once.

    Python's `and`, `or` and `not` force their operands to a truth value, which
    an array does not have, and a chained comparison is `and` underneath. The
    array equivalents are `&`, `|` and `~`.

    The rewrite is done on the parse tree rather than on the text because `&`
    binds tighter than a comparison does: `a == 1 & b == 2` parses as
    `a == (1 & b) == 2`, so substituting the characters would silently change
    what the condition means.
    """

    def visit_BoolOp(self, node):
        self.generic_visit(node)
        operator = _BOOL_OPS[type(node.op)]()
        folded = node.values[0]
        for value in node.values[1:]:
            folded = ast.BinOp(left=folded, op=operator, right=value)
        return ast.copy_location(folded, node)

    def visit_UnaryOp(self, node):
        self.generic_visit(node)
        if isinstance(node.op, ast.Not):
            return ast.copy_location(
                ast.UnaryOp(op=ast.Invert(), operand=node.operand), node
            )
        return node

    def visit_Compare(self, node):
        self.generic_visit(node)
        if len(node.ops) == 1 and not isinstance(node.ops[0], (ast.In, ast.NotIn)):
            return node

        operands = [node.left] + node.comparators
        parts = []
        for i, op in enumerate(node.ops):
            left, right = operands[i], operands[i + 1]
            if isinstance(op, (ast.In, ast.NotIn)):
                # `x in [...]` is a membership test per element
                part = ast.Call(
                    func=ast.Name(id="_isin", ctx=ast.Load()),
                    args=[left, right],
                    keywords=[],
                )
                if isinstance(op, ast.NotIn):
                    part = ast.UnaryOp(op=ast.Invert(), operand=part)
            else:
                part = ast.Compare(left=left, ops=[op], comparators=[right])
            parts.append(part)

        folded = parts[0]
        for part in parts[1:]:
            folded = ast.BinOp(left=folded, op=ast.BitAnd(), right=part)
        return ast.copy_location(ast.fix_missing_locations(folded), node)


def _by_distinct_value(values, predicate):
    """Apply a scalar `predicate` once per distinct value of a column.

    A chromosome column has a handful of distinct values and millions of rows,
    so the wildcard and regex helpers cost a few calls rather than a few
    million.
    """
    codes, uniques = pd.factorize(values)
    answers = np.array([bool(predicate(value)) for value in uniques], dtype=bool)
    return answers[codes]


def _vectorized_helpers():
    """Column-wise versions of the helpers a CONDITION can call.

    Each falls back to pairtools' own scalar implementation when handed a
    scalar, so a condition mixing the two still means the same thing.
    """

    def lift(scalar_version):
        def helper(x, pattern):
            if not isinstance(x, np.ndarray):
                return scalar_version(x, pattern)
            return _by_distinct_value(x, lambda value: scalar_version(value, pattern))

        return helper

    def region_match(chrom, pos, region_chrom, region_start=-1, region_end=-1):
        if not isinstance(chrom, np.ndarray):
            return pairtools_select.region_match(
                chrom, pos, region_chrom, region_start, region_end
            )
        within = chrom == region_chrom
        within &= pos >= region_start
        if region_end != -1:
            within &= pos <= region_end
        return within

    return {
        "csv_match": lambda x, csv: (
            np.isin(x, csv.split(","))
            if isinstance(x, np.ndarray)
            else pairtools_select.csv_match(x, csv)
        ),
        "wildcard_match": lift(pairtools_select.wildcard_match),
        "regex_match": lift(pairtools_select.regex_match),
        "region_match": region_match,
        "_isin": lambda x, values: np.isin(x, list(values)),
    }


def compile_vectorized(condition):
    """Compile `condition` for whole-column evaluation, or None if unparseable."""
    try:
        tree = ast.parse(condition.strip(), mode="eval")
    except SyntaxError:
        return None
    tree = ast.fix_missing_locations(_ArrayRewriter().visit(tree))
    try:
        return compile(tree, "<condition>", "eval")
    except (SyntaxError, ValueError):
        return None


def evaluate_vectorized(columns, code, n_rows):
    """Evaluate a compiled condition over whole columns, or None if it cannot be.

    Returns None rather than raising for anything the rewrite does not cover --
    a helper that only works on one value at a time, an expression that
    collapses to a scalar -- so the caller can fall back to pairtools' row-wise
    evaluation and get the same answer, only slower.
    """
    namespace = dict(vars(pairtools_select))
    namespace.update(_vectorized_helpers())
    namespace.update(columns)
    try:
        result = eval(code, namespace)
    except Exception:
        return None

    if not isinstance(result, np.ndarray) or result.shape != (n_rows,):
        return None
    if result.dtype != bool:
        return None
    return result


def evaluate_batch(batch, condition, type_cast=(), startup_code=None, code=None):
    """Return a boolean mask of the rows of `batch` satisfying `condition`.

    Tries whole-column evaluation first and falls back to pairtools'
    ``evaluate_df``, which evaluates the condition once per row. The fallback
    is the definition of what the condition means; the fast path has to agree
    with it, which `tests/test_select.py` checks condition by condition.
    """
    if startup_code is not None:
        exec(startup_code, vars(pairtools_select))

    if code is not None:
        columns = {
            name: batch.column(name).to_numpy(zero_copy_only=False)
            for name in batch.schema.names
        }
        mask = evaluate_vectorized(columns, code, batch.num_rows)
        if mask is not None:
            return mask

    # engine="python" is the one that matches `pairtools select`: the pandas
    # engine goes through DataFrame.eval, which cannot call the helper
    # functions and rejects `and`/`or`.
    mask = evaluate_df(
        batch.to_pandas(), condition, type_cast, startup_code, engine="python"
    )
    return np.asarray(mask, dtype=bool)


def select_pairs(
    input_path,
    output,
    condition,
    output_rest=None,
    chrom_subset=None,
    startup_code=None,
    type_cast=(),
    remove_columns="",
    compress_program="auto",
    row_group_size=None,
    util_name=UTIL_NAME,
    **kwargs
):
    """Write the pairs of `input_path` satisfying `condition` to `output`.

    Parameters
    ----------
    input_path, output : str
        Paths; the format of each is taken from its extension.
    condition : str
        A pairtools CONDITION expression.
    output_rest : str, optional
        Where to write the pairs that did not satisfy the condition.
    chrom_subset : str, optional
        Path to a chromsizes file; additionally requires both sides of a pair
        to be on those chromosomes, and subsets the header accordingly.
    startup_code : str, optional
        Executed before filtering, to define helpers used in `condition`.
    type_cast : tuple of (str, str)
        Extra column-to-type casts for the evaluation.
    remove_columns : str
        Comma-separated columns to drop from the output.
    """
    header, reader = open_pairs(
        input_path,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )

    new_header = header_update(header, util_name, remove_columns, chrom_subset)
    keep_columns = headerops.extract_column_names(new_header)
    # As an Arrow value set rather than a Python set: `pc.is_in` does the
    # membership test in C++, where `np.isin` over a column materialized as
    # Python strings costs seconds per batch on a large file.
    chroms = (
        pa.array(sorted(set(read_chrom_subset(chrom_subset))), type=pa.string())
        if chrom_subset
        else None
    )

    writer_kwargs = dict(
        compress_program=compress_program,
        row_group_size=row_group_size,
        nproc_out=kwargs.get("nproc_out", 8),
    )
    selected_schema = pa.schema([reader.schema.field(c) for c in keep_columns])
    # --type-cast changes what a column *is* for the evaluation, and NumPy
    # compares mismatched types elementwise-false rather than raising, so the
    # fast path could differ without ever failing. Rare enough to just skip.
    code = None if type_cast else compile_vectorized(condition)

    with PairsWriter(
        output, new_header, schema=selected_schema, **writer_kwargs
    ) as selected_writer:
        rest_writer = (
            PairsWriter(output_rest, new_header, schema=selected_schema,
                        **writer_kwargs)
            if output_rest
            else None
        )
        try:
            for batch in reader:
                mask = evaluate_batch(batch, condition, type_cast, startup_code, code)
                if chroms is not None:
                    for side in ("chrom1", "chrom2"):
                        in_subset = pc.is_in(
                            batch.column(side).cast(pa.string()), value_set=chroms
                        )
                        mask &= in_subset.to_numpy(zero_copy_only=False).astype(bool)

                projected = batch.select(keep_columns)
                selected_writer.write(projected.filter(pa.array(mask)))
                if rest_writer is not None:
                    rest_writer.write(projected.filter(pa.array(~mask)))
        finally:
            if rest_writer is not None:
                rest_writer.close()

"""Assigning pairs to restriction fragments.

`pairtools restrict` calls ``find_rfrag`` once per side per row, each call doing
its own ``searchsorted``. Here the rows of a batch are grouped by chromosome and
looked up in one vectorized ``searchsorted`` per chromosome, which is the same
computation with the per-row overhead removed.

The fragment boundaries are built exactly as upstream builds them: per
chromosome, ``[0] + [end + 1 for each fragment]``, sorted.
"""

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pacsv
from pairtools.lib import headerops, pairsam_format

from .._logging import get_logger
from .arrowio import PairsWriter, open_pairs

UTIL_NAME = "pairtools_parquet_restrict"

logger = get_logger()

RFRAG_COLUMNS = [
    "rfrag1",
    "rfrag_start1",
    "rfrag_end1",
    "rfrag2",
    "rfrag_start2",
    "rfrag_end2",
]


def load_rfrags(frags_path):
    """Read a restriction-fragment BED into per-chromosome boundary arrays.

    Produces exactly what upstream's ``np.genfromtxt`` route produces -- sorted
    by (chrom, start, end), then ``[0] + ends + 1`` per chromosome -- but reads
    the file with pyarrow. On a genome-scale fragment file (~1M fragments)
    genfromtxt takes about 5s and dominates the whole command; this is roughly
    two orders of magnitude faster and is why `restrict` is faster here at all.
    """
    table = pacsv.read_csv(
        frags_path,
        read_options=pacsv.ReadOptions(autogenerate_column_names=True),
        parse_options=pacsv.ParseOptions(delimiter="\t"),
    )
    # BED files may carry extra columns; only the first three are used.
    table = table.select(table.column_names[:3]).rename_columns(
        ["chrom", "start", "end"]
    )
    # `comments="#"` upstream: drop commented and track lines.
    keep = pa.compute.invert(
        pa.compute.starts_with(table.column("chrom").cast(pa.string()), "#")
    )
    table = table.filter(keep)
    table = table.sort_by(
        [("chrom", "ascending"), ("start", "ascending"), ("end", "ascending")]
    )

    if table.num_rows == 0:
        return {}

    # The BED can be millions of rows, so the chromosome column is compared as
    # integer codes rather than as Python strings -- the table is already
    # sorted, so run boundaries are just where the code changes.
    codes, names = chrom_codes(table.column("chrom"))
    ends = table.column("end").to_numpy(zero_copy_only=False)

    borders = np.r_[0, 1 + np.where(codes[:-1] != codes[1:])[0], len(codes)]
    return {
        names[codes[i]]: np.concatenate([[0], ends[i:j] + 1])
        for i, j in zip(borders[:-1], borders[1:])
    }


def chrom_codes(column):
    """Dictionary-encode a chromosome column into (integer codes, names).

    Everything downstream then works on small integers: the distinct values are
    the dictionary, so there is no ``np.unique`` over 5.6M strings, and
    selecting one chromosome's rows is an integer comparison rather than an
    elementwise comparison of Python string objects.
    """
    encoded = pc.dictionary_encode(column)
    if isinstance(encoded, pa.ChunkedArray):
        encoded = encoded.combine_chunks()
    codes = (
        pc.fill_null(encoded.indices, -1)
        .to_numpy(zero_copy_only=False)
        .astype(np.int64)
    )
    return codes, encoded.dictionary.to_pylist()


def annotate_side(column, positions, rfrags, warned):
    """Return (index, start, end) arrays for one side of every pair in a batch.

    Unmapped sides, and sides on a chromosome with no annotated fragments, get
    the unannotated sentinels -- as `find_rfrag` is written to do.
    """
    codes, names = chrom_codes(column)
    n = len(codes)
    index = np.full(n, pairsam_format.UNANNOTATED_RFRAG, dtype=np.int64)
    start = np.full(n, pairsam_format.UNMAPPED_POS, dtype=np.int64)
    end = np.full(n, pairsam_format.UNMAPPED_POS, dtype=np.int64)

    for code, chrom in enumerate(names):
        if chrom == pairsam_format.UNMAPPED_CHROM:
            continue
        if chrom not in rfrags:
            if chrom not in warned:
                logger.warning(
                    "Chromosome %s does not have annotated restriction "
                    "fragments, return empty.",
                    chrom,
                )
                warned.add(chrom)
            continue

        sites = rfrags[chrom]
        mask = codes == code
        idx = np.minimum(
            np.maximum(0, np.searchsorted(sites, positions[mask], "right") - 1),
            len(sites) - 2,
        )
        index[mask] = idx
        start[mask] = sites[idx]
        end[mask] = sites[idx + 1]

    return index, start, end


def restrict_pairs(
    input_path,
    output,
    frags,
    compress_program="auto",
    row_group_size=None,
    util_name=UTIL_NAME,
    **kwargs
):
    """Annotate every pair of `input_path` with its restriction fragments."""
    rfrags = load_rfrags(frags)

    header, reader = open_pairs(
        input_path,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )
    new_header = headerops.append_new_pg(
        list(header), ID=util_name, PN=util_name
    )
    new_header = headerops.append_columns(new_header, RFRAG_COLUMNS)

    warned = set()

    with PairsWriter(
        output,
        new_header,
        compress_program=compress_program,
        row_group_size=row_group_size,
        nproc_out=kwargs.get("nproc_out", 8),
    ) as writer:
        for batch in reader:
            table = pa.Table.from_batches([batch])
            annotations = []
            for chrom_col, pos_col in (("chrom1", "pos1"), ("chrom2", "pos2")):
                positions = table.column(pos_col).to_numpy(zero_copy_only=False)
                annotations.extend(
                    annotate_side(
                        table.column(chrom_col), positions, rfrags, warned
                    )
                )

            writer.write(
                pa.Table.from_arrays(
                    list(table.columns) + [pa.array(a) for a in annotations],
                    names=table.column_names + RFRAG_COLUMNS,
                )
            )

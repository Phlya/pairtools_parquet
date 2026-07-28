"""Extracting pairs from aligned SAM/BAM data.

The parser is pairtools' ``streaming_classify`` -- 1500 lines of walk rescue,
orientation and position reporting, SAM tag extraction -- called unchanged.
The only thing that changes is where its output goes.

``streaming_classify`` emits one formatted, tab-separated row at a time through
``out_file.write(line)``. Rather than copy ``write_pairsam``, which knows how to
pack SAM records, drop sequences and append tag columns, this module hands it a
sink that keeps the lines and converts them in bulk. The conversion is the same
``pyarrow.csv`` reader every .pairs file goes through, so ``parse -o x.parquet``
gives exactly what ``parse -o x.pairs`` followed by ``csv-to-parquet`` would.

What this buys is the end of the pipeline's first serialization round trip: the
parsed pairs go to Parquet without ever being written as text and read back.
"""

import sys

import pyarrow as pa
import pyarrow.csv as pacsv
from pairtools.lib import fileio, headerops, pairsam_format
from pairtools.lib.parse import streaming_classify
from pairtools.lib.stats import PairCounter

from .arrowio import PairsWriter, csv_convert_options, csv_parse_options
from .schema import schema_from_columns

UTIL_NAME = "pairtools_parquet_parse"
UTIL_NAME_PARSE2 = "pairtools_parquet_parse2"

#: Rows buffered before a batch is converted and written.
DEFAULT_BATCH_SIZE = 1 << 16

#: The header line `--output-parsed-alignments` starts with.
ALIGNMENTS_HEADER = (
    "readID\tside\tchrom\tpos\tstrand\tmapq\tcigar\tdist_5_lo\tdist_5_hi"
    "\tmatched_bp\n"
)


class PairsBatchSink:
    """A write-only file object that turns .pairsam rows into Arrow batches.

    ``streaming_classify`` writes rows with ``out_file.write(line)``; this
    collects them and, every `batch_size` rows, parses the accumulated block
    with the same reader that reads a .pairs file from disk and hands the table
    to a :class:`PairsWriter`.

    Going through the CSV reader rather than splitting each line in Python is
    both faster and, more usefully, identical: whatever a value would become on
    the way in from text, it becomes here.
    """

    def __init__(self, writer, columns, batch_size=DEFAULT_BATCH_SIZE):
        self.writer = writer
        self.columns = list(columns)
        self.batch_size = batch_size
        self._schema = schema_from_columns(self.columns)
        self._lines = []

    def write(self, line):
        self._lines.append(line.encode("utf-8") if isinstance(line, str) else line)
        if len(self._lines) >= self.batch_size:
            self.flush()

    def writelines(self, lines):
        for line in lines:
            self.write(line)

    def flush(self):
        if not self._lines:
            return
        block = b"".join(self._lines)
        self._lines = []
        self.writer.write(
            pacsv.read_csv(
                pa.BufferReader(block),
                read_options=pacsv.ReadOptions(column_names=self.columns),
                parse_options=csv_parse_options(),
                convert_options=csv_convert_options(self._schema),
            )
        )

    def close(self):
        self.flush()


def output_columns(add_columns="", drop_sam=True, add_pair_index=False):
    """The columns `parse` writes, in order.

    Restates the column bookkeeping of ``pairtools/cli/parse.py:parse_py``,
    which is inline there rather than a function; UPSTREAM.md records that.
    """
    add_columns = [col for col in str(add_columns).split(",") if col]
    for col in add_columns:
        if not (
            col in pairsam_format.EXTRA_COLUMNS or (len(col) == 2 and col.isupper())
        ):
            raise ValueError("{} is not a valid extra column".format(col))

    columns = pairsam_format.COLUMNS + [
        col + side for col in add_columns for side in ["1", "2"]
    ]

    if drop_sam:
        columns.pop(columns.index("sam1"))
        columns.pop(columns.index("sam2"))
    if not add_pair_index:
        columns.pop(columns.index("walk_pair_index"))
        columns.pop(columns.index("walk_pair_type"))

    return columns


def parse_pairs(
    sam_path,
    chroms_path,
    output,
    output_parsed_alignments="",
    output_stats="",
    parse2=False,
    compress_program="auto",
    row_group_size=None,
    batch_size=DEFAULT_BATCH_SIZE,
    util_name=None,
    **kwargs
):
    """Parse `sam_path` into pairs, writing them to `output`.

    Parameters
    ----------
    sam_path : str
        A .sam/.bam file, or empty for stdin.
    chroms_path : str
        A chrom.sizes file giving the chromosome order.
    output : str
        Where to write the pairs; the extension picks .pairs, .pairs.gz or
        .parquet.
    parse2 : bool
        Use the complex-walk parser, i.e. ``pairtools parse2``.

    Every other keyword goes straight to ``streaming_classify``.
    """
    if util_name is None:
        util_name = UTIL_NAME_PARSE2 if parse2 else UTIL_NAME
    if not output:
        raise ValueError(
            "An output path is required: Parquet cannot be written to stdout"
        )

    # Imported here so pysam is only needed by the commands that read SAM.
    from pairtools.lib.parse_pysam import AlignmentFilePairtoolized

    input_sam = AlignmentFilePairtoolized(
        sam_path if sam_path else "-", "r", threads=kwargs.get("nproc_in", 3)
    )

    samheader = input_sam.header
    if not samheader:
        raise ValueError(
            "The input sam is missing a header! If reading a bam file, please "
            "use `samtools view -h` to include the header."
        )

    columns = output_columns(
        add_columns=kwargs.get("add_columns", ""),
        drop_sam=kwargs.get("drop_sam", True),
        add_pair_index=kwargs.get("add_pair_index", False),
    )

    sam_chromsizes = headerops.get_chromsizes_from_pysam_header(samheader)
    chromosomes = headerops.get_chrom_order(chroms_path, list(sam_chromsizes.keys()))

    header = headerops.make_standard_pairsheader(
        # `--assembly` defaults to None, which make_standard_pairsheader turns
        # into "unknown" -- passing "" instead would write an empty value.
        assembly=kwargs.get("assembly", ""),
        chromsizes=[(chrom, sam_chromsizes[chrom]) for chrom in chromosomes],
        columns=columns,
        shape="whole matrix" if not kwargs["flip"] else "upper triangle",
    )
    header = headerops.insert_samheader_pysam(header, samheader)
    header = headerops.append_new_pg(header, ID=util_name, PN=util_name)

    out_alignments_stream = _open_side_output(
        output_parsed_alignments, kwargs
    )
    if out_alignments_stream:
        out_alignments_stream.write(ALIGNMENTS_HEADER)

    out_stats_stream = _open_side_output(output_stats, kwargs)
    out_stat = PairCounter() if output_stats else None

    with PairsWriter(
        output,
        header,
        schema=schema_from_columns(columns),
        compress_program=compress_program,
        row_group_size=row_group_size,
        nproc_out=kwargs.get("nproc_out", 8),
    ) as writer:
        sink = PairsBatchSink(writer, columns, batch_size=batch_size)
        try:
            streaming_classify(
                input_sam,
                sink,
                chromosomes,
                out_alignments_stream,
                out_stat,
                parse2=parse2,
                **kwargs
            )
            sink.flush()
        finally:
            if out_alignments_stream and out_alignments_stream is not sys.stdout:
                out_alignments_stream.close()

    if out_stat is not None:
        out_stat.save(out_stats_stream)
    if out_stats_stream and out_stats_stream is not sys.stdout:
        out_stats_stream.close()


def _open_side_output(path, kwargs):
    """Open one of the plain-text side outputs, or return None."""
    if not path:
        return None
    return fileio.auto_open(
        path,
        mode="w",
        nproc=kwargs.get("nproc_out", 8),
        command=kwargs.get("cmd_out", None),
    )

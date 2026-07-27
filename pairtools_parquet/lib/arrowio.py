"""Format-agnostic streaming I/O for .pairs data as Arrow record batches.

Every tool reads through :func:`open_pairs` and writes through
:class:`PairsWriter`, so no tool needs to know whether it is looking at a
``.pairs`` text file, a ``.pairs.gz``, or a ``.parquet``. Both directions are
streaming: a file is never held in memory in full.

The module deliberately depends only on ``pyarrow`` and ``pairtools``, not on
anything else in this package, so it can move into pairtools unchanged.

Known limitation: pyarrow's CSV writer refuses to emit a value containing a
double quote when quoting is disabled, and .pairs is an unquoted format, so a
field with a ``"`` in it can be read but not written back to text. Such a value
round-trips through Parquet without complaint, and pairtools itself writes it
happily; only our text output rejects it, loudly rather than silently.
"""

import contextlib
import subprocess
import sys

import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq
from pairtools.lib import fileio, headerops

from . import headerio
from .schema import schema_from_columns

#: Rows per batch when reading Parquet.
DEFAULT_BATCH_SIZE = 1 << 17

#: Bytes per block when reading text; pyarrow sizes text batches by bytes.
DEFAULT_BLOCK_SIZE = 1 << 24

PARQUET_SUFFIXES = (".parquet", ".pq")


def is_parquet(path):
    """Whether `path` names a Parquet file, by extension."""
    return str(path).lower().endswith(PARQUET_SUFFIXES)


def read_header(path, nproc_in=3, cmd_in=None):
    """Read just the .pairs header of a file, in either format."""
    if is_parquet(path):
        return headerio.read_header(path)

    instream = fileio.auto_open(path, mode="r", nproc=nproc_in, command=cmd_in)
    try:
        header, _ = headerops.get_header(instream)
    finally:
        if instream is not sys.stdin:
            instream.close()
    return header


def _binary_stream(stream):
    """Return the underlying binary stream of a possibly-text stream.

    ``headerops.get_header`` consumes header lines through ``stream.buffer``
    when the stream is text, so the binary layer -- not the text wrapper -- is
    the one positioned at the first data line.
    """
    return stream.buffer if hasattr(stream, "buffer") else stream


def _csv_read_options(columns, block_size):
    return pacsv.ReadOptions(column_names=list(columns), block_size=block_size)


def _csv_parse_options():
    # .pairs is strictly tab-separated with no quoting: a readID containing a
    # double quote is data, not a quoted field.
    return pacsv.ParseOptions(
        delimiter="\t",
        quote_char=False,
        double_quote=False,
        escape_char=False,
        newlines_in_values=False,
    )


def _csv_convert_options(schema, columns=None):
    # Nothing in .pairs is null: unmapped reads use the '!' and 0 sentinels, and
    # a readID of "NA" or "null" must survive as that string.
    return pacsv.ConvertOptions(
        column_types=schema,
        null_values=[],
        strings_can_be_null=False,
        quoted_strings_can_be_null=False,
        include_columns=list(columns) if columns else None,
    )


def open_pairs(
    path,
    columns=None,
    batch_size=DEFAULT_BATCH_SIZE,
    block_size=DEFAULT_BLOCK_SIZE,
    nproc_in=3,
    cmd_in=None,
):
    """Open a .pairs or .parquet file as a stream of Arrow record batches.

    Parameters
    ----------
    path : str
        Path to a ``.pairs``/``.pairs.gz``/``.parquet`` file, or ``-`` for stdin
        (text only -- Parquet cannot be streamed from a pipe).
    columns : list of str, optional
        Read only these columns. Column projection is where Parquet pays off, so
        pass this whenever a tool touches a subset.
    batch_size : int
        Rows per batch, for Parquet input.
    block_size : int
        Bytes per block, for text input. pyarrow sizes text batches by bytes.
    nproc_in, cmd_in
        Passed to ``pairtools.lib.fileio.auto_open`` for text input.

    Returns
    -------
    header : list of str
        The .pairs header lines.
    reader : pyarrow.RecordBatchReader
        Iterable of record batches. DuckDB can scan one of these directly, so a
        tool can hand it to SQL without materializing the file.
    """
    if is_parquet(path):
        header, batches, fallback = _open_parquet(path, columns, batch_size)
    else:
        header, batches, fallback = _open_text(
            path, columns, block_size, nproc_in, cmd_in
        )
    return header, _reader_from_batches(batches, fallback)


def _reader_from_batches(batches, fallback_schema):
    """Wrap a batch generator in a RecordBatchReader.

    The schema is taken from the first batch rather than declared up front, so
    a file written with types we did not choose (int64 positions, a dictionary
    column) streams through untouched instead of failing a schema check. An
    empty input has no first batch, and falls back to the header's schema.
    """
    iterator = iter(batches)
    try:
        first = next(iterator)
    except StopIteration:
        return pa.RecordBatchReader.from_batches(fallback_schema, iter(()))

    def chained():
        yield first
        for batch in iterator:
            yield batch

    return pa.RecordBatchReader.from_batches(first.schema, chained())


def _open_parquet(path, columns, batch_size):
    parquet_file = pq.ParquetFile(path)
    header = headerio.metadata_to_header(parquet_file.schema_arrow.metadata)

    def batches():
        try:
            for batch in parquet_file.iter_batches(
                batch_size=batch_size, columns=columns
            ):
                yield batch
        finally:
            parquet_file.close()

    fallback = parquet_file.schema_arrow.remove_metadata()
    if columns:
        fallback = pa.schema([fallback.field(name) for name in columns])
    return header, batches(), fallback


def _open_text(path, columns, block_size, nproc_in, cmd_in):
    instream = fileio.auto_open(path, mode="r", nproc=nproc_in, command=cmd_in)
    header, instream = headerops.get_header(instream)
    column_names = headerops.extract_column_names(header)
    if not column_names:
        raise ValueError(
            "{} has no '#columns:' header line, so its columns are unknown; "
            "add one with `pairtools header generate`".format(path)
        )

    schema = schema_from_columns(column_names)
    if columns:
        projected = pa.schema([schema.field(name) for name in columns])
    else:
        projected = schema

    try:
        reader = pacsv.open_csv(
            _binary_stream(instream),
            read_options=_csv_read_options(column_names, block_size),
            parse_options=_csv_parse_options(),
            convert_options=_csv_convert_options(schema, columns),
        )
    except pa.ArrowInvalid as error:
        # A header with no data rows is a valid .pairs file -- `select` that
        # matched nothing produces one -- but pyarrow rejects the empty body.
        if "Empty CSV file" not in str(error):
            raise
        if instream is not sys.stdin:
            instream.close()
        return header, iter(()), projected

    def batches():
        try:
            for batch in reader:
                yield batch
        finally:
            reader.close()
            if instream is not sys.stdin:
                instream.close()

    return header, batches(), projected


class PairsWriter:
    """Write Arrow record batches out as .pairs, .pairs.gz or .parquet.

    The output format is chosen by the extension of `path`. Use as a context
    manager, or call :meth:`close` when done -- Parquet in particular is invalid
    without its footer, which is only written on close.

    Batches whose schema differs from `schema` are cast on the way out, so a
    caller can hand over whatever its engine produced without matching types
    exactly.
    """

    def __init__(
        self,
        path,
        header,
        schema=None,
        compression="snappy",
        row_group_size=None,
        compress_program="auto",
        nproc_out=8,
    ):
        self.path = str(path)
        self.header = list(header)
        if schema is None:
            columns = headerops.extract_column_names(self.header)
            schema = schema_from_columns(columns)
        self.schema = schema

        self._is_parquet = is_parquet(self.path)
        self._row_group_size = row_group_size
        self._proc = None
        self._closed = False
        self._stack = contextlib.ExitStack()

        try:
            if self._is_parquet:
                self._open_parquet(compression)
            else:
                self._open_text(compress_program, nproc_out)
        except Exception:
            self._stack.close()
            raise

    def _open_parquet(self, compression):
        self._writer = self._stack.enter_context(
            pq.ParquetWriter(
                self.path,
                headerio.apply_to_schema(self.schema, self.header),
                compression=compression,
            )
        )

    def _open_text(self, compress_program, nproc_out):
        # Imported here so this module's import graph stays free of the DuckDB
        # side of the package.
        from .csv_parquet_converter import choose_compressor

        _, command = choose_compressor(compress_program, threads=nproc_out)

        output_file = self._stack.enter_context(open(self.path, "wb"))
        if command:
            # Registered before the process is entered, so it runs after
            # Popen.__exit__ has closed stdin and waited: the exit code is only
            # meaningful once the compressor has seen EOF.
            self._stack.callback(self._check_compressor)
            self._proc = self._stack.enter_context(
                subprocess.Popen(command, stdin=subprocess.PIPE, stdout=output_file)
            )
            sink_file = self._proc.stdin
        else:
            sink_file = output_file

        sink = pa.output_stream(sink_file)
        sink.write(
            "".join(line.rstrip() + "\n" for line in self.header).encode("utf-8")
        )

        self._writer = self._stack.enter_context(
            pacsv.CSVWriter(
                sink,
                self.schema,
                write_options=pacsv.WriteOptions(
                    include_header=False, delimiter="\t", quoting_style="none"
                ),
            )
        )

    def _check_compressor(self):
        if self._proc.returncode not in (0, None):
            raise RuntimeError(
                "compression command failed with exit code {}".format(
                    self._proc.returncode
                )
            )

    def write(self, data):
        """Write a ``RecordBatch`` or a ``Table``."""
        if data.num_rows == 0:
            return
        if isinstance(data, pa.RecordBatch):
            data = pa.Table.from_batches([data])
        if not data.schema.equals(self.schema, check_metadata=False):
            data = data.cast(self.schema)
        if self._is_parquet and self._row_group_size:
            self._writer.write_table(data, row_group_size=self._row_group_size)
        else:
            self._writer.write_table(data)

    def write_all(self, batches):
        """Write every batch of an iterable."""
        for batch in batches:
            self.write(batch)

    def close(self):
        if self._closed:
            return
        self._closed = True
        self._stack.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()
        return False

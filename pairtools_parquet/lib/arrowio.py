"""Format-agnostic streaming I/O for .pairs data as Arrow record batches.

Every tool reads through :func:`open_pairs` and writes through
:class:`PairsWriter`, so no tool needs to know whether it is looking at a
``.pairs`` text file, a ``.pairs.gz``, or a ``.parquet``. Both directions are
streaming: a file is never held in memory in full.

The module deliberately depends only on ``pyarrow``, ``numpy`` and
``pairtools`` -- all of which pairtools itself already requires -- and on
nothing else in this package, so it can move into pairtools unchanged.

Known limitation: pyarrow's CSV writer refuses to emit a value containing a
double quote when quoting is disabled, and .pairs is an unquoted format, so a
field with a ``"`` in it can be read but not written back to text. Such a value
round-trips through Parquet without complaint, and pairtools itself writes it
happily; only our text output rejects it, loudly rather than silently.
"""

import contextlib
import io
import subprocess
import sys

import numpy as np
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

#: Arrow IPC *stream* format, which -- unlike Parquet, whose footer is written
#: last -- can be read as it arrives, so it is the one binary format that can
#: go through a pipe.
ARROW_SUFFIXES = (".arrow", ".ipc")

#: First bytes of an Arrow IPC stream: the continuation token that precedes
#: every message. .pairs text always starts with '#', so the two cannot be
#: confused, and an input stream can be recognised rather than declared.
ARROW_STREAM_MAGIC = b"\xff\xff\xff\xff"

#: The path that means stdin on the way in and stdout on the way out.
STDIO_PATH = "-"

#: stdin can only be read once, so the header and the stream positioned just
#: past it are remembered. A tool that reads the header to decide what to do --
#: `sort` resolving its keys, `dedup` its columns -- and then opens the input
#: for real would otherwise find the header already consumed and the columns
#: unknowable.
_stdin_state = None

#: The same, for an Arrow stream: opening it reads the schema -- and so the
#: header -- without touching the batches, so a tool can ask for the header
#: first and still get every row afterwards.
_stdin_arrow = None

#: Extensions pairtools' auto_open treats as compressed.
COMPRESSED_SUFFIXES = (".gz", ".bz2", ".lz4")


def is_parquet(path):
    """Whether `path` names a Parquet file, by extension."""
    return str(path).lower().endswith(PARQUET_SUFFIXES)


def is_arrow(path):
    """Whether `path` names an Arrow IPC stream, by extension."""
    return str(path).lower().endswith(ARROW_SUFFIXES)


def is_stdio(path):
    """Whether `path` names the standard stream rather than a file.

    ``-`` is the .pairs stream both tools agree on, and so are an empty path
    and no path at all -- `pairtools` reads ``not path or path == "-"``, which
    is what makes an omitted ``-o`` mean stdout. ``-.arrow`` is the same stream
    carrying Arrow IPC instead: a stream has no name to hang an extension on,
    and giving one to ``-`` keeps the format in the path where every other
    format already lives, rather than in a flag that would have to be threaded
    through every tool.
    """
    if not path:
        return True
    path = str(path)
    return path == STDIO_PATH or path.startswith(STDIO_PATH + ".")


def compresses_by_extension(path):
    """Whether text written to `path` should be compressed, by extension."""
    return str(path).lower().endswith(COMPRESSED_SUFFIXES)


def read_header(path, nproc_in=3, cmd_in=None):
    """Read just the .pairs header of a file, in either format."""
    if is_parquet(path):
        return headerio.read_header(path)

    if is_arrow(path) and not is_stdio(path):
        with pa.ipc.open_stream(path) as reader:
            return headerio.metadata_to_header(reader.schema.metadata)

    if is_stdio(path):
        if _stdin_holds_arrow():
            return headerio.metadata_to_header(_open_stdin_arrow().schema.metadata)
        return list(_open_stdin(nproc_in, cmd_in)[0])

    instream = fileio.auto_open(path, mode="r", nproc=nproc_in, command=cmd_in)
    try:
        header, _ = headerops.get_header(instream)
    finally:
        if instream is not sys.stdin:
            instream.close()
    return header


def _open_stdin(nproc_in, cmd_in, ignore_warning=False):
    """Read stdin's header once, and keep the stream positioned after it."""
    global _stdin_state
    if _stdin_state is None:
        instream = fileio.auto_open(
            STDIO_PATH, mode="r", nproc=nproc_in, command=cmd_in
        )
        _stdin_state = headerops.get_header(
            instream, ignore_warning=ignore_warning
        )
    return _stdin_state


def _binary_stream(stream):
    """Return the underlying binary stream of a possibly-text stream.

    ``headerops.get_header`` consumes header lines through ``stream.buffer``
    when the stream is text, so the binary layer -- not the text wrapper -- is
    the one positioned at the first data line.
    """
    return stream.buffer if hasattr(stream, "buffer") else stream


def _csv_read_options(columns, block_size):
    return pacsv.ReadOptions(column_names=list(columns), block_size=block_size)


def csv_parse_options():
    # .pairs is strictly tab-separated with no quoting: a readID containing a
    # double quote is data, not a quoted field.
    return pacsv.ParseOptions(
        delimiter="\t",
        quote_char=False,
        double_quote=False,
        escape_char=False,
        newlines_in_values=False,
    )


def csv_convert_options(schema, columns=None):
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
    column_names=None,
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
    column_names : list of str, optional
        Read the columns positionally under these names, ignoring whatever the
        file calls them -- and, for text, whether it names them at all. This is
        what lets ``header generate`` put a header on a headerless file, and
        ``header set-columns`` rename the columns of one that has them. Column
        types follow the given names, since that is all a .pairs field's type
        is ever derived from.

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
        if column_names is not None:
            batches, fallback = _rename(batches, fallback, column_names, path)
    elif is_arrow(path) or (is_stdio(path) and _stdin_holds_arrow()):
        header, batches, fallback = _open_arrow(path, columns)
        if column_names is not None:
            batches, fallback = _rename(batches, fallback, column_names, path)
    else:
        header, batches, fallback = _open_text(
            path, columns, block_size, nproc_in, cmd_in, column_names
        )
    return header, _reader_from_batches(batches, fallback)


def _stdin_holds_arrow():
    """Whether what is waiting on stdin is an Arrow stream rather than text.

    Peeked, not consumed, so the reader that follows still sees the first
    bytes. A stream already claimed by either reader is not re-examined.
    """
    if _stdin_arrow is not None:
        return True
    if _stdin_state is not None:
        return False
    try:
        return sys.stdin.buffer.peek(4)[:4] == ARROW_STREAM_MAGIC
    except (AttributeError, ValueError, OSError):
        return False


def _open_stdin_arrow():
    """Open stdin as an Arrow stream once, and keep the reader."""
    global _stdin_arrow
    if _stdin_arrow is None:
        _stdin_arrow = pa.ipc.open_stream(sys.stdin.buffer)
    return _stdin_arrow


def _open_arrow(path, columns):
    """Read an Arrow IPC stream, from a file or from stdin."""
    reader = _open_stdin_arrow() if is_stdio(path) else pa.ipc.open_stream(path)
    header = headerio.metadata_to_header(reader.schema.metadata)

    if not columns:
        return header, reader, reader.schema

    projected = pa.schema([reader.schema.field(name) for name in columns])
    batches = (
        pa.RecordBatch.from_arrays(
            [batch.column(name) for name in columns], names=list(columns)
        )
        for batch in reader
    )
    return header, batches, projected


class _keep_open(io.RawIOBase):
    """A writable view of a stream that flushes, but never closes, the original.

    pyarrow closes whatever it was handed, and stdout is the caller's, not
    ours: a tool that writes several outputs must not lose stdout when the
    first of them finishes.
    """

    def __init__(self, stream):
        self._stream = stream

    def writable(self):
        return True

    def write(self, data):
        return self._stream.write(data)

    def flush(self):
        if not self.closed:
            self._stream.flush()

    def close(self):
        if not self.closed:
            self._stream.flush()
        # Marks this wrapper closed; the wrapped stream stays open.
        super().close()


def with_row_ids(reader, columns=None, counter=None, name="rid"):
    """Wrap `reader` in one that prepends a file-order row id.

    DuckDB has no equivalent for text input: ``row_number() OVER ()`` counts in
    the order rows reach the window operator, and its CSV scanner is parallel,
    so that is not the order they are in the file. ``read_parquet`` does expose
    ``file_row_number``, but only for Parquet. Counting here is exact for every
    input format, because an Arrow reader is sequential.

    `counter` is a one-element list, so a caller that needs the total row count
    can read it back once the reader has been drained.

    Returns a ``RecordBatchReader``, so the result can be registered with
    DuckDB and scanned like a table.
    """
    columns = list(columns) if columns is not None else reader.schema.names
    if counter is None:
        counter = [0]

    def batches():
        for batch in reader:
            rid = pa.array(
                np.arange(counter[0], counter[0] + batch.num_rows, dtype=np.int64)
            )
            counter[0] += batch.num_rows
            yield pa.RecordBatch.from_arrays(
                [rid] + [batch.column(column) for column in columns],
                names=[name] + columns,
            )

    schema = pa.schema(
        [pa.field(name, pa.int64())]
        + [reader.schema.field(column) for column in columns]
    )
    return pa.RecordBatchReader.from_batches(schema, batches())


def _rename(batches, fallback, column_names, path):
    """Rename the columns of every batch, positionally."""
    column_names = list(column_names)
    if len(fallback) != len(column_names):
        raise ValueError(
            "{} has {} columns, but {} names were given: {}".format(
                path, len(fallback), len(column_names), ", ".join(column_names)
            )
        )
    renamed = pa.schema(
        [field.with_name(name) for field, name in zip(fallback, column_names)]
    )
    return (batch.rename_columns(column_names) for batch in batches), renamed


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


def _open_text(path, columns, block_size, nproc_in, cmd_in, column_names=None):
    # A caller that supplied the column names is not relying on the file to
    # declare them, so a headerless input is not worth warning about.
    if is_stdio(path):
        header, instream = _open_stdin(
            nproc_in, cmd_in, ignore_warning=column_names is not None
        )
    else:
        instream = fileio.auto_open(
            path, mode="r", nproc=nproc_in, command=cmd_in
        )
        header, instream = headerops.get_header(
            instream, ignore_warning=column_names is not None
        )
    if column_names is None:
        column_names = headerops.extract_column_names(header)
    if not column_names:
        raise ValueError(
            "{} has no '#columns:' header line, so its columns are unknown; "
            "add one with `pairtools_parquet header generate`".format(path)
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
            parse_options=csv_parse_options(),
            convert_options=csv_convert_options(schema, columns),
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
        write_header=True,
    ):
        self.path = str(path)
        self.header = list(header)
        # `--send-header-to` exists so text outputs can be concatenated, which
        # is a text-only concern. Parquet keeps its header either way: there it
        # is key-value metadata rather than leading lines, and a file without it
        # cannot be read back as pairs at all.
        self._write_header = bool(write_header)
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
            elif is_arrow(self.path):
                self._open_arrow()
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

    def _open_arrow(self):
        """Write an Arrow IPC stream, to a file or to stdout.

        The header rides in the schema metadata, exactly as it does in
        Parquet, so a stream carries everything a file would.
        """
        if is_stdio(self.path):
            sink = self._stack.enter_context(_keep_open(sys.stdout.buffer))
        else:
            sink = self._stack.enter_context(open(self.path, "wb"))
        self._writer = self._stack.enter_context(
            pa.ipc.new_stream(sink, headerio.apply_to_schema(self.schema, self.header))
        )

    def _open_text(self, compress_program, nproc_out):
        # Imported here so this module's import graph stays free of the DuckDB
        # side of the package.
        from .csv_parquet_converter import choose_compressor

        if is_stdio(self.path):
            # `-` is a stream, so there is no extension to decide compression
            # from and nothing to reopen -- pipe through a compressor yourself
            # if you want one, exactly as `pairtools <tool> -o -` requires.
            sink_file = _keep_open(sys.stdout.buffer)
        else:
            # Whether to compress is decided by the extension, as pairtools'
            # auto_open does; --compress-program only picks which compressor. A
            # plain `.pairs` output must be plain text, however that flag is set.
            if compresses_by_extension(self.path):
                _, command = choose_compressor(compress_program, threads=nproc_out)
            else:
                command = []

            output_file = self._stack.enter_context(open(self.path, "wb"))
            if command:
                # Registered before the process is entered, so it runs after
                # Popen.__exit__ has closed stdin and waited: the exit code is
                # only meaningful once the compressor has seen EOF.
                self._stack.callback(self._check_compressor)
                self._proc = self._stack.enter_context(
                    subprocess.Popen(
                        command, stdin=subprocess.PIPE, stdout=output_file
                    )
                )
                sink_file = self._proc.stdin
            else:
                sink_file = output_file

        sink = pa.output_stream(sink_file)
        # Only the line terminator is stripped, not trailing whitespace:
        # `headerops.get_header` keeps trailing spaces, and pairtools emits at
        # least one line that ends in one (`#genome_assembly: ` with no
        # assembly given), so stripping them would break byte parity.
        if self._write_header:
            sink.write(
                "".join(
                    line.rstrip("\n") + "\n" for line in self.header
                ).encode("utf-8")
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

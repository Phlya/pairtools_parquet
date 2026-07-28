"""Re-cutting a stream of DataFrames into fixed-size chunks.

Arrow batch boundaries are set by the reader -- by row-group size for Parquet,
by block size in bytes for text -- and have nothing to do with the ``--chunksize``
a pairtools command was given. Where a tool's result depends on which rows share
a chunk, matching pairtools means matching its chunk boundaries exactly, so the
batches have to be re-cut before the kernel sees them.

`dedup` needs this because a chunk is the window within which duplicates are
searched. `scaling` needs it because per-chunk results are combined with
``DataFrame.add``, which promotes integer counts to floats -- so whether a file
is read in one chunk or several is visible in the output.
"""

import pandas as pd


def rechunk(frames, chunksize):
    """Re-cut an iterable of DataFrames into frames of exactly `chunksize` rows.

    The last frame is whatever is left over.
    """
    buffered = []
    buffered_rows = 0

    for frame in frames:
        buffered.append(frame)
        buffered_rows += len(frame)

        while buffered_rows >= chunksize:
            combined = (
                buffered[0] if len(buffered) == 1 else pd.concat(buffered, axis=0)
            )
            combined = combined.reset_index(drop=True)
            yield combined.iloc[:chunksize].reset_index(drop=True)

            leftover = combined.iloc[chunksize:].reset_index(drop=True)
            buffered = [leftover] if len(leftover) else []
            buffered_rows = len(leftover)

    if buffered_rows:
        combined = buffered[0] if len(buffered) == 1 else pd.concat(buffered, axis=0)
        yield combined.reset_index(drop=True)

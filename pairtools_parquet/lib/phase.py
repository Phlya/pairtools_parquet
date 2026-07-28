"""Assigning pairs to parental haplotypes.

The phasing decision for one side of one pair is pairtools'
``phase_side_XB``/``phase_side_XA``, called unchanged. Those inspect the
aligner's XB/XA, AS, XS and NM tags, so this tool is row-wise both here and
upstream; what it saves is the text parsing and re-serialization around it.

The input must have been parsed with the alignment tags kept
(``--add-columns XB,AS,XS`` or ``--add-columns XA,NM,AS,XS``), as upstream
requires.
"""

import pyarrow as pa
from pairtools.lib import headerops, pairsam_format
from pairtools.lib.phase import phase_side_XA, phase_side_XB

from .arrowio import PairsWriter, open_pairs

UTIL_NAME = "pairtools_parquet_phase"

#: Extra columns each tag mode reports with --report-scores.
SCORE_COLUMNS = {
    "XB": ["S1_1", "S1_2", "S2_1", "S2_2", "S3_1", "S3_2"],
    "XA": ["M1_1", "M1_2", "M2_1", "M2_2", "M3_1", "M3_2"],
}

#: Columns each tag mode needs in the input.
REQUIRED_COLUMNS = {
    "XB": ["XB1", "XB2", "AS1", "AS2", "XS1", "XS2"],
    "XA": ["XA1", "XA2", "NM1", "NM2", "AS1", "AS2", "XS1", "XS2"],
}


def output_columns(input_columns, tag_mode, clean_output, report_scores):
    """The columns of the phased output, in order."""
    if clean_output:
        columns = [c for c in input_columns if c in pairsam_format.COLUMNS]
    else:
        columns = list(input_columns)

    columns += ["phase1", "phase2"]
    if report_scores:
        columns += SCORE_COLUMNS[tag_mode]
    return columns


def check_required_columns(input_columns, tag_mode):
    missing = [c for c in REQUIRED_COLUMNS[tag_mode] if c not in input_columns]
    if missing:
        raise ValueError(
            "The input pairs file must be parsed with the flag "
            "--add-columns {} --min-mapq 0 (missing: {})".format(
                "XB,AS,XS" if tag_mode == "XB" else "XA,NM,AS,XS",
                ",".join(missing),
            )
        )


def phase_side(row, side, tag_mode, phase_suffixes):
    """Phase one side of one pair, via pairtools' own per-side function."""
    chrom = row["chrom{}".format(side)]
    if tag_mode == "XB":
        return phase_side_XB(
            chrom,
            row["XB{}".format(side)],
            int(row["AS{}".format(side)]),
            int(row["XS{}".format(side)]),
            phase_suffixes,
        )
    return phase_side_XA(
        chrom,
        row["XA{}".format(side)],
        int(row["AS{}".format(side)]),
        int(row["XS{}".format(side)]),
        int(row["NM{}".format(side)]),
        phase_suffixes,
    )


def phase_row(row, tag_mode, phase_suffixes, report_scores):
    """Return the phased version of one row, as a dict."""
    out = dict(row)
    out["phase1"] = "!"
    out["phase2"] = "!"
    if report_scores:
        for name in SCORE_COLUMNS[tag_mode]:
            out[name] = -1

    pair_type = row.get("pair_type", "")

    for side in (1, 2):
        if row["chrom{}".format(side)] == pairsam_format.UNMAPPED_CHROM:
            continue

        phase, chrom_base, s1, s2, s3 = phase_side(
            row, side, tag_mode, phase_suffixes
        )
        out["phase{}".format(side)] = phase
        if report_scores:
            names = SCORE_COLUMNS[tag_mode]
            # the six score columns are interleaved: S1_1, S1_2, S2_1, ...
            out[names[0 + (side - 1)]] = s1
            out[names[2 + (side - 1)]] = s2
            out[names[4 + (side - 1)]] = s3

        out["chrom{}".format(side)] = chrom_base
        if chrom_base == "!":
            out["chrom{}".format(side)] = pairsam_format.UNMAPPED_CHROM
            out["pos{}".format(side)] = pairsam_format.UNMAPPED_POS
            out["strand{}".format(side)] = pairsam_format.UNMAPPED_STRAND
            # the unmappable side becomes 'M' in the pair type
            pair_type = (
                "M" + pair_type[1:] if side == 1 else pair_type[:1] + "M"
            )

    if "pair_type" in out:
        out["pair_type"] = pair_type
    return out


def phase_pairs(
    input_path,
    output,
    phase_suffixes,
    tag_mode="XB",
    clean_output=False,
    report_scores=False,
    compress_program="auto",
    row_group_size=None,
    util_name=UTIL_NAME,
    **kwargs
):
    """Phase the pairs of `input_path` onto parental haplotypes."""
    if len(phase_suffixes) != 2:
        raise ValueError("Exactly two phase suffixes are expected")

    header, reader = open_pairs(
        input_path,
        nproc_in=kwargs.get("nproc_in", 3),
        cmd_in=kwargs.get("cmd_in", None),
    )
    input_columns = headerops.extract_column_names(header)
    check_required_columns(input_columns, tag_mode)

    new_header = headerops.append_new_pg(
        list(header), ID=util_name, PN=util_name
    )
    columns = output_columns(input_columns, tag_mode, clean_output, report_scores)
    new_header = headerops._update_header_entry(
        new_header, "columns", " ".join(columns)
    )

    with PairsWriter(
        output,
        new_header,
        compress_program=compress_program,
        row_group_size=row_group_size,
        nproc_out=kwargs.get("nproc_out", 8),
    ) as writer:
        for batch in reader:
            rows = [
                phase_row(row, tag_mode, phase_suffixes, report_scores)
                for row in batch.to_pylist()
            ]
            writer.write(
                pa.Table.from_pylist([{c: row[c] for c in columns} for row in rows])
            )

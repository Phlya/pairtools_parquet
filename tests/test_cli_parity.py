"""Every pairtools command and option must exist here under the same name.

The point of this package is that `pairtools_parquet <tool>` can stand in for
`pairtools <tool>`, so an option we renamed or never added is a break even when
the underlying feature works. Nothing else catches that: the rest of the suite
invokes options by whichever spelling we chose, so a renamed option still
passes. `flip`'s chromosome-order file was called `--chrom-subset` here and
`--chroms-path` upstream for exactly that reason -- every test used `-c`.

Where we deliberately do not accept an option, it is listed below with the
reason, so the omission is a decision on the record rather than an oversight.
"""

import pytest

pairtools_cli = pytest.importorskip("pairtools.cli")

from pairtools_parquet import cli as our_cli


#: pairtools options we do not accept, and why. Anything not listed here has to
#: exist, with the same name.
DOCUMENTED_OMISSIONS = {
    ("dedup", "--sep"): "the .pairs spec fixes the separator at a tab, and it "
    "is meaningless for parquet input; arrowio reads tab-separated text only",
    ("filterbycov", "--sep"): "as dedup --sep",
    ("filterbycov", "--comment-char"): "arrowio delegates header parsing to "
    "pairtools' headerops, which hardcodes '#'",
}

SHARED_COMMANDS = sorted(set(pairtools_cli.cli.commands) & set(our_cli.cli.commands))


def option_names(command):
    return {
        opt
        for param in command.params
        for opt in param.opts
        if opt.startswith("--")
    }


def test_every_pairtools_command_exists():
    missing = set(pairtools_cli.cli.commands) - set(our_cli.cli.commands)
    assert not missing, "commands with no pairtools_parquet equivalent: {}".format(
        sorted(missing)
    )


@pytest.mark.parametrize("name", SHARED_COMMANDS)
def test_every_pairtools_option_exists(name):
    theirs = option_names(pairtools_cli.cli.commands[name])
    ours = option_names(our_cli.cli.commands[name])

    missing = {
        opt for opt in theirs - ours if (name, opt) not in DOCUMENTED_OMISSIONS
    }
    assert not missing, (
        "pairtools_parquet {} does not accept {}. Add it, or record why not in "
        "DOCUMENTED_OMISSIONS.".format(name, ", ".join(sorted(missing)))
    )


@pytest.mark.parametrize("name", SHARED_COMMANDS)
def test_documented_omissions_are_real(name):
    """An omission that upstream does not have either is stale bookkeeping."""
    theirs = option_names(pairtools_cli.cli.commands[name])
    ours = option_names(our_cli.cli.commands[name])
    for command, opt in DOCUMENTED_OMISSIONS:
        if command != name:
            continue
        assert opt in theirs, "{} {} is not a pairtools option".format(name, opt)
        assert opt not in ours, "{} {} is listed as omitted but we accept it".format(
            name, opt
        )


@pytest.mark.parametrize("name", SHARED_COMMANDS)
def test_short_flags_keep_their_long_names(name):
    """A short flag must not be paired with a different long name than upstream.

    This is the specific shape of the `flip -c` bug: the short flag matched, so
    every call site kept working while the long name said something else.
    """
    theirs = {
        short: sorted(o for o in param.opts if o.startswith("--"))
        for param in pairtools_cli.cli.commands[name].params
        for short in param.opts
        if short.startswith("-") and not short.startswith("--")
    }
    ours = {
        short: sorted(o for o in param.opts if o.startswith("--"))
        for param in our_cli.cli.commands[name].params
        for short in param.opts
        if short.startswith("-") and not short.startswith("--")
    }

    for short, long_names in theirs.items():
        if short in ours and long_names and ours[short]:
            assert ours[short] == long_names, (
                "{} {} is {} upstream but {} here".format(
                    name, short, long_names, ours[short]
                )
            )

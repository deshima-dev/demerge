__all__ = []


# standard library
from pathlib import Path
from subprocess import run
from sys import argv


def demerge_sh() -> None:
    """Run cli.sh (formerly run.sh)."""

    cmd = Path(__file__).parent / "cli.sh"
    run([cmd, *argv[1:]])

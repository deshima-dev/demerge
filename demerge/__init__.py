__all__ = ["merge_function", "merge_to_dems"]
__version__ = "2.13.0"


# submodules
from . import merge_function
from . import merge_to_dems


def main() -> None:
    """Run run.sh (works on Linux-like systems)."""
    from pathlib import Path
    from subprocess import run
    from sys import argv

    cmd = Path(__file__).parent / "run.sh"
    run([cmd, *argv[1:]])

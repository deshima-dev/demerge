__all__ = [
    "calc_resonance_params",
    "demerge",
    "make_divided_data",
    "make_reduced_fits",
    "merge_function",
    "merge_to_dems",
    "plot",
]
__version__ = "2.13.0"


# submodules
from . import calc_resonance_params
from . import demerge
from . import make_divided_data
from . import make_reduced_fits
from . import merge_function
from . import merge_to_dems
from . import plot


def main() -> None:
    """Run run.sh (works on Linux-like systems)."""
    from pathlib import Path
    from subprocess import run
    from sys import argv

    cmd = Path(__file__).parent / "run.sh"
    run([cmd, *argv[1:]])

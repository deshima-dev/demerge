__all__ = ["analyze"]


# standard library
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory


# dependencies
from fire import Fire


# constants
SCRIPTS = Path(__file__).parent / "utils" / "scripts" / "aste"


def analyze(data_dir: str, output_dir: str) -> None:
    """Analyze a DESHIMA TOD FITS and create a reduced FITS.

    Args:
        data_dir: Path of raw data directory (e.g. ``cosmos_20240101000000``).
        output_dir: Path of output directory (e.g. ``analysis_20240101000000``).

    """
    data_dir = Path(data_dir).resolve()
    output_dir = Path(output_dir).resolve()

    with TemporaryDirectory() as work_dir:
        run(["python", SCRIPTS / "Configure.py", data_dir, output_dir], cwd=work_dir)
        run(["python", SCRIPTS / "FitSweep.py"], cwd=work_dir)
        run(["python", SCRIPTS / "SaveFits.py"], cwd=work_dir)


def cli() -> None:
    """Command line interface of the analyze function."""
    Fire(analyze)

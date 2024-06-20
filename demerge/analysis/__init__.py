__all__ = ["analyze"]


# standard library
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory


# dependencies
from fire import Fire


# constants
SCRIPTS = Path(__file__).parent / "utils" / "scripts" / "aste"


def analyze(data_dir: Path, output_dir: Path) -> Path:
    """Analyze a DESHIMA TOD FITS and create a reduced FITS.

    Args:
        data_dir: Path of raw data directory (e.g. ``cosmos_20240101000000``).
        output_dir: Path of output directory (e.g. ``analysis_20240101000000``).

    Returns:
        Path of the created reduced FITS.

    """
    data_dir = Path(data_dir).resolve()
    output_dir = Path(output_dir).resolve()

    with TemporaryDirectory() as work_dir:
        run(
            ["python", SCRIPTS / "Configure.py", data_dir, output_dir],
            capture_output=True,
            cwd=work_dir,
        )
        run(
            ["python", SCRIPTS / "FitSweep.py"],
            capture_output=True,
            cwd=work_dir,
        )
        run(
            ["python", SCRIPTS / "SaveFits.py"],
            capture_output=True,
            cwd=work_dir,
        )

    return list(output_dir.glob("reduced_*.fits"))[0]


def cli() -> None:
    """Command line interface of the analyze function."""
    Fire(analyze)

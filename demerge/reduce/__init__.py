__all__ = ["reduce"]


# standard library
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory


# dependencies
from fire import Fire


# constants
SCRIPTS = Path(__file__).parent / "utils" / "scripts" / "aste"


def reduce(data_dir: Path, output_dir: Path, /) -> Path:
    """Reduce raw data of KID measurements into a single "reduced" FITS.

    Args:
        data_dir: Path of raw data directory (e.g. ``cosmos_YYYYmmddHHMMSS``).
        output_dir: Path of output directory (e.g. ``output_YYYYmmddHHMMSS``).

    Returns:
        Path of the created reduced FITS (in the output directory).

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
    """Command line interface of the reduce function."""
    Fire(reduce)

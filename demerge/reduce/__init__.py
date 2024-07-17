__all__ = ["reduce"]


# standard library
from contextlib import contextmanager
from logging import DEBUG, basicConfig, getLogger
from pathlib import Path
from shutil import rmtree
from subprocess import run
from tempfile import TemporaryDirectory
from typing import Union


# dependencies
from fire import Fire


# type hints
PathLike = Union[Path, str]


# constants
LOGGER = getLogger(__name__)
SCRIPTS = Path(__file__).parent / "utils" / "scripts" / "aste"


@contextmanager
def set_logger(debug: bool):
    level = LOGGER.level

    if debug:
        LOGGER.setLevel(DEBUG)

    try:
        yield
    finally:
        LOGGER.setLevel(level)


def reduce(
    *,
    data_dir: PathLike,
    reduced_dir: PathLike,
    overwrite: bool = False,
    debug: bool = False,
) -> Path:
    """Reduce raw data of KID measurements into a single "reduced" FITS.

    Args:
        data_dir: Path of raw data directory (e.g. ``cosmos_YYYYmmddHHMMSS``).
        reduced_dir: Path of reduced data directory (e.g. ``reduced_YYYYmmddHHMMSS``).
        overwrite: If True, ``reduced_dir`` will be overwritten even if it exists.
        debug: If True, detailed logs for debugging will be printed.

    Returns:
        Path of the created reduced FITS (in the output directory).

    Raises:
        FileNotFoundError: Raised if ``data_dir`` does not exist.
        FileExistsError: Raised if ``reduced_dir`` exists.

    """
    with set_logger(debug):
        for key, val in locals().items():
            LOGGER.debug(f"{key}: {val!r}")

    # Resolve paths (must be done before changing working directory)
    if not (data_dir := Path(data_dir).resolve()).exists():
        raise FileNotFoundError(data_dir)

    if (reduced_dir := Path(reduced_dir).resolve()).exists() and not overwrite:
        raise FileExistsError(reduced_dir)

    if overwrite:
        rmtree(reduced_dir, ignore_errors=True)

    # Run scripts in a temporary directory (to isolate intermediate files)
    with TemporaryDirectory() as work_dir:
        run(
            ["python", SCRIPTS / "Configure.py", data_dir, reduced_dir],
            check=True,
            cwd=work_dir,
            # False if logging is implemented
            capture_output=True,
        )
        run(
            ["python", SCRIPTS / "FitSweep.py"],
            check=True,
            cwd=work_dir,
            # False if logging is implemented
            capture_output=True,
        )
        run(
            ["python", SCRIPTS / "SaveFits.py"],
            check=True,
            cwd=work_dir,
            # False if logging is implemented
            capture_output=True,
        )

    return list(reduced_dir.glob("reduced_*.fits"))[0]


def cli() -> None:
    """Command line interface of the reduce function."""
    basicConfig(
        datefmt="%Y-%m-%d %H:%M:%S",
        format="[%(asctime)s %(name)s %(funcName)s %(levelname)s] %(message)s",
    )

    Fire(reduce)

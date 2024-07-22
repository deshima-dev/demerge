__all__ = ["reduce"]


# standard library
from collections.abc import Iterator
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


def set_dir(dir: PathLike, /) -> Path:
    """Resolve a directory."""
    return Path(dir).expanduser().resolve()


@contextmanager
def set_logger(debug: bool, /) -> Iterator[None]:
    """Temporarily set the level of the module logger."""
    level = LOGGER.level

    if debug:
        LOGGER.setLevel(DEBUG)

    try:
        yield
    finally:
        LOGGER.setLevel(level)


def reduce(
    *,
    data_pack: PathLike,
    reduced_pack: PathLike,
    overwrite: bool = False,
    debug: bool = False,
) -> Path:
    """Reduce KID measurements into a single "reduced" FITS.

    Args:
        data_pack: Path of data package (e.g. ``cosmos_YYYYmmddHHMMSS``).
        reduced_pack: Path of reduced package (e.g. ``reduced_YYYYmmddHHMMSS``).
        overwrite: If True, ``reduced_pack`` will be overwritten even if it exists.
        debug: If True, detailed logs for debugging will be printed.

    Returns:
        Path of the created reduced FITS (in the reduced package).

    Raises:
        FileNotFoundError: Raised if ``data_pack`` does not exist.
        FileExistsError: Raised if ``reduced_pack`` exists and overwrite is False.

    """
    with set_logger(debug):
        for key, val in locals().items():
            LOGGER.debug(f"{key}: {val!r}")

    # Resolve paths (must be done before changing working directory)
    if not (data_pack := set_dir(data_pack)).exists():
        raise FileNotFoundError(data_pack)

    if (reduced_pack := set_dir(reduced_pack)).exists() and not overwrite:
        raise FileExistsError(reduced_pack)

    if overwrite:
        rmtree(reduced_pack, ignore_errors=True)

    # Run scripts in a temporary directory (to isolate intermediate files)
    with TemporaryDirectory() as work_dir:
        run(
            ["python", SCRIPTS / "Configure.py", data_pack, reduced_pack],
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

    return list(reduced_pack.glob("*.fits"))[0]


def reduce_cli() -> None:
    """Command line interface of the reduce function."""
    basicConfig(
        datefmt="%Y-%m-%d %H:%M:%S",
        format="[%(asctime)s %(name)s %(funcName)s %(levelname)s] %(message)s",
    )

    Fire(reduce)

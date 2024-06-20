__all__ = ["reduce"]


# standard library
from logging import DEBUG, basicConfig, getLogger
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory


# dependencies
from fire import Fire


# constants
LOGGER = getLogger(__name__)
SCRIPTS = Path(__file__).parent / "utils" / "scripts" / "aste"


def reduce(
    data_dir: Path,
    output_dir: Path,
    /,
    *,
    debug: bool = False,
) -> Path:
    """Reduce raw data of KID measurements into a single "reduced" FITS.

    Args:
        data_dir: Path of raw data directory (e.g. ``cosmos_YYYYmmddHHMMSS``).
        output_dir: Path of output directory (e.g. ``output_YYYYmmddHHMMSS``).
        debug: If True, detailed logs for debugging will be printed.

    Returns:
        Path of the created reduced FITS (in the output directory).

    Raises:
        FileNotFoundError: Raised if ``data_dir`` does not exist.
        FileExistsError: Raised if ``output_dir`` exists.

    """
    if debug:
        LOGGER.setLevel(DEBUG)

    basicConfig(
        datefmt="%Y-%m-%d %H:%M:%S",
        format="[%(asctime)s %(name)s %(levelname)s] %(message)s",
    )

    for key, val in locals().items():
        LOGGER.debug(f"{key}: {val!r}")

    # Resolve paths (must be done before changing working directory)
    if not (data_dir := Path(data_dir).resolve()).exists():
        raise FileNotFoundError(data_dir)

    if (output_dir := Path(output_dir).resolve()).exists():
        raise FileExistsError(output_dir)

    # Run scripts in a temporary directory (to isolate intermediate files)
    with TemporaryDirectory() as work_dir:
        run(
            ["python", SCRIPTS / "Configure.py", data_dir, output_dir],
            cwd=work_dir,
            # False if logging is implemented
            capture_output=True,
        )
        run(
            ["python", SCRIPTS / "FitSweep.py"],
            cwd=work_dir,
            # False if logging is implemented
            capture_output=True,
        )
        run(
            ["python", SCRIPTS / "SaveFits.py"],
            cwd=work_dir,
            # False if logging is implemented
            capture_output=True,
        )

    return list(output_dir.glob("reduced_*.fits"))[0]


def cli() -> None:
    """Command line interface of the reduce function."""
    Fire(reduce)

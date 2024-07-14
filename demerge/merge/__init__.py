__all__ = ["merge"]


# standard library
from contextlib import contextmanager
from logging import DEBUG, basicConfig, getLogger
from pathlib import Path


# dependencies
from fire import Fire
from .utils import create_dems


# constants
LOGGER = getLogger(__name__)


@contextmanager
def set_logger(debug: bool):
    level = LOGGER.level

    if debug:
        LOGGER.setLevel(DEBUG)

    try:
        yield
    finally:
        LOGGER.setLevel(level)


def merge(
    dems: Path,
    /,
    *,
    ddb: str,
    corresp: str,
    obsinst: str,
    antenna: str,
    readout: str,
    skychop: str,
    weather: str,
    misti: str = "",
    cabin: str = "",
    coordinate: str = "azel",
    measure: str = "df/f",
    overwrite: bool = False,
    debug: bool = False,
    offset_time_antenna: int = 0,
) -> Path:
    """Merge datasets of an observation into a single DEMS file.

    Args:
        dems: Path of the output DEMS file (.zarr.zip).
        ddb: Path of the DDB file (.fits or .fits.gz).
        corresp: Path of the Master-to-KID ID correspondence file (.json).
        obsinst: Path of the observation instruction file (.obs).
        antenna: Path of the antenna log file (.ant).
        readout: Path of the reduced FITS file (.fits).
        skychop: Path of the Sky chopper file (.skychop).
        weather: Path of the weather log file (.weather).
        misti: Path of the MiSTI log file (.misti).
        cabin: Path of the cabin log file (.cabin).
        coordinate: Coordinate system of the output data (azel or radec).
        measure: Output data format (df/f or brightness).
        overwrite: If True, ``dems`` will be overwritten even if it exists.
        debug: If True, detailed logs for debugging will be printed.
        offset_time_antenna: Time diff (ms) between reduced FITS and antenna log.

    Returns:
        Path of the merged DEMS file.

    Raises:
        FileExistsError: Raised if ``dems`` exists.

    """
    with set_logger(debug):
        for key, val in locals().items():
            LOGGER.debug(f"{key}: {val!r}")

    da = create_dems(
        ddbfits_path=ddb,
        corresp_path=corresp,
        obsinst_path=obsinst,
        antenna_path=antenna,
        readout_path=readout,
        skychop_path=skychop,
        weather_path=weather,
        misti_path=misti,
        cabin_path=cabin,
        coordinate=coordinate,
        measure=measure,
        offset_time_antenna=offset_time_antenna,
    )

    if (dems := Path(dems)).exists() and not overwrite:
        raise FileExistsError(dems)

    if overwrite:
        dems.unlink(missing_ok=True)

    dems.parent.mkdir(exist_ok=True, parents=True)
    da.to_zarr(dems, mode="w")
    return dems.resolve()


def cli() -> None:
    """Command line interface of the merge function."""
    basicConfig(
        datefmt="%Y-%m-%d %H:%M:%S",
        format="[%(asctime)s %(name)s %(funcName)s %(levelname)s] %(message)s",
    )

    Fire(merge)

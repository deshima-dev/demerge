__all__ = ["merge"]


# standard library
from contextlib import contextmanager
from logging import DEBUG, basicConfig, getLogger
from pathlib import Path
from typing import Literal, Optional, Union


# dependencies
from fire import Fire
from .utils import to_brightness, to_dems


# type hints
PathLike = Union[Path, str]


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
    # required datasets
    corresp: PathLike,
    ddb: PathLike,
    obsinst: PathLike,
    readout: PathLike,
    # optional datasets
    antenna: Optional[PathLike] = None,
    cabin: Optional[PathLike] = None,
    misti: Optional[PathLike] = None,
    skychop: Optional[PathLike] = None,
    weather: Optional[PathLike] = None,
    # optional time offsets
    dt_antenna: Union[int, str] = "0 ms",
    dt_cabin: Union[int, str] = "0 ms",
    dt_misti: Union[int, str] = "0 ms",
    dt_skychop: Union[int, str] = "0 ms",
    dt_weather: Union[int, str] = "0 ms",
    # merge options
    measure: Literal["df/f", "brightness"] = "df/f",
    overwrite: bool = False,
    debug: bool = False,
) -> Path:
    """Merge observation datasets into a single DEMS.

    Args:
        dems: Path of the merged DEMS.
        corresp: Path of the KID correspondence.
        ddb: Path of DDB FITS.
        obsinst: Path of the observation instruction.
        readout: Path of the reduced readout FITS.
        antenna: Path of the antenna log.
        cabin: Path of the cabin log.
        misti: Path of the MiSTI log.
        skychop: Path of the sky chopper log.
        weather: Path of the weather log.
        dt_antenna: Time offset of the antenna log with explicit
            unit such that (dt_antenna = t_antenna - t_readout).
        dt_cabin: Time offset of the cabin log with explicit
            unit such that (dt_cabin = t_cabin - t_readout).
        dt_misti: Time offset of the MiSTI log with explicit
            unit such that (dt_misti = t_misti - t_readout).
        dt_skychop: Time offset of the sky chopper log with explicit
            unit such that (dt_skychop = t_skychop - t_readout).
        dt_weather: Time offset of the weather log with explicit
            unit such that (dt_weather = t_weather - t_readout).
        measure: Measure of the DEMS (either df/f or brightness).
        overwrite: If True, ``dems`` will be overwritten even if it exists.
        debug: If True, detailed logs for debugging will be printed.

    Returns:
        Path of the merged DEMS.

    Raises:
        FileExistsError: Raised if ``dems`` exists and ``overwrite`` is False.

    """
    with set_logger(debug):
        for key, val in locals().items():
            LOGGER.debug(f"{key}: {val!r}")

    da = to_dems(
        # required datasets
        corresp=corresp,
        ddb=ddb,
        obsinst=obsinst,
        readout=readout,
        # optional datasets
        antenna=antenna,
        cabin=cabin,
        misti=misti,
        skychop=skychop,
        weather=weather,
        # optional time offsets
        dt_antenna=dt_antenna,
        dt_cabin=dt_cabin,
        dt_misti=dt_misti,
        dt_skychop=dt_skychop,
        dt_weather=dt_weather,
    )

    if measure == "brightness":
        da = to_brightness(da)

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

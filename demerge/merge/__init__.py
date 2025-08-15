__all__ = ["merge", "utils"]


# standard library
from collections.abc import Iterator
from contextlib import contextmanager
from logging import DEBUG, basicConfig, getLogger
from os import PathLike
from pathlib import Path
from typing import Any, Optional, Union


# dependencies
from . import utils
from fire import Fire
from .utils import to_dems


# type hints
StrPath = Union[PathLike[str], str]


# constants
LOGGER = getLogger(__name__)


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


def merge(
    dems: StrPath,
    /,
    *,
    # required datasets
    cdb: StrPath,
    ddb: StrPath,
    obsinst: StrPath,
    readout: StrPath,
    # optional datasets
    antenna: Optional[StrPath] = None,
    cabin: Optional[StrPath] = None,
    misti: Optional[StrPath] = None,
    skychop: Optional[StrPath] = None,
    weather: Optional[StrPath] = None,
    # optional time offsets
    dt_antenna: Union[int, str] = "0 ms",
    dt_cabin: Union[int, str] = "0 ms",
    dt_misti: Union[int, str] = "0 ms",
    dt_skychop: Union[int, str] = "9 ms",
    dt_weather: Union[int, str] = "0 ms",
    # optional merge strategies
    include_disabled_mkids: bool = False,
    include_filterless_mkids: bool = False,
    overwrite: bool = False,
    debug: bool = False,
    **_: Any,
) -> Path:
    """Merge observation datasets into a single DEMS.

    Args:
        dems: Path of the merged DEMS.
        cdb: Path of CDB (KID correspondence database) file.
        ddb: Path of DDB (DESHIMA database) file.
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
            Defaults to 9 ms (for DESHIMA campaign in 2024).
        dt_weather: Time offset of the weather log with explicit
            unit such that (dt_weather = t_weather - t_readout).
        include_disabled_mkids: Whether to include disabled
            (e.g. fit-failed) MKID responses in the merged DEMS.
            Note that such data will be all filled with NaN.
        include_filterless_mkids: Whether to include wideband and/or
            no-filter-information MKID responses in the merged DEMS.
            Note that such data will be all filled with NaN.
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
        cdb=cdb,
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
        # optional merge strategies
        include_disabled_mkids=include_disabled_mkids,
        include_filterless_mkids=include_filterless_mkids,
    )

    if (dems := Path(dems)).exists() and not overwrite:
        raise FileExistsError(dems)

    if overwrite:
        dems.unlink(missing_ok=True)

    dems.parent.mkdir(exist_ok=True, parents=True)
    da.to_zarr(dems, mode="w")
    return dems.resolve()


def merge_cli() -> None:
    """Command line interface of the merge function."""
    basicConfig(
        datefmt="%Y-%m-%d %H:%M:%S",
        format="[%(asctime)s %(name)s %(funcName)s %(levelname)s] %(message)s",
    )

    Fire(merge)

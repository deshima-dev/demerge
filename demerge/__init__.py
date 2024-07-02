__all__ = ["demerge", "merge", "reduce"]
__version__ = "3.0.1"


# standard library
from logging import DEBUG, basicConfig, getLogger
from pathlib import Path
from typing import Any


# dependencies
from fire import Fire
from . import merge, reduce


# constants
DEFAULT_DDB_FILE = Path(__file__).parent / "data" / "ddb_20231123.fits.gz"
DEFAULT_OFFSET_TIME_ANTENNA = 20  # ms
LOGGER = getLogger(__name__)


def demerge(
    obsid: str,
    /,
    *,
    data_dir: Path = Path(),
    dems_dir: Path = Path(),
    reduced_dir: Path = Path(),
    ddb: Path = DEFAULT_DDB_FILE,
    debug: bool = False,
    **merge_options: Any,
) -> Path:
    """Run reduce and merge commands to create a single DEMS file.

    Args:
        obsid: Observation ID (YYYYmmddHHMMSS).
        data_dir: Path where raw data directory is placed,
            i.e. expecting ``${data_dir}/cosmos_YYYYmmddHHMMSS``.
        dems_dir: Path where merged DEMS file will be placed,
            i.e. expecting ``${dems_dir}/dems_YYYYmmddHHMMSS.zarr.zip``.
        reduced_dir: Path where reduced data directory will be placed,
            i.e. expecting ``${reduced_dir}/reduced_YYYYmmddHHMMSS``.
        ddb: Path of DDB (DESHIMA database) file.
        debug: If True, detailed logs for debugging will be printed.
        **merge_options: Other merge options for the merge command.

    Returns:
        Path of the merged DEMS file.

    """
    if debug:
        LOGGER.setLevel(DEBUG)

    basicConfig(
        datefmt="%Y-%m-%d %H:%M:%S",
        format="[%(asctime)s %(name)s %(levelname)s] %(message)s",
    )

    for key, val in locals().items():
        LOGGER.debug(f"{key}: {val!r}")

    data_dir_ = Path(data_dir).resolve() / f"cosmos_{obsid}"
    reduced_dir_ = Path(reduced_dir).resolve() / f"reduced_{obsid}"
    dems_dir_ = Path(dems_dir).resolve()
    ddb = Path(ddb).resolve()

    # Run reduce function
    readout = reduce.reduce(data_dir_, reduced_dir_, debug=debug)

    # Run merge function
    if (dems := dems_dir_ / f"dems_{obsid}.zarr.zip").exists():
        raise FileExistsError(dems)

    if not ddb.exists():
        raise FileNotFoundError(ddb)

    if not (corresp := data_dir_ / "kid_corresp.json").exists():
        raise FileNotFoundError(corresp)

    if not (obsinst := data_dir_ / f"{obsid}.obs").exists():
        raise FileNotFoundError(obsinst)

    if not (antenna := data_dir_ / f"{obsid}.ant").exists():
        raise FileNotFoundError(antenna)

    if not (skychop := data_dir_ / f"{obsid}.skychopper.dat.xz").exists():
        raise FileNotFoundError(skychop)

    if not (weather := data_dir_ / f"{obsid}.wea").exists():
        raise FileNotFoundError(weather)

    if not (misti := data_dir_ / f"{obsid}.misti").exists():
        misti = None

    if not (cabin := data_dir_ / f"{obsid}.cabin").exists():
        cabin = None

    merge_options = {
        "offset_time_antenna": DEFAULT_OFFSET_TIME_ANTENNA,
        **merge_options,
    }

    merge.merge(
        dems,
        ddb=ddb,
        corresp=corresp,
        readout=readout,
        obsinst=obsinst,
        antenna=antenna,
        skychop=skychop,
        weather=weather,
        misti=misti,
        cabin=cabin,
        debug=debug,
        **merge_options,
    )


def cli() -> None:
    """Command line interface of the demerge function."""
    Fire(demerge)

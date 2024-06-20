__all__ = ["demerge", "merge", "reduce"]
__version__ = "2.13.0"


# standard library
from pathlib import Path
from typing import Any


# dependencies
from fire import Fire
from . import merge, reduce


# constants
DEFAULT_CACHE_DIR = Path("cache").resolve()
DEFAULT_DATA_DIR = Path("data").resolve()
DEFAULT_DEMS_DIR = Path("dems").resolve()
DEFAULT_DDB_FILE = Path(__file__).parent / "data" / "ddb_20231123.fits.gz"
DEFAULT_OFFSET_TIME_ANTENNA = 20  # ms


def demerge(
    obsid: str,
    /,
    *,
    cache_dir: Path = DEFAULT_CACHE_DIR,
    data_dir: Path = DEFAULT_DATA_DIR,
    dems_dir: Path = DEFAULT_DEMS_DIR,
    ddb_file: Path = DEFAULT_DDB_FILE,
    **merge_options: Any,
) -> Path:
    """Run reduce and merge commands to create a single DEMS file.

    Args:
        obsid: Observation ID (e.g. YYYYmmddHHMMSS).
        cache_dir: Path of cache directory (e.g. ./cache).
        data_dir: Path of data directory (e.g. ./data).
        dems_dir: Path of DEMS directory (e.g. ./dems).
        ddb_file: Path of DDB (DESHIMA database) file.
            Defaults to the one shipped with the de:merge package.
        **merge_options: Other merge options for the merge command.

    Returns:
        Path of the merged DEMS file.

    """
    cache_dir_ = Path(cache_dir) / str(obsid)
    data_dir_ = Path(data_dir) / f"cosmos_{obsid}"
    dems_dir_ = Path(dems_dir) / str(obsid)

    # Run reduce function
    readout = reduce.reduce(data_dir_, cache_dir_)

    # Run merge function
    if (dems := dems_dir_ / f"dems_{obsid}.zarr.zip").exists():
        raise FileExistsError(dems)

    if not ddb_file.exists():
        raise FileNotFoundError(ddb_file)

    if not (corresp := data_dir_ / "kid_corresp.json").exists():
        raise FileNotFoundError(corresp)

    if not (obs := data_dir_ / f"{obsid}.obs").exists():
        raise FileNotFoundError(obs)

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
        ddb=ddb_file,
        corresp=corresp,
        readout=readout,
        obs=obs,
        antenna=antenna,
        skychop=skychop,
        weather=weather,
        misti=misti,
        cabin=cabin,
        **merge_options,
    )


def cli() -> None:
    """Command line interface of the demerge function."""
    Fire(demerge)

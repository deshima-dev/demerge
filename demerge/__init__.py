__all__ = ["demerge", "merge", "reduce"]
__version__ = "2024.7.0"


# standard library
from contextlib import contextmanager
from logging import DEBUG, basicConfig, getLogger
from pathlib import Path
from typing import Any, Literal


# dependencies
from fire import Fire
from . import merge, reduce


# constants
LOGGER = getLogger(__name__)
PACKAGE_DATA = Path(__file__).parent / "data"


@contextmanager
def set_logger(debug: bool):
    level = LOGGER.level

    if debug:
        LOGGER.setLevel(DEBUG)

    try:
        yield
    finally:
        LOGGER.setLevel(level)


def demerge(
    obsid: str,
    /,
    *,
    # data paths
    data_dir: Path = Path(),
    dems_dir: Path = Path(),
    reduced_dir: Path = Path(),
    ddb: Path = PACKAGE_DATA / "ddb_20240713.fits.gz",
    # merge options
    measure: Literal["df/f", "brightness"] = "df/f",
    overwrite: bool = False,
    debug: bool = False,
    **options: Any,
) -> Path:
    """Run reduce and merge commands to create a single DEMS.

    Args:
        obsid: Observation ID (YYYYmmddHHMMSS).
        data_dir: Path where raw data directory is placed,
            i.e. expecting ``${data_dir}/cosmos_YYYYmmddHHMMSS``.
        dems_dir: Path where merged DEMS file will be placed,
            i.e. expecting ``${dems_dir}/dems_YYYYmmddHHMMSS.zarr.zip``.
        reduced_dir: Path where reduced data directory will be placed,
            i.e. expecting ``${reduced_dir}/reduced_YYYYmmddHHMMSS``.
        ddb: Path of DDB (DESHIMA database) file.
        measure: Measure of the DEMS (either df/f or brightness).
        overwrite: If True, reduced data directory and merged DEMS file
            will be overwritten even if they exist.
        debug: If True, detailed logs for debugging will be printed.
        **options: Other merge options for the merge command.

    Returns:
        Path of the merged DEMS.

    """
    with set_logger(debug):
        for key, val in locals().items():
            LOGGER.debug(f"{key}: {val!r}")

    data_dir_ = Path(data_dir).resolve() / f"cosmos_{obsid}"
    reduced_dir_ = Path(reduced_dir).resolve() / f"reduced_{obsid}"
    dems_dir_ = Path(dems_dir).resolve()
    ddb = Path(ddb).resolve()

    # Run reduce function
    readout = reduce.reduce(
        data_dir=data_dir_,
        reduced_dir=reduced_dir_,
        overwrite=overwrite,
        debug=debug,
    )

    # Run merge function
    if (dems := dems_dir_ / f"dems_{obsid}.zarr.zip").exists():
        raise FileExistsError(dems)

    if not (corresp := data_dir_ / "kid_corresp.json").exists():
        raise FileNotFoundError(corresp)

    if not ddb.exists():
        raise FileNotFoundError(ddb)

    if not (obsinst := data_dir_ / f"{obsid}.obs").exists():
        raise FileNotFoundError(obsinst)

    if not (antenna := data_dir_ / f"{obsid}.ant").exists():
        antenna = None

    if not (cabin := data_dir_ / f"{obsid}.cabin").exists():
        cabin = None

    if not (misti := data_dir_ / f"{obsid}.misti").exists():
        misti = None

    if not (skychop := data_dir_ / f"{obsid}.skychopper.dat.xz").exists():
        skychop = None

    if not (weather := data_dir_ / f"{obsid}.wea").exists():
        weather = None

    return merge.merge(
        dems,
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
        # merge options
        measure=measure,
        overwrite=overwrite,
        debug=debug,
        **options,
    )


def cli() -> None:
    """Command line interface of the demerge function."""
    basicConfig(
        datefmt="%Y-%m-%d %H:%M:%S",
        format="[%(asctime)s %(name)s %(funcName)s %(levelname)s] %(message)s",
    )

    Fire(demerge)

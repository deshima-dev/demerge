__all__ = ["data", "demerge", "merge", "reduce"]
__version__ = "2024.9.0"


# standard library
from collections.abc import Iterator
from contextlib import contextmanager
from logging import DEBUG, basicConfig, getLogger
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Literal, Optional, Union


# dependencies
from fire import Fire
from . import data, merge, reduce


# type hints
PathLike = Union[Path, str]


# constants
LOGGER = getLogger(__name__)
PACKAGE_DATA = Path(__file__).parent / "data"


@contextmanager
def set_dir(dir: Optional[PathLike] = None, /) -> Iterator[Path]:
    """Resolve a directory or set a temporary directory."""
    if dir is None:
        with TemporaryDirectory() as temp_dir:
            yield Path(temp_dir)
    else:
        yield Path(dir).expanduser().resolve()


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


def demerge(
    obsid: str,
    /,
    *,
    # data paths
    data_dir: PathLike = Path(),
    dems_dir: PathLike = Path(),
    reduced_dir: Optional[Path] = None,
    ddb: PathLike = PACKAGE_DATA / "ddb_20240713.fits.gz",
    # merge options
    measure: Literal["df/f", "brightness"] = "df/f",
    overwrite: bool = False,
    debug: bool = False,
    **options: Any,
) -> Path:
    """Run reduce and merge commands to create a single DEMS.

    Args:
        obsid: Observation ID (YYYYmmddHHMMSS).
        data_dir: Path of directory where data packages are placed,
            i.e. expecting ``${data_dir}/cosmos_YYYYmmddHHMMSS``.
        dems_dir: Path of directory where merged DEMS will be placed,
            i.e. expecting ``${dems_dir}/dems_YYYYmmddHHMMSS.zarr.zip``.
        reduced_dir: Path of directory where reduced packages are placed,
            i.e. expecting ``${reduced_dir}/reduced_YYYYmmddHHMMSS``.
            If not specified, a temporary directory will be used.
        ddb: Path of DDB (DESHIMA database) file.
        measure: Measure of the DEMS (either df/f or brightness).
        overwrite: If True, the reduced package and the merged DEMS file
            will be overwritten even if they exist.
        debug: If True, detailed logs for debugging will be printed.
        **options: Other merge options for the merge command.

    Returns:
        Path of the merged DEMS.

    """
    with set_logger(debug):
        for key, val in locals().items():
            LOGGER.debug(f"{key}: {val!r}")

    with (
        set_dir(data_dir) as data_dir,
        set_dir(dems_dir) as dems_dir,
        set_dir(reduced_dir) as reduced_dir,
    ):
        data_pack = Path(data_dir).resolve() / f"cosmos_{obsid}"
        reduced_pack = reduced_dir / f"reduced_{obsid}"
        data_pack_ = data.parse(data_pack)

        # Run reduce function
        readout = reduce.reduce(
            data_pack=data_pack,
            reduced_pack=reduced_pack,
            overwrite=overwrite,
            debug=debug,
        )

        # Run merge function
        return merge.merge(
            dems_dir / f"dems_{obsid}.zarr.zip",
            # required datasets
            corresp=data_pack_.corresp,
            ddb=ddb,
            obsinst=data_pack_.obsinst,
            readout=readout,
            # optional datasets
            antenna=data_pack_.antenna,
            cabin=data_pack_.cabin,
            misti=data_pack_.misti,
            skychop=data_pack_.skychop,
            weather=data_pack_.weather,
            # merge options
            measure=measure,
            overwrite=overwrite,
            debug=debug,
            **options,
        )


def demerge_cli() -> None:
    """Command line interface of the demerge function."""
    basicConfig(
        datefmt="%Y-%m-%d %H:%M:%S",
        format="[%(asctime)s %(name)s %(funcName)s %(levelname)s] %(message)s",
    )

    Fire(demerge)

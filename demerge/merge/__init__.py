"""demsオブジェクトを生成する

python 3.9
dems   0.8.0

(C) 2023 内藤システムズ
"""

__all__ = ["merge"]


# standard library
from logging import DEBUG, basicConfig, getLogger
from pathlib import Path


# dependencies
from fire import Fire
from .utils import create_dems


# constants
LOGGER = getLogger(__name__)


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
    findR: bool = False,
    ch: int = 0,
    Rth: float = 280.0,
    skyth: float = 150.0,
    cutnum: int = 1,
    still: bool = False,
    period: int = 2,
    shuttle: bool = False,
    lon_min_off: float = 0.0,
    lon_max_off: float = 0.0,
    lon_min_on: float = 0.0,
    lon_max_on: float = 0.0,
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
        findR: If specified, execute FindR and Sky.
        ch: Channel to use for FindR as an integer.
        Rth: R threshold for FindR as a float.
        skyth: Sky threshold for FindR as a float.
        cutnum: Cut number for FindR as an integer.
        still: If specified, perform analysis for still observation.
        period: Half period (seconds) for still observation as an integer.
        shuttle: If specified, perform analysis for shuttle observation.
        lon_min_off: Minimum longitude of OFF during shuttle observation.
        lon_max_off: Maximum longitude of OFF during shuttle observation.
        lon_min_on: Minimum longitude of ON during shuttle observation.
        lon_max_on: Maximum longitude of ON during shuttle observation.
        debug: If True, detailed logs for debugging will be printed.
        offset_time_antenna: Time diff (ms) between reduced FITS and antenna log.

    Returns:
        Path of the merged DEMS file.

    Raises:
        FileExistsError: Raised if ``dems`` exists.

    """
    # ロガーの設定
    if debug:
        LOGGER.setLevel(DEBUG)

    basicConfig(
        datefmt="%Y-%m-%d %H:%M:%S",
        format="[%(asctime)s %(name)s %(levelname)s] %(message)s",
    )

    # 引数と値をロガーに記録
    for key, val in locals().items():
        LOGGER.debug(f"{key}: {val!r}")

    # マージの実行
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
        findR=findR,
        ch=ch,
        Rth=Rth,
        skyth=skyth,
        cutnum=cutnum,
        still=still,
        period=period,
        shuttle=shuttle,
        lon_min_off=lon_min_off,
        lon_max_off=lon_max_off,
        lon_min_on=lon_min_on,
        lon_max_on=lon_max_on,
        offset_time_antenna=offset_time_antenna,
    )

    if (dems := Path(dems)).exists():
        raise FileExistsError(dems)

    dems.parent.mkdir(exist_ok=True, parents=True)
    da.to_zarr(dems, mode="w")
    return dems.resolve()


def cli() -> None:
    """Command line interface of the merge function."""
    Fire(merge)

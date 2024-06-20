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
    filename: str,
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
    loadtype: str = "fshift",
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
        filename: 出力ファイルへのパス (.zarr.zip)
        ddb: DDBファイルへのパス (.fits or .fits.gz)
        corresp: Master-to-KID ID対応ファイルへのパス (.json)
        obsinst: 指示書ファイルへのパス (.obs)
        antenna: アンテナログファイルへのパス (.ant)
        readout: Reduced FITSファイルへのパス (.fits)
        skychop: Sky chopperファイルへのパス (.skychop)
        weather: 気象ログファイルへのパス (.weather)
        misti: MiSTIログファイルへのパス (.misti)
        cabin: キャビンログファイルへのパス (.cabin)
        coordinate: 出力データの座標系 (azel or radec)
        loadtype: 出力データ形式 (fshift or Tsignal)
        findR: 指定するとFindR, Skyを実行します
        ch: FindRの利用するチャネルを整数で指定します
        Rth: FindRのR閾値を実数で指定します
        skyth: FindRのsky閾値を実数で指定します
        cutnum: FindRのカット数を整数で指定します
        still: 指定するとstill観測用の解析を行います
        period: Still観測の1/2周期(秒)を整数で指定します
        shuttle: 指定するとshuttle観測用の解析を行います
        lon_min_off: Shuttle観測時のOFFにするlongitudeの最小値
        lon_max_off: Shuttle観測時のOFFにするlongitudeの最大値
        lon_min_on: Shuttle観測時のONにするlongitudeの最小値
        lon_max_on: Shuttle観測時のONにするlongitudeの最大値
        debug: 指定すると全ての引数の値をログとして表示します
        offset_time_antenna: Reduced FITSとアンテナログの時刻のずれの補正値 (ms)

    Returns:
        Path of the merged DEMS file.

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
    dems = create_dems(
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
        loadtype=loadtype,
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

    # DEMSをZarr形式で保存
    dems.to_zarr(filename, mode="w")
    return Path(filename).resolve()


def cli() -> None:
    """Command line interface of the merge function."""
    Fire(merge)

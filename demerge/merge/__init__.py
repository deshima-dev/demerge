"""demsオブジェクトを生成する

python 3.9
dems   0.8.0

(C) 2023 内藤システムズ
"""

__all__ = ["merge"]


# standard library
import argparse
from logging import DEBUG, basicConfig, getLogger


# dependencies
from .utils import create_dems


# module logger
logger = getLogger(__name__)


def cli() -> None:
    """Demsオブジェクトを作成する"""
    parser = argparse.ArgumentParser()

    # 必須引数
    parser.add_argument(
        "filename",
        type=str,
        help="出力ファイルへのパスを指定して下さい(.zarr.zip)",
    )
    parser.add_argument(
        "--ddb",
        type=str,
        required=True,
        help="DDBファイルへのパスを指定して下さい(.fits.gz)",
    )
    parser.add_argument(
        "--corresp",
        type=str,
        required=True,
        help="Master-to-KID ID対応ファイルへのパスを指定して下さい(.json)",
    )
    parser.add_argument(
        "--obs",
        type=str,
        required=True,
        help="obsファイルへのパスを指定して下さい(.obs)",
    )
    parser.add_argument(
        "--antenna",
        type=str,
        required=True,
        help="antennaファイルへのパスを指定して下さい(.antenna)",
    )
    parser.add_argument(
        "--readout",
        type=str,
        required=True,
        help="reduced readoutファイルへのパスを指定して下さい(.fits)",
    )
    parser.add_argument(
        "--skychop",
        type=str,
        required=True,
        help="skychopファイルへのパスを指定して下さい(.skychop)",
    )
    parser.add_argument(
        "--weather",
        type=str,
        required=True,
        help="weatherファイルへのパスを指定して下さい(.weather)",
    )

    # オプション引数
    parser.add_argument(
        "--misti",
        type=str,
        default="",
        help="mistiファイルへのパスを指定して下さい(.misti)",
    )
    parser.add_argument(
        "--cabin",
        type=str,
        default="",
        help="cabinファイルへのパスを指定して下さい(.cabin)",
    )
    parser.add_argument(
        "--coordinate",
        type=str,
        default="azel",
        help="座標系(azel/radec)を文字列で指定します",
    )
    parser.add_argument(
        "--loadtype",
        type=str,
        default="fshift",
        help="読み込むデータを文字列で指定します(既定値: fshift, fshiftかTsignalを指定できます)",
    )
    parser.add_argument(
        "--findR",
        action="store_true",
        help="指定するとFindR, Skyを実行します",
    )
    parser.add_argument(
        "--ch",
        type=int,
        default=0,
        help="findRに利用するチャネルを整数で指定します",
    )
    parser.add_argument(
        "--Rth",
        type=float,
        default=280.0,
        help="R閾値を実数で指定します",
    )
    parser.add_argument(
        "--skyth",
        type=float,
        default=150.0,
        help="sky閾値を実数で指定します",
    )
    parser.add_argument(
        "--cutnum",
        type=int,
        default=1,
        help="findRでのカット数を整数で指定します",
    )
    parser.add_argument(
        "--still",
        action="store_true",
        help="指定するとstill観測用の解析を行います",
    )
    parser.add_argument(
        "--period",
        type=int,
        default=2,
        help="still観測の1/2周期(秒)を整数で指定します",
    )
    parser.add_argument(
        "--shuttle",
        action="store_true",
        help="指定するとshuttle観測用の解析を行います",
    )
    parser.add_argument(
        "--lon_min_off",
        type=float,
        default=0.0,
        help="shuttle観測時のOFFにするlongitudeの最小値を実数で指定します",
    )
    parser.add_argument(
        "--lon_max_off",
        type=float,
        default=0.0,
        help="shuttle観測時のOFFにするlongitudeの最大値を実数で指定します",
    )
    parser.add_argument(
        "--lon_min_on",
        type=float,
        default=0.0,
        help="shuttle観測時のONにするlongitudeの最小値を実数で指定します",
    )
    parser.add_argument(
        "--lon_max_on",
        type=float,
        default=0.0,
        help="shuttle観測時のONにするlongitudeの最大値を実数で指定します",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="指定すると全ての引数の値をログとして表示します",
    )
    parser.add_argument(
        "--offset_time_antenna",
        type=int,
        default=0,
        help="TODとAntennaログの時刻のずれの補正値(ms)",
    )

    # 引数の読み取り
    a = parser.parse_args()

    # ロガーの設定
    if a.debug:
        logger.setLevel(DEBUG)

    basicConfig(
        datefmt="%Y-%m-%d %H:%M:%S",
        format="[%(asctime)s %(name)s %(levelname)s] %(message)s",
    )

    # 引数と値をロガーに記録
    for key, val in vars(a).items():
        logger.debug(f"{key}: {val!r}")

    # マージの実行
    dems = create_dems(
        ddbfits_path=a.ddb,
        corresp_path=a.corresp,
        obsinst_path=a.obs,
        antenna_path=a.antenna,
        readout_path=a.readout,
        skychop_path=a.skychop,
        weather_path=a.weather,
        misti_path=a.misti,
        cabin_path=a.cabin,
        coordinate=a.coordinate,
        loadtype=a.loadtype,
        findR=a.findR,
        ch=a.ch,
        Rth=a.Rth,
        skyth=a.skyth,
        cutnum=a.cutnum,
        still=a.still,
        period=a.period,
        shuttle=a.shuttle,
        lon_min_off=a.lon_min_off,
        lon_max_off=a.lon_max_off,
        lon_min_on=a.lon_min_on,
        lon_max_on=a.lon_max_on,
        offset_time_antenna=a.offset_time_antenna,
    )

    dems.to_zarr(a.filename, mode="w")

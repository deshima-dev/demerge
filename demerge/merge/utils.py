__all__ = ["create_dems"]


# standard library
import json
import lzma
import re
from datetime import datetime as dt
from functools import partial, reduce
from pathlib import Path
from typing import Any, Literal


# dependencies
import numpy as np
import pandas as pd
import xarray as xr
from astropy.io import ascii, fits
from dems.d2 import MS
from numpy.typing import NDArray
from .. import __version__ as DEMERGE_VERSION


# constants
COLUMN_NAMES_ANTENNA = (
    "time",
    "ra-prg",  # deg
    "dec-prg",  # deg
    "az-prg",  # deg
    "el-prg",  # deg
    "az-real",  # deg
    "el-real",  # deg
    "x",  # mm
    "y",  # mm
    "z",  # mm
    "xt",  # deg
    "yt",  # deg
    "zt",  # deg
    "lst",
    "az-prg(no-cor)",  # deg
    "el-prog(no-cor)",  # deg
    "az-prog(center)",  # deg
    "el-prog(center)",  # deg
    "type",
)
COLUMN_NAMES_CABIN = (
    "date",
    "time",
    "upper_cabin_temperature",  # degC
    "main_cabin_temperature",  # degC
    *[f"_{i}" for i in range(18)],
)
COLUMN_NAMES_MISTI = (
    "date",
    "time",
    "unix_time",
    "az",  # deg
    "el",  # deg
    "pwv",  # um
    "ground_temperature",  # K
)
COLUMN_NAMES_SKYCHOP = (
    "time",
    "is_blocking",
)
COLUMN_NAMES_WEATHER = (
    "time",
    "temperature",  # degC
    "pressure",  # hPa
    "humidity",  # %
    "wind_speed",  # m/s
    "wind_direction",  # deg
    "_",
)
DATE_PARSER_ANTENNA = lambda s: dt.strptime(s, "%Y%m%d%H%M%S.%f")
DATE_PARSER_CABIN = lambda s: dt.strptime(s, "%Y/%m/%d %H:%M")
DATE_PARSER_MISTI = lambda s: dt.strptime(s, "%Y/%m/%d %H:%M:%S.%f")
DATE_PARSER_SKYCHOP = lambda s: dt.fromtimestamp(float(s))
DATE_PARSER_WEATHER = lambda s: dt.strptime(s, "%Y%m%d%H%M%S")


def get_antenna(antenna: Path, /) -> xr.Dataset:
    """Load an antenna log as xarray Dataset."""
    return pd.read_csv(
        antenna,
        # read settings
        names=COLUMN_NAMES_ANTENNA,
        delimiter="\s+",
        comment="#",
        # index settings
        index_col=0,
        parse_dates=[0],
        date_parser=DATE_PARSER_ANTENNA,
    ).to_xarray()


def get_cabin(cabin: Path, /) -> xr.Dataset:
    """Load a cabin log as xarray Dataset."""
    return (
        pd.read_csv(
            cabin,
            # read settings
            names=COLUMN_NAMES_CABIN,
            delimiter="\s+",
            comment="#",
            # index settings
            index_col=[0],
            parse_dates=[[0, 1]],
            date_parser=DATE_PARSER_CABIN,
        )
        .rename_axis("time")
        .to_xarray()
    )


def get_misti(misti: Path, /) -> xr.Dataset:
    """Load a MiSTI log as xarray Dataset."""
    return (
        pd.read_csv(
            misti,
            # read settings
            names=COLUMN_NAMES_MISTI,
            delimiter="\s+",
            comment="#",
            # index settings
            index_col=[0],
            parse_dates=[[0, 1]],
            date_parser=DATE_PARSER_MISTI,
        )
        .rename_axis("time")
        .to_xarray()
    )


def get_obstable(obstable: Path, /) -> dict[str, str]:
    """Load an observation table to get parameters."""
    with open(obstable) as f:
        lines = f.read()

    def search(pattern: str) -> str:
        if (match := re.search(pattern, lines)) is None:
            return ""
        else:
            return match[1]

    return {
        # DES
        "obs_file": search("SET DES OBS_FILE\s*'(.*)'"),
        "obs_user": search("SET DES OBS_USER\s*'(.*)'"),
        "project": search("SET DES PROJECT\s*'(.*)'"),
        # ANTENNA_G
        "epoch": search("SET ANTENNA_G EPOCH\s*'(.*)'"),
        "scan_cood": search("SET ANTENNA_G SCAN_COOD\s*'(.*)'"),
        "src_name": search("SET ANTENNA_G SRC_NAME\s*'(.*)'"),
        "src_pos": search("SET ANTENNA_G SRC_POS\s*\((.*)\)"),
    }


def get_readout(readout: Path, /) -> xr.DataArray:
    """Load a reduced readout FITS as xarray DataArray."""
    with fits.open(readout) as hdus:
        readout_data = hdus["READOUT"].data
        kidcols = readout_data.columns[2:].names
        time = pd.to_datetime(readout_data["timestamp"], unit="s")
        linph = np.array([readout_data[col] for col in kidcols]).T[1]

        kidsinfo_data = hdus["KIDSINFO"].data
        linyfc = kidsinfo_data["yfc, linyfc"].T[1]
        Qr = kidsinfo_data["Qr, dQr (Sky)"].T[0]
        fr = kidsinfo_data["fr, dfr (Sky)"].T[0]
        fr_room = kidsinfo_data["fr, dfr (Room)"].T[0]

    if np.isnan(fr_room).all():
        dfof = (linph - linyfc) / (4.0 * Qr)
    else:
        dfof = (linph - linyfc) / (4.0 * Qr) - (fr - fr_room) / fr

    return xr.DataArray(
        dfof,
        name="df/f",
        dims=("time", "kidid"),
        coords={
            "time": time,
            "kidid": np.arange(dfof.shape[1]),
        },
    )


def get_skychop(skychop: Path, /) -> xr.Dataset:
    """Load a sky chopper log as xarray Dataset."""
    return pd.read_csv(
        skychop,
        # read settings
        names=COLUMN_NAMES_SKYCHOP,
        delimiter="\s+",
        comment="#",
        # index settings
        index_col=0,
        parse_dates=[0],
        date_parser=DATE_PARSER_SKYCHOP,
    ).to_xarray()


def get_weather(weather: Path, /) -> xr.Dataset:
    """Load a weather log as xarray Dataset."""
    return pd.read_csv(
        weather,
        # read settings
        names=COLUMN_NAMES_WEATHER,
        delimiter="\s+",
        comment="#",
        # index settings
        index_col=0,
        parse_dates=[0],
        date_parser=DATE_PARSER_WEATHER,
    ).to_xarray()


def load_obsinst(obsinst: Path) -> dict[str, Any]:
    """Get various values from an obsinst file."""

    # default values
    epoch = 2000.0
    obs = ""
    obs_user = ""
    project = ""
    src_name = ""
    src_pos = ["0.0", "0.0"]
    trk_type = ""

    with open(obsinst, "r") as f:
        for line in f:
            if "SET ANTENNA_G EPOCH" in line:
                epoch = line.split()[-1].strip("'JB")
            elif "SET DES OBS_FILE" in line:
                obs_file = line.split()[-1].strip("'")
            elif "SET DES OBS_USER" in line:
                obs_user = line.split()[-1].strip("'")
            elif "SET DES PROJECT" in line:
                project = line.split()[-1].strip("'")
            elif "SET ANTENNA_G SRC_NAME" in line:
                src_name = line.split()[-1].strip("'")
            elif "SET ANTENNA_G SRC_POS" in line:
                src_pos = line.split()[-1].strip("()").split(",")
            elif "SET ANTENNA_G TRK_TYPE" in line:
                trk_type = line.split()[-1].strip("'")

    if trk_type == "RADEC":
        ra, dec = src_pos
    else:
        ra, dec = 0.0, 0.0

    return {
        "observation": obs_file,
        "observer": obs_user,
        "project": project,
        "object": src_name,
        "equinox": float(epoch),
        "ra": float(ra),
        "dec": float(dec),
    }


def get_corresp_frame(ddb: fits.HDUList, corresp_file: str) -> pd.DataFrame:
    """Get correspondence between KID ID and each KID attribute.

    Args:
        ddb: DESHIMA database (DDB) as an HDU list.
        corresp_file: Master-to-KID ID correspondence as JSON file.

    Returns:
        DataFrame of correspondence between KID ID and each KID attribute.

    """

    def native(array: NDArray[Any]) -> NDArray[Any]:
        """Convert the byte order of an array to native."""
        return array.astype(array.dtype.type)

    frames: list[pd.DataFrame] = []

    # DataFrame of KIDDES HDU
    data = ddb["KIDDES"].data
    frame = pd.DataFrame(
        index=pd.Index(
            native(data["masterid"]),
            name="masterid",
        ),
        data={
            "kidtype": native(data["attribute"]),
        },
    )
    frames.append(frame)

    # DataFrame of KIDFILT HDU
    data = ddb["KIDFILT"].data
    frame = pd.DataFrame(
        index=pd.Index(
            data=native(data["masterid"]),
            name="masterid",
        ),
        data={
            "kidfreq": native(data["F_filter, df_filter"][:, 0]),
            "kidQ": native(data["Q_filter, dQ_filter"][:, 0]),
        },
    )
    frame["kidfreq"] *= 1e9
    frames.append(frame)

    # DataFrame of KIDRESP HDU
    if "KIDRESP" in ddb:
        data = ddb["KIDRESP"].data
        frame = pd.DataFrame(
            index=pd.Index(
                data=native(data["masterid"]),
                name="masterid",
            ),
            data={
                "p0": native(data["cal params"][:, 0]),
                "etaf": native(data["cal params"][:, 1]),
                "T0": native(data["cal params"][:, 2]),
            },
        )
        frames.append(frame)

    # Outer-join DataFrames
    join = partial(pd.DataFrame.join, how="outer")
    frame = reduce(join, frames)

    # Assign KID ID as index
    with open(corresp_file, mode="r") as f:
        corresp = json.load(f)

    index = pd.Index(
        data=[corresp.get(str(i), -1) for i in frame.index],
        name="kidid",
    )
    frame = frame.reset_index().set_index(index)

    # Drop rows with invalid master or KID IDs (-1)
    frame = frame[frame.index != -1]
    frame = frame[frame.masterid != -1]
    return frame


def convert_readout(
    readout: fits.HDUList,
    corresp: pd.DataFrame,
    to: Literal["brightness", "df/f"],
    T_room: float,
    T_amb: float,
):
    """Reduced readoutの値をDEMS出力形式に変換（校正）する

    Args:
        readout: Reduced readout FITSのHDUListオブジェクト
        corresp: KID IDと各KIDの測定値を対応づけるDataFrame
        to: 変換形式（brightness or df/f)
        T_room: キャビン温度(K)
        T_amb: 外気温(K)

    """
    kidcols = readout["READOUT"].data.columns[2:].names
    linph = np.array([readout["READOUT"].data[n] for n in kidcols]).T[1]
    linyfc = np.array(readout["KIDSINFO"].data["yfc, linyfc"]).T[1]
    Qr = np.array(readout["KIDSINFO"].data["Qr, dQr (Sky)"]).T[0]
    fr = np.array(readout["KIDSINFO"].data["fr, dfr (Sky)"]).T[0]
    fr_room = np.array(readout["KIDSINFO"].data["fr, dfr (Room)"]).T[0]

    if np.isnan(fr_room).all():
        fshift = (linph - linyfc) / (4.0 * Qr)
    else:
        fshift = (linph - linyfc) / (4.0 * Qr) - (fr - fr_room) / fr

    if to == "df/f":
        return fshift[:, corresp.index.values]

    if to == "brightness":
        return Tlos_model(
            dx=fshift[:, corresp.index.values],
            p0=corresp.p0.values,
            etaf=corresp.etaf.values,
            T0=corresp.T0.values,
            Troom=T_room,
            Tamb=T_amb,
        )

    raise ValueError(f"Invalid output type: {to}")


def Tlos_model(dx, p0, etaf, T0, Troom, Tamb):
    """Calibrate 'amplitude' and 'phase' to 'power'"""
    return (
        (dx + p0 * np.sqrt(Troom + T0)) ** 2 / (p0**2 * etaf)
        - T0 / etaf
        - (1 - etaf) / etaf * Tamb
    )


def convert_asciitime(asciitime, form_fitstime):
    """Ascii time"""
    asciitime = [dt.strptime("%14.6f" % t, "%Y%m%d%H%M%S.%f") for t in asciitime]
    asciitime = [dt.strftime(t, form_fitstime) for t in asciitime]
    return np.array(asciitime)


def convert_timestamp(timestamp):
    """Timestamp"""
    timestamp = [dt.utcfromtimestamp(t) for t in timestamp]
    timestamp = [dt.strftime(t, "%Y-%m-%dT%H:%M:%S.%f") for t in timestamp]
    return np.array(timestamp)


def retrieve_cabin_temps(filename=None):
    """キャビン内温度を取得する

    Args:
        str ファイル名

    Keyword Args:
        なし

    Returns:
        tuple (timestames, upperCabinTemps, lowerCabinTemps)
               tupleの各要素はnumpy.array。要素数は同じ。
               また、ファイル名が空の場合はデフォルト値が入った配列が返される。
    """
    if filename == "" or filename is None:
        return (
            np.array(["1970-01-01"]).astype("datetime64[ns]"),
            np.array([20.0]).astype(np.float64),
            np.array([20.0]).astype(np.float64),
        )

    table = ascii.read(filename, format="no_header")

    # 日付と時刻を取得して文字列でタイムスタンプを作成しそれをnumpy.datetime64へ変換する
    # テーブルの1列目と2列目がそれぞれ日付と時刻
    datetimes = []
    for date, time in zip(table["col1"], table["col2"]):
        s = "{}T{}".format(date, time)
        s = s.replace("/", "-")
        datetimes.append(s)

    datetimes = np.array(datetimes).astype("datetime64[ns]")
    upper_cabin_temps = np.array(table["col3"]).astype(np.float64)
    lower_cabin_temps = np.array(table["col4"]).astype(np.float64)

    return (datetimes, upper_cabin_temps, lower_cabin_temps)


def retrieve_skychop_states(
    skychop: Path,
) -> tuple[NDArray[np.float64], NDArray[np.int8]]:
    """skychopファイル(text file)からskychopの時系列状態を取得する

    Args:
        str ファイル名(xzで圧縮されていても読める)

    Returns:
        tuple (timestames, states)
          tupleの各要素はnumpy.array。要素数は同じ。

    説明:
        時刻について
        ============
        skychopファイルに記録されている時刻はUNIX時間。

        ファイル形式
        ============
        1列目 UNIX時刻
        2列目 0/1による状態
        "#"から始まるコメントがファイル冒頭に数行ある。
    """
    if Path(skychop).suffix == ".xz":
        with lzma.open(skychop, "rt") as f:
            data = f.read()
    else:
        with open(skychop, "rt") as f:
            data = f.read()

    table = ascii.read(
        data,
        guess=False,
        format="no_header",
        delimiter=" ",
        names=["datetime", "state"],
    )
    datetimes = np.array(table["datetime"]).astype(np.float64)
    states = np.array(table["state"]).astype(np.int8)
    return (datetimes, states)


def retrieve_misti_log(filename):
    """mistiファイルからの時系列データを取得する

    Args:
        str ファイル名

    Returns:
        tuple (timestames, az, el, pwv)
            tupleの各要素はnumpy.array。要素数は同じ。
            また、ファイル名が空の場合はNaNが1つだけ入った配列を返す。


    説明:
        ファイル形式:
        1列目 年/月/日
        2列目 時:分:6列目 秒(小数点以下2桁も含む)
        3列目 UNIX時間
        4列目 az(deg)
        5列目 el(deg)
        6列目 pwv(um)
        7列目 Tground(K)
        "#"から始まるコメントがファイル冒頭に数行ある。

    """
    if filename == "" or filename is None:
        return (
            np.array([np.nan]).astype("datetime64[ns]"),
            np.array([np.nan]).astype(np.float64),
            np.array([np.nan]).astype(np.float64),
            np.array([np.nan]).astype(np.float64),
        )

    column_names = [
        "date",
        "time",
        "unixtime",
        "az",
        "el",
        "pwv",
        "Tround",
    ]
    table = ascii.read(
        filename,
        guess=False,
        format="no_header",
        delimiter=" ",
        names=column_names,
    )

    az = np.array(table["az"]).astype(np.float64)
    el = np.array(table["el"]).astype(np.float64)
    pwv = np.array(table["pwv"]).astype(np.float64) / 1000.0  # umからmmへ変換

    datetimes = []
    for row in table:
        datetimes.append(
            dt.strptime(
                "{} {}".format(row["date"], row["time"]),
                format="%Y/%m/%d %H:%M:%S.%f",
            )
        )

    return np.array(datetimes).astype("datetime64[ns]"), az, el, pwv


def create_dems(
    ddbfits_path: str,
    corresp_path: str,
    obsinst_path: str,
    antenna_path: str,
    readout_path: str,
    skychop_path: str,
    weather_path: str,
    misti_path: str,
    cabin_path: str,
    coordinate: str,
    measure: str,
    offset_time_antenna: int,
):
    # 時刻と各種データを読み込む
    readout_hdul = fits.open(readout_path, mode="readonly")
    ddbfits_hdul = fits.open(ddbfits_path, mode="readonly")
    weather_table = ascii.read(weather_path)
    weather_table["tmperature"] += 273.15
    # 最後の1行は終端を表す意味のないデータが入っているため無視する
    antenna_table = ascii.read(antenna_path)[:-1]
    # 観測スクリプトに含まれているパラメタを抽出する
    obsinst_params = load_obsinst(obsinst_path)

    # 必要に応じて時刻はnp.datetime64[ns]へ変換する
    times = convert_timestamp(readout_hdul["READOUT"].data["timestamp"])
    times = np.array(times).astype("datetime64[ns]")

    times_misti, az_misti, el_misti, pwv_misti = retrieve_misti_log(misti_path)
    times_cabin, _, lower_cabin_temp = retrieve_cabin_temps(cabin_path)
    lower_cabin_temp = lower_cabin_temp + 273.15  # 度CからKへ変換

    times_weather = convert_asciitime(
        asciitime=weather_table["time"],
        form_fitstime="%Y-%m-%dT%H:%M:%S.%f",
    )
    times_weather = np.array(times_weather).astype("datetime64[ns]")

    times_skychop, states_skychop = retrieve_skychop_states(skychop_path)
    times_skychop = convert_timestamp(times_skychop)
    times_skychop = np.array(times_skychop).astype("datetime64[ns]")

    times_antenna = convert_asciitime(
        asciitime=antenna_table["time"],
        form_fitstime="%Y-%m-%dT%H:%M:%S.%f",
    )
    # fmt: off
    times_antenna = (
        np.array(times_antenna).astype("datetime64[ns]")
        + np.timedelta64(offset_time_antenna, "ms")
    )
    # fmt: on

    ddb_version = ddbfits_hdul["PRIMARY"].header["DDB_ID"]

    corresp = get_corresp_frame(ddbfits_hdul, corresp_path)
    response = convert_readout(
        readout=readout_hdul,
        corresp=corresp,
        to=measure,
        T_room=lower_cabin_temp[0],
        T_amb=np.nanmean(weather_table["tmperature"]),
    )

    if measure == "brightness":
        long_name = "Brightness"
        units = "K"
    elif measure == "df/f":
        long_name = "df/f"
        units = "dimensionless"
    else:
        raise KeyError("Invalid measure: {}".format(measure))

    ddbfits_hdul.close()
    readout_hdul.close()

    # モードに応じて経度(lon)と緯度(lat)を選択(azelかradecか)する
    if coordinate == "azel":
        if "az-prg(no-cor)" in antenna_table.colnames:
            az_prg = antenna_table["az-prg(no-cor)"]
            el_prog = antenna_table["el-prog(no-cor)"]
        elif "az-prg(no-col)" in antenna_table.colnames:
            az_prg = antenna_table["az-prg(no-col)"]
            el_prog = antenna_table["el-prog(no-col)"]
        else:
            raise KeyError(
                str(antenna_path)
                + "ファイルにaz-prg(no-cor)列またはaz-prg(no-col)列がありません。"
            )

        lon = az_prg + antenna_table["az-real"] - antenna_table["az-prg"]
        lat = el_prog + antenna_table["el-real"] - antenna_table["el-prg"]
        lon_origin = antenna_table["az-prog(center)"]
        lat_origin = antenna_table["el-prog(center)"]
    elif coordinate == "radec":
        lon = antenna_table["ra-prg"]
        lat = antenna_table["dec-prg"]
        # 観測スクリプトに設定されているRA,DEC
        lon_origin = np.full_like(lon, obsinst_params["ra"])
        lat_origin = np.full_like(lat, obsinst_params["dec"])
    else:
        raise KeyError("Invalid coodinate type: {}".format(coordinate))

    # 補間関数で扱うためにSCANTYPE(文字列)を適当な整数に対応させる
    states = np.array(antenna_table["type"])
    state_types = {state_type: i for i, state_type in enumerate(np.unique(states))}
    state_type_numbers = np.zeros(states.shape[0], dtype=int)
    for state_type, i in state_types.items():
        state_type_numbers[states == state_type] = i

    # 補間のためにDataArrayへ格納する
    response_xr = xr.DataArray(
        data=response,
        dims=["time", "chan"],
        coords=[times, corresp.index],
    )
    lon_xr = xr.DataArray(
        data=lon,
        coords={"time": times_antenna},
    )
    lat_xr = xr.DataArray(
        data=lat,
        coords={"time": times_antenna},
    )
    lon_origin_xr = xr.DataArray(
        data=lon_origin,
        coords={"time": times_antenna},
    )
    lat_origin_xr = xr.DataArray(
        data=lat_origin,
        coords={"time": times_antenna},
    )
    temperature_xr = xr.DataArray(
        data=weather_table["tmperature"],
        coords={"time": times_weather},
    )
    humidity_xr = xr.DataArray(
        data=weather_table["vapor-pressure"],
        coords={"time": times_weather},
    )
    pressure_xr = xr.DataArray(
        data=weather_table["presure"],
        coords={"time": times_weather},
    )
    wind_speed_xr = xr.DataArray(
        data=weather_table["aux1"],
        coords={"time": times_weather},
    )
    wind_direction_xr = xr.DataArray(
        data=weather_table["aux2"],
        coords={"time": times_weather},
    )
    skychop_state_xr = xr.DataArray(
        data=states_skychop,
        coords={"time": times_skychop},
    )
    aste_cabin_temperature_xr = xr.DataArray(
        data=lower_cabin_temp,
        coords={"time": times_cabin},
    )
    aste_subref_x_xr = xr.DataArray(
        data=antenna_table["x"],
        coords={"time": times_antenna},
    )
    aste_subref_y_xr = xr.DataArray(
        data=antenna_table["y"],
        coords={"time": times_antenna},
    )
    aste_subref_z_xr = xr.DataArray(
        data=antenna_table["z"],
        coords={"time": times_antenna},
    )
    aste_subref_xt_xr = xr.DataArray(
        data=antenna_table["xt"],
        coords={"time": times_antenna},
    )
    aste_subref_yt_xr = xr.DataArray(
        data=antenna_table["yt"],
        coords={"time": times_antenna},
    )
    aste_subref_zt_xr = xr.DataArray(
        data=antenna_table["zt"],
        coords={"time": times_antenna},
    )
    aste_misti_lon_xr = xr.DataArray(
        data=az_misti,
        coords={"time": times_misti},
    )
    aste_misti_lat_xr = xr.DataArray(
        data=el_misti,
        coords={"time": times_misti},
    )
    aste_misti_pwv_xr = xr.DataArray(
        data=pwv_misti,
        coords={"time": times_misti},
    )
    state_type_numbers_xr = xr.DataArray(
        data=state_type_numbers,
        coords={"time": times_antenna},
    )

    # Tsignalsの時刻に合わせて補間する
    lon = lon_xr.interp_like(response_xr)
    lat = lat_xr.interp_like(response_xr)
    lon_origin = lon_origin_xr.interp_like(response_xr)
    lat_origin = lat_origin_xr.interp_like(response_xr)
    temperature = temperature_xr.interp_like(response_xr)
    humidity = humidity_xr.interp_like(response_xr)
    pressure = pressure_xr.interp_like(response_xr)
    wind_speed = wind_speed_xr.interp_like(response_xr)
    wind_direction = wind_direction_xr.interp_like(response_xr)
    aste_subref_x = aste_subref_x_xr.interp_like(response_xr)
    aste_subref_y = aste_subref_y_xr.interp_like(response_xr)
    aste_subref_z = aste_subref_z_xr.interp_like(response_xr)
    aste_subref_xt = aste_subref_xt_xr.interp_like(response_xr)
    aste_subref_yt = aste_subref_yt_xr.interp_like(response_xr)
    aste_subref_zt = aste_subref_zt_xr.interp_like(response_xr)
    skychop_state = skychop_state_xr.interp_like(
        response_xr,
        method="nearest",
    )
    state_type_numbers = state_type_numbers_xr.interp_like(
        response_xr,
        method="nearest",
    )

    aste_cabin_temperature = np.nan
    aste_misti_lon = np.nan
    aste_misti_lat = np.nan
    aste_misti_pwv = np.nan

    if cabin_path != "" and cabin_path is not None:
        aste_cabin_temperature = aste_cabin_temperature_xr.interp_like(response_xr)

    if misti_path != "" and misti_path is not None:
        aste_misti_lon = aste_misti_lon_xr.interp_like(response_xr)
        aste_misti_lat = aste_misti_lat_xr.interp_like(response_xr)
        aste_misti_pwv = aste_misti_pwv_xr.interp_like(response_xr)

    # 補間後のSTATETYPEを文字列に戻す
    state = np.full_like(state_type_numbers, "GRAD", dtype="<U8")
    for state_type, i in state_types.items():
        state[state_type_numbers == i] = state_type

    # Sky chopperの状態からビームラベルを割り当てる(1 -> B, 0 -> A)
    beam = np.where(skychop_state, "B", "A")

    return MS.new(
        data=response,
        long_name=long_name,
        units=units,
        time=times,
        chan=corresp.masterid,
        beam=beam,
        state=state,
        lon=lon,
        lat=lat,
        lon_origin=lon_origin,
        lat_origin=lat_origin,
        temperature=temperature,
        pressure=pressure,
        humidity=humidity,
        wind_speed=wind_speed,
        wind_direction=wind_direction,
        frequency=corresp.kidfreq,
        aste_cabin_temperature=aste_cabin_temperature,
        aste_subref_x=aste_subref_x,
        aste_subref_y=aste_subref_y,
        aste_subref_z=aste_subref_z,
        aste_subref_xt=aste_subref_xt,
        aste_subref_yt=aste_subref_yt,
        aste_subref_zt=aste_subref_zt,
        aste_misti_lon=aste_misti_lon,
        aste_misti_lat=aste_misti_lat,
        aste_misti_pwv=aste_misti_pwv,
        d2_mkid_id=corresp.index,
        d2_mkid_type=corresp.kidtype,
        d2_mkid_frequency=corresp.kidfreq,
        d2_skychopper_isblocking=skychop_state,
        d2_demerge_version=DEMERGE_VERSION,
        d2_ddb_version=ddb_version,
        # 18 arcsec MergeToDfits()でも固定値が指定されていた
        beam_major=0.005,
        # 18 arcsec MergeToDfits()でも固定値が指定されていた
        beam_minor=0.005,
        # 18 arcsec MergeToDfits()でも固定値が指定されていた
        beam_pa=0.005,
        # MergeToDfits()でも固定値が指定されていた
        exposure=1.0 / 196,
        # MergeToDfits()でも固定値が指定されていた
        interval=1.0 / 196,
        observation=obsinst_params["observation"],
        observer=obsinst_params["observer"],
        object=obsinst_params["object"],
    )

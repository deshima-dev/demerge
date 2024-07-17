__all__ = ["to_brightness", "to_dems"]


# standard library
import json
import re
from datetime import datetime as dt
from pathlib import Path
from typing import Any, Optional, Union
from warnings import catch_warnings, simplefilter


# dependencies
import numpy as np
import pandas as pd
import xarray as xr
from astropy.io import fits
from astropy.units import Quantity
from dems.d2 import MS
from numpy.typing import NDArray
from .. import __version__ as DEMERGE_VERSION


# type hints
PathLike = Union[Path, str]


# constants
COLUMN_NAMES_ANTENNA = (
    "time",  # %Y%m%d%H%M%S.%f
    "ra_prog",  # deg
    "dec_prog",  # deg
    "az_prog",  # deg
    "el_prog",  # deg
    "az_real",  # deg
    "el_real",  # deg
    "x",  # mm
    "y",  # mm
    "z",  # mm
    "xt",  # deg
    "yt",  # deg
    "zt",  # deg
    "lst",  # unknown
    "az_prog_no_cor",  # deg
    "el_prog_no_cor",  # deg
    "az_prog_center",  # deg
    "el_prog_center",  # deg
    "scan_type",
)
COLUMN_NAMES_CABIN = (
    "date",  # %Y/%m/%d
    "time",  # %H:%M
    "upper_temperature",  # degC
    "main_temperature",  # degC
    *[f"_{i}" for i in range(18)],  # unknown
)
COLUMN_NAMES_MISTI = (
    "date",  # %Y/%m/%d
    "time",  # %H:%M:%S.%f
    "unix_time",
    "az",  # deg
    "el",  # deg
    "pwv",  # um
    "ground_temperature",  # K
)
COLUMN_NAMES_SKYCHOP = (
    "time",  # unix
    "is_blocking",  # bool
)
COLUMN_NAMES_WEATHER = (
    "time",  # %Y%m%d%H%M%S
    "temperature",  # degC
    "pressure",  # hPa
    "humidity",  # %
    "wind_speed",  # m/s
    "wind_direction",  # deg
    "_",  # unknown
)
DATE_PARSER_ANTENNA = lambda s: dt.strptime(s, "%Y%m%d%H%M%S.%f")
DATE_PARSER_CABIN = lambda s: dt.strptime(s, "%Y/%m/%d %H:%M")
DATE_PARSER_MISTI = lambda s: dt.strptime(s, "%Y/%m/%d %H:%M:%S.%f")
DATE_PARSER_SKYCHOP = lambda s: dt.fromtimestamp(float(s))
DATE_PARSER_WEATHER = lambda s: dt.strptime(s, "%Y%m%d%H%M%S")
PACKAGE_DATA = Path(__file__).parents[1] / "data"


def get_antenna(antenna: PathLike, /) -> xr.Dataset:
    """Load an antenna log as xarray Dataset."""
    return pd.read_csv(  # type: ignore
        antenna,
        # read settings
        names=COLUMN_NAMES_ANTENNA,
        delimiter=r"\s+",
        comment="#",
        # index settings
        index_col=0,
        parse_dates=[0],
        date_parser=DATE_PARSER_ANTENNA,
    ).to_xarray()


def get_cabin(cabin: PathLike, /) -> xr.Dataset:
    """Load a cabin log as xarray Dataset."""
    return (
        pd.read_csv(  # type: ignore
            cabin,
            # read settings
            names=COLUMN_NAMES_CABIN,
            delimiter=r"\s+",
            comment="#",
            # index settings
            index_col=[0],
            parse_dates=[[0, 1]],
            date_parser=DATE_PARSER_CABIN,
        )
        .rename_axis("time")
        .to_xarray()
    )


def get_corresp(corresp: PathLike, /) -> xr.DataArray:
    """Load a KID correspondence as xarray DataArray."""
    with open(corresp) as f:
        masterid, kidid = zip(*json.load(f).items())

    return xr.DataArray(
        np.array(masterid, np.int64),
        name="masterid",
        dims=("kidid",),
        coords={"kidid": np.array(kidid, np.int64)},
    )


def get_ddb(ddb: PathLike, /) -> xr.Dataset:
    """Load a DDB FITS as xarray Dataset."""
    dim = "masterid"

    with fits.open(ddb) as hdus:
        # read from PRIMARY HDU
        version = hdus["PRIMARY"].header["DDB_ID"]

        # read from KIDDES HDU
        ds_kiddes = xr.Dataset(
            coords={
                dim: to_native((data := hdus["KIDDES"].data)[dim]),
            },
            data_vars={
                "type": (dim, to_native(data["attribute"])),
            },
        ).drop_duplicates(dim)

        # read from KIDFILT HDU
        ds_kidfilt = xr.Dataset(
            coords={
                dim: to_native((data := hdus["KIDFILT"].data)[dim]),
            },
            data_vars={
                "F": (dim, to_native(data["F_filter, df_filter"].T[0])),
                "Q": (dim, to_native(data["Q_filter, dQ_filter"].T[0])),
            },
        ).drop_duplicates(dim)

        # read from KIDRESP HDU
        ds_kidresp = xr.Dataset(
            coords={
                dim: to_native((data := hdus["KIDRESP"].data)[dim]),
            },
            data_vars={
                "p0": (dim, to_native(data["cal params"].T[0])),
                "fwd": (dim, to_native(data["cal params"].T[1])),
                "T0": (dim, to_native(data["cal params"].T[2])),
            },
        ).drop_duplicates(dim)

    ds = xr.merge([ds_kiddes, ds_kidfilt, ds_kidresp])
    ds = ds.where(ds.masterid >= 0, drop=True)
    return ds.assign_attrs(version=version)


def get_misti(misti: PathLike, /) -> xr.Dataset:
    """Load a MiSTI log as xarray Dataset."""
    return (
        pd.read_csv(  # type: ignore
            misti,
            # read settings
            names=COLUMN_NAMES_MISTI,
            delimiter=r"\s+",
            comment="#",
            # index settings
            index_col=[0],
            parse_dates=[[0, 1]],
            date_parser=DATE_PARSER_MISTI,
        )
        .rename_axis("time")
        .to_xarray()
    )


def get_obsinst(obsinst: PathLike, /) -> dict[str, str]:
    """Load an observation instruction to get parameters."""
    with open(obsinst) as f:
        lines = f.read()

    if match := re.search(r"(\d{14})", Path(obsinst).name):
        obs_id = match[1]
    else:
        obs_id = ""

    def search(pattern: str) -> str:
        if match := re.search(pattern, lines):
            return match[1]
        else:
            return ""

    return {
        # read from DES section
        "group": search(r"SET DES GROUP\s*'(.*)'"),
        "obs_file": search(r"SET DES OBS_FILE\s*'(.*)'"),
        "obs_user": search(r"SET DES OBS_USER\s*'(.*)'"),
        "project": search(r"SET DES PROJECT\s*'(.*)'"),
        # read from ANTENNA_G section
        "scan_cood": search(r"SET ANTENNA_G SCAN_COOD\s*'(.*)'"),
        "src_name": search(r"SET ANTENNA_G SRC_NAME\s*'(.*)'"),
        "src_pos": search(r"SET ANTENNA_G SRC_POS\s*\((.*)\)"),
        # read from file name
        "obs_id": obs_id,
    }


def get_readout(readout: PathLike, /) -> xr.DataArray:
    """Load a reduced readout FITS as xarray DataArray."""
    with fits.open(readout) as hdus:
        kidsinfo = hdus["KIDSINFO"].data
        readout_ = hdus["READOUT"].data

        # read from KIDSINFO HDU
        Qr = kidsinfo["Qr, dQr (Sky)"].T[0]
        fr = kidsinfo["fr, dfr (Sky)"].T[0]
        fr_room = kidsinfo["fr, dfr (Room)"].T[0]
        linyfc = kidsinfo["yfc, linyfc"].T[1]

        # read from READOUT HDU
        cols = readout_.columns[2:].names
        time = pd.to_datetime(readout_["timestamp"], unit="s")
        linph = np.array([readout_[col] for col in cols]).T[1]

    # calculate df/f (or fshift, dx)
    if np.isnan(fr_room).all():
        dfof = (linph - linyfc) / (4.0 * Qr)
    else:
        dfof = (linph - linyfc) / (4.0 * Qr) - (fr - fr_room) / fr

    return xr.DataArray(
        dfof,
        name="df/f",
        dims=("time", "kidid"),
        coords={"time": time, "kidid": np.arange(dfof.shape[1])},
    )


def get_skychop(skychop: PathLike, /) -> xr.Dataset:
    """Load a sky chopper log as xarray Dataset."""
    return pd.read_csv(  # type: ignore
        skychop,
        # read settings
        names=COLUMN_NAMES_SKYCHOP,
        delimiter=r"\s+",
        comment="#",
        # index settings
        index_col=0,
        parse_dates=[0],
        date_parser=DATE_PARSER_SKYCHOP,
    ).to_xarray()


def get_weather(weather: PathLike, /) -> xr.Dataset:
    """Load a weather log as xarray Dataset."""
    return pd.read_csv(  # type: ignore
        weather,
        # read settings
        names=COLUMN_NAMES_WEATHER,
        delimiter=r"\s+",
        comment="#",
        # index settings
        index_col=0,
        parse_dates=[0],
        date_parser=DATE_PARSER_WEATHER,
    ).to_xarray()


def to_brightness(dfof: xr.DataArray, /) -> xr.DataArray:
    """Convert a DEMS of df/f to that of brightness."""
    if np.isnan(T_room := dfof.aste_cabin_temperature.mean().data):
        T_room = 293.0

    if np.isnan(T_amb := dfof.temperature.mean().data):
        T_amb = 273.0

    fwd = dfof.d2_resp_fwd.data
    p0 = dfof.d2_resp_p0.data
    T0 = dfof.d2_resp_t0.data

    return (
        dfof.copy(
            deep=True,
            data=(dfof.data + p0 * np.sqrt(T_room + T0)) ** 2 / (p0**2 * fwd)
            - T0 / fwd
            - (1 - fwd) / fwd * T_amb,
        )
        .rename("Brightness")
        .assign_attrs(long_name="Brightness", units="K")
    )


def to_dems(
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
) -> xr.DataArray:
    """Merge observation datasets into a single DEMS of df/f.

    Args:
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

    Returns:
        Merged DEMS of df/f as xarray DataArray.

    """
    # load required datasets
    corresp_ = get_corresp(corresp)
    ddb_ = get_ddb(ddb)
    readout_ = get_readout(readout)
    obsinst_ = get_obsinst(obsinst)

    # load optional datasets
    if antenna is None:
        antenna = PACKAGE_DATA / "missing.ant"

    if cabin is None:
        cabin = PACKAGE_DATA / "missing.cabin"

    if misti is None:
        misti = PACKAGE_DATA / "missing.misti"

    if skychop is None:
        skychop = PACKAGE_DATA / "missing.skychop"

    if weather is None:
        weather = PACKAGE_DATA / "missing.wea"

    with catch_warnings():
        simplefilter("ignore")
        antenna_ = get_antenna(antenna)
        cabin_ = get_cabin(cabin)
        misti_ = get_misti(misti)
        skychop_ = get_skychop(skychop)
        weather_ = get_weather(weather)

    # merge datasets
    mkid = xr.merge([corresp_, readout_], join="left")
    mkid = mkid.swap_dims({"kidid": "masterid"})
    mkid = xr.merge([mkid, ddb_], join="left")

    # correct for time offset and sampling
    # fmt: off
    antenna_ = (
        antenna_
        .assign_coords(time=antenna_.time + to_timedelta(dt_antenna))
        .interp_like(readout_, kwargs={"fill_value": "extrapolate"})
    )
    cabin_ = (
        cabin_
        .assign_coords(time=cabin_.time + to_timedelta(dt_cabin))
        .interp_like(readout_, kwargs={"fill_value": "extrapolate"})
    )
    misti_ = (
        misti_
        .assign_coords(time=misti_.time + to_timedelta(dt_misti))
        .interp_like(readout_, kwargs={"fill_value": "extrapolate"})
    )
    skychop_ = (
        skychop_
        .assign_coords(time=skychop_.time + to_timedelta(dt_skychop))
        .interp_like(readout_, kwargs={"fill_value": "extrapolate"})
    )
    weather_ = (
        weather_
        .assign_coords(time=weather_.time + to_timedelta(dt_weather))
        .interp_like(readout_, kwargs={"fill_value": "extrapolate"})
    )
    # fmt: on

    # calculate coordinates
    if obsinst_["scan_cood"] == "RAZEL":
        lon = antenna_.az_prog_no_cor + (antenna_.az_real - antenna_.az_prog)
        lat = antenna_.el_prog_no_cor + (antenna_.el_real - antenna_.el_prog)
        lon_origin = antenna_.az_prog_center
        lat_origin = antenna_.el_prog_center
        frame = "altaz"
    else:  # == "RRADEC"
        lon = antenna_.ra_prog
        lat = antenna_.dec_prog
        lon_origin = np.full_like(lon, obsinst_["src_pos"].split(",")[0], float)
        lat_origin = np.full_like(lat, obsinst_["src_pos"].split(",")[1], float)
        frame = "fk5"

    return MS.new(
        # data
        data=mkid["df/f"].data,
        long_name="df/f",
        units="dimensionless",
        name="df/f",
        # dimensions
        time=mkid.time.data,
        chan=mkid.masterid.data,
        # labels
        beam=np.where(skychop_.is_blocking.data, "B", "A"),
        state=antenna_.scan_type.data,
        # telescope pointing
        lon=lon.data,
        lat=lat.data,
        lon_origin=lon_origin.data,
        lat_origin=lat_origin.data,
        frame=frame,
        # weather information
        temperature=weather_.temperature.data + 273.15,  # degC -> K
        pressure=weather_.pressure.data * 100,  # Pa -> hPa
        humidity=weather_.humidity.data,
        wind_speed=weather_.wind_speed.data,
        wind_direction=weather_.wind_direction.data,
        # data information
        frequency=mkid.F.data * 1e9,  # GHz -> Hz
        exposure=1 / 160,
        interval=1 / 160,
        # observation information
        observation=obsinst_["obs_file"],
        observer=obsinst_["obs_user"],
        project=f"{obsinst_['group']}/{obsinst_['project']}",
        object=obsinst_["src_name"],
        telescope_name="ASTE",
        telescope_diameter=10.0,
        telescope_coordinates=(
            +2230817.2140945992,
            -5440188.022176585,
            -2475718.801708271,
        ),
        # aste specific
        aste_cabin_temperature=cabin_.main_temperature.data + 273.15,  # degC -> K
        aste_obs_id=obsinst_["obs_id"],
        aste_obs_group=obsinst_["group"],
        aste_obs_project=obsinst_["project"],
        aste_obs_file=obsinst_["obs_file"],
        aste_obs_user=obsinst_["obs_user"],
        aste_subref_x=antenna_.x.data,
        aste_subref_y=antenna_.y.data,
        aste_subref_z=antenna_.z.data,
        aste_subref_xt=antenna_.xt.data,
        aste_subref_yt=antenna_.yt.data,
        aste_subref_zt=antenna_.zt.data,
        aste_misti_lon=misti_.az.data,
        aste_misti_lat=misti_.el.data,
        aste_misti_pwv=misti_.pwv.data * 1e-3,  # um -> mm
        aste_misti_frame="altaz",
        # deshima 2.0 specific
        d2_mkid_id=mkid.masterid.data,
        d2_mkid_type=mkid.type.data,
        d2_mkid_frequency=mkid.F.data * 1e9,  # GHz -> Hz
        d2_mkid_q=mkid.Q.data,
        d2_resp_fwd=mkid.fwd.data,
        d2_resp_p0=mkid.p0.data,
        d2_resp_t0=mkid.T0.data,
        d2_skychopper_isblocking=skychop_.is_blocking.data,
        d2_ddb_version=ddb_.version,
        d2_demerge_version=DEMERGE_VERSION,
    )


def to_native(array: NDArray[Any], /) -> NDArray[Any]:
    """Convert the byte order of an array to native."""
    return array.astype(array.dtype.type)


def to_timedelta(dt: Union[int, str], unit: str = "ms", /) -> np.timedelta64:
    """Convert a time offset to NumPy timedelta (float will be rounded)."""
    return np.timedelta64(int(Quantity(dt, unit=unit).to(unit).value), unit)

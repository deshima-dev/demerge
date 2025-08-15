__all__ = [
    "open_antenna",
    "open_cabin",
    "open_cdb",
    "open_ddb",
    "open_misti",
    "open_obsinst",
    "open_readout",
    "open_skychop",
    "open_weather",
    "to_dems",
]


# standard library
import re
from datetime import datetime as dt
from os import PathLike
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
from metpy.calc import saturation_vapor_pressure
from metpy.units import units
from numpy.typing import NDArray
from .. import __version__ as DEMERGE_VERSION


# type hints
StrPath = Union[PathLike[str], str]


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
    "vapor_pressure",  # hPa
    "wind_speed",  # m/s
    "wind_direction",  # deg
    "_",  # unknown
)
DATE_PARSER_ANTENNA = lambda s: dt.strptime(s, "%Y%m%d%H%M%S.%f")
DATE_PARSER_CABIN = lambda s: dt.strptime(s, "%Y/%m/%d %H:%M")
DATE_PARSER_MISTI = lambda s: dt.strptime(s, "%Y/%m/%d %H:%M:%S.%f")
DATE_PARSER_OBSID = lambda s: dt.strptime(s, "%Y%m%d%H%M%S")
DATE_PARSER_SKYCHOP = lambda s: dt.utcfromtimestamp(float(s))
DATE_PARSER_WEATHER = lambda s: dt.strptime(s, "%Y%m%d%H%M%S")
MASTERID_MISSING = -1
PACKAGE_DATA = Path(__file__).parents[1] / "data"


def open_antenna(antenna: StrPath, /) -> xr.Dataset:
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


def open_cabin(cabin: StrPath, /) -> xr.Dataset:
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


def open_cdb(cdb: StrPath, /) -> xr.DataArray:
    """Load a CDB Zarr as xarray DataArray."""
    return xr.open_dataarray(cdb, engine="zarr").compute()


def open_ddb(ddb: StrPath, /) -> xr.Dataset:
    """Load a DDB FITS as xarray Dataset."""
    dim = "masterid"

    with fits.open(ddb) as hdus:  # type: ignore
        # read from PRIMARY HDU
        version = hdus["PRIMARY"].header["DDB_ID"]

        # read from KIDDES HDU
        kiddes = xr.Dataset(
            coords={
                dim: (data := to_native(hdus["KIDDES"].data))[dim],
            },
            data_vars={
                "type": (dim, np.array(data["attribute"])),
            },
        )

        # read from KIDFILT HDU
        kidfilt = xr.Dataset(
            coords={
                dim: (data := to_native(hdus["KIDFILT"].data))[dim],
            },
            data_vars={
                "F": (dim, data["F_filter, df_filter"].T[0]),
                "Q": (dim, data["Q_filter, dQ_filter"].T[0]),
            },
        )

        # read from KIDRESP HDU
        kidresp = xr.Dataset(
            coords={
                dim: (data := to_native(hdus["KIDRESP"].data))[dim],
            },
            data_vars={
                "p0": (dim, data["cal params"].T[0]),
                "fwd": (dim, data["cal params"].T[1]),
                "T0": (dim, data["cal params"].T[2]),
            },
        )

    return xr.merge(
        [
            kiddes.where(kiddes.masterid != MASTERID_MISSING, drop=True),
            kidfilt.where(kidfilt.masterid != MASTERID_MISSING, drop=True),
            kidresp.where(kidresp.masterid != MASTERID_MISSING, drop=True),
        ]
    ).assign_attrs(version=version)


def open_misti(misti: StrPath, /) -> xr.Dataset:
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


def open_obsinst(obsinst: StrPath, /) -> dict[str, str]:
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


def open_readout(readout: StrPath, /) -> xr.DataArray:
    """Load a reduced readout FITS as xarray DataArray."""
    with fits.open(readout) as hdus:  # type:ignore
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
        coords={
            "time": time,
            "kidid": np.arange(dfof.shape[1]),
            "enabled": ("kidid", ~np.isnan(fr)),
        },
    )


def open_skychop(skychop: StrPath, /) -> xr.Dataset:
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


def open_weather(weather: StrPath, /) -> xr.Dataset:
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


def to_dems(
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
) -> xr.DataArray:
    """Merge observation datasets into a single DEMS of df/f.

    Args:
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

    Returns:
        Merged DEMS of df/f as xarray DataArray.

    """
    # record merge options first
    d2_merge_options = to_merge_options(locals())

    # load required datasets
    cdb_ = open_cdb(cdb)
    ddb_ = open_ddb(ddb)
    readout_ = open_readout(readout)
    obsinst_ = open_obsinst(obsinst)

    if not include_filterless_mkids:
        ddb_ = ddb_.dropna("masterid")

    if not include_disabled_mkids:
        readout_ = readout_.where(readout_.enabled, drop=True)

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
        antenna_ = open_antenna(antenna)
        cabin_ = open_cabin(cabin)
        misti_ = open_misti(misti)
        skychop_ = open_skychop(skychop)
        weather_ = open_weather(weather)

    # select KID correspondence
    cdb_ = cdb_.reindex(
        time=[DATE_PARSER_OBSID(obsinst_["obs_id"])],
        method="ffill",
        fill_value=MASTERID_MISSING,  # type: ignore
    )
    cdb_ = cdb_.where(cdb_ != MASTERID_MISSING, drop=True)

    # merge datasets
    merged_ = xr.merge([cdb_[0], readout_], join="inner")
    merged_ = merged_.swap_dims({"kidid": "masterid"})
    merged_ = xr.merge([ddb_, merged_], join="inner")

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
        data=merged_["df/f"].data,
        long_name="df/f",
        units="dimensionless",
        name=obsinst_["obs_id"],
        # dimensions
        time=merged_.time.data,
        chan=merged_.masterid.data,
        # labels
        observation=obsinst_["obs_id"],
        scan=(scan := to_phase(antenna_.scan_type)).data,
        subscan=skychop_.is_blocking.astype(bool).groupby(scan).apply(to_phase).data,
        state=antenna_.scan_type.data,
        beam=np.where(skychop_.is_blocking.data, "B", "A"),
        # telescope pointing
        lon=lon.data,
        lat=lat.data,
        lon_origin=lon_origin.data,
        lat_origin=lat_origin.data,
        frame=frame,
        # weather information
        temperature=weather_.temperature.data + 273.15,  # degC -> K
        pressure=weather_.pressure.data * 100,  # Pa -> hPa
        humidity=to_humidity(weather_.vapor_pressure.data, weather_.temperature.data),
        pwv=misti_.pwv.data * 1e-3,  # um -> mm
        wind_speed=weather_.wind_speed.data,
        wind_direction=weather_.wind_direction.data,
        # data information
        frequency=merged_.F.data * 1e9,  # GHz -> Hz
        exposure=1 / 160,
        interval=1 / 160,
        # observation information
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
        d2_mkid_id=merged_.masterid.data,
        d2_mkid_type=merged_.type.data,
        d2_mkid_frequency=merged_.F.data * 1e9,  # GHz -> Hz
        d2_mkid_q=merged_.Q.data,
        d2_resp_fwd=merged_.fwd.data,
        d2_resp_p0=merged_.p0.data,
        d2_resp_t0=merged_.T0.data,
        d2_skychopper_isblocking=skychop_.is_blocking.data,
        d2_ddb_version=ddb_.version,
        d2_merge_options=d2_merge_options,
        d2_demerge_version=DEMERGE_VERSION,
    )


def to_humidity(
    vapor_pressure: NDArray[Any],
    temperature: NDArray[Any],
    /,
) -> NDArray[Any]:
    """Convert water vapor pressure in hPa to relative humidity in %.

    Args:
        vapor_pressure: Water vaper pressure (array) in hPa.
        temperature: Atmospheric temperautre (array) in degC.

    Returns:
        Relative humidity (array) in %.

    """
    return (
        (
            (vapor_pressure * units.hPa)
            / saturation_vapor_pressure(temperature * units.degC)
        )
        .to("%")
        .magnitude
    )


def to_merge_options(locals: dict[str, Any], /) -> dict[str, Any]:
    """Convert local dictionary to merge-option dictionary."""
    merge_options = {}

    for key, val in locals.items():
        if val is None:
            continue

        if isinstance(val, Path):
            merge_options[key] = str(val)
        else:
            merge_options[key] = val

    return merge_options


def to_native(array: NDArray[Any], /) -> NDArray[Any]:
    """Convert the byte order of an array to native."""
    return array.astype(array.dtype.newbyteorder())


def to_phase(array: xr.DataArray, /) -> xr.DataArray:
    """Assign a phase to each value in a 1D DataArray."""
    if array.ndim != 1:
        raise ValueError("Input array must be 1D.")

    is_transision = xr.zeros_like(array, bool)
    is_transision.data[1:] = array.data[1:] != array.data[:-1]
    return is_transision.cumsum()


def to_timedelta(dt: Union[int, str], unit: str = "ms", /) -> np.timedelta64:
    """Convert a time offset to NumPy timedelta (float will be rounded)."""
    return np.timedelta64(int(Quantity(dt, unit=unit).to(unit).value), unit)

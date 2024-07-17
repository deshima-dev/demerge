# standard library
import json
import re
from datetime import datetime as dt
from pathlib import Path
from typing import Any, Union


# dependencies
import numpy as np
import pandas as pd
import xarray as xr
from astropy.io import fits
from numpy.typing import NDArray


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
    "upper_cabin_temperature",  # degC
    "main_cabin_temperature",  # degC
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


def get_antenna(antenna: PathLike, /) -> xr.Dataset:
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


def get_cabin(cabin: PathLike, /) -> xr.Dataset:
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


def get_obsinst(obsinst: PathLike, /) -> dict[str, str]:
    """Load an observation instruction to get parameters."""
    with open(obsinst) as f:
        lines = f.read()

    def search(pattern: str) -> str:
        if (match := re.search(pattern, lines)) is None:
            return ""
        else:
            return match[1]

    return {
        # DES
        "group": search("SET DES GROUP\s*'(.*)'"),
        "obs_file": search("SET DES OBS_FILE\s*'(.*)'"),
        "obs_user": search("SET DES OBS_USER\s*'(.*)'"),
        "project": search("SET DES PROJECT\s*'(.*)'"),
        # ANTENNA_G
        "scan_cood": search("SET ANTENNA_G SCAN_COOD\s*'(.*)'"),
        "src_name": search("SET ANTENNA_G SRC_NAME\s*'(.*)'"),
        "src_pos": search("SET ANTENNA_G SRC_POS\s*\((.*)\)"),
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


def get_weather(weather: PathLike, /) -> xr.Dataset:
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


def to_native(array: NDArray[Any], /) -> NDArray[Any]:
    """Convert the byte order of an array to native."""
    return array.astype(array.dtype.type)

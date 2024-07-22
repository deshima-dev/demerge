__all__ = ["Data", "parse_data"]


# standard library
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union


# type hints
PathLike = Union[Path, str]


# constants
PACKAGE_DATA = Path(__file__).parent


@dataclass
class Data:
    """Parsed data package structure."""

    antenna: Path
    """Path of the antenna log."""

    cabin: Path
    """Path of the cabin log."""

    corresp: Path
    """Path of the KID correspondence."""

    misti: Path
    """Path of the MiSTI log."""

    obsinst: Path
    """Path of the observation instruction."""

    readout: Path
    """Path of the KID readout FITS."""

    skychop: Path
    """Path of the sky chopper log."""

    weather: Path
    """Path of the weather log."""


def first(glob_results: Iterator[Path]) -> Optional[Path]:
    """Return the first Path.glob results if exists."""
    for path in glob_results:
        return path


def parse_data(data: PathLike, /) -> Data:
    """Parse a data package (data directory)."""
    if not (data := Path(data)).exists():
        raise FileNotFoundError(data)

    if (antenna := first(data.glob("*.ant"))) is None:
        antenna = PACKAGE_DATA / "missing.ant"

    if (cabin := first(data.glob("*.cabin"))) is None:
        cabin = PACKAGE_DATA / "missing.cabin"

    if (corresp := first(data.glob("*.json"))) is None:
        raise FileNotFoundError("KID correspondence (*.json).")

    if (misti := first(data.glob("*.misti"))) is None:
        misti = PACKAGE_DATA / "missing.misti"

    if (obsinst := first(data.glob("*.obs"))) is None:
        raise FileNotFoundError(f"Observation instruction (*.obs).")

    if (readout := first(data.glob("*.fits*"))) is None:
        raise FileNotFoundError(f"KID readout FITS (*.fits).")

    if (skychop := first(data.glob("*.skychopper*"))) is None:
        skychop = PACKAGE_DATA / "missing.skychop"

    if (weather := first(data.glob("*.wea"))) is None:
        weather = PACKAGE_DATA / "missing.wea"

    return Data(
        antenna=antenna,
        cabin=cabin,
        corresp=corresp,
        misti=misti,
        obsinst=obsinst,
        readout=readout,
        skychop=skychop,
        weather=weather,
    )

__all__ = ["DataPackage", "parse"]


# standard library
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union


# type hints
PathLike = Union[Path, str]


@dataclass
class DataPackage:
    """Parsed data package structure."""

    antenna: Optional[Path]
    """Path of the antenna log (optional)."""

    cabin: Optional[Path]
    """Path of the cabin log (optional)."""

    corresp: Path
    """Path of the KID correspondence (required)."""

    misti: Optional[Path]
    """Path of the MiSTI log (optional)."""

    obsinst: Path
    """Path of the observation instruction (required)."""

    readout: Path
    """Path of the KID readout FITS (required)."""

    skychop: Optional[Path]
    """Path of the sky chopper log (optional)."""

    weather: Optional[Path]
    """Path of the weather log (optional)."""


def first(glob_results: Iterator[Path]) -> Optional[Path]:
    """Return the first Path.glob results if exists."""
    for path in glob_results:
        return path


def parse(data_pack: PathLike, /) -> DataPackage:
    """Parse a data package (data directory)."""
    if not (data_pack := Path(data_pack)).exists():
        raise FileNotFoundError(data_pack)

    if (corresp := first(data_pack.glob("*.json"))) is None:
        raise FileNotFoundError("KID correspondence (*.json).")

    if (obsinst := first(data_pack.glob("*.obs"))) is None:
        raise FileNotFoundError(f"Observation instruction (*.obs).")

    if (readout := first(data_pack.glob("*.fits*"))) is None:
        raise FileNotFoundError(f"KID readout FITS (*.fits).")

    return DataPackage(
        antenna=first(data_pack.glob("*.ant")),
        cabin=first(data_pack.glob("*.cabin")),
        corresp=corresp,
        misti=first(data_pack.glob("*.misti")),
        obsinst=obsinst,
        readout=readout,
        skychop=first(data_pack.glob("*.skychopper*")),
        weather=first(data_pack.glob("*.wea")),
    )

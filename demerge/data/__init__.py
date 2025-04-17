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


def last(glob_results: Iterator[Path], /) -> Optional[Path]:
    """Return the last Path.glob result if exists."""
    for path in reversed(sorted(glob_results)):
        return path


def parse(data_pack: PathLike, /) -> DataPackage:
    """Parse a data package (data directory)."""
    if not (data_pack := Path(data_pack)).exists():
        raise FileNotFoundError(data_pack)

    if (corresp := last(data_pack.glob("*.json"))) is None:
        raise FileNotFoundError("KID correspondence (*.json).")

    if (obsinst := last(data_pack.glob("*.obs"))) is None:
        raise FileNotFoundError("Observation instruction (*.obs).")

    if (readout := last(data_pack.glob("*.fits*"))) is None:
        raise FileNotFoundError("KID readout FITS (*.fits).")

    return DataPackage(
        antenna=last(data_pack.glob("*.ant")),
        cabin=last(data_pack.glob("*.cabin")),
        corresp=corresp,
        misti=last(data_pack.glob("*.misti")),
        obsinst=obsinst,
        readout=readout,
        skychop=last(data_pack.glob("*.skychopper*")),
        weather=last(data_pack.glob("*.wea")),
    )

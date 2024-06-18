__all__ = ["analyze"]


# dependencies
from fire import Fire


def analyze() -> None:
    """Analyze a TOD FITS and create a reduced FITS."""
    ...


def cli() -> None:
    """Command line interface of the analyze function."""
    Fire(analyze)

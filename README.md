# de:merge

[![Release](https://img.shields.io/pypi/v/demerge?label=Release&color=cornflowerblue&style=flat-square)](https://pypi.org/project/demerge/)
[![Python](https://img.shields.io/pypi/pyversions/demerge?label=Python&color=cornflowerblue&style=flat-square)](https://pypi.org/project/demerge/)
[![Downloads](https://img.shields.io/pypi/dm/demerge?label=Downloads&color=cornflowerblue&style=flat-square)](https://pepy.tech/project/demerge)
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.10015892-cornflowerblue?style=flat-square)](https://doi.org/10.5281/zenodo.10015892)
[![Tests](https://img.shields.io/github/actions/workflow/status/deshima-dev/demerge/tests.yaml?label=Tests&style=flat-square)](https://github.com/deshima-dev/demerge/actions)

DESHIMA merge code for observed datasets

## Installation

```shell
pip install demerge==2025.8.4
```

## Command line interface

de:merge ships with three commands, `demerge`, `reduce`, and `merge`.
`demerge` is the main command that (1) reduces raw data of KID measurements into a single "reduced" FITS and (2) merge datasets of a DESHIMA observation (including the reduced FITS) into a single DEMS file.
`reduce` and `merge` are the subcommands that perform (1) and (2), respectively, and they should be used only for debugging.

Here is (part of) the `demerge` command help:
```plaintext
SYNOPSIS
    demerge OBSID <flags>

DESCRIPTION
    Run reduce and merge commands to create a single DEMS.

POSITIONAL ARGUMENTS
    OBSID
        Type: str
        Observation ID (YYYYMMDDHHMMSS).

FLAGS
    --data_dir=DATA_DIR
        Type: Path
        Default: PosixPath('.')
        Path of directory where data packages are placed,
        i.e. expecting ``${data_dir}/cosmos_YYYYMMDDHHMMSS``.
    --dems_dir=DEMS_DIR
        Type: Path
        Default: PosixPath('.')
        Path of directory where merged DEMS will be placed,
        i.e. expecting ``${dems_dir}/dems_YYYYMMDDHHMMSS.zarr.zip``.
    --reduced_dir=REDUCED_DIR
        Type: Path
        Default: PosixPath('.')
        Path of directory where reduced packages are placed,
        i.e. expecting ``${reduced_dir}/reduced_YYYYMMDDHHMMSS``.
        If not specified, a temporary directory will be used.
    --cdb=CDB
        Type: Path
        Default: PosixPath('/path/to/demerge/data/cdb_20250528.zarr.zip')
        Path of CDB (KID correspondence database) file.
    --ddb=DDB
        Type: Path
        Default: PosixPath('/path/to/demerge/data/ddb_20250819.fits.gz')
        Path of DDB (DESHIMA database) file.
    --overwrite=OVERWRITE
        Type: bool
        Default: False
        If True, the reduced package and the merged DEMS file
        will be overwritten even if they exist.
    --debug=DEBUG
        Type: bool
        Default: False
        If True, detailed logs for debugging will be printed.
```

A typical command run is as follows:
```shell
demerge YYYYMMDDHHMMSS --data_dir /path/to/data --reduced_dir .cache
```
where the raw data directory, `/path/to/data/cosmos_YYYYMMDDHHMMSS`, will be selected, then an intermediate reduced data directory will be created as `./.cache/reduced_YYYYMMDDHHMMSS`, and finally the merged DEMS file will be placed as `./dems_YYYYMMDDHHMMSS.zarr.zip`.

## Raw data directory

The name and the file structure of a raw data directory must be as follows:

```plaintext
cosmos_YYYYMMDDHHMMSS/
├── YYYYMMDDHHMMSS.ant
├── YYYYMMDDHHMMSS.cabin
├── YYYYMMDDHHMMSS.fits
├── YYYYMMDDHHMMSS.misti
├── YYYYMMDDHHMMSS.obs
├── YYYYMMDDHHMMSS.skychopper.dat.xz
├── YYYYMMDDHHMMSS.wea
├── YYYYMMDDHHMMSS_info.txt
├── kid_corresp.json
├── kids.list
├── localsweep.sweep
├── localsweep_info.txt
└── pretune
    ├── kids.list
    ├── localsweep.sweep
    └── localsweep_info.txt
```

If multiple files with the same extension exist, then the last one sorted by name will be selected.

## Use as a Python package

All (sub)commands are also available as Python functions of the same names in the `demerge` package.
For example, the following two runs are equivalent:

```shell
# run by a shell script
demerge YYYYMMDDHHMMSS --data_dir /path/to/data --reduced_dir .cache
```

```python
# run by a Python script
from demerge import demerge
demerge("YYYYMMDDHHMMSS", data_dir="/path/to/data", reduced_dir=".cache")
```

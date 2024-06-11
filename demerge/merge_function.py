"""merge_to_dfits.py: Read logging data and merge them into a FITS object

 Author : Tetsutaro Ueda, Junya Suzuki, Kenichi Karatsu, Tatsuya Takekoshi
 Created: 2017/11/02
 Revision History:
     2018/02/08 - KK - rewrite using class.
     2018/06/08 - TT - apply new calibration method.
     2021         NAITO systems modfied.
     2023         NAITO systems modfied.
"""
__all__ = [
    'FORM_FITSTIME',
    'FORM_FITSTIME_P',
    'DEFAULT_ROOM_T',
    'create_bintablehdu',
    'load_obsinst',
    'get_maskid_corresp'
    'Tlos_model',
    'convert_readout',
    'convert_asciitime',
    'convert_timestamp',
    'update_corresp',
]


# standard library
import json
import lzma
from datetime import datetime
from functools import partial, reduce
from typing import Any, Literal


# dependencies
import numpy as np
import pandas as pd
from astropy.io import fits, ascii
from astropy.io.fits import BinTableHDU
from numpy.typing import NDArray


# constants
FORM_FITSTIME   = '%Y-%m-%dT%H:%M:%S'                          # YYYY-mm-ddTHH:MM:SS
FORM_FITSTIME_P = '%Y-%m-%dT%H:%M:%S.%f'                       # YYYY-mm-ddTHH:MM:SS.ss

CABIN_Q_MARGIN  = 5*60 # seconds. Margin for cabin data query.
DEFAULT_ROOM_T  = 17. + 273. # Kelvin
DEFAULT_AMB_T   = 0.  + 273. # Kelvin


# constants (master-to-KID correspondence)
CORRESP_IGNORES = "pixelid", "runid, framelen, refid"
CORRESP_NOTFOUND = -1
KIDFILT = "KIDFILT"
KIDID = "kidid"
MASTERID = "masterid"


def create_bintablehdu(hd):
    """Create Binary Table HDU from 'hdu_dict'"""
    header = fits.Header()
    for (i, j) in zip(hd['hdr_vals'].items(), hd['hdr_coms'].items()):
        header[i[0]] = i[1], j[1]
    columns = [
        fits.Column(name=i[0], format=j[1], array=i[1], unit=k[1])
        for (i, j, k) in zip(
            hd['col_vals'].items(),
            hd['col_form'].items(),
            hd['col_unit'].items()
        )
    ]
    hdu = fits.BinTableHDU.from_columns(columns, header)
    for i in hd['hdr_coms'].items():
        hdu.header.comments[i[0]] = i[1]
    return hdu

def load_obsinst(obsinst):
    """Get data for 'OBSINFO'"""
    if not '.obs' in obsinst:
        raise ValueError('The input file must be an observational instruction!!')

    with open(obsinst, 'r') as f:
        equinox = 2000  # Default parameter
        for line in f:
            if 'SET ANTENNA_G TRK_TYPE' in line:
                trktype = line.split()[-1].strip('\'')
            elif 'SET ANTENNA_G SRC_NAME' in line:
                obs_object = line.split()[-1].strip('\'')
            elif 'SET ANTENNA_G SRC_POS' in line:
                srcpos = [float(c) for c in line.split()[-1].strip('()').split(',')]
            elif 'SET ANTENNA_G EPOCH' in line:
                equinox = line.split()[-1].strip('\'JB')
            elif 'SET DES OBS_USER' in line:
                observer = line.split()[-1].strip('\'')
            elif 'SET DES PROJECT' in line:
                project = line.split()[-1].strip('\'')
            elif 'SET DES PROJECT' in line:
                project = line.split()[-1].strip('\'')
            elif '% OBS=' in line:
                observation = line.split('=')[-1].strip()
    if trktype == 'RADEC':
        ra  = srcpos[0]
        dec = srcpos[1]
    else:
        ra  = 0
        dec = 0
    return {'observer': observer, 'obs_object': obs_object,  'ra': ra, 'dec': dec, 'equinox': equinox, 'project': project, 'observation': observation}


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
    data = ddb['KIDDES'].data
    frame = pd.DataFrame(
        index=pd.Index(
            native(data['masterid']),
            name='masterid',
        ),
        data = {
            'kidtype': native(data['attribute']),
        },
    )
    frames.append(frame)

    # DataFrame of KIDFILT HDU
    data = ddb['KIDFILT'].data
    frame = pd.DataFrame(
        index=pd.Index(
            native(data['masterid']),
            name='masterid'
        ),
        data={
            'kidfreq': native(data["F_filter, df_filter"][:, 0]),
            'kidQ': native(data['Q_filter, dQ_filter'][:, 0]),
        },
    )
    frame['kidfreq'] *= 1e9
    frames.append(frame)

    # DataFrame of KIDRESP HDU
    if 'KIDRESP' in ddb:
        data = ddb['KIDRESP'].data
        frame = pd.DataFrame(
            index=pd.Index(
                native(data['masterid']),
                name='masterid',
            ),
            data={
                "p0": native(data['cal params'][:, 0]),
                "etaf": native(data['cal params'][:, 1]),
                "T0": native(data['cal params'][:, 2]),
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
        [corresp.get(str(i), -1) for i in frame.index],
        name='kidid'
    )
    frame = frame.reset_index().set_index(index)

    # Drop rows with invalid master or KID IDs (-1)
    frame = frame[frame.index != -1]
    frame = frame[frame.masterid != -1]
    return frame


def convert_readout(
    readout: fits.HDUList,
    corresp: pd.DataFrame,
    to: Literal['Tsignal', 'fshift'],
    T_room: float,
    T_amb: float,
):
    """Reduced readoutの値をDEMS出力形式に変換（校正）する

    Args:
        readout: Reduced readout FITSのHDUListオブジェクト
        corresp: KID IDと各KIDの測定値を対応づけるDataFrame
        to: 変換形式（Tsignal→Tsky, fshift→df/f)
        T_room: キャビン温度(K)
        T_amb: 外気温(K)

    """
    kidcols = readout['READOUT'].data.columns[2:].names
    linph = np.array([readout['READOUT'].data[n] for n in kidcols]).T[2]
    linyfc = np.array(readout['KIDSINFO'].data['yfc, linyfc']).T[1]
    Qr = np.array(readout['KIDSINFO'].data['Qr, dQr (300K)']).T[0]
    fshift = (linph - linyfc) / (4.0 * Qr)

    if to == 'fshift':
        return fshift[:, corresp.index.values]

    if to == 'Tsignal':
        return Tlos_model(
            dx=fshift[:, corresp.index.values],
            p0=corresp.p0.values,
            etaf=corresp.etaf.values,
            T0=corresp.T0.values,
            Troom=T_room,
            Tamb=T_amb,
        )

    raise ValueError(f'Invalid output type: {to}')


def Tlos_model(dx, p0, etaf, T0, Troom, Tamb):
    """Calibrate 'amplitude' and 'phase' to 'power'"""
    return (dx + p0*np.sqrt(Troom+T0))**2 / (p0**2 * etaf) - T0/etaf - (1-etaf)/etaf*Tamb

def convert_asciitime(asciitime, form_fitstime):
    """Ascii time"""
    asciitime = [datetime.strptime('%14.6f' %t, '%Y%m%d%H%M%S.%f') for t in asciitime]
    asciitime = [datetime.strftime(t, form_fitstime) for t in asciitime]
    return np.array(asciitime)

def convert_timestamp(timestamp):
    """Timestamp"""
    timestamp = [datetime.utcfromtimestamp(t) for t in timestamp]
    timestamp = [datetime.strftime(t, FORM_FITSTIME_P) for t in timestamp]
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
    if filename=='' or filename==None:
        return (
            np.array(["1970-01-01"]).astype('datetime64[ns]'),
            np.array([20.0]).astype(np.float64),
            np.array([20.0]).astype(np.float64),
        )

    table = ascii.read(filename, format='no_header')

    # 日付と時刻を取得して文字列でタイムスタンプを作成しそれをnumpy.datetime64へ変換する
    # テーブルの1列目と2列目がそれぞれ日付と時刻
    datetimes = []
    for date, time in zip(table['col1'], table['col2']):
        s = '{}T{}'.format(date, time)
        s = s.replace('/', '-')
        datetimes.append(s)

    datetimes         = np.array(datetimes).astype('datetime64[ns]')
    upper_cabin_temps = np.array(table['col3']).astype(np.float64)
    lower_cabin_temps = np.array(table['col4']).astype(np.float64)

    return (datetimes, upper_cabin_temps, lower_cabin_temps)

def retrieve_skychop_states(filename):
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
    data = None
    if filename.endswith('.xz'):
        with lzma.open(filename, 'rt') as f:
            data = f.read()
    else:
        with open(filename, 'rt') as f:
            data = f.read()

    table = ascii.read(data, guess=False, format='no_header', delimiter=' ', names=['datetime', 'state'])
    datetimes = np.array(table['datetime']).astype(np.float64)
    states    = np.array(table['state']).astype(np.int8)
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
        3列目 az(deg)
        4列目 el(deg)
        5列目 pwv(um)
        6列目 Tground(K)
        "#"から始まるコメントがファイル冒頭に数行ある。

    """
    if filename=='' or filename==None:
        return (np.array([np.nan]).astype('datetime64[ns]'),
                np.array([np.nan]).astype(np.float64),
                np.array([np.nan]).astype(np.float64),
                np.array([np.nan]).astype(np.float64))

    column_names = [
        'date',
        'time',
        'az',
        'el',
        'pwv',
        'Tround',
    ]
    table = ascii.read(filename, guess=False, format='no_header', delimiter=' ', names=column_names)

    az  = np.array(table['az']).astype(np.float64)
    el  = np.array(table['el']).astype(np.float64)
    pwv = np.array(table['pwv']).astype(np.float64)/1000.0 # umからmmへ変換

    datetimes = []
    for row in table:
        datetimes.append(datetime.strptime('{} {}'.format(row['date'], row['time']), '%Y/%m/%d %H:%M:%S.%f'))

    return (np.array(datetimes).astype('datetime64[ns]'), az, el, pwv)

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
    'convert_timestamp'
]


# standard library
import lzma
from datetime import datetime
from typing import Literal


# dependencies
import numpy as np
from astropy.io import fits, ascii
from numpy.typing import NDArray


#-------------------------------- CONSTANTS
FORM_FITSTIME   = '%Y-%m-%dT%H:%M:%S'                          # YYYY-mm-ddTHH:MM:SS
FORM_FITSTIME_P = '%Y-%m-%dT%H:%M:%S.%f'                       # YYYY-mm-ddTHH:MM:SS.ss

CABIN_Q_MARGIN  = 5*60 # seconds. Margin for cabin data query.
DEFAULT_ROOM_T  = 17. + 273. # Kelvin
DEFAULT_AMB_T   = 0.  + 273. # Kelvin

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

def get_maskid_corresp(ddb: fits.HDUList):
    """Get Correspondance of 'master' and 'kid'"""
    kidnames = dict(
        zip(
            ddb['KIDDES'].data['masterid'],
            ddb['KIDDES'].data['attribute'],
        )
    )

    masterids: list[int]  = []
    kidids: list[int]     = []
    kidtypes: list[str]   = []
    kidfreqs: list[float] = []
    kidQs: list[float]    = []

    for i in range(len(ddb['KIDFILT'].data)):
        kidfilt_i = ddb['KIDFILT'].data[i]

        masterids.append(kidfilt_i['masterid'])
        kidids.append(kidfilt_i['kidid'])
        kidfreqs.append(kidfilt_i['F_filter, dF_filter'][0] * 1e9)
        kidQs.append(kidfilt_i['Q_filter, dQ_filter'][0])

        if (masterid := kidfilt_i['masterid']) < 0:
            kidtypes.append("unknown")
        else:
            kidtypes.append(kidnames[masterid])

    return masterids, kidids, kidtypes, kidfreqs, kidQs

def convert_readout(
    ro: fits.HDUList,
    ddb: fits.HDUList,
    to: Literal['Tsignal', 'fshift'],
    T_room: float,
    T_amb: float,
) -> NDArray:
    """Reduced readoutの値をDEMS出力形式に変換（校正）する

    Args:
        ro: Reduced readout FITSのHDUListオブジェクト
        ddb: DDB FITSのHDUListオブジェクト
        to: 変換形式（Tsignal→Tsky, fshift→df/f)
        T_room: キャビン温度(K)
        T_amb: 外気温(K)

    """
    kidcols = ro['READOUT'].data.columns[2:]
    linph   = np.array([ro['READOUT'].data[n] for n in kidcols.names]).T[2]
    linyfc  = np.array(ro['KIDSINFO'].data['yfc, linyfc']).T[1]
    Qr      = np.array(ro['KIDSINFO'].data['Qr, dQr (300K)']).T[0]
    fshift  = (linph - linyfc) / (4.0 * Qr)
    # ここまではKID ID = indexの対応となっている

    n_time = len(fshift)
    n_chan = len(ddb['KIDFILT'].data)
    output = np.zeros([n_time, n_chan], np.float32)

    for i in range(n_chan):
        kidid = ddb['KIDFILT'].data[i]['kidid']
        masterid = ddb['KIDFILT'].data[i]['masterid']
        fshift_id = fshift[:, kidid]

        if masterid < 0:
            output[:, i] = np.nan
        elif to == "fshift":
            output[:, i] = fshift_id
        elif to == "Tsignal":
            p0, etaf, T0 = ddb['KIDRESP'].data[i]['cal params']
            output[:, i] = Tlos_model(fshift_id, p0, etaf, T0, T_room, T_amb)
        else:
            raise ValueError(f'Invalid output type: {to}')

    return np.array(output)

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

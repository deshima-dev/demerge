"""merge_to_dfits.py: Read logging data and merge them into a FITS object

 Author : Tetsutaro Ueda, Junya Suzuki, Kenichi Karatsu, Tatsuya Takekoshi
 Created: 2017/11/02
 Revision History:
     2018/02/08 - KK - rewrite using class.
     2018/06/08 - TT - apply new calibration method.
     2021         NAITO systems modfied.
     2023         NAITO systems modfied.
"""
from __future__ import print_function

__all__ = [
    'FORM_FITSTIME',
    'FORM_FITSTIME_P',
    'DEFAULT_ROOM_T',
    'create_bintablehdu',
    'load_obsinst',
    'get_maskid_corresp'
    'Tlos_model',
    'calibrate_to_power',
    'convert_asciitime',
    'convert_timestamp'
]

from datetime import datetime
from datetime import timedelta
from calendar import timegm
from astropy.io import fits, ascii
import numpy as np
import scipy.interpolate
import sys

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

def get_maskid_corresp(pixelid, ddb):
    """Get Correspondance of 'master' and 'kid'"""
    nkid = ddb['KIDFILT'].header['NKID%d' %pixelid]
    kiddict, kidfilt = {}, {}
    for (i, j, k, l) in zip(ddb['KIDFILT'].data['kidid'],
                            ddb['KIDFILT'].data['masterid'],
                            ddb['KIDFILT'].data['F_filter, dF_filter'],
                            ddb['KIDFILT'].data['Q_filter, dQ_filter']):
        kiddict[i] = j
        kidfilt[i] = (k[0], l[0])
    kidname = {}
    for (i, j) in zip(ddb['KIDDES'].data['masterid'], ddb['KIDDES'].data['attribute']):
        kidname[i] = j

    masterids, kidids, kidtypes, kidfreqs, kidQs = [], [], [], [], []
    for i in range(nkid):
        masterid = kiddict[i]
        if masterid < 0:
            kind = 'unknown'
        else:
            kind = kidname[masterid]

        masterids.append( masterid )
        kidids.append( i )
        kidtypes.append( kind )
        kidfreqs.append( kidfilt[i][0] * 1e9 )
        kidQs.append( kidfilt[i][1] )
    return masterids, kidids, kidtypes, kidfreqs, kidQs

def Tlos_model(dx, p0, etaf, T0, Troom, Tamb):
    """Calibrate 'amplitude' and 'phase' to 'power'"""
    return (dx + p0*np.sqrt(Troom+T0))**2 / (p0**2 * etaf) - T0/etaf - (1-etaf)/etaf*Tamb

def calibrate_to_power(pixelid, Troom, Tamb, rhdus, ddb):
    nkid = rhdus['READOUT'].header['NKID%d' %pixelid]
    kiddict = {}
    for (i, j) in zip(ddb['KIDFILT'].data['kidid'], ddb['KIDFILT'].data['masterid']):
        kiddict[i] = j

    linphase = np.transpose([rhdus['READOUT'].data['Amp, Ph, linPh %d' %i].T[2] for i in range(nkid)])
    linyfc   = rhdus['KIDSINFO'].data['yfc, linyfc'].T[1]
    Qr       = rhdus['KIDSINFO'].data['Qr, dQr (300K)'].T[0]
    fshift = np.array((linphase - linyfc)/ (4.*Qr)).T
    fshift_err = np.zeros( len(fshift) )
    #---- Responsivity curve
    (p0, etaf, T0) = ddb['KIDRESP'].data['cal params'].T
    Tsignal = []
    for i in range(nkid):
        masterid = kiddict[i]
        if masterid<0:
            Tsignal.append( [np.nan for j in range( len(fshift[i]) )] )
            continue
        #---- Convert to power
        Tsignal.append(Tlos_model(fshift[i], p0[i], etaf[i], T0[i], Troom, Tamb))
    return np.array(Tsignal).T

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

def retrieve_cabin_temps(filename):
    """キャビン内温度を取得する
    引数
    ====
    str ファイル名

    戻り値
    ======
    tuple (timestames, upperCabinTemps, lowerCabinTemps)
      tupleの各要素はnumpy.array。要素数は同じ。
    """
    table = ascii.read(filename, format='no_header')

    # 日付と時刻を取得して文字列でタイムスタンプを作成しそれをnumpy.datetime64へ変換する
    # テーブルの1列目と2列目がそれぞれ日付と時刻
    datetimes = []
    for date, time in zip(table['col1'], table['col2']):
        s = '{}T{}'.format(date, time)
        s = s.replace('/', '-')
        datetimes.append(s)
    datetimes       = np.array(datetimes).astype('datetime64[ns]')
    upper_cabin_temps = np.array(table['col3']).astype(np.float64)
    lower_cabin_temps = np.array(table['col4']).astype(np.float64)

    return (datetimes, upper_cabin_temps, lower_cabin_temps)

def retrieve_skychop_states(filename):
    """skychopファイルからskychopの時系列状態を取得する
    引数
    ====
    str ファイル名

    戻り値
    ======
    tuple (timestames, states)
      tupleの各要素はnumpy.array。要素数は同じ。

    時刻について
    ============
    skychopファイルに記録されている時刻はUNIX時間。

    ファイル形式
    ============
    1列目 UNIX時刻
    2列目 0/1による状態
    "#"から始まるコメントがファイル冒頭に数行ある。
    """
    table = ascii.read(filename, guess=False, format='no_header', delimiter=' ', names=['datetime', 'state'])

    datetimes = np.array(table['datetime']).astype(np.float64)
    states    = np.array(table['state']).astype(np.int8)
    return (datetimes, states)

def retrieve_misti_log(filename):
    """mistiファイルからの時系列データを取得する
    引数
    ====
    str ファイル名

    戻り値
    ======
    tuple (timestames, az, el, pwv)
      tupleの各要素はnumpy.array。要素数は同じ。

    ファイル形式
    ============
    1列目 年/月/日
    2列目 時:分:6列目 秒(小数点以下2桁も含む)
    3列目 az(deg)
    4列目 el(deg)
    5列目 pwv(um)
    6列目 Tground(K)
    "#"から始まるコメントがファイル冒頭に数行ある。
    """
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

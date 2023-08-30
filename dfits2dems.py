"""dfitsをdemsへ変換するプログラム

File name: dfits2dems.py
Python 3.9
(C) 2023 内藤システムズ
"""

import numpy as np
#import ascii

from astropy.io        import fits, ascii
from dems.d2           import MS
from scipy.interpolate import interp1d

def convert_dfits_to_dems(filename, **kwargs):
    # Optionalな引数の処理
    coodinate           = kwargs.pop('coodinate', 'azel')
    loadmode            = kwargs.pop('loadmode', 0)
    # still data
    still_period        = kwargs.pop('still_period', None) # 秒(整数)
    # shuttle observation
    shuttle_min_lon_off = kwargs.pop('shuttle_min_lon_off', None)
    shuttle_max_lon_off = kwargs.pop('shuttle_max_lon_off', None)
    shuttle_min_lon_on  = kwargs.pop('shuttle_min_lon_on', None)
    shuttle_max_lon_on  = kwargs.pop('shuttle_max_lon_on', None)
    # find R
    ch                  = kwargs.pop('ch', 0)      # チャネル
    Rth                 = kwargs.pop('Rth', 280)   # R閾値
    skyth               = kwargs.pop('skyth', 150) # sky閾値
    cutnum              = kwargs.pop('cutnum', 1)  # カット番号

    # shuttle観測に関する引数が1つでも与えられたら与えられなかった変数に規定値(0.0)を設定する。
    # 1つも引数が与えられなかった場合はshuttle観測に関するマスクは無効とする。
    if not (shuttle_min_lon_off == None and shuttle_max_lon_off == None and shuttle_min_lon_on  == None and shuttle_max_lon_on  == None):
        if shuttle_min_lon_off == None:
            shuttle_min_lon_off = 0.0
        if shuttle_max_lon_off == None:
            shuttle_max_lon_off = 0.0
        if shuttle_min_lon_on  == None:
            shuttle_min_lon_on  = 0.0
        if shuttle_max_lon_on  == None:
            shuttle_max_lon_on  = 0.0

    with fits.open(filename) as hdul:
        readout   = hdul['READOUT'].data
        obsinfo   = hdul['OBSINFO'].data
        obsheader = hdul['OBSINFO'].header
        antenna   = hdul['ANTENNA'].data
        weather   = hdul['WEATHER'].data
        cabin     = hdul['CABIN_T'].data

    time         = np.array(readout['starttime']).astype(np.datetime64)
    time_antenna = np.array(antenna['time']).astype(np.datetime64)
    time_weather = np.array(weather['time']).astype(np.datetime64)
    time_cabin   = np.array(cabin['time']).astype(np.datetime64)

    # 補間のために時刻(年月日時分秒)を時間(秒)に変更する。READOUTの最初の時刻を基準とする。
    seconds         = (time         - time[0])/np.timedelta64(1, 's')
    seconds_antenna = (time_antenna - time[0])/np.timedelta64(1, 's')
    seconds_weather = (time_weather - time[0])/np.timedelta64(1, 's')
    seconds_cabin   = (time_cabin   - time[0])/np.timedelta64(1, 's')

    # モードに応じて経度(lon)と緯度(lat)を選択(azelかradecか)する
    if coodinate == 'azel':
        lon        = antenna['az']
        lat        = antenna['el']
        lon_origin = np.median(antenna['az_center'])
        lat_origin = np.median(antenna['el_center'])
        if loadmode in [0, 1]:
            lon -= antenna['az_center']
            lat -= antenna['el_center']
            if loadmode == 0:
                lon *= np.cos(np.deg2rad(antenna['el']))
    elif coodinate == 'radec':
        lon        = antenna['ra']
        lat        = antenna['dec']
        lon_origin = obsheader['RA']
        lat_origin = obsheader['DEC']
    else:
        raise KeyError('Invalid coodinate type: {}'.format(coodinate))

    scan = np.array(antenna['scantype'])

    # 補間関数で扱うためにSCANTYPE(文字列)を整数に置き換える
    scan_types = {scantype:i for i, scantype in enumerate(np.unique(scan))}
    scan_type_numbers = np.zeros(scan.shape[0], dtype=int)
    for scantype, i in scan_types.items():
        scan_type_numbers[scan == scantype] = i

    # READOUTの時間に合わせるためにSCANTYPEを補間する
    f_scantype = interp1d(seconds_antenna, scan_type_numbers, kind='nearest', bounds_error=False, fill_value=(scan_type_numbers[0], scan_type_numbers[-1]))
    scan_type_numbers = f_scantype(seconds)

    # 補間後のSCANTYPEを文字列に戻す
    scan = np.full_like(scan_type_numbers, 'GRAD', dtype='<U8')
    for scantype, i in scan_types.items():
        scan[scan == i] = scantype

    # 座標、気象情報、キャビン情報もREADOUTの時間に合わせて補間する
    lon                    = np.interp(seconds, seconds_antenna, lon)
    lat                    = np.interp(seconds, seconds_antenna, lat)
    temperature            = np.interp(seconds, seconds_weather, weather['temperature'])
    pressure               = np.interp(seconds, seconds_weather, weather['pressure'])
    humidity               = np.interp(seconds, seconds_weather, weather['vapor-pressure'])
    wind_speed             = np.interp(seconds, seconds_weather, weather['windspd'])
    wind_direction         = np.interp(seconds, seconds_weather, weather['winddir'])
    aste_cabin_temperature = np.interp(seconds, seconds_cabin,   cabin['main_cabin'])

    # 静止データの周期に応じてOFFマスクとSCANマスクを設定する
    if still_period != None:
        for i in range(int(seconds[-1]) // still_period + 1):
            off_mask = (still_period*(2*i)     <= seconds) & (seconds < still_period*(2*i + 1))
            on_mask  = (still_period*(2*i + 1) <= seconds) & (seconds < still_period*(2*i + 2))
            scan[off_mask] = 'OFF'
            scan[on_mask]  = 'SCAN'

    # shuttle観測のマスクを設定する
    if shuttle_min_lon_off != None:
        #print(shuttle_min_lon_on < lon)
        off_mask = (shuttle_min_lon_off < lon) & (lon < shuttle_max_lon_off)
        on_mask  = (shuttle_min_lon_on  < lon) & (lon < shuttle_max_lon_on)
        scan[off_mask]                 = 'OFF'
        scan[on_mask]                  = 'SCAN'
        scan[(~off_mask) & (~on_mask)] = 'JUNK'

    response = readout['Tsignal']
        
    # findR
    if True:
        mask = np.where(response[:, ch] >= Rth)
        scan[mask] = 'R'
            
    ms = MS.new(
        data=response,
        time=time,
        chan=obsinfo['kidids'][0].astype(np.int64),
        scan=scan,
        lon=lon,
        lat=lat,
        lon_origin=lon_origin,
        lat_origin=lat_origin,
        temperature=temperature,
        pressure=pressure,
        humidity=humidity,
        wind_speed=wind_speed,
        wind_direction=wind_direction,
        aste_cabin_temperature=aste_cabin_temperature,

        d2_mkid_id=obsinfo['kidids'][0].astype(np.int64),
        d2_mkid_type=obsinfo['kidtypes'][0],
        d2_mkid_frequency=obsinfo['kidfreqs'][0].astype(np.float64),

        exposure=obsinfo['integtime'][0],
        interval=obsinfo['interval'][0],
        beam_major=obsinfo['beamsize'][0],
        beam_minor=obsinfo['beamsize'][0],
        beam_pa=obsinfo['beamsize'][0],

        observer=obsheader['OBSERVER'],
        object=obsheader['OBJECT'],
        telescope_name=obsheader['TELESCOP'],
    )
    return ms

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
    table = ascii.read(filename)

    # 日付と時刻を取得して文字列でタイムスタンプを作成しそれをnumpy.datetime64へ変換する
    # テーブルの1列目と2列目がそれぞれ日付と時刻
    datetimes = []
    for date, time in zip(table['col1'], table['col2']):
        s = '{}T{}'.format(date, time)
        s = s.replace('/', '-')
        datetimes.append(s)
    datetimes       = np.array(datetimes).astype(np.datetime64)
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
    """
    return

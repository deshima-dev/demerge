"""dfitsをdemsへ変換するプログラム

File name: dfits2dems.py
Python 3.9
(C) 2023 内藤システムズ
"""

import numpy as np

from astropy.io        import fits
from dems.d2           import MS
from scipy.interpolate import interp1d

def convert_dfits_to_dems(filename):
    coodinate='azel'
    loadmode=0
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

    chan         = obsinfo['kidids'][0].astype(np.int64)
    scan         = np.array(antenna['scantype'])

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

    # 補間のために時刻(年月日時分秒)を時間(秒)に変更する。READOUTの最初の時刻を基準とする。
    seconds         = (time         - time[0])/np.timedelta64(1, 's')
    seconds_antenna = (time_antenna - time[0])/np.timedelta64(1, 's')
    seconds_weather = (time_weather - time[0])/np.timedelta64(1, 's')
    seconds_cabin   = (time_cabin   - time[0])/np.timedelta64(1, 's')

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
    
    ms = MS.new(
        data=readout['Tsignal'],
        time=time,
        chan=chan,
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

"""demsオブジェクトを生成する

(C) 2023 内藤システムズ
"""
import numpy          as np
import xarray         as xa
import merge_function as mf
import dfits2dems

from scipy.interpolate import interp1d
from astropy.io        import fits, ascii
from dems.d2           import MS

def merge_to_dems(
        ddbfits_path='',
        obsinst_path='',
        antenna_path='',
        readout_path='',
        skychop_path='',
        weather_path='',
        misti_path='',
        cabin_path='',
        pixel_id=0,
        coordinate='azel',
        loadmode=0):

    readout_hdul   = fits.open(readout_path)
    ddbfits_hdul   = fits.open(ddbfits_path)
    weather_table  = ascii.read(weather_path)
    antenna_table  = ascii.read(antenna_path)[:-1] # 最後の1行は終端を表す意味のないデータが入っているため無視する
    obsinst_params = mf.load_obsinst(obsinst_path) # 観測スクリプトに含まれているパラメタを抽出する
    times_misti, az_misti, el_misti, pwv_misti = dfits2dems.retrieve_misti_log(misti_path)

    times_cabin, upper_cabin_temp, lower_cabin_temp = dfits2dems.retrieve_cabin_temps(cabin_path)
    lower_cabin_temp = lower_cabin_temp + 273.15 # 度CからKへ変換

    times_weather = mf.convert_asciitime(weather_table['time'], '%Y-%m-%dT%H:%M:%S.%f')
    times_weather = np.array(times_weather).astype(np.datetime64)

    times_skychop, states_skychop = dfits2dems.retrieve_skychop_states(skychop_path)
    times_skychop = mf.convert_timestamp(times_skychop)
    times_skychop = np.array(times_skychop).astype(np.datetime64)

    # モードに応じて経度(lon)と緯度(lat)を選択(azelかradecか)する
    if coordinate == 'azel':
        if   'az-prg(no-cor)' in antenna_table.colnames:
            az_prg  = antenna_table['az-prg(no-cor)']
            el_prog = antenna_table['el-prog(no-cor)']
        elif 'az-prg(no-col)' in antenna_table.colnames:
            az_prg  = antenna_table['az-prg(no-col)']
            el_prog = antenna_table['el-prog(no-col)']
        else:
            raise KeyError('{}ファイルにaz-prg(no-cor)列とaz-prg(no-col)列がありません。どちらか片方が存在する必要があります。ANTENNAファイルの形式を確認してください。'.format(antenna_path))
        
        lon        = az_prg  + antenna_table['az-real'] - antenna_table['az-prg']
        lat        = el_prog + antenna_table['el-real'] - antenna_table['el-prg']
        lon_origin = np.median(antenna_table['az-prog(center)'])
        lat_origin = np.median(antenna_table['el-prog(center)'])
        if loadmode in [0, 1]:
            lon -= antenna_table['az-prog(center)']
            lat -= antenna_table['el-prog(center)']
            if loadmode == 0:
                lon *= np.cos(np.deg2rad(lat))
    elif coordinate == 'radec':
        lon        = antenna_table['ra-prg']
        lat        = antenna_table['dec-prg']
        lon_origin = obsinst_params['ra'] # 観測スクリプトに設定されているRA,DEC
        lat_origin = obsinst_params['dec']
    else:
        raise KeyError('Invalid coodinate type: {}'.format(coodinate))
    
    times_antenna = mf.convert_asciitime(antenna_table['time'], '%Y-%m-%dT%H:%M:%S.%f')
    times_antenna = np.array(times_antenna).astype(np.datetime64)

    times = mf.convert_timestamp(readout_hdul['READOUT'].data['timestamp'])
    times = np.array(times).astype(np.datetime64)

    # 補間のために時刻(年月日時分秒)を時間(秒)に変更する。READOUTの最初の時刻を基準とする。
    seconds         = (times         - times[0])/np.timedelta64(1, 's')
    seconds_cabin   = (times_cabin   - times[0])/np.timedelta64(1, 's')
    seconds_antenna = (times_antenna - times[0])/np.timedelta64(1, 's')
    seconds_weather = (times_weather - times[0])/np.timedelta64(1, 's')
    seconds_misti   = (times_misti   - times[0])/np.timedelta64(1, 's')
    seconds_skychop = (times_skychop - times[0])/np.timedelta64(1, 's')

    lon                    = np.interp(seconds, seconds_antenna, lon)
    lat                    = np.interp(seconds, seconds_antenna, lat)
    temperature            = np.interp(seconds, seconds_weather, weather_table['tmperature'])
    humidity               = np.interp(seconds, seconds_weather, weather_table['vapor-pressure'])
    pressure               = np.interp(seconds, seconds_weather, weather_table['presure'])
    wind_speed             = np.interp(seconds, seconds_weather, weather_table['aux1'])
    wind_direction         = np.interp(seconds, seconds_weather, weather_table['aux2'])
    skychop_state          = np.interp(seconds, seconds_skychop, states_skychop)
    aste_cabin_temperature = np.interp(seconds, seconds_cabin,   lower_cabin_temp)
    aste_misti_lon         = np.interp(seconds, seconds_misti,   az_misti)
    aste_misti_lat         = np.interp(seconds, seconds_misti,   el_misti)
    aste_misti_pwv         = np.interp(seconds, seconds_misti,   pwv_misti)

    states = np.array(antenna_table['type'])

    # 補間関数で扱うためにSCANTYPE(文字列)を整数に置き換える
    state_types = {state_type:i for i, state_type in enumerate(np.unique(states))}
    state_type_numbers = np.zeros(states.shape[0], dtype=int)
    for state_type, i in state_types.items():
        state_type_numbers[states == state_type] = i

    # READOUTの時間に合わせるためにSTATETYPEを補間する
    f_statetype = interp1d(seconds_antenna, state_type_numbers, kind='nearest', bounds_error=False, fill_value=(state_type_numbers[0], state_type_numbers[-1]))
    state_type_numbers = f_statetype(seconds)

    # 補間後のSTATETYPEを文字列に戻す
    state = np.full_like(state_type_numbers, 'GRAD', dtype='<U8')
    for state_type, i in state_types.items():
        state[state_type_numbers == i] = state_type

    T_amb     = np.nanmean(weather_table['tmperature']) + 273.15
    n_kid     = readout_hdul['KIDSINFO'].header['NKID{}'.format(pixel_id)]
    T_signals = mf.calibrate_to_power(pixel_id, lower_cabin_temp[0], T_amb, readout_hdul, ddbfits_hdul)

    master_id, kid_id, kid_type, kid_freq, kid_Q = mf.get_maskid_corresp(pixel_id, ddbfits_hdul)

    ddbfits_hdul.close()
    readout_hdul.close()

    return MS.new(
        data=T_signals,
        time=times,
        chan=kid_id,
        state=state,
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
        aste_misti_lon=aste_misti_lon,
        aste_misti_lat=aste_misti_lat,
        aste_misti_pwv=aste_misti_pwv,
        d2_mkid_id=kid_id,
        d2_mkid_type=kid_type,
        d2_mkid_frequency=kid_freq,
        d2_skychopper_isblocking=skychop_state,
        beam_major=0.005, # 18 arcsec MergeToDfits()でも固定値が指定されていた
        beam_minor=0.005, # 18 arcsec MergeToDfits()でも固定値が指定されていた
        beam_pa=0.005,    # 18 arcsec MergeToDfits()でも固定値が指定されていた
        exposure=1./196,  #           MergeToDfits()でも固定値が指定されていた
        interval=1./196,  #           MergeToDfits()でも固定値が指定されていた
        observation=obsinst_params['observation'],
        observer=obsinst_params['observer'],
        object=obsinst_params['obs_object'],
    )

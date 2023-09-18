"""demsオブジェクトを生成する

python 3.9
dems   0.4.0

(C) 2023 内藤システムズ
"""
import numpy          as np
import xarray         as xr
import merge_function as mf

from astropy.io import fits, ascii
from dems.d2    import MS

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
        loadmode=0,
        **kwargs
        ):
    # kwargsの処理と既定値の設定

    # 時刻と各種データを読み込む(必要に応じて時刻はnp.datetime64[ns]へ変換する)
    readout_hdul   = fits.open(readout_path)
    ddbfits_hdul   = fits.open(ddbfits_path)
    weather_table  = ascii.read(weather_path)
    antenna_table  = ascii.read(antenna_path)[:-1] # 最後の1行は終端を表す意味のないデータが入っているため無視する
    obsinst_params = mf.load_obsinst(obsinst_path) # 観測スクリプトに含まれているパラメタを抽出する

    times = mf.convert_timestamp(readout_hdul['READOUT'].data['timestamp'])
    times = np.array(times).astype('datetime64[ns]')

    times_misti, az_misti, el_misti, pwv_misti = mf.retrieve_misti_log(misti_path)

    times_cabin, upper_cabin_temp, lower_cabin_temp = mf.retrieve_cabin_temps(cabin_path)
    lower_cabin_temp = lower_cabin_temp + 273.15 # 度CからKへ変換

    times_weather = mf.convert_asciitime(weather_table['time'], '%Y-%m-%dT%H:%M:%S.%f')
    times_weather = np.array(times_weather).astype('datetime64[ns]')

    times_skychop, states_skychop = mf.retrieve_skychop_states(skychop_path)
    times_skychop = mf.convert_timestamp(times_skychop)
    times_skychop = np.array(times_skychop).astype('datetime64[ns]')

    times_antenna = mf.convert_asciitime(antenna_table['time'], '%Y-%m-%dT%H:%M:%S.%f')
    times_antenna = np.array(times_antenna).astype('datetime64[ns]')

    # T_signalsを計算する
    master_id, kid_id, kid_type, kid_freq, kid_Q = mf.get_maskid_corresp(pixel_id, ddbfits_hdul)
    T_amb     = np.nanmean(weather_table['tmperature']) + 273.15 # 度CからKへ変換
    T_signals = mf.calibrate_to_power(pixel_id, lower_cabin_temp[0], T_amb, readout_hdul, ddbfits_hdul)

    ddbfits_hdul.close()
    readout_hdul.close()

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

    # 補間関数で扱うためにSCANTYPE(文字列)を適当な整数に対応させる
    states = np.array(antenna_table['type'])
    state_types = {state_type:i for i, state_type in enumerate(np.unique(states))}
    state_type_numbers = np.zeros(states.shape[0], dtype=int)
    for state_type, i in state_types.items():
        state_type_numbers[states == state_type] = i

    # 補間のためにDataArrayへ格納する
    T_signals_xr              = xr.DataArray(data=T_signals, dims=['time', 'chan'], coords=[times, kid_id])
    lon_xr                    = xr.DataArray(data=lon,                              coords={'time': times_antenna})
    lat_xr                    = xr.DataArray(data=lat,                              coords={'time': times_antenna})
    temperature_xr            = xr.DataArray(data=weather_table['tmperature'],      coords={'time': times_weather})
    humidity_xr               = xr.DataArray(data=weather_table['vapor-pressure'],  coords={'time': times_weather})
    pressure_xr               = xr.DataArray(data=weather_table['presure'],         coords={'time': times_weather})
    wind_speed_xr             = xr.DataArray(data=weather_table['aux1'],            coords={'time': times_weather})
    wind_direction_xr         = xr.DataArray(data=weather_table['aux2'],            coords={'time': times_weather})
    skychop_state_xr          = xr.DataArray(data=states_skychop,                   coords={'time': times_skychop})
    aste_cabin_temperature_xr = xr.DataArray(data=lower_cabin_temp,                 coords={'time': times_cabin})
    aste_misti_lon_xr         = xr.DataArray(data=az_misti,                         coords={'time': times_misti})
    aste_misti_lat_xr         = xr.DataArray(data=el_misti,                         coords={'time': times_misti})
    aste_misti_pwv_xr         = xr.DataArray(data=pwv_misti,                        coords={'time': times_misti})
    state_type_numbers_xr     = xr.DataArray(data=state_type_numbers,               coords={'time': times_antenna})

    # Tsignalsの時刻に合わせて補間する
    lon                    =                    lon_xr.interp_like(T_signals_xr, kwargs={'fill_value': 'extrapolate'})
    lat                    =                    lat_xr.interp_like(T_signals_xr, kwargs={'fill_value': 'extrapolate'})
    temperature            =            temperature_xr.interp_like(T_signals_xr, kwargs={'fill_value': 'extrapolate'})
    humidity               =               humidity_xr.interp_like(T_signals_xr, kwargs={'fill_value': 'extrapolate'})
    pressure               =               pressure_xr.interp_like(T_signals_xr, kwargs={'fill_value': 'extrapolate'})
    wind_speed             =             wind_speed_xr.interp_like(T_signals_xr, kwargs={'fill_value': 'extrapolate'})
    wind_direction         =         wind_direction_xr.interp_like(T_signals_xr, kwargs={'fill_value': 'extrapolate'})
    skychop_state          =          skychop_state_xr.interp_like(T_signals_xr, kwargs={'fill_value': 'extrapolate'})
    aste_cabin_temperature = aste_cabin_temperature_xr.interp_like(T_signals_xr, kwargs={'fill_value': 'extrapolate'})
    aste_misti_lon         =         aste_misti_lon_xr.interp_like(T_signals_xr, kwargs={'fill_value': 'extrapolate'})
    aste_misti_lat         =         aste_misti_lat_xr.interp_like(T_signals_xr, kwargs={'fill_value': 'extrapolate'})
    aste_misti_pwv         =         aste_misti_pwv_xr.interp_like(T_signals_xr, kwargs={'fill_value': 'extrapolate'})
    state_type_numbers     =     state_type_numbers_xr.interp_like(T_signals_xr, kwargs={'fill_value': 'extrapolate'}, method='nearest')

    # 補間後のSTATETYPEを文字列に戻す
    state = np.full_like(state_type_numbers, 'GRAD', dtype='<U8')
    for state_type, i in state_types.items():
        state[state_type_numbers == i] = state_type

    # find R
    findR = False
    if findR:
        None

    return MS.new(
        data                    =T_signals,
        time                    =times,
        chan                    =kid_id,
        state                   =state,
        lon                     =lon,
        lat                     =lat,
        lon_origin              =lon_origin,
        lat_origin              =lat_origin,
        temperature             =temperature,
        pressure                =pressure,
        humidity                =humidity,
        wind_speed              =wind_speed,
        wind_direction          =wind_direction,
        aste_cabin_temperature  =aste_cabin_temperature,
        aste_misti_lon          =aste_misti_lon,
        aste_misti_lat          =aste_misti_lat,
        aste_misti_pwv          =aste_misti_pwv,
        d2_mkid_id              =kid_id,
        d2_mkid_type            =kid_type,
        d2_mkid_frequency       =kid_freq,
        d2_skychopper_isblocking=skychop_state,
        beam_major              =0.005,  # 18 arcsec MergeToDfits()でも固定値が指定されていた
        beam_minor              =0.005,  # 18 arcsec MergeToDfits()でも固定値が指定されていた
        beam_pa                 =0.005,  # 18 arcsec MergeToDfits()でも固定値が指定されていた
        exposure                =1./196, #           MergeToDfits()でも固定値が指定されていた
        interval                =1./196, #           MergeToDfits()でも固定値が指定されていた
        observation             =obsinst_params['observation'],
        observer                =obsinst_params['observer'],
        object                  =obsinst_params['obs_object'],
    )

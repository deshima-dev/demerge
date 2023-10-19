"""demsオブジェクトを生成する

python 3.9
dems   0.4.0

(C) 2023 内藤システムズ
"""
import argparse
import numpy          as np
import xarray         as xr
from . import merge_function as mf

from astropy.io import fits, ascii
from dems.d2    import MS


# constants
from . import __version__ as DEMERGE_VERSION


def merge_to_dems(
        ddbfits_path='',
        obsinst_path='',
        antenna_path='',
        readout_path='',
        skychop_path='',
        weather_path='',
        misti_path='',
        cabin_path='',
        **kwargs
        ):
    # その他の引数の処理と既定値の設定
    pixel_id   = kwargs.pop('pixel_id',   0)
    coordinate = kwargs.pop('coordinate', 'azel')
    mode       = kwargs.pop('mode',   0)
    loadtype   = kwargs.pop('loadtype',   'Tsignal')
    # find R, sky
    findR  = kwargs.pop("findR",  False)
    ch     = kwargs.pop("ch",     0)
    Rth    = kwargs.pop("Rth",    280)
    skyth  = kwargs.pop("skyth",  150)
    cutnum = kwargs.pop("cutnum", 1)
    # still
    still  = kwargs.pop("still",  False)
    period = kwargs.pop("period", 2) # 秒
    # shuttle
    shuttle     = kwargs.pop("shuttle", False)
    lon_min_off = kwargs.pop("lon_min_off", 0)
    lon_max_off = kwargs.pop("lon_max_off", 0)
    lon_min_on  = kwargs.pop("lon_min_on",  0)
    lon_max_on  = kwargs.pop("lon_max_on",  0)

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

    response = None
    if loadtype == 'Tsignal':
        # T_signalsを計算する
        master_id, kid_id, kid_type, kid_freq, kid_Q = mf.get_maskid_corresp(pixel_id, ddbfits_hdul)
        T_amb     = np.nanmean(weather_table['tmperature']) + 273.15 # 度CからKへ変換
        T_signals = mf.calibrate_to_power(pixel_id, lower_cabin_temp[0], T_amb, readout_hdul, ddbfits_hdul)
        response  = T_signals
    else:
        raise KeyError('Invalid loadtype: {}'.format(loadtype))

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
        if mode in [0, 1]:
            lon -= antenna_table['az-prog(center)']
            lat -= antenna_table['el-prog(center)']
            if mode == 0:
                lon *= np.cos(np.deg2rad(lat))
    elif coordinate == 'radec':
        lon        = antenna_table['ra-prg']
        lat        = antenna_table['dec-prg']
        lon_origin = obsinst_params['ra'] # 観測スクリプトに設定されているRA,DEC
        lat_origin = obsinst_params['dec']
    else:
        raise KeyError('Invalid coodinate type: {}'.format(coordinate))

    # 補間関数で扱うためにSCANTYPE(文字列)を適当な整数に対応させる
    states = np.array(antenna_table['type'])
    state_types = {state_type:i for i, state_type in enumerate(np.unique(states))}
    state_type_numbers = np.zeros(states.shape[0], dtype=int)
    for state_type, i in state_types.items():
        state_type_numbers[states == state_type] = i

    # 補間のためにDataArrayへ格納する
    response_xr               = xr.DataArray(data=response, dims=['time', 'chan'], coords=[times, kid_id])
    lon_xr                    = xr.DataArray(data=lon,                             coords={'time': times_antenna})
    lat_xr                    = xr.DataArray(data=lat,                             coords={'time': times_antenna})
    temperature_xr            = xr.DataArray(data=weather_table['tmperature'],     coords={'time': times_weather})
    humidity_xr               = xr.DataArray(data=weather_table['vapor-pressure'], coords={'time': times_weather})
    pressure_xr               = xr.DataArray(data=weather_table['presure'],        coords={'time': times_weather})
    wind_speed_xr             = xr.DataArray(data=weather_table['aux1'],           coords={'time': times_weather})
    wind_direction_xr         = xr.DataArray(data=weather_table['aux2'],           coords={'time': times_weather})
    skychop_state_xr          = xr.DataArray(data=states_skychop,                  coords={'time': times_skychop})
    aste_cabin_temperature_xr = xr.DataArray(data=lower_cabin_temp,                coords={'time': times_cabin})
    aste_misti_lon_xr         = xr.DataArray(data=az_misti,                        coords={'time': times_misti})
    aste_misti_lat_xr         = xr.DataArray(data=el_misti,                        coords={'time': times_misti})
    aste_misti_pwv_xr         = xr.DataArray(data=pwv_misti,                       coords={'time': times_misti})
    state_type_numbers_xr     = xr.DataArray(data=state_type_numbers,              coords={'time': times_antenna})

    # Tsignalsの時刻に合わせて補間する
    lon                    =                    lon_xr.interp_like(response_xr)
    lat                    =                    lat_xr.interp_like(response_xr)
    temperature            =            temperature_xr.interp_like(response_xr)
    humidity               =               humidity_xr.interp_like(response_xr)
    pressure               =               pressure_xr.interp_like(response_xr)
    wind_speed             =             wind_speed_xr.interp_like(response_xr)
    wind_direction         =         wind_direction_xr.interp_like(response_xr)
    aste_cabin_temperature = aste_cabin_temperature_xr.interp_like(response_xr)
    aste_misti_lon         =         aste_misti_lon_xr.interp_like(response_xr)
    aste_misti_lat         =         aste_misti_lat_xr.interp_like(response_xr)
    aste_misti_pwv         =         aste_misti_pwv_xr.interp_like(response_xr)
    skychop_state          =          skychop_state_xr.interp_like(response_xr, method='nearest')
    state_type_numbers     =     state_type_numbers_xr.interp_like(response_xr, method='nearest')

    # 補間後のSTATETYPEを文字列に戻す
    state = np.full_like(state_type_numbers, 'GRAD', dtype='<U8')
    for state_type, i in state_types.items():
        state[state_type_numbers == i] = state_type

    # 静止データの周期に応じてOFFマスクとSCANマスクを設定する
    if still:
        seconds = (times - times[0])/np.timedelta64(1, 's')
        for i in range(int(seconds[-1]) // period + 1):
            off_mask = (period*(2*i)     <= seconds) & (seconds < period*(2*i + 1))
            on_mask  = (period*(2*i + 1) <= seconds) & (seconds < period*(2*i + 2))
            state[off_mask] = 'OFF'
            state[on_mask]  = 'SCAN'

    # shuttle観測のマスクを設定する
    if shuttle:
        mask_off = (lon_min_off < lon) & (lon < lon_max_off)
        mask_on  = (lon_min_on  < lon) & (lon < lon_max_on)
        state[mask_off]                 = 'OFF'
        state[mask_on]                  = 'SCAN'
        state[(~mask_off) & (~mask_on)] = 'JUNK'

    # Rとskyの部分を探し、その変化点も含めてJUNKな部分を調べる。
    if findR:
        # Rの部分とその変化の部分を探す
        indices = np.where(response[:, ch] >= Rth)
        state[indices] = 'R'

        # cutnum個だけ左右を切り取った配列を作り、互いに異なる部分を探す。そこはおおよそ変化が起きている部分と考えられる。
        #
        # cutnum = 2 の例
        #
        # ** 状態Rへ変化していく場合 **
        # state                   = 0 0 0 0 0 1 1 1 1 1
        # [cutnum: ]              = 0 0 0 1 1 1 1 1
        # [:-cutnum]              = 0 0 0 0 0 1 1 1
        # [cutnum:] != [:-cutnum] = 0 0 0 1 1 0 0 0
        # state_right_shift       = 0 0 0 0 0 1 1 0 0 0
        # state_left_shift        = 0 0 0 1 1 0 0 0 0 0
        # state_R                 = 0 0 0 0 0 1 1 1 1 1
        # mask_moving             = 0 0 0 0 0 1 1 0 0 0
        #
        # ** 状態Rから別の状態へ変化していく場合 **
        # state                   = 1 1 1 1 1 0 0 0 0 0
        # [cutnum: ]              = 1 1 1 0 0 0 0 0
        # [:-cutnum]              = 1 1 1 1 1 0 0 0
        # [cutnum:] != [:-cutnum] = 0 0 0 1 1 0 0 0
        # state_right_shift       = 0 0 0 0 0 1 1 0 0 0
        # state_left_shift        = 0 0 0 1 1 0 0 0 0 0
        # state_R                 = 1 1 1 1 1 0 0 0 0 0
        # mask_moving             = 0 0 0 1 1 1 1 0 0 0
        #
        # 状態がRへ変化する場合と、状態Rから別の状態へ変化する場合でmask_movingのでき方が違う。
        #
        state_cut = state[cutnum:] != state[:-cutnum]

        state_right_shift = np.hstack( [[False]*cutnum, state_cut] ) # 左側をFalseで埋めて右にずらす
        state_left_shift  = np.hstack( [state_cut, [False]*cutnum] ) # 右側をFalseで埋めて左にずらす
        state_R           =          ( state == 'R' )

        mask_moving = state_R & state_left_shift | state_right_shift
        state[mask_moving] = 'JUNK'

        indices = (response[:, ch] >  skyth) & (state != 'R')
        state[indices] = 'JUNK'

        indices = (response[:, ch] <= skyth) & (state == 'R')
        state[indices] = 'JUNK'

        # SKYの部分とその変化の部分を探す
        indices      = np.where(response[:, ch] <= skyth)
        tmp          = state.copy() # 最終的にSKYを残さないためにコピーを扱う
        tmp[indices] = 'SKY'        # 一時的にSKYをマークする

        tmp_cut = tmp[cutnum:] != tmp[:-cutnum] # cutnum個だけ左右にずらした配列を作り、変化を探す。

        tmp_right_shift = np.hstack( [[False]*cutnum, tmp_cut] )
        tmp_left_shift  = np.hstack( [tmp_cut, [False]*cutnum] )
        tmp_sky         =          ( tmp == 'SKY' )

        mask_moving = tmp_sky & tmp_left_shift | tmp_right_shift
        state[mask_moving] = 'JUNK' # 変化の部分はJUNKに置き換える(Rとは違いSKYは残らない)

    return MS.new(
        data                    =response,
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
        d2_demerge_version       =DEMERGE_VERSION,
        beam_major              =0.005,  # 18 arcsec MergeToDfits()でも固定値が指定されていた
        beam_minor              =0.005,  # 18 arcsec MergeToDfits()でも固定値が指定されていた
        beam_pa                 =0.005,  # 18 arcsec MergeToDfits()でも固定値が指定されていた
        exposure                =1./196, #           MergeToDfits()でも固定値が指定されていた
        interval                =1./196, #           MergeToDfits()でも固定値が指定されていた
        observation             =obsinst_params['observation'],
        observer                =obsinst_params['observer'],
        object                  =obsinst_params['obs_object'],
    )

def main() -> None:
    """Demsオブジェクトを作成する"""
    parser = argparse.ArgumentParser()

    # 必須引数
    parser.add_argument('filename',  type=str, help='出力ファイルへのパスを指定して下さい(.nc)')
    parser.add_argument('--ddb',     type=str, required=True, help='DDBファイルへのパスを指定して下さい(.fits.gz)')
    parser.add_argument('--obs',     type=str, required=True, help='obsファイルへのパスを指定して下さい(.obs)')
    parser.add_argument('--antenna', type=str, required=True, help='antennaファイルへのパスを指定して下さい(.antenna)')
    parser.add_argument('--readout', type=str, required=True, help='reduced readoutファイルへのパスを指定して下さい(.fits)')
    parser.add_argument('--skychop', type=str, required=True, help='skychopファイルへのパスを指定して下さい(.skychop)')
    parser.add_argument('--weather', type=str, required=True, help='weatherファイルへのパスを指定して下さい(.weather)')
    parser.add_argument('--misti',   type=str, required=True, help='mistiファイルへのパスを指定して下さい(.misti)')
    parser.add_argument('--cabin',   type=str, required=True, help='cabinファイルへのパスを指定して下さい(.cabin)')

    # オプション引数
    parser.add_argument('--pixel_id',    type=int,   default=0,         help='pixel_idを整数で指定します')
    parser.add_argument('--coordinate',  type=str,   default='azel',    help='座標系(azel/radec)を文字列で指定します')
    parser.add_argument('--mode',        type=int,   default=0,         help='座標の読み込み方法を整数で指定します(0:相対座標cos射影あり(おすすめ既定値), 1:相対座標cos射影なし, 2:絶対座標)')
    parser.add_argument('--loadtype',    type=str,   default='Tsignal', help='読み込むデータを文字列で指定します(既定値: Tsignal, 現在はTsignalの読み込みしかできません)')
    parser.add_argument('--findR',       type=bool,  default=False,     help='findRにTrueを指定するとFindR, Skyを実行します')
    parser.add_argument('--ch',          type=int,   default=0,         help='findRに利用するチャネルを整数で指定します')
    parser.add_argument('--Rth',         type=float, default=280.0,     help='R閾値を実数で指定します')
    parser.add_argument('--skyth',       type=float, default=150.0,     help='sky閾値を実数で指定します')
    parser.add_argument('--cutnum',      type=int,   default=1,         help='findRでのカット数を整数で指定します')
    parser.add_argument('--still',       type=bool,  default=False,     help='Trueに設定するとstill観測用の解析を行います')
    parser.add_argument('--period',      type=int,   default=2,         help='still観測の1/2周期(秒)を整数で指定します')
    parser.add_argument('--shuttle',     type=bool,  default=False,     help='Trueを指定するとshuttle観測用の解析を行います')
    parser.add_argument('--lon_min_off', type=float, default=0.0,       help='shuttle観測時のOFFにするlongitudeの最小値を実数で指定します')
    parser.add_argument('--lon_max_off', type=float, default=0.0,       help='shuttle観測時のOFFにするlongitudeの最大値を実数で指定します')
    parser.add_argument('--lon_min_on',  type=float, default=0.0,       help='shuttle観測時のONにするlongitudeの最小値を実数で指定します')
    parser.add_argument('--lon_max_on',  type=float, default=0.0,       help='shuttle観測時のONにするlongitudeの最大値を実数で指定します')

    a = parser.parse_args()

    dems = merge_to_dems(
        ddbfits_path=a.ddb,
        obsinst_path=a.obs,
        antenna_path=a.antenna,
        readout_path=a.readout,
        skychop_path=a.skychop,
        weather_path=a.weather,
        misti_path  =a.misti,
        cabin_path  =a.cabin,

        pixel_id   =a.pixel_id,
        coordinate =a.coordinate,
        mode       =a.mode,
        loadtype   =a.loadtype,
        findR      =a.findR,
        ch         =a.ch,
        Rth        =a.Rth,
        skyth      =a.skyth,
        cutnum     =a.cutnum,
        still      =a.still,
        period     =a.period,
        shuttle    =a.shuttle,
        lon_min_off=a.lon_min_off,
        lon_max_off=a.lon_max_off,
        lon_min_on =a.lon_min_on,
        lon_max_on =a.lon_max_on,
    )

    dems.to_netcdf(a.filename)

if __name__ == '__main__':
    main()

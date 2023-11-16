"""merge_to_dems()関数をテストする

Python 3.9
dems   0.4.0

テスト用のデータについて
========================
同ディレクトリにあるmktd.shを実行すると本テストに必要なダミーデータが生成される。
 $ ./mktd.sh

(C) 2023 内藤システムズ
"""
import unittest
import numpy as np
import xarray as xa

from demerge import merge_function as mf
from demerge import merge_to_dems as mtd

from astropy.io import fits, ascii
from datetime import datetime


# constants
from dems import __version__ as DEMS_VERSION
from demerge import __version__ as DEMERGE_VERSION


class MergeToDemsTestDrive(unittest.TestCase):
    """merge_to_dems.pyモジュールの単体テスト"""
    def setUp(self):
        self.obsid = '20171103184436'
        self.reduced_fits = 'cache/{0}/reduced_{0}.fits'.format(self.obsid)
        self.dfits = 'cache/{0}/dfits_{0}.fits.gz'.format(self.obsid)
        return

    def test_short_measure(self):
        """READOUT時間が環境観測時間よりも長い場合に補間に与える影響を確認する"""
        prefix = 'testdata_short_measure'
        dems = mtd.merge_to_dems(
            ddbfits_path='{}_DDB.fits.gz'.format(prefix),
            obsinst_path='../data/cosmos/cosmos_{0}/{0}.obs'.format(self.obsid),
            antenna_path='{}.ant'.format(prefix),
            readout_path='{}_reduced_readout.fits'.format(prefix),
            skychop_path='{}.skychop'.format(prefix),
            weather_path='{}.wea'.format(prefix),
            misti_path='{}.misti'.format(prefix),
            cabin_path='{}.cabin'.format(prefix),
            loadtype='Tsignal'
        )
        times_skychop, states_skychop = mf.retrieve_skychop_states('{}.skychop'.format(prefix))
        times_skychop = mf.convert_timestamp(times_skychop)
        times_skychop = np.array(times_skychop).astype('datetime64[ns]')

        antenna_table = ascii.read('{}.ant'.format(prefix))[:-1] # 最後の1行は終端を表す意味のないデータが入っているため無視する
        times_antenna = mf.convert_asciitime(antenna_table['time'], '%Y-%m-%dT%H:%M:%S.%f')
        times_antenna = np.array(times_antenna).astype('datetime64[ns]')

        # d2_skychopper_isblockingの外挿部分の論理値の確認
        #
        # readout time   |---------------------------------------------------------|
        #
        # skychop time   |--------------------------------------|
        #
        # skychop state  11111111111111111110000000000000000000011111111111111111111
        #                                この境界を確認する--> ~~
        #                                                       |<~~ exporalate ~~>|
        #
        # d2_skychopper_isblockingが外挿されるとそれはnumpy.nanになる。これがdemsのbool型に格納されるとTrueに変換される。
        # ダミーデータtestdata_short_measure.skychopは上図のようにreadout timeよりskychop timeの方が短くなっている。
        # またskychopは半分の時間で0に切り替わるようになっている。このため外挿が始まると0から1への切り替えが発生して、外挿部分がわかる。
        # まず、skychopの時刻とdems全体の時刻から外挿が発生すべき時刻のインデックスを調べる。
        # 次に、skychopがFalseになっている部分のインデックスを調べ、外挿が発生したインデックスの値がFalseのインデックス+1であることを確認する。
        # これで外挿された時の振る舞いが分かる。
        #
        self.assertEqual(len(dems.time.values), len(dems.d2_skychopper_isblocking.values), 'd2_skychopper_isblockingの要素数の確認')

        indices_extrapolate = np.where(dems.time > times_skychop[-1])[0]
        self.assertTrue((dems.d2_skychopper_isblocking.values[indices_extrapolate] == True).all(), 'd2_skychopper_isblockingの外挿部分がTrueになることを確認')

        indices_false = np.where(dems.d2_skychopper_isblocking == False)[0]
        self.assertEqual(indices_false[-1] + 1, indices_extrapolate[0], '外挿部分と内挿部分(Falseになっている部分)がk隣り合っていることを確認')

        # stateの外挿部分の値の確認(確認の方法はd2_skychopper_isblockingと同じ)
        self.assertEqual(len(dems.time.values), len(dems.state.values), 'stateの要素数の確認')

        indices_extrapolate = np.where(dems.time > times_antenna[-1])[0]
        self.assertTrue((dems.state[indices_extrapolate] == 'GRAD').all(), 'stateの外挿部分がGRADになることを確認')

        indices_on = np.where(dems.state == 'ON')[0]
        self.assertEqual(indices_on[-1] + 1, indices_extrapolate[0], 'stateの外挿部分と内挿部分(Falseになっている部分)が隣り合っていることを確認')

        return

    def test_shuttle(self):
        lon_min_off = 20;
        lon_max_off = 80;
        lon_min_on  = 100;
        lon_max_on  = 150;
        prefix = 'testdata_linear_antenna'
        dems = mtd.merge_to_dems(
            ddbfits_path='{}_DDB.fits.gz'.format(prefix),
            obsinst_path='../data/cosmos/cosmos_{0}/{0}.obs'.format(self.obsid),
            antenna_path='{}.ant'.format(prefix),
            readout_path='{}_reduced_readout.fits'.format(prefix),
            skychop_path='{}.skychop'.format(prefix),
            weather_path='{}.wea'.format(prefix),
            misti_path='{}.misti'.format(prefix),
            cabin_path='{}.cabin'.format(prefix),
            shuttle=True,
            lon_min_off=lon_min_off,
            lon_max_off=lon_max_off,
            lon_min_on=lon_min_on,
            lon_max_on=lon_max_on,
            loadtype='Tsignal'
        )
        state = dems.state
        time  = dems.time
        lon   = dems.lon
        indices = np.where((lon_min_off < lon) & (lon < lon_max_off))
        self.assertTrue(np.array(state[indices] == 'OFF').all(), 'OFFの区間を確認')

        indices = np.where((lon_min_on < lon) & (lon < lon_max_on))
        self.assertTrue(np.array(state[indices] == 'SCAN').all(), 'ONの区間を確認')

        indices = np.where(lon <= lon_min_off)
        self.assertTrue(np.array(state[indices] == 'JUNK').all(), 'JUNKの区間を確認(1)')

        indices = np.where((lon_max_off <= lon) & (lon <= lon_min_on))
        self.assertTrue(np.array(state[indices] == 'JUNK').all(), 'JUNKの区間を確認(2)')

        indices = np.where(lon_max_on <= lon)
        self.assertTrue(np.array(state[indices] == 'JUNK').all(), 'JUNKの区間を確認(3)')

        return

    def test_still(self):
        period = 2 # 秒
        prefix = 'testdata'
        dems = mtd.merge_to_dems(
            ddbfits_path='{}_DDB.fits.gz'.format(prefix),
            obsinst_path='../data/cosmos/cosmos_{0}/{0}.obs'.format(self.obsid),
            antenna_path='{}.ant'.format(prefix),
            readout_path='{}_reduced_readout.fits'.format(prefix),
            skychop_path='{}.skychop'.format(prefix),
            weather_path='{}.wea'.format(prefix),
            misti_path='{}.misti'.format(prefix),
            cabin_path='{}.cabin'.format(prefix),
            still=True,
            period=period,
            loadtype='Tsignal'
        )
        state = dems.state
        time  = dems.time

        # period毎にOFF/SCANを繰り返す
        # GRADはなくなる
        self.assertTrue(np.array(state != 'GRAD').all(), 'GRADは無いことを確認')

        seconds = (time - time[0])/np.timedelta64(1, 's')
        off = (seconds < period)
        on  = (period <= seconds) & (seconds < 2*period)
        self.assertTrue(np.array(state[off] == 'OFF').all(),  'OFFであることを確認(1)')
        self.assertTrue(np.array(state[on]  == 'SCAN').all(), 'SCANであることを確認(1)')
        off = (2*period <= seconds) & (seconds < 3*period)
        on  = (3*period <= seconds) & (seconds < 4*period)
        self.assertTrue(np.array(state[off] == 'OFF').all(),  'OFFであることを確認(2)')
        self.assertTrue(np.array(state[on]  == 'SCAN').all(), 'SCANであることを確認(2)')
        return

    def test_find_R_linear_inc(self):
        cutnum = 2
        prefix = 'testdata_linear_inc'
        dems = mtd.merge_to_dems(
            ddbfits_path='{}_DDB.fits.gz'.format(prefix),
            obsinst_path='../data/cosmos/cosmos_{0}/{0}.obs'.format(self.obsid),
            antenna_path='{}.ant'.format(prefix),
            readout_path='{}_reduced_readout.fits'.format(prefix),
            skychop_path='{}.skychop'.format(prefix),
            weather_path='{}.wea'.format(prefix),
            misti_path='{}.misti'.format(prefix),
            cabin_path='{}.cabin'.format(prefix),
            findR=True,
            cutnum=cutnum,
            loadtype='Tsignal',
        )
        response = dems.data.T[0]
        state    = dems.state

        # skyth以下の部分を探す
        skyth = 150
        indices = np.where(response <= skyth)[0]
        i_sky = indices[-1]

        # Rthを越える部分を探す
        Rth = 280
        indices = np.where(response > Rth)[0]
        i_R = indices[0]

        self.assertTrue(np.array(state[i_R + 2*cutnum:         ] == 'R').all(),    'Rの部分を確認')
        self.assertTrue(np.array(state[i_sky         :i_R      ] == 'JUNK').all(), 'JUNKの部分を確認')
        self.assertTrue(np.array(state[              :i_sky - 1] == 'GRAD').all(), 'GRADの部分を確認')

        self.assertTrue(np.array(state != 'SKY').all(), 'SKYの部分は無いことを確認')

        return

    def test_find_R_linear_dec(self):
        cutnum = 2
        prefix = 'testdata_linear_dec'
        dems = mtd.merge_to_dems(
            ddbfits_path='{}_DDB.fits.gz'.format(prefix),
            obsinst_path='../data/cosmos/cosmos_{0}/{0}.obs'.format(self.obsid),
            antenna_path='{}.ant'.format(prefix),
            readout_path='{}_reduced_readout.fits'.format(prefix),
            skychop_path='{}.skychop'.format(prefix),
            weather_path='{}.wea'.format(prefix),
            misti_path='{}.misti'.format(prefix),
            cabin_path='{}.cabin'.format(prefix),
            findR=True,
            cutnum=cutnum,
            loadtype='Tsignal',
        )
        response = dems.data.T[0]
        state    = dems.state

        # skyth以下の部分を探す
        skyth = 150
        indices = np.where(response <= skyth)[0]
        i_sky = indices[0]

        # Rthを越える部分を探す
        Rth = 280
        indices = np.where(response > Rth)[0]
        i_R = indices[-1]

        self.assertTrue(np.array(state[              :i_R - cutnum] == 'R').all(),    'Rの部分を確認')
        self.assertTrue(np.array(state[i_R           :i_sky       ] == 'JUNK').all(), 'JUNKの部分を確認')
        self.assertTrue(np.array(state[i_sky + cutnum:            ] == 'GRAD').all(), 'GRADの部分を確認')

        self.assertTrue(np.array(state != 'SKY').all(), 'SKYの部分は無いことを確認')

        return

    def test_merge_to_dems(self):
        prefix = 'testdata'

        n_time = None
        n_kid  = None
        times  = None
        with fits.open('testdata_reduced_readout.fits') as hdul:
            n_time = hdul['READOUT'].header['NAXIS2']
            n_kid  = hdul['KIDSINFO'].header['NAXIS2']
            times  = np.array(mf.convert_timestamp(hdul['READOUT'].data['timestamp'])).astype(np.datetime64)
            kidids = np.array(hdul['KIDSINFO'].data['kidid']).astype(np.int64)

        dems = mtd.merge_to_dems(
            ddbfits_path='{}_DDB.fits.gz'.format(prefix),
            obsinst_path='../data/cosmos/cosmos_{0}/{0}.obs'.format(self.obsid),
            antenna_path='{}.ant'.format(prefix),
            readout_path='{}_reduced_readout.fits'.format(prefix),
            skychop_path='{}.skychop'.format(prefix),
            weather_path='{}.wea'.format(prefix),
            misti_path='{}.misti'.format(prefix),
            cabin_path='{}.cabin'.format(prefix),
            loadtype='Tsignal',
        )
        self.assertEqual(xa.DataArray, type(dems), 'demsオブジェクトの型の確認')

        #
        # MS::data(Tsignal)の計算の試験
        #
        # merge_function.pyに定義されているcalibrate_to_power()関数で計算される。
        #
        # f_shift = (linPh - linyfc)/(4*Qr)
        #
        # Tlos_model = ( f_shift + p0*sqrt(Troom + T0) )**2 / (p0**2 * etaf) - T0/etaf - ((1 - etaf)/etaf)*Tamb
        #
        # ここでtdmaker.pyのTestDataMaker()クラスで単純な以下のようなテストデータ与えられる。これらの値は0による割り算が発生しないように適当に選ばれている。
        # linPh  = 1.0
        # linyfc = 0.25
        # Qr     = 1.1
        # p0     = 1.0
        # etaf   = 0.5
        # T0     = 1.0
        # Troom  = 273.15 + 15 ( = T_lower_cabin)
        # Tamb   = 273.15 + 15 ( = weather temperature)
        #
        # 上記の値で計算すると以下のようになる。
        # f_shift    =  15/88
        # Tlos_model = (15/88 + sqrt(289.15))**2/0.5 - 2 - 288.15
        #
        # この結果をfloat32表現するとdataの値に一致する。
        #
        expected = np.array([((15/88 + np.sqrt(289.15))**2)/0.5 - 2 - 288.15]).astype(np.float32) # 注意 float32
        self.assertEqual(round(expected[0], 4), round(dems.data[0][0], 4), 'MS::dataの計算値が正しいことを確認')
        self.assertEqual(n_time, len(dems.data),    'MS::dataの打刻数の確認')
        self.assertEqual(n_kid,  len(dems.data[0]), 'MS::dataのチャネル数の確認')

        # MS::mask (既定値)
        self.assertFalse(dems.mask[0 ][0 ], 'maskに既定値が格納されていることを確認[ 0, 0]')
        self.assertFalse(dems.mask[0 ][-1], 'maskに既定値が格納されていることを確認[ 0,-1]')
        self.assertFalse(dems.mask[-1][-1], 'maskに既定値が格納されていることを確認[-1,-1]')
        self.assertFalse(dems.mask[-1][0 ], 'maskに既定値が格納されていることを確認[-1, 0]')
        self.assertTrue(np.array(dems.mask == False).all(), '全maskに既定値が格納されていることを確認')

        # MS::weight (既定値)
        self.assertEqual(1.0, dems.weight[0 ][0 ], 'weightに既定値が格納されていることを確認[ 0, 0]')
        self.assertEqual(1.0, dems.weight[0 ][-1], 'weightに既定値が格納されていることを確認[ 0,-1]')
        self.assertEqual(1.0, dems.weight[-1][-1], 'weightに既定値が格納されていることを確認[-1,-1]')
        self.assertEqual(1.0, dems.weight[-1][0 ], 'weightに既定値が格納されていることを確認[-1, 0]')
        self.assertTrue(np.array(dems.weight == 1.0).all(), '全weightに既定値が格納されていることを確認')

        # MS::long_name (既定値)
        self.assertEqual('Brightness', dems.attrs['long_name'], 'long_nameに既定値が確認されていることを確認')

        # MS::units (既定値)
        self.assertEqual('K', dems.attrs['units'], 'unitsに既定値が確認されていることを確認')

        # MS::time
        self.assertEqual(n_time, len(dems.time),            'MS::timeの打点数の確認')
        self.assertTrue(np.array(times == dems.time).all(), '時刻がREADOUTの打刻と一致することを確認')
        self.assertEqual(times[0],  dems.time.values[0],    'demsの開始時刻が元のデータに一致することを確認')
        self.assertEqual(times[-1], dems.time.values[-1],   'demsの終了時刻が元のデータに一致することを確認')

        # MS::chan
        self.assertEqual(n_kid, len(dems.chan),              'チャネル数の確認')
        self.assertTrue(np.array(kidids == dems.chan).all(), 'KIDIDがチャネル番号と一致することを確認')

        # MS::beam (既定値)
        self.assertEqual(n_time, len(dems.beam),         'MS::beamの打点数が打刻数に一致することを確認')
        self.assertTrue(np.array(dems.beam != '').all(), 'MS::beamに既定値が格納されていることを確認')
        self.assertEqual(0.5, round(np.count_nonzero(dems.beam == 'A')/n_time, 1), 'MS::beamのおよそ半数がAであることを確認')
        self.assertEqual(0.5, round(np.count_nonzero(dems.beam == 'B')/n_time, 1), 'MS::beamのおよそ半数がBであることを確認')
        self.assertEqual('B', dems.beam[0],  'MS::beamの最初の要素がGRADであることを確認')
        self.assertEqual('A', dems.beam[-1], 'MS::beamの最後の要素がONであることを確認')

        # MS::scan (既定値)
        self.assertEqual(n_time, len(dems.scan),         'MS::scanの打点数が打刻数に一致することを確認')
        self.assertTrue(np.array(dems.scan == '').all(), 'MS::scanに既定値が格納されていることを確認')

        # MS::state
        #
        # stateはTestDataMaker()によって前半にGRAD、後半にONが格納されるようになっている。
        #
        self.assertEqual(n_time, len(dems.state), 'MS::stateの打点数が打刻数に一致することを確認')
        self.assertEqual(0.5, round(np.count_nonzero(dems.state == 'GRAD')/n_time, 1), 'MS::stateのおよそ半数がGRADであることを確認')
        self.assertEqual(0.5, round(np.count_nonzero(dems.state == 'ON')/n_time, 1),   'MS::stateのおよそ半数がONであることを確認')
        self.assertEqual('GRAD', dems.state[0],  'MS::stateの最初の要素がGRADであることを確認')
        self.assertEqual('ON',   dems.state[-1], 'MS::stateの最後の要素がONであることを確認')

        # この操作は非常に時間がかかる
        # n_state_change = 0
        # prev = dems.state[0]
        # for state in dems.state:
        #     if (prev != state):
        #         n_state_change += 1
        #     prev = state
        # self.assertEqual(1, n_state_change, 'MS::stateの切り替えが1回だけ発生していることを確認')

        # Telescope Pointing
        self.assertTrue((dems.lon        == -1.1).all(), 'MS::lonの全ての値が正しいことを確認')
        self.assertTrue((dems.lat        == -1.1).all(), 'MS::latの全ての値が正しいことを確認')
        self.assertTrue((dems.lon_origin ==  1.1).all(), 'MS::lon_originの全ての値が正しいことを確認')
        self.assertTrue((dems.lat_origin ==  1.1).all(), 'MS::lat_originの全ての値が正しいことを確認')
        self.assertEqual('altaz', dems.frame, 'MS::frameが既定値であることを確認')

        # Weather Informations
        self.assertTrue((dems.temperature    == 15.0).all(), 'MS::temperatureが既定値でないことを確認')
        self.assertTrue((dems.pressure       == 15.0).all(), 'MS::pressureが既定値でないことを確認')
        self.assertTrue((dems.humidity       == 15.0).all(), 'MS::humidityが既定値でないことを確認')
        self.assertTrue((dems.wind_speed     == 15.0).all(), 'MS::wind_speedが既定値でないことを確認')
        self.assertTrue((dems.wind_direction == 15.0).all(), 'MS::wind_directionが既定値でないことを確認')

        # Data Informations
        #
        # 注意
        # beam_major, beam_minor, beam_pa, exposure, intervalは固定値。
        # MergeToDfits()クラスでも固定値が指定されていた。
        #
        self.assertTrue(np.array(dems.bandwidth  == 0.0  ).all(), 'MS::bandwidthが既定値であることを確認')
        self.assertTrue(np.array(dems.frequency  == 1.5e9).all(), 'MS::frequencyが既定値であることを確認')
        self.assertTrue(np.array(dems.beam_major == 0.005).all(), 'MS::beam_majorが既定値で無いことを確認')
        self.assertTrue(np.array(dems.beam_minor == 0.005).all(), 'MS::beam_minorが既定値で無いことを確認')
        self.assertTrue(np.array(dems.beam_pa    == 0.005).all(), 'MS::beam_paが既定値で無いことを確認')
        self.assertTrue(np.array(dems.exposure   == 1/196).all(), 'MS::exposureが既定値で無いことを確認')
        self.assertTrue(np.array(dems.interval   == 1/196).all(), 'MS::intervalが既定値で無いことを確認')

        # Observation Informations
        self.assertEqual(dems.observation,    'saturn_zscan_v05_2_rtn_pattern6', 'MS::observer')
        self.assertEqual(dems.observer,       'clumsy', 'MS::observer')
        self.assertEqual(dems.object,         'SATURN', 'MS::object')
        self.assertEqual(dems.telescope_name, 'ASTE',   'MS::telescope_name')
        self.assertEqual(dems.telescope_diameter, 10.0, 'MS::telescope_diameterが既定値であることを確認')

        expected = (+2230817.2140945992, -5440188.022176585, -2475718.801708271)
        self.assertEqual(dems.telescope_coordinates, expected, 'MS::telescope_coordinatesが既定値であることを確認')

        # ASTE Specific
        self.assertTrue((dems.aste_cabin_temperature == 15.0 + 273.15).all(), 'MS::aste_cabin_temperatureが既定値でないことを確認')
        self.assertTrue((dems.aste_subref_x.values   == 1.1          ).all(), 'MS::aste_subref_xが既定値で無いことを確認')
        self.assertTrue((dems.aste_subref_y.values   == 1.1          ).all(), 'MS::aste_subref_yが既定値で無いことを確認')
        self.assertTrue((dems.aste_subref_z.values   == 1.1          ).all(), 'MS::aste_subref_zが既定値で無いことを確認')
        self.assertTrue((dems.aste_subref_xt.values  == 1.1          ).all(), 'MS::aste_subref_xtが既定値で無いことを確認')
        self.assertTrue((dems.aste_subref_yt.values  == 1.1          ).all(), 'MS::aste_subref_ytが既定値で無いことを確認')
        self.assertTrue((dems.aste_subref_zt.values  == 1.1          ).all(), 'MS::aste_subref_ztが既定値で無いことを確認')
        self.assertTrue((dems.aste_misti_lon.values  == 180          ).all(), 'MS::aste_misti_lonが既定値で無いことを確認')
        self.assertTrue((dems.aste_misti_lat.values  == 90           ).all(), 'MS::aste_misti_latが既定値で無いことを確認')
        self.assertTrue((dems.aste_misti_pwv.values  == 0.61         ).all(), 'MS::aste_misti_pwvが既定値で無いことを確認')
        self.assertEqual('altaz', dems.aste_misti_frame,                      'MS::aste_misti_frameが既定値であることを確認')

        # DESHIMA 2.0 specific
        self.assertTrue(np.array(dems.d2_mkid_id                != 0).any())
        self.assertTrue(np.array(dems.d2_mkid_type              != '').all())
        self.assertTrue(np.array(dems.d2_mkid_frequency         == 1.5e9).all(), 'd2_mkid_frequencyの値を確認(DDB.KIDDES.F_filter)')
        self.assertTrue(np.array(dems.d2_roomchopper_isblocking == False).all(), 'd2_roomchopper_isblockingの値が既定値であることを確認')
        self.assertTrue(np.array(dems.d2_skychopper_isblocking  == False).any(), 'd2_skychopper_isblockingの値が既定値で無いことを確認')
        self.assertTrue(np.array(dems.d2_skychopper_isblocking  == True).any(),  'd2_skychopper_isblockingの値が既定値で無いことを確認')
        self.assertEqual(0.5, round(np.count_nonzero(dems.d2_skychopper_isblocking == False)/n_time, 1), 'MS::d2_skychopper_isblockingのおよそ半数がFalseであることを確認')
        self.assertEqual(0.5, round(np.count_nonzero(dems.d2_skychopper_isblocking == True)/n_time, 1),  'MS::d2_skychopper_isblockingのおよそ半数がTrueであることを確認')
        self.assertEqual(DEMS_VERSION,    dems.d2_dems_version)
        self.assertEqual(DEMERGE_VERSION, dems.d2_demerge_version)
        return

    def test_no_cabin_file(self):
        prefix = 'testdata'
        dems = mtd.merge_to_dems(
            ddbfits_path='{}_DDB.fits.gz'.format(prefix),
            obsinst_path='../data/cosmos/cosmos_{0}/{0}.obs'.format(self.obsid),
            antenna_path='{}.ant'.format(prefix),
            readout_path='{}_reduced_readout.fits'.format(prefix),
            skychop_path='{}.skychop'.format(prefix),
            weather_path='{}.wea'.format(prefix),
            misti_path='{}.misti'.format(prefix),
            cabin_path='',
        )
        self.assertTrue((np.isnan(dems.aste_cabin_temperature)).all(), 'MS::aste_cabin_temperatureがNaNであることを確認')
        return

    def test_no_misti_file(self):
        prefix = 'testdata'
        dems = mtd.merge_to_dems(
            ddbfits_path='{}_DDB.fits.gz'.format(prefix),
            obsinst_path='../data/cosmos/cosmos_{0}/{0}.obs'.format(self.obsid),
            antenna_path='{}.ant'.format(prefix),
            readout_path='{}_reduced_readout.fits'.format(prefix),
            skychop_path='{}.skychop'.format(prefix),
            weather_path='{}.wea'.format(prefix),
            misti_path='',
            cabin_path='',
        )
        self.assertTrue((np.isnan(dems.aste_cabin_temperature)).all(), 'MS::aste_cabin_temperatureがNaNであることを確認')
        self.assertTrue((np.isnan(dems.aste_misti_lon)).all(), 'MS::aste_misti_lonがNaNであることを確認')
        self.assertTrue((np.isnan(dems.aste_misti_lat)).all(), 'MS::aste_misti_latがNaNであることを確認')
        self.assertTrue((np.isnan(dems.aste_misti_pwv)).all(), 'MS::aste_misti_pwvがNaNであることを確認')
        return

    def test_retrieve_cabin_temps(self):
        """cabinの温度をロードする"""
        datetimes, upper, lower = mf.retrieve_cabin_temps('../data/cosmos/cosmos_20171103184436/20171103184436.cabin')
        self.assertEqual(4, len(datetimes), '時刻の行数を確認')
        self.assertEqual(4, len(upper),     'upper_cabin_tempの行数を確認')
        self.assertEqual(4, len(lower),     'lower_cabin_tempの行数を確認')
        expected = np.array([datetime(2017, 11, 3, hour=18, minute=44)]).astype('datetime64[ns]')
        self.assertEqual(expected[0], datetimes[0])
        self.assertEqual(13.2,        upper[0])
        self.assertEqual(16.6,        lower[0])
        expected = np.array([datetime(2017, 11, 3, hour=18, minute=47)]).astype('datetime64[ns]')
        self.assertEqual(expected[0], datetimes[3])
        self.assertEqual(13.2,        upper[3])
        self.assertEqual(16.6,        lower[3])

        datetimes, upper, lower = mf.retrieve_cabin_temps()
        self.assertTrue(np.isnan(datetimes[0]), '時刻がNaT(filename=None)')
        self.assertTrue(np.isnan(upper[0]),     '温度がNaT(upper)(filename=None)')
        self.assertTrue(np.isnan(lower[0]),     '温度がNaT(lower)(filename=None)')

        datetimes, upper, lower = mf.retrieve_cabin_temps('')
        self.assertTrue(np.isnan(datetimes[0]), '時刻がNaT')
        self.assertTrue(np.isnan(upper[0]),     '温度がNaT(upper)')
        self.assertTrue(np.isnan(lower[0]),     '温度がNaT(lower)')
        return

    def test_retrieve_skychop_states(self):
        datetimes, states = mf.retrieve_skychop_states('../data/cosmos/cosmos_20171103184436/20171103184436.skychop')
        self.assertEqual(192201,            len(datetimes))
        self.assertEqual(192201,            len(states))
        self.assertEqual(1509734678.900000, datetimes[0])
        self.assertEqual(1,                 states[0])
        self.assertEqual(1,                 states[1])
        self.assertEqual(1,                 states[2])

        datetimes, states = mf.retrieve_skychop_states('testdata.skychop.dat.xz')
        self.assertEqual(240000,       len(datetimes))
        self.assertEqual(240000,       len(states))
        self.assertEqual(1695373200.0, datetimes[0])
        self.assertEqual(1,            states[0])
        return

    def test_retrieve_misti_log(self):
        time, az, el, pwv = mf.retrieve_misti_log('../data/cosmos/cosmos_20171103184436/20171103184436.misti')
        self.assertEqual(1848, len(time), 'mistiファイルの行数を確認')
        expected = np.array([datetime(2017, 11, 3, hour=18, minute=44, second=38, microsecond=900000)]).astype('datetime64[ns]')
        self.assertEqual(expected[0], time[0], '時刻が正しいことを確認')
        self.assertEqual(180.0,       az[0],   'azが正しいことを確認')
        self.assertEqual(90.0,        el[0],   'elが正しいことを確認')
        self.assertEqual(0.610,       pwv[0],  'pwvが正しいことを確認')
        return

if __name__=='__main__':
    unittest.main()

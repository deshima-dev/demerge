"""dfitsをdemsへ変換するプログラムのテスト

File name: test_dfits2dems.py
Python 3.9
(C) 2023 内藤システムズ
"""
import unittest
import dfits2dems as dd
import numpy      as np
import pandas     as pd
import math
import os

from astropy.io     import fits
from dems.d2        import MS
from merge_to_dfits import MergeToDfits
from datetime       import datetime

import merge_function as fc

class Dfits2demsTestDrive(unittest.TestCase):
    """dfits2dems.pyモジュールの単体テスト"""
    def setUp(self):
        # testdata = 'dfits_testdata.fits.gz'
        # if not os.path.exists(testdata):
        #     obsid = '20171110114116'
        #     path  = 'data/deshima2.0/cosmos_{}'.format(obsid)
        #     mtd = MergeToDfits(ddbfits    = 'DDB_20180619.fits.gz',
        #                        dfitsdict  = 'dfits_dict.yaml',
        #                        obsinst    = '{}/{}.obs'.format(path, obsid),
        #                        antennalog = '{}/{}.ant'.format(path, obsid),
        #                        weatherlog = '{}/{}.wea'.format(path, obsid),
        #                        cabinlog   = '{}/{}.cabin'.format(path, obsid),
        #                        skychoplog = 'skychop_testdata2.skychop',
        #                        rout_data  = 'cache/{}/reduced_{}.fits'.format(obsid, obsid))
            
        #     dfits_hdus = mtd.dfits
        #     dfits_hdus.writeto(testdata, overwrite=True)
        #     mtd.kidsinfo_hdus.close()

        #self.filename = '../dmerge/cache/20171110114116/dfits_20171110114116.fits.gz'
        self.filename = 'dfits_dummy.fits.gz'
        with fits.open(self.filename) as hdul:
            self.readout = hdul['READOUT'].data
            self.obsinfo = hdul['OBSINFO'].data
            self.antenna = hdul['ANTENNA'].data
            self.weather = hdul['WEATHER'].data
            self.cabin   = hdul['CABIN_T'].data
            self.skychop = hdul['SKYCHOP'].data
            self.misti   = hdul['MISTI'].data
        
        return

    def test_dfits2dems(self):
        """dfits2dems関数のテスト"""
        ms       = dd.convert_dfits_to_dems(self.filename, still_period=2, shuttle_min_lon_on=-0.0001, shuttle_max_lon_on=0.1)
        expected = self.readout['Tsignal'].astype(np.float64)

        # (NaN==NaN)はall()による判定がFalseになるので一時的に-1へ変換する
        expected = np.where(np.isnan(expected), -1, expected)
        result   = np.where(np.isnan(ms.data),  -1, ms.data)

        self.assertTrue((result == expected).all(), 'MS::data')

        expected = np.array(self.readout['starttime']).astype(np.datetime64)
        result   = ms.time

        self.assertTrue((result == expected).all(), 'MS::time')

        expected = self.obsinfo['kidids'][0].astype(np.int64)
        result   = ms.chan

        self.assertTrue((result == expected).all(), 'MS::chan')

        result = np.array(np.where(ms.scan != ''))
        self.assertTrue(result.any(), 'MS::scanに規定値以外がセットされたことを確認する')
        
        result = np.array(np.where(ms.temperature > 0.0))
        self.assertTrue(result.any(), 'MS::temperatureに規定値以外がセットされたことを確認する')

        result = np.array(np.where(ms.pressure > 0.0))
        self.assertTrue(result.any(), 'MS::pressureに規定値以外がセットされたことを確認する')

        result = np.array(np.where(ms.humidity > 0.0))
        self.assertTrue(result.any(), 'MS::humidityに規定値以外がセットされたことを確認する')

        result = np.array(np.where(ms.wind_speed > 0.0))
        self.assertTrue(result.any(), 'MS::wind_speedに規定値以外がセットされたことを確認する')

        result = np.array(np.where(ms.wind_direction > 0.0))
        self.assertTrue(result.any(), 'MS::wind_directionに規定値以外がセットされたことを確認する')

        result = np.array(np.where(ms.aste_cabin_temperature > 0.0))
        self.assertTrue(result.any(), 'MS::aste_cabin_temperatureに規定値以外がセットされたことを確認する')

        result = np.array(np.where(ms.d2_mkid_id > 0))
        self.assertTrue(result.any(), 'MS::d2_mkid_idに規定値以外がセットされたことを確認する')

        result = np.array(np.where(ms.d2_mkid_type != ''))
        self.assertTrue(result.any(), 'MS::d2_mkid_typeに規定値以外がセットされたことを確認する')

        result = np.array(np.where(ms.d2_mkid_frequency > 0.0))
        self.assertTrue(result.any(), 'MS::d2_mkid_frequencyに規定値以外がセットされたことを確認する')

        result = np.array(np.where(ms.beam_major > 0.0))
        self.assertTrue(result.any(), 'MS::beam_majorに規定値以外がセットされたことを確認する')

        result = np.array(np.where(ms.beam_minor > 0.0))
        self.assertTrue(result.any(), 'MS::beam_minorに規定値以外がセットされたことを確認する')

        result = np.array(np.where(ms.beam_pa > 0.0))
        self.assertTrue(result.any(), 'MS::beam_paに規定値以外がセットされたことを確認する')

        result = np.array(np.where(ms.lon != 0.0))
        self.assertTrue(result.any(), 'MS::lonに規定値以外がセットされたことを確認する')

        result = np.array(np.where(ms.lat != 0.0))
        self.assertTrue(result.any(), 'MS::latに規定値以外がセットされたことを確認する')

        self.assertTrue(ms.lon_origin > 0.0, 'MS::lon_originに規定値以外がセットされたことを確認する')
        self.assertTrue(ms.lat_origin > 0.0, 'MS::lat_originに規定値以外がセットされたことを確認する')
        
        self.assertEqual(ms.observer,       'dummy_observer',  'MS::observer')
        self.assertEqual(ms.object,         'dummy_object',    'MS::object')
        self.assertEqual(ms.telescope_name, 'dummy_telescope', 'MS::telescope_name')

        self.assertEqual(ms.exposure, self.obsinfo['integtime'][0], 'MS::exposure')
        self.assertEqual(ms.interval, self.obsinfo['interval'][0],  'MS::interval')

        # skychopの状態が変化する時点を検索
        indexes = np.where(self.skychop['state'] == 1)
        # そこインデックスに対する時刻を調べる
        change_time1 = datetime.fromtimestamp(self.skychop['time'][indexes[0][0]])
        change_time1 = np.datetime64(change_time1)

        # 同じことをms.d2_skychopper_isblockingにも行う
        indexes = np.where(ms.d2_skychopper_isblocking == 1)
        change_time2 = ms.time[indexes[0][0]].values

        e = 100 # ns
        self.assertTrue(e > change_time2 - change_time1, 'skychopの状態が変化した時刻がdfitsとdemsでだいたい同じであることを確認する')
        
        return

    def test_retrieve_cabin_temps(self):
        """cabinの温度をロードする"""
        datetimes, upper, lower = dd.retrieve_cabin_temps('data/deshima2.0/cosmos_20171110114116/20171110114116.cabin')
        self.assertEqual(8, len(datetimes))
        self.assertEqual(8, len(upper))
        self.assertEqual(8, len(lower))
        return

    def test_retrieve_skychop_states(self):
        datetimes, states = dd.retrieve_skychop_states('skychop_testdata3.skychop')
        self.assertEqual(991, len(datetimes))
        self.assertEqual(991, len(states))
        self.assertEqual(2, states[0])
        self.assertEqual(1, states[1])
        self.assertEqual(0, states[2])
        self.assertEqual(1, states[3])
        #print(np.array([datetime.fromtimestamp(t) for t in datetimes]).astype(np.datetime64))
        return

    def test_MergeToDfits(self):
        """MergeToDfitsクラスのテスト
        このテストを行うにはreduced_XXX.fitsがあらかじめ作成されている必要がある。
        """
        obsid = '20171110114116'
        path  = 'data/deshima2.0/cosmos_{}'.format(obsid)
        mtd = MergeToDfits(ddbfits    = 'DDB_20180619.fits.gz',
                           dfitsdict  = 'dfits_dict.yaml',
                           obsinst    = '{}/{}.obs'.format(path, obsid),
                           antennalog = '{}/{}.ant'.format(path, obsid),
                           weatherlog = '{}/{}.wea'.format(path, obsid),
                           cabinlog   = '{}/{}.cabin'.format(path, obsid),
                           #skychoplog = '{}/{}.skychop'.format(path, obsid),
                           skychoplog = 'skychop_testdata.skychop',
                           rout_data  = 'cache/{}/reduced_{}.fits'.format(obsid, obsid))
        
        dfits_hdus = mtd.dfits
        mtd.kidsinfo_hdus.close()

        # 20171110114116.cabinには8行のデータがある
        self.assertEqual(len(dfits_hdus['cabin_t'].data['time']),        8)
        self.assertEqual(len(dfits_hdus['cabin_t'].data['upper_cabin']), 8)
        self.assertEqual(len(dfits_hdus['cabin_t'].data['main_cabin']),  8)

        self.assertTrue('skychop' in dfits_hdus)
        #self.assertEqual(len(dfits_hdus['skychop'].data['time']), 428201)
        #self.assertEqual(len(dfits_hdus['skychop'].data['state']), 428201)
        self.assertEqual(len(dfits_hdus['skychop'].data['time']), 146)
        self.assertEqual(len(dfits_hdus['skychop'].data['state']), 146)
        self.assertEqual(2, dfits_hdus['skychop'].data['state'][0])
        self.assertEqual(1, dfits_hdus['skychop'].data['state'][1])
        self.assertEqual(0, dfits_hdus['skychop'].data['state'][2])
        self.assertEqual(1, dfits_hdus['skychop'].data['state'][3])
        return

    def test_dfits2dems_dummy_dfits(self):
        """dummy_dfitsを使ったdfits2dems関数のテスト"""
        ms = dd.convert_dfits_to_dems('dfits_dummy.fits.gz', still_period=2, shuttle_min_lon_on=-0.0001, shuttle_max_lon_on=0.1)
        #print(ms.d2_skychopper_isblocking)
        
if __name__=='__main__':
    unittest.main()

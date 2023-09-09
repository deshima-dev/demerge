"""merge_to_dems()関数をテストする

Python 3.9

(C) 2023 内藤システムズ
"""
import unittest
import dfits2dems as dd
import numpy      as np
import xarray     as xa
import pandas     as pd
import math
import os

import merge_function as mf
import merge_to_dems  as mtd

from astropy.io     import fits
from dems.d2        import MS
from merge_to_dfits import MergeToDfits
from datetime       import datetime

class MergeToDemsTestDrive(unittest.TestCase):
    """merge_to_dems.pyモジュールの単体テスト"""
    def setUp(self):
        self.obsid = '20171103184436'
        self.reduced_fits = 'cache/{0}/reduced_{0}.fits'.format(self.obsid)
        self.dfits = 'cache/{0}/dfits_{0}.fits.gz'.format(self.obsid)
        return

    def test_merge_to_dems(self):
        dems = mtd.merge_to_dems(
            ddbfits_path='DDB_20180619.fits.gz',
            obsinst_path='data/deshima2.0/cosmos_{0}/{0}.obs'.format(self.obsid),
            antenna_path='data/deshima2.0/cosmos_{0}/{0}.ant'.format(self.obsid),
            readout_path=self.reduced_fits,
            skychop_path='',
            weather_path='data/deshima2.0/cosmos_{0}/{0}.wea'.format(self.obsid),
            misti_path='data/deshima2.0/cosmos_{0}/{0}.misti'.format(self.obsid),
            cabin_path='data/deshima2.0/cosmos_{0}/{0}.cabin'.format(self.obsid),
        )
        self.assertEqual(xa.DataArray, type(dems))

        with fits.open(self.reduced_fits) as hdul:
            n_time = hdul['READOUT'].header['NAXIS2']
            n_kid  = hdul['KIDSINFO'].header['NAXIS2']
            self.assertEqual(n_time, len(dems.time))

            time = np.array(mf.convert_timestamp(hdul['READOUT'].data['timestamp'])).astype(np.datetime64)
            self.assertEqual(time[0],  dems.time.values[0])
            self.assertEqual(time[-1], dems.time.values[-1])

            self.assertEqual(n_time, len(dems.data), '時間軸のデータ数の確認')
            self.assertEqual(n_kid,  len(dems.data[0]), 'KID軸のデータ数の確認')

        with fits.open(self.dfits) as hdul:
            expected = hdul['READOUT'].data['Tsignal'].astype(np.float64)

            # (NaN==NaN)はall()による判定がFalseになるので一時的に-1へ変換する
            expected = np.where(np.isnan(expected), -1, expected)
            result   = np.where(np.isnan(dems.data),  -1, dems.data)
            self.assertTrue((result == expected).all(), 'MS::data')

        self.assertTrue(np.array(dems.aste_cabin_temperature.values > 0).all())

        self.assertTrue(np.array(dems.lon != 0).all())
        self.assertTrue(np.array(dems.lat != 0).all())
        self.assertTrue(np.array(dems.lon_origin > 0).all())
        self.assertTrue(np.array(dems.lat_origin > 0).all())

        self.assertTrue(np.array(dems.temperature != 0).all())
        self.assertTrue(np.array(dems.pressure != 0).all())
        self.assertTrue(np.array(dems.humidity != 0).all())
        self.assertTrue(np.array(dems.wind_speed != 0).all())
        self.assertTrue(np.array(dems.wind_direction != 0).all())
        self.assertTrue(np.array(dems.d2_mkid_id != 0).any())
        self.assertTrue(np.array(dems.d2_mkid_type != '').all())

        result = np.where(dems.d2_mkid_frequency > 0, dems.d2_mkid_frequency, 1) # NaNを一時的に1に置き換える
        self.assertTrue(np.array(result).all())

        self.assertTrue(np.array(dems.beam_major > 0).all())
        self.assertTrue(np.array(dems.beam_minor > 0).all())
        self.assertTrue(np.array(dems.beam_pa > 0).all())

        self.assertEqual(dems.observation,  'saturn_zscan_v05_2_rtn_pattern6', 'MS::observer')
        self.assertEqual(dems.observer,     'clumsy', 'MS::observer')
        self.assertEqual(dems.object,       'SATURN', 'MS::object')
        self.assertEqual(dems.telescope_name, 'ASTE',   'MS::telescope_name')
        
        self.assertEqual(1./196, float(dems.exposure), 'MS::exposure')
        self.assertEqual(1./196, float(dems.interval), 'MS::interval')

        self.assertEqual('altaz', dems.aste_misti_frame.values)
        self.assertTrue((dems.aste_misti_lon.values != 0).all())
        self.assertTrue((dems.aste_misti_lat.values != 0).all())

        self.assertTrue((dems.scan != '').all())

        return

if __name__=='__main__':
    unittest.main()

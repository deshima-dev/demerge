"""dfitsをdemsへ変換するプログラムのテスト

File name: test_dfits2dems.py
Python 3.9
(C) 2023 内藤システムズ
"""
import unittest
import dfits2dems as dd
import numpy      as np
import math

from astropy.io import fits
from dems.d2    import MS

class Dfits2demsTestDrive(unittest.TestCase):
    """dfits2dems.pyモジュールの単体テスト"""
    def setUp(self):
        filename = '../dmerge/cache/20171110114116/dfits_20171110114116.fits.gz'
        with fits.open(filename) as hdul:
            self.readout = hdul['READOUT'].data
            self.obsinfo = hdul['OBSINFO'].data
            self.antenna = hdul['ANTENNA'].data
            self.weather = hdul['WEATHER'].data
            self.cabin   = hdul['CABIN_T'].data
        
        return

    def test_dfits2dems(self):
        """dfits2dems関数のテスト"""
        ms       = dd.dfits2dems()
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
        
        self.assertEqual(ms.observer,       'clumsy', 'MS::observer')
        self.assertEqual(ms.object,         'MARS',   'MS::object')
        self.assertEqual(ms.telescope_name, 'ASTE',   'MS::telescope_name')

        self.assertEqual(ms.exposure, self.obsinfo['integtime'][0], 'MS::exposure')
        self.assertEqual(ms.interval, self.obsinfo['interval'][0],  'MS::interval')

        return 

if __name__=='__main__':
    unittest.main()

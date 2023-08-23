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

        expected = len(self.readout['starttime'])

        self.assertEqual(len(ms.scan), expected, 'MS::scanの要素数')
        self.assertEqual(len(ms.temperature), expected, 'MS::temperatureの要素数')
        self.assertEqual(len(ms.pressure), expected, 'MS::pressureの要素数')
        self.assertEqual(len(ms.humidity), expected, 'MS::humidityの要素数')
        self.assertEqual(len(ms.wind_speed), expected, 'MS::wind_speedの要素数')
        self.assertEqual(len(ms.wind_direction), expected, 'MS::wind_directionの要素数')

        return 

if __name__=='__main__':
    unittest.main()

"""merge_functionモジュールをテストする

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
from datetime   import datetime

# constants
from dems    import __version__ as DEMS_VERSION
from demerge import __version__ as DEMERGE_VERSION

class MergeFunctionTestDrive(unittest.TestCase):
    """merge_function.pyモジュールの単体テスト"""
    def setUp(self):
        
        return

    def test_fshift(self):
        readout_hdul = fits.open('testdata_reduced_readout.fits')
        pixelid      = 0

        result   = mf.fshift(readout_hdul, pixelid)
        expected = np.array([(1.0 - 0.25)/(4*1.1)]).astype('float32')
        self.assertEqual(result[0][0],   expected[0],                                    'fshiftの計算値')
        self.assertEqual(len(result),    63,                                             'fshiftのKID数')
        self.assertEqual(len(result[0]), len(readout_hdul['READOUT'].data['timestamp']), 'fshiftのデータ点数')

        readout_hdul.close()
        return

    def test_calibrate_to_power(self):
        
        return

if __name__=='__main__':
    unittest.main()

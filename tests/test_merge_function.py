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
        self.assertEqual(len(result),    readout_hdul['READOUT'].header['NKID0'],        'fshiftのKID数')
        self.assertEqual(len(result[0]), len(readout_hdul['READOUT'].data['timestamp']), 'fshiftのデータ数')

        readout_hdul.close()
        return

    def test_calibrate_to_power(self):
        readout_hdul = fits.open('testdata_reduced_readout.fits')
        ddbfits_hdul = fits.open('testdata_DDB.fits.gz')
        pixelid      = 0

        linPh   = 1.0
        linyfc  = 0.25
        Qr      = 1.1
        fshift  = np.array([(linPh - linyfc)/(4*Qr)]).astype('float32')[0]

        T_amb   =  0  + 273.15
        T_cabin = 15  + 273.15
        p0      = 1.0
        etaf    = 0.5
        T0      = 1.0
        Tlos_model = (fshift + p0*np.sqrt(T_cabin + T0))**2/(p0*p0*etaf) - (T0/etaf) - ((1 - etaf)/etaf)*T_amb
        Tlos_model = np.array([Tlos_model]).astype('float32')[0]
        
        fshift = mf.fshift(readout_hdul, pixelid)
        result = mf.calibrate_to_power(T_cabin, T_amb, fshift, ddbfits_hdul)
        self.assertEqual(result[0][0],   Tlos_model,                                     'calibrate_to_powerの計算値')
        self.assertEqual(len(result),    len(readout_hdul['READOUT'].data['timestamp']), 'データ数')
        self.assertEqual(len(result[0]), readout_hdul['READOUT'].header['NKID0'],        'KID数')

        ddbfits_hdul.close()
        readout_hdul.close()
        return

if __name__=='__main__':
    unittest.main()

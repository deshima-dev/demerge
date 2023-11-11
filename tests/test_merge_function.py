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
        
        return

if __name__=='__main__':
    unittest.main()

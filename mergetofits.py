"""フィットの結果をfitsファイルへ統合する

Filename: mergetofits.py
(C) 2021 内藤システムズ

使い方
------
例)
  $ CACHE_DIR="./cache"
  $ DATA_DIR="../raw_dataset/obs"
  $ OBSID="20171103184836"
  $ python mergetofits.py \
       "${CACHE_DIR}/${OBSID}/reduced_${OBSID}.fits" \
       "${CACHE_DIR}/${OBSID}/dfits_${OBSID}.fits.gz" \
       "${DATA_DIR}/cosmos_${OBSID}/${OBSID}.obs" \
       "${DATA_DIR}/cosmos_${OBSID}/${OBSID}.ant" \
       "${DATA_DIR}/cosmos_${OBSID}/${OBSID}.wea" \
       ../caldb/DDB_20180619.fits.gz \
       dfits_dict.yaml \
       cabin.db
"""
import pickle
import os
import sys
import shutil
import unittest
import pickle
import numpy as np
import lmfit
import glob

from astropy.io import fits
from scipy import interpolate
import dmerge
from merge_to_dfits import MergeToDfits

if __name__=='__main__':
    """
    コマンドライン引数
    ------------------
    args[1] string reduced fitsファイルへの相対パス
    args[2] string 作成するdfitsファイルへの相対パス
    args[3] string obsファイルへの相対パス
    args[4] string antファイルへの相対パス
    args[5] string weaファイルへの相対パス
    args[6] string caldb fitsファイルへの相対パス
    args[7] string yamlファイルへの相対パス
    args[8] string cabin.dbファイルへの相対パス
    """
    args = sys.argv
    reducedfitsfile = args[1]
    dfitsfile = args[2]
    obsinst = args[3]
    antennalog = args[4]
    weatherlog = args[5]
    ddb_fits = args[6]
    dfitsdict = args[7]
    cabinlog = args[8]

    if os.path.exists(dfitsfile):
        os.remove(dfitsfile)

    mtd = MergeToDfits(ddbfits=ddb_fits,
                       dfitsdict=dfitsdict,
                       obsinst=obsinst,
                       antennalog=antennalog,
                       rout_data=reducedfitsfile,
                       weatherlog=weatherlog,
                       cabinlog=cabinlog)
    dfits_hdus = mtd.dfits
    dfits_hdus.writeto(dfitsfile)
    mtd.kidsinfo_hdus.close()

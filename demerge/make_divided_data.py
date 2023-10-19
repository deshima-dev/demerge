"""SweepデータやTODをKID毎に分割されたファイル(*.pkl)にする

Filename: make_divided_data.py
(C) 2021 内藤システムズ

コマンド引数
------------
kids listファイルの名前
local sweepファイルの名前
観測データfitsファイルの名前
一時ファイルを格納するディレクトリ名

使い方
------
コマンドラインでkids.list、Sweepデータ、TODのファイル(fits)名を指定する。
KID毎に分割されたファイル(*.pkl)はcacheディレクトリ内に作成される。

例)
  $ python make_divided_data.py ../raw_dataset/obs/cosmos_20171103184836/kids.list \
                                ../raw_dataset/obs/cosmos_20171103184836/localsweep.sweep \
                                ../raw_dataset/obs/cosmos_20171103184836/20171103184836.fits \
                                cache/20171103184836
"""
import os
import sys
import shutil
import pickle
from . import demerge

def main() -> None:
    """
    コマンドライン引数
    ------------------
    args[1] string kidsリストファイルへのパス(テキストファイル)
    args[2] string ローカルsweepのデータファイルへのパス(テキストファイル)
    args[3] string TODファイルへのパス(fitsファイル)
    args[4] string 分割されたファイルの保管先ディレクトリ
    """
    args = sys.argv
    kidslist = demerge.load_kidslist(args[1])
    localsweeps = demerge.load_localsweep(args[2], framelen=kidslist[0]['framelen'])
    tods = demerge.TODs(args[3])
    kids = demerge.divide_data(kidslist, localsweeps, tods)
    tods.close()
    if os.path.exists(args[4]):
        shutil.rmtree(args[4])
    os.mkdir(args[4])
    for kid in kids:
        with open(args[4] + '/kid{:05}.pkl'.format(int(kid['kidid'])), 'wb') as f:
            pickle.dump(kid, f)

if __name__=='__main__':
    main()

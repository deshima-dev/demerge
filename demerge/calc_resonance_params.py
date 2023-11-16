"""1つのKIDについて共振パラメータを計算する

Filename: calc_resonance_params.py
(C) 2021 内藤システムズ

使い方
------
コマンドラインで第一引数に一時ファイル(*.pkl)の名前を指定する。

例)
  $ python calc_resonance_params.py cache/kid00001.pkl
"""
import sys
import pickle
from . import demerge

def calc_resonance_params(kid, nfwhm=5, minq=100, maxratio=1.0):
    """1つのKIDに対するresonanceパラメータを求める
    引数
    ----
    dict  KIDに関する測定値などをまとめた辞書(divide_data()関数によって作成される)
    float FWHMを単位としたフィット範囲[オプション]
    float qの最小値(フィットの良し悪しを判定する基準として利用する:demerge.find_peaks()関数で利用)[オプション]
    float フィットの良し悪しを判定する基準((bg-d)/bgの最大値:demerge.find_peaks()関数で利用)[オプション]

    戻り値
    ------
    dict フィットの結果を引数で与えられたkidに追加した辞書
    """
    fc = kid['tod'].frequency
    peaks = demerge.find_peaks(kid['localsweep'].x, kid['localsweep'].amplitude, fc, minq=minq, maxratio=maxratio)
    if len(peaks) == 0:
        print('calc_resonance_params(): no peaks')
        kid['resonance_params'] = None
        kid['fitrange']         = None
        return kid

    result, fitrange = demerge.fit_onepeak(kid['localsweep'], peaks, nfwhm=nfwhm)

    fr  = result.params['fr'].value # GHz
    dfr = result.params['fr'].stderr # GHz
    Qr  = result.params['Qr'].value
    dQr = result.params['Qr'].stderr
    Qc  = result.params['Qc'].value
    dQc = result.params['Qc'].stderr
    Qi  = result.params['Qi'].value
    dQi = result.params['Qi'].stderr

    if dfr is None or dQr is None or dQc is None or dQi is None:
        kid['enabled'] = False
    elif dfr/fr<0. or dQr/Qr<0. or dQc/Qc<0. or dQi/Qi<0.:
        kid['enabled'] = False
    elif fr!=fr or dfr!=dfr or Qr!=Qr or dQr!=dQr: # reject Nan
        kid['enabled'] = False
    elif abs(Qr)==float('inf') or abs(dQr)==float('inf'): # reject +/-inf
        kid['enabled'] = False

    if kid['enabled']:
        kid['resonance_params'] = result.params
        kid['fitrange'] = fitrange
    else:
        kid['resonance_params'] = None
        kid['fitrange'] = None

    return kid

def main() -> None:
    """
    コマンドライン引数
    ------------------
    args[1] string divide_data()関数で生成したpklファイルへのパス
    """
    args = sys.argv
    filename = args[1]

    with open(filename, 'rb') as f:
        kid = pickle.load(f)

    kid = calc_resonance_params(kid)

    with open(filename, 'wb') as f:
        pickle.dump(kid, f)

if __name__=='__main__':
    main()

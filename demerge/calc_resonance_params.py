"""1つのKIDについて共振パラメータを計算する

Filename: calc_resonance_params.py
(C) 2021 内藤システムズ

使い方
------
コマンドラインで第一引数に一時ファイル(*.pkl)の名前を指定する。

例)
  $ python calc_resonance_params.py cache/kid00001.pkl
"""
# standard library
import sys
import pickle
from logging import basicConfig, getLogger


# dependencies
from .demerge import find_peaks, fit_onepeak


# module logger
logger = getLogger(__name__)


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
    def on_no_peaks(message: str = None) -> dict:
        if message is not None:
            kidid = kid['kidid']
            logger.warning(f'KID ID {kidid}: {message}')

        kid['fitrange'] = None
        kid['resonance_params'] = None
        return kid

    def on_fit_failed(message: str = None) -> dict:
        if message is not None:
            kidid = kid['kidid']
            logger.warning(f'KID ID {kidid}: {message}')

        kid['enabled'] = False
        kid['fitrange'] = None
        kid['resonance_params'] = None
        return kid

    peaks = find_peaks(
        freq=kid['localsweep'].x,
        ampl=kid['localsweep'].amplitude,
        fc=kid['tod'].frequency,
        minq=minq,
        maxratio=maxratio,
    )

    if not peaks:
        return on_no_peaks('Found no peaks')

    try:
        result, fitrange = fit_onepeak(kid['localsweep'], peaks, nfwhm)
    except Exception as error:
        return on_fit_failed(error)

    fr  = result.params['fr'].value # GHz
    dfr = result.params['fr'].stderr # GHz
    Qr  = result.params['Qr'].value
    dQr = result.params['Qr'].stderr
    Qc  = result.params['Qc'].value
    dQc = result.params['Qc'].stderr
    Qi  = result.params['Qi'].value
    dQi = result.params['Qi'].stderr

    if dfr is None or dQr is None or dQc is None or dQi is None:
        return on_fit_failed('Got invalid fit parameter(s)')

    if dfr/fr<0. or dQr/Qr<0. or dQc/Qc<0. or dQi/Qi<0.:
        return on_fit_failed('Got invalid fit parameter(s)')

    if fr!=fr or dfr!=dfr or Qr!=Qr or dQr!=dQr: # reject Nan
        return on_fit_failed('Got invalid fit parameter(s)')

    if abs(Qr)==float('inf') or abs(dQr)==float('inf'): # reject +/-inf
        return on_fit_failed('Got invalid fit parameter(s)')

    kid['fitrange'] = fitrange
    kid['resonance_params'] = result.params
    return kid


def main() -> None:
    """
    コマンドライン引数
    ------------------
    args[1] string divide_data()関数で生成したpklファイルへのパス
    """
    basicConfig(
        datefmt='%Y-%m-%d %H:%M:%S',
        format='[%(asctime)s %(name)s %(levelname)s] %(message)s',
    )

    filename = sys.argv[1]

    with open(filename, 'rb') as f:
        kid = pickle.load(f)

    kid = calc_resonance_params(kid)

    with open(filename, 'wb') as f:
        pickle.dump(kid, f)


if __name__=='__main__':
    main()

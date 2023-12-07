"""DESHIMA解析プログラム

File name: demerge.py
Python 3.7
(C) 2021 内藤システムズ
"""
# standard library
from collections.abc import Mapping
from logging import getLogger
from warnings import catch_warnings, simplefilter


# dependencies
import lmfit
import scipy
import sympy
import numpy as np
from astropy.io import fits


# module logger
logger = getLogger(__name__)


#-------------------------
# For make_divided_data.py
def divide_data(kidslist, localsweep, tods):
    """KID毎にSweepデータやTODを分割して辞書に格納する
    引数
    ----
    list demerge.load_kidslist()関数で得られるタプル(info, kid_bins, blind_bins, powers)
    list SweepDataオブジェクトのリスト
    TODs TODsオブジェクト

    戻り値
    ------
    リスト
    """
    result = []
    for i, (nbin, l, r) in enumerate(find_nearest_blind_tones(kidslist)):
        result.append({
            'pixelid': 0,
            'kidid': i,
            'bin': nbin,
            'fftgain': tods.fftgain,
            'framert': tods.framert,
            'nbins': tods.nbins,
            'npoints': tods.npoints,
            'lofreq': tods.lofreq,
            'localsweep': localsweep[nbin],
            'readpower': kidslist[3][nbin],
            'blind_tone_left': tods[l], #FixedData class
            'blind_tone_right': tods[r], #FixedData class
            'tod': tods[nbin], #FixedData class
            'tod_timestamp': tods.timestamp,
            'enabled': True,
        })
    return result

def find_nearest_blind_tones(kidslist):
    """ KID信号の左右にある最も近いblind信号を探す
    引数
    ----
    kidslist(demerge.read_kidslist()関数の戻り値)

    戻り値
    ------
    list [[KID信号, 左blind tone, 右blind tone], ...]

    戻り値の例
    ----------
    [[-196483 -201203 -192949]
     [-189634 -192949 -184695]
     [-188948 -192949 -184695]
    ...
     [ 154774  153170  161962]
     [ 158692  153170  161962]
     [ 169653  161962  171230]]

    詳細
    ----
    この関数は引数に与えられたkidslistに含まれるkidsとblindsを利用する。
    kidsとblindsを組み合わせて各kid信号の両側にある最も近いblind信号を探す。
    見つかった左右のblind信号とkid信号をリストにまとめて返す。戻り値はリストの入れ子構造となる。
    blind信号が見つからない場合はNoneが格納される。
    戻り値はkidのbin番号で昇順になっている。
    この関数はbbsweeplib.TODcalibに定義されていたkids_with_both_blinds()関数を単機能にしたもの。

    探索の方法
    ----------
    kid信号とblind信号を合わせたリストを作成する。このリストはbin番号で昇順に並べ替えられている。
    あるkid信号のindexを取得し、そこを起点に左右に探索する(indexを減少または増加させる)。起点は含まない。
    見つかった信号のbin番号がblind信号のリストに含まれる場合は、それをそのkid信号のblind信号とする。
    """
    info, kids, blinds, powers = kidslist
    all_tones = sorted(list(kids) + list(blinds))
    result = []
    for kid in sorted(kids):
        left_tone = None
        i = all_tones.index(kid) - 1
        while i >= 0:
            tone = all_tones[i]
            if tone in blinds:
                left_tone = tone
                break
            i -= 1
        right_tone = None
        i = all_tones.index(kid) + 1
        while i < len(all_tones):
            tone = all_tones[i]
            if tone in blinds:
                right_tone = tone
                break
            i += 1
        result.append([kid, left_tone, right_tone])
    return result

def load_kidslist(filename):
    """kids.listファイルを読み込みBIN番号やローカル発振器の周波数などを取得する
    引数
    ----
    kids.listファイルへのパス

    戻り値
    ------
    dict ローカル発振器の周波数とフレーム長を辞書にまとめたもの
    list KID信号のBIN番号のリスト
    list blind信号のBIN番号のリスト
    dict BIN番号をキー(Int型)にしたpower(float型)の辞書

    詳細
    ----
    kids.listファイルの最初の2行はヘッダーとして処理する。
    1行目はローカル発振器の周波数、2行目はフレーム長である信じて読み込む。
    3行目以降は「#KIDs」か「#blinds」かBIN番号とpowerのペアが記されているものとする。
    「#KIDs」があった場合はそれ以降のBIN番号はkid_bins変数のリストに追加する。
    「#blinds」があった場合はそれ以降のBIN番号はblind_bins変数のリストに追加する。
    bins変数はkid_binsリストへの参照かblind_binsリストへの参照かになる。
    bins変数によってmode変数やフラグ変数による判定が必要なくなる。

    読み込むkids.listファイルの形式
    -------------------------------
    #LO:5600.0000
    #framelen:19
    #KIDs
    -196483 -75
    -189634 -75
    -188948 -75
        :
        :
    154774  -75
    158692  -78
    169653  -78
    #blinds
    -232027 -75
    -224371 -75
    -216715 -75
        :
        :

    注意
    ----
    ファイルから読み取った1行(line)は必ずstrip()で処理しないと正しく読み込めない。
    bins変数とnbin変数を混同しないように注意されたい。
    bins変数はkid_binsリストかblind_binsリストへの参照であり、nbin変数はBIN番号を一時的に格納する変数である。
    """
    MHz = 1e6
    info = {}
    powers = {}
    kid_bins = []
    blind_bins = []
    with open(filename, 'rt') as f:
        # 最初の2行はヘッダーとして処理する
        line = f.readline().strip() # 行末の空白や改行などを削除する
        name, lofreq = line.split(':')
        info[name[1:]] = float(lofreq)*MHz # name変数に格納されている文字列から最初の"#"を除いた部分をキーとする

        line = f.readline().strip()
        name, framelen = line.split(':')
        info[name[1:]] = int(framelen)

        bins = None # kid_binsリストかblind_binsリストかへの参照になる
        for line in f:
            line = line.strip()
            if line[1:] == 'KIDs':
                bins = kid_bins
                continue
            elif line[1:] == 'blinds':
                bins = blind_bins
                continue
            nbin, power = line.split()
            bins.append(int(nbin))
            powers[int(nbin)] = float(power)
    return info, kid_bins, blind_bins, powers

def load_localsweep(filename, framelen):
    """ローカルsweepデータを読み込みBIN番号とSweepDataオブジェクトを対応づけた辞書として返す
    引数
    ----
    string  ローカルsweepファイルへのパス
    integer フレーム長

    戻り値
    ------
    dict BIN番号とSweepDataオブジェクト対応づけた辞書

    説明
    ----
    元のプログラムではread_localsweep()関数と_read_sweep()関数に分かれていた処理を一緒にした。
    read_localsweep()関数ではOrderedDict()を用いてSweepDataをBIN番号順に並べ替えていたが、
    この後の計算や処理でBIN番号順に並んでいる必要がないので、単純な辞書にSweepDataを格納することとした。
    並べ替えなどの処理を省いた。
    またSweepDataの辞書を作る際に辞書包括表記を用いた。

    参考文献
    --------
    spaceKIDS Readout System Manula(SK-SRON-018-V0.3) page.7
    上記のページにローカルsweepデータの読み方と各種計算方法が記されている。
    """
    MHz = 1e6
    sample_rate = 2e9 # Sample/sec
    rawdata = np.loadtxt(filename)
    bins    = map(int, rawdata[0, 1::3]) #リスト内包表記よりもmap()関数の方が速い
    idata   = rawdata[:, 2::3].T
    qdata   = rawdata[:, 3::3].T
    lofreqs = rawdata[:, 0]*MHz
    dfreq   = sample_rate/(2**framelen) # frequency resolution
    return {b: SweepData(i, q, lofreqs + b*dfreq) for b, i, q in zip(bins, idata, qdata)}

#-----------------------------
# For calc_resonance_params.py
def fit_onepeak(sweepdata, peaks, nfwhm=5):
    """Sweepデータをgaolinbgでフィットする
    引数
    ----
    SweepData Sweepデータオブジェクト
    list      find_peaks()関数で得られる各ピークのパラメタのリスト
    integer   nfwhm

    戻り値
    ------
    MinimizerResult gaolinbgによるフィットの結果
    slice           フィット範囲(インデックス)
    """
    if len(peaks) == 0:
        return (None, None)

    initial_params = gaolinbg_guess(sweepdata, peaks) #params = guess(data)
    gaolinbg_param_names = initial_params.keys()
    s = slice(None)
    if nfwhm > 0:
        if len(peaks) == 1:
            pdict = peaks[0]
        else:
            raise Exception("number of peak not 1: %d" % len(peaks))
        s = adjust_fitrange(nfwhm, len(sweepdata), len(gaolinbg_param_names), pdict)

    # gaolinbg関数をパラメタで全微分した関数を取得
    Dfun = gaolinbg_param_gradient_function()

    def wrapped_Dfun(x, arga, absa, tau, fr, Qr, Qc, phi0, c):
        """全微分の導関数の結果を2次元ベクトルに変換する
        説明
        ----
        lmfit.Minimize()関数では複素関数を扱えない為、このような変換を行う。
        """
        return complex_to_cartesian2darray(Dfun(x, arga, absa, tau, fr, Qr, Qc, phi0, c))

    def Dfun_pars(p, x, dat=None):
        v = dict((k, v) for (k, v) in p.valuesdict().items() if k in gaolinbg_param_names)
        return wrapped_Dfun(x, **v)

    kws = {}
    kws['Dfun'] = Dfun_pars
    kws['col_deriv'] = 1

    gaolinbg_func = gaolinbg_function()
    ys = complex_to_cartesian2darray(sweepdata.iq[s])
    def residue(p, x):
        """データとの差を計算する
        引数
        ----
        lmfit.Parameters
        numpy.array 周波数

        戻り値
        ------
        lmfit.MinimizerResult フィットの結果

        注意
        ----
        lmfit.minimize()は複素数の最小化はできない。
        なので、複素数を一旦2次元ベクトルに変換して計算する。
        """
        v = {k: v for (k, v) in p.valuesdict().items() if k in gaolinbg_param_names}
        return complex_to_cartesian2darray(gaolinbg_func(x, **v)) - ys

    with catch_warnings():
        simplefilter("ignore")
        params = lmfit.Parameters()
        params.add('arga', initial_params['arga'])
        params.add('absa', initial_params['absa'])
        params.add('tau',  initial_params['tau'])
        params.add('fr',   initial_params['fr'])
        params.add('Qr',   initial_params['Qr'])
        params.add('Qc',   initial_params['Qc'])
        params.add('phi0', initial_params['phi0'])
        params.add('c',    initial_params['c'])
        params.add('Qi',   expr='1/(1/Qr - 1/Qc*cos(phi0))')
        minimizer_result = lmfit.minimize(residue, params, args=(sweepdata.x[s],), **kws)
        return minimizer_result, s # fit結果, フィット範囲

def search_peak(sweepdata, peaks):
    """複数のピークが見つかった場合は中心周波数に一番近いピークを採用する
    引数
    ----
    SweepData Sweepデータオブジェクト
    list      find_peaks()関数によって得られる各ピークのパラメタをまとめたリスト

    戻り値
    ------
    dict 選択されたピークのパラメタをまとめた辞書
    """
    center  = (sweepdata.x[0] + sweepdata.x[-1])/2.0
    minind  = 0
    mindist = abs(peaks[0]['f0'] - center)
    for i in range(1, len(peaks)):
        if mindist > abs(peaks[i]['f0'] - center):
            minind  = i
            mindist = abs(peaks[i]['f0'] - center)
    return peaks[minind]

def find_peaks(freq, ampl, fc, smooth=None, threshold=None, minq=10000, maxratio=0.5):
    """Search peaks in data.
    引数
    ----
    freq      周波数の配列
    ampl      amplitudeの配列
    fc        TODのキャリア周波数
    smooth    平滑化のwindow長(正の整数)[オプション]
    threshold ピークの判定基準[オプション]
    minq      Qの最小値(ピークの良し悪しの判定に利用する)[オプション]
    maxratio  (bg - d)/bgの最小値(ピークの良し悪しの判定に利用する)[オプション]

    戻り値
    ------
    list Lorentzianのフィット共にgaolinbgのパラメタの初期値を決定するのに必要なパラメタを辞書にまとめたリスト
    """
    dfreq = freq[1] - freq[0] # frequency resolution
    if smooth is None:
        smooth = 15

    # 平滑化
    deriv2 = scipy.signal.savgol_filter(ampl, oddify(smooth), 3, 2, delta=dfreq)
    ampl_s = scipy.signal.savgol_filter(ampl, oddify(smooth), 3, 0, delta=dfreq)

    # ピーク探索
    if threshold is None:
        d2sigma   = np.std(deriv2)
        d2sigma   = np.std(deriv2[np.abs(deriv2)<2*d2sigma])
        threshold = 3*d2sigma
    kid_indices = []
    for i in np.where(deriv2 > threshold)[0]:
        if i < len(deriv2) - 1 and deriv2[i - 1] <= deriv2[i] and deriv2[i] >= deriv2[i + 1]:
            kid_indices.append(i)

    kids = [] # 各KIDをLorentzianでフィットした結果(Lorentzianのパラメタ)
    nbadq = 0 # Qが条件を満たさなかったKIDの数
    nbadd = 0 # maxratioを満たさなかったKIDの数
    if not kid_indices:
        return []

    for i in kid_indices:
        #
        # 詳細
        # ----
        # ピークの左右でamplitudeの2階の微係数が0になる部分(変曲点)を利用してフィット範囲を設定する。
        # フィットの結果Qと(bg-d)/bgが条件を満たさない場合は、そのピークは無視される。
        #
        l, r = i, i
        while l > 0           and deriv2[l] > 0: l -= 1
        while r < len(deriv2) and deriv2[r] > 0: r += 1
        w = (r - l + 1)*dfreq
        w = w*6.0/np.sqrt(3) # convert to FWHM
        l = int((l - i)*6.0/np.sqrt(3) + i)
        r = int((r - i)*6.0/np.sqrt(3) + i)
        if l < 0         : l = 0
        if r >= len(freq): r = len(freq) - 1

        q0 = freq[i]/w
        f, q, d, bg = fitLorentzian(freq[l:r], ampl[l:r], freq[i], q0)
        ##### refitting by extending fit range
        if (bg - d)/bg > maxratio:
            n = 10
            if l - n >= 0 and r + n < len(freq):
                f, q, d, bg = fitLorentzian(freq[l-n:r+n], ampl[l-n:r+n], freq[i], q0)
                logger.debug('Rough refitting with Lorentzian')
        if q < minq:
            nbadq += 1
            continue
        if (bg - d)/bg > maxratio:
            nbadd += 1
            continue
        kids.append((f, q, d, bg))
    del l, r, f, q, d, bg
    if nbadq > 0:
        logger.debug(f'Removed {nbadq} peaks with bad Q')
    if nbadd > 0:
        logger.debug(f'Removed {nbadd} peaks with bad S21min')
    kids.sort() #sort by frequency

    # pick up a peak which is closest to the carrier frequency (when fc in freq range)
    if len(kids) > 0 and fc >= freq[0] and fc <= freq[-1]:
        i = np.argmin(abs(np.array(kids).T[0] - fc))
        kids = [kids[i]]

    # gaolinbgのパラメタの初期値を決めるためのbackgroundなどを計算する
    for i, (f, q, depth, bg) in enumerate(kids):
        f0ind = np.argmin(abs(freq - f))
        w     = f/q
        dl    = w/2.0/dfreq
        dr    = w/2.0/dfreq
        bg_l  = ampl_s[max(int(f0ind - 3*dl), 0)]
        if int(f0ind - 3*dl)<0:
            bg_r = ampl_s[0]
        else:
            bg_r  = ampl_s[min(int(f0ind - 3*dl), len(freq) - 1)]
        a_off = (bg_l + bg_r)/2.0
        a_on  = ampl_s[f0ind]
        kids[i] = {
            'Q': q,
            'f0': f,
            'f0ind': f0ind,
            'dl': int(dl),
            'dr': int(dr),
            'a_off': a_off,
            'a_on': a_on,
            'ampl': depth,
            'bg': bg
        }
    return kids

def oddify(n):
    """n/2を切り上げそれを2倍して1足し奇数にして返す
    引数
    ----
    integer

    戻り値
    ------
    integer 奇数
    """
    return int(np.ceil(n/2.0)*2+1)

def fitLorentzian(freq, ampl, f0, q0):
    """Fit data with Lorentzian curve.
    引数
    ----
    freq a 1-D array of frequency
    ampl a 1-D array of amplitude
    f0   initial parameter for center frequency
    q0   initial parameter for quality factor

    戻り値
    ------
    (fc, q, d, bg) フィットの結果をまとめたタプル

    パラメタの説明
    --------------
    f : center frequency
    q : quality factor
    d : amplitude for Lorentzian curve
    bg: constant background level
    """
    def f(x):
        (a, b, c, d) = x # background, amplitude, 2/FWHM, freq. center
        y = a + b / (((freq-d)*c)**2 + 1)
        return y - ampl
    a = np.median(ampl)
    b = -0.8 * a
    c = 2.0 * q0 / f0
    d = f0
    x0 = np.array([a, b, c, d])

    with catch_warnings():
        simplefilter("ignore")
        x1 = scipy.optimize.leastsq(f, x0)[0]

    (a, b, c, d) = x1
    fc = d
    q = abs(c * d / 2.0)
    return (fc, q, -b, a)

def adjust_fitrange(nfwhm, ndata, nparams, peakparams, factor=1):
    """Adjust fit range to make sure :math:`N_{free} \geq nparam`
    引数
    ----
    float   nfwhm: minimum fit length in unit of FWHM
    integer ndata: length of data
    integer nparams: number of parameter
    dict    peakparams: peak information
    integer factor

    戻り値
    ------
    tuple (rbegin, rend) フィット範囲の始まりと終わりのインデックス

    factorについて
    --------------
    Specify 2 as 'factor' if you fit complex variable(as 2d-vector).

    利用するピークパラメタ
    ----------------------
     - f0ind: a integer, index of peak center
     - dl   : a integer, index count to reach the left half-maximum
     - dr   : a integer, index count to reach the right half-maximum
    """
    f0ind = peakparams['f0ind']
    l, c, r = f0ind-peakparams['dl'], f0ind, f0ind+peakparams['dr']
    if (r - l + 1)*factor >= nparams:
        rbegin, rend = int(c-nfwhm*(c-l)), int(c+nfwhm*(r-c))
    else:
        n = nparams/(float(r - l + 1)*factor)
        rbegin, rend = int(f0ind-nfwhm*n*(c-l)), int(c+nfwhm*n*(r-c))
    if rbegin < 0:
        if rend + (-rbegin) >= ndata:
            rend = ndata-1
            rbegin = 0
        else:
            rend = rend + (-rbegin)
            rbegin = 0
    if rend >= ndata:
        if rbegin - (rend - ndata) < 0:
           rbegin = 0
           rend = ndata-1
        else:
           rbegin = rbegin - (rend - ndata)
           rend = ndata-1
    if (rend - rbegin + 1)*factor < 11:
        raise Exception("Fit range guess error")
    return slice(rbegin, rend)

#-------------------------
# For make_reduced_fits.py
def readout_dict():
    hdr_key_lis = ['EXTNAME', 'FILENAME', 'FRAMERT', 'FRAMELEN', 'DSAMPLE',]
    hdr_val_lis = ['READOUT', None, None, None, None]
    hdr_com_lis = ['name of binary data',
                   'input filename',
                   'sampling rate',
                   '2-log of frame length',
                   'number of down sampling',
                   'label for field 0', 'data format of field 0',
                   'label for field 1', 'data format of field 1',]
    cols_key_lis = ['timestamp', 'pixelid',]
    cols_data_lis = []
    tform = ['D', 'I',]
    tunit = [None, None,]

    r_dict = {'hdr_key_lis': hdr_key_lis,
              'hdr_val_lis': hdr_val_lis,
              'hdr_com_lis': hdr_com_lis,
              'cols_key_lis': cols_key_lis,
              'cols_data_lis': cols_data_lis,
              'tform': tform,
              'tunit': tunit}
    return r_dict

def kids_dict():
    hdr_key_lis = ['EXTNAME', 'FILENAME',]
    hdr_val_lis = ['KIDSINFO', None,]
    hdr_com_lis = ['name of binary data',
                   'localsweep filename',
                   'label for field 0', 'data format of field 0',
                   'label for field 1', 'data format of field 1',
                   'label for field 2', 'data format of field 2', 'data unit of field 2',
                   'label for field 3', 'data format of field 3', 'data unit of field 3',
                   'label for field 4', 'data format of field 4',
                   'label for field 5', 'data format of field 5', 'data unit of field 5',
                   'label for field 6', 'data format of field 6',
                   'label for field 7', 'data format of field 7',
                   'label for field 8', 'data format of field 8']
    cols_key_lis = ['pixelid', 'kidid', 'Pread', 'fc', 'yfc, linyfc',
                    'fr, dfr (300K)', 'Qr, dQr (300K)', 'Qc, dQc (300K)', 'Qi, dQi (300K)']
    cols_data_lis = []
    tform = ['I', 'I', 'E', 'E', '2E', '2E', '2E', '2E', '2E']
    tunit = [None, None, 'dBm', 'GHz', None, 'GHz', None, None, None]

    k_dict = {'hdr_key_lis': hdr_key_lis,
              'hdr_val_lis': hdr_val_lis,
              'hdr_com_lis': hdr_com_lis,
              'cols_key_lis': cols_key_lis,
              'cols_data_lis': cols_data_lis,
              'tform': tform,
              'tunit': tunit}
    return k_dict

def rebin_array(inarray, rebin):
    return inarray[:(inarray.size//rebin)*rebin].reshape(-1, rebin).mean(axis=1)

def find_glitch(yss,
                baseline_thresh = 6.0, glitch_thresh = 5.0, clusterize_thresh = 2, offset = 0,
                opt_diff2=True, smrange=1, returnval=False):
    """
    find glitches common to `yss`, assuming glitch exists at the same
    time of all data in `yss`. Too close glitches or broad glitch are
    treated as one glitch.

    引数
    ----
    yss              : an array of 1-D arrays of data to deglitch
    baseline_thresh  : threshold to decide baseline
    glitch_thresh    : threshold to decide as glitch
    clusterize_thresh: if gap between glitches are less than or equal to this, treat them as one glitch.
    offset           : an integer offset added to returned bad array
    opt_diff2        : bool to select whether to use 2nd derivative for glitch identification
    smrange          : an integer to set range of rolling mean (used for glitch identification based on std)
    returnval        : bool to choose whether to return threshold value or not

    戻り値(returnvalの値によって変化する)
    -------------------------------------
    an 1-D array of boolean, that can be used as index of ys (elements of yss).
    (option) threshold value for glitch identification (return when returnval=True)
    """
    ave = np.average(yss, axis=0)
    xs  = np.arange(len(yss[0]))
    dx  = xs[1] - xs[0]
    if opt_diff2:
        diff2 = np.array( _numdif_2(ave, dx) )
        sigma = np.std(diff2)
        good  = (np.abs(diff2) < (baseline_thresh*sigma))
        sigma = np.std(diff2[good])
        bad   = (np.abs(diff2) >= (glitch_thresh*sigma))
        if returnval:
            thre = glitch_thresh*sigma/np.sqrt(6)*(dx**2)
    else:
        ##### rolling mean
        weights = np.ones(smrange)/smrange
        assert( len(ave)>len(weights) )
        smooth  = np.convolve(ave, weights, 'same')

        invalid = int( (len(weights)-1)/2 )
        smooth[:invalid] = np.mean(ave[:invalid]) * np.ones( len(ave[:invalid]) )
        smooth[-invalid-1:] = np.mean(ave[-invalid-1:]) * np.ones( len(ave[-invalid-1:]) )

        ave_sm = ave - smooth

        mean  = np.mean(ave_sm)
        sigma = np.std(ave_sm)
        good  = (np.abs(ave_sm-mean) < (baseline_thresh*sigma))
        mean  = np.mean(ave_sm[good])
        sigma = np.std(ave_sm[good])
        bad   = (np.abs(ave_sm-mean) >= (glitch_thresh*sigma))
        if returnval:
            thre = glitch_thresh*sigma
    ## treat broad glitch (or too close glitches) as one glitch
    bad = _clusterize_indices(bad, clusterize_thresh)
    if opt_diff2:
        bad = np.concatenate( ([bad[0]], bad, [bad[-1]]) )

    bad_ = np.copy(bad)
    if offset>0:
        for i,b in enumerate(bad):
            if b:
                bad_[i] = True
                for j in range(offset):
                    if (i+j+1)<len(bad): bad_[i+j+1] = True
    if returnval:
        return bad_, thre
    else:
        return bad_

def _numdif_2(y, dx):
    return (y[2:]+y[:-2]-2*y[1:-1])/dx**2

def _clusterize_indices(indices, threshold):
    """Fill small gap (within `threshold`) in indices.
    引数
    ----
    indices:   an 1-D array of Boolean
    threshold: an interger, allowed gap between True's in indices

    戻り値
    ------
    list 真理値のリスト

    実行結果の例
    ------------
     If threshold == 1,
      [True, False, False, True, True] => [True, False, False, True, True].
     If threshold == 2,
      [True, False, False, True, True] => [True, True, True, True, True].
    """

    results = np.copy(indices)
    prev  = 0
    first = True
    for i, x in enumerate(results):
        if x:
            if (not first) and i - prev <= threshold + 1:
                for j in range(prev, i):
                    results[j] = True
            prev  = i
            first = False
    return results

def calibrate_with_blind_tones(signal, blind_left, blind_right):
    """do blind tone calibration.
    引数
    ----
    signal     : signal tone TOD.
    blind_left : left blind tone TOD.
    blind_right: right blind tone TOD.

    戻り値
    ------
    a 1-D array of complex value (calibrated IQ)

    詳細
    ----
    each parameter is either a FixedData or a tuple (freq, IQ)
    where freq is a frequency [Hz],
    IQ is a 1-D array of IQ data as complex value.
    """
    def get_f_iq(fx):
        if isinstance(blind_left, FixedData):
            return fx._frequency, fx.iq
        else:
            return fx
    f_s, iq_s = get_f_iq(signal)
    f_l, iq_l = get_f_iq(blind_left)
    f_r, iq_r = get_f_iq(blind_right)

    calib_l = iq_l / np.average(iq_l)
    calib_r = iq_r / np.average(iq_r)
    calib_s = calib_l + (calib_r - calib_l)/(f_r - f_l) * (f_s - f_l) # interpolate

    calibrated = (iq_s / calib_s)
    return calibrated

def interpolate_bad(ys, bad):
    """Linear interpolate `ys` in region `bad`.
    引数
    ----
    ys : an 1-D array
    bad: an 1-D array of boolean, with True indicating bad region in `ys`.

    戻り値
    ------
    an 1-D array, with bad region interpolated.
    """
    xs  = np.arange(len(ys))
    deglitched = np.interp(xs, xs[~bad], ys[~bad])
    return deglitched

def createBinTableHDU(data_dict):
    #-------- Set Header and Comments
    header = fits.Header()
    for (i, j, k) in zip(data_dict['hdr_key_lis'], data_dict['hdr_val_lis'], data_dict['hdr_com_lis']):
        header[i] = j, k
    #-------- Create Collumns of the Binary Table
    columns = []
    for i in range(len(data_dict['cols_key_lis'])):
        columns.append(fits.Column(name=data_dict['cols_key_lis'][i],
                                    format=data_dict['tform'][i],
                                    array=data_dict['cols_data_lis'][i],
                                    unit=data_dict['tunit'][i]))
    hdu = fits.BinTableHDU.from_columns(columns, header)
    #-------- Add comments
    for i in range(-(len(data_dict['hdr_com_lis'])-data_dict['hdr_com_lis'].index('label for field 0')), 0):
        addHeaderComments(hdu.header, list(hdu.header)[i], data_dict['hdr_com_lis'][i])
    return hdu

def addHeaderComments(hdr, key, com):
    """Add Comments to Header"""
    hdr.comments[key] = com

def f_Lorentzian(params, freq, data=None):
    """Lorengzianモデル関数
    mkid_data/peak_search.pyで使われている関数と同じ。
    引数のdataが指定された場合はモデル関数の計算値とdataの差を返す。

    引数
    ----
    numpy.Parameters
    numpy.array 周波数
    numpy.array amplitude

    戻り値
    ------
    モデルの計算値|計算値とdataの差
    """
    params = params.valuesdict()
    b      = params['bg']
    a      = params['amp']
    fc     = params['f0']
    fwhm   = params['fwhm']
    model = b + a/(((freq - fc)*(2.0/fwhm))**2 + 1)
    if data is None:
        return model
    return model - data

def gaolinbg_guess(sweepdata, peaks):
    gaolinbg_paramnames = 'arga absa tau fr Qr Qc phi0 c'.split()
    x      = sweepdata.x
    deltax = x[1] - x[0]
    pdict  = search_peak(sweepdata, peaks)
    y0     = sweepdata.iq[pdict['f0ind']]
    FWHM   = pdict['f0']/pdict['Q']
    ddeg   = sweepdata.deg[1:] - sweepdata.deg[:-1]
    tau    = -np.average(ddeg[abs(ddeg)<180])*np.pi/180.0/deltax/2/np.pi
    f0     = pdict['f0']
    theta  = np.angle(y0)
    arga   = np.angle(y0*np.exp(1j*2*np.pi*tau*f0))
    absa   = pdict['a_off']
    fr     = f0
    Qr     = f0/FWHM
    Qc     = Qr
    phi0   = 0
    c      = 0
    return dict(zip(gaolinbg_paramnames, (arga, absa, tau, fr, Qr, Qc, phi0, c)))

def gaolinbg_rewind(x, y, arga, absa, tau, fr, Qr, Qc, phi0, c):
    tmp = y/absa/np.exp(-1j*(2*np.pi*x*tau - arga)) - c*(x - fr)
    return (tmp - 1)*Qc/Qr/np.exp(1j*phi0) + 0.5

def gaolinbg_expression():
    arg_symbols = sympy.symbols(r'(\arg{a}) |a| tau f_r Q_r Q_c phi_0 c')
    arg_names = 'arga absa tau fr Qr Qc phi0 c'.split()
    arga, absa, tau, fr, Qr, Qc, phi0, c = arg_symbols
    x = sympy.symbols('x')
    from sympy import exp, I, pi
    expr = (absa*exp(-I*(2*pi*x*tau - arga))*(1 + c*(x - fr) - Qr/Qc*exp(I*phi0)/(1 + 2*I*Qr*((x - fr)/fr))))
    return expr, arg_symbols, arg_names

def gaolinbg_bg_expression():
    arg_symbols = sympy.symbols(r'(\arg{a}) |a| tau f_r Q_r Q_c phi_0 c')
    arg_names = 'arga absa tau fr Qr Qc phi0 c'.split()
    arga, absa, tau, fr, Qr, Qc, phi0, c = arg_symbols
    x = sympy.symbols('x')
    from sympy import exp, I, pi
    expr = (absa*exp(-I*(2*pi*x*tau - arga))*(1 + c*(x - fr)))
    return expr, arg_symbols, arg_names

def gaolinbg_function():
    expr, arg_symbols, arg_names = gaolinbg_expression()
    expr_ = expr.subs(zip(arg_symbols, arg_names))
    args = sympy.symbols(['x'] + arg_names)
    return sympy.lambdify(args, expr_, 'numpy', dummify=False)

def gaolinbg_bg_function():
    expr, arg_symbols, arg_names = gaolinbg_bg_expression()
    expr_ = expr.subs(zip(arg_symbols, arg_names))
    args = sympy.symbols(['x'] + arg_names)
    return sympy.lambdify(args, expr_, 'numpy', dummify=False)

def gaolinbg_param_gradient_function():
    expr, arg_symbols, arg_names = gaolinbg_expression()
    funcargs = ['x'] + arg_names
    arga, absa, tau, fr, Qr, Qc, phi0, c = arg_symbols
    gradfuncs = [ #各パラメータに関する導関数を求め、関数化(Pythonで使える)する。
        sympy.lambdify(funcargs, sympy.diff(expr, arga).subs(zip(arg_symbols, arg_names)), 'numpy', dummify=False),
        sympy.lambdify(funcargs, sympy.diff(expr, absa).subs(zip(arg_symbols, arg_names)), 'numpy', dummify=False),
        sympy.lambdify(funcargs, sympy.diff(expr, tau ).subs(zip(arg_symbols, arg_names)), 'numpy', dummify=False),
        sympy.lambdify(funcargs, sympy.diff(expr, fr  ).subs(zip(arg_symbols, arg_names)), 'numpy', dummify=False),
        sympy.lambdify(funcargs, sympy.diff(expr, Qr  ).subs(zip(arg_symbols, arg_names)), 'numpy', dummify=False),
        sympy.lambdify(funcargs, sympy.diff(expr, Qc  ).subs(zip(arg_symbols, arg_names)), 'numpy', dummify=False),
        sympy.lambdify(funcargs, sympy.diff(expr, phi0).subs(zip(arg_symbols, arg_names)), 'numpy', dummify=False),
        sympy.lambdify(funcargs, sympy.diff(expr, c   ).subs(zip(arg_symbols, arg_names)), 'numpy', dummify=False),
    ]
    def total_derivative(x, arga, absa, tau, fr, Qr, Qc, phi0, c):
        return np.array(np.broadcast_arrays(
            gradfuncs[0](x, arga, absa, tau, fr, Qr, Qc, phi0, c),
            gradfuncs[1](x, arga, absa, tau, fr, Qr, Qc, phi0, c),
            gradfuncs[2](x, arga, absa, tau, fr, Qr, Qc, phi0, c),
            gradfuncs[3](x, arga, absa, tau, fr, Qr, Qc, phi0, c),
            gradfuncs[4](x, arga, absa, tau, fr, Qr, Qc, phi0, c),
            gradfuncs[5](x, arga, absa, tau, fr, Qr, Qc, phi0, c),
            gradfuncs[6](x, arga, absa, tau, fr, Qr, Qc, phi0, c),
            gradfuncs[7](x, arga, absa, tau, fr, Qr, Qc, phi0, c)))

    return total_derivative

def complex_to_cartesian2darray(x):
    x = np.atleast_1d(x)
    shape = x.shape
    return np.concatenate([np.real(x), np.imag(x)], axis=len(shape)-1)

def cartesian2darray_to_complex(x):
    assert x.shape[-1] % 2 == 0
    size = x.shape[-1] / 2
    return x[...,:size] + 1j*x[...,size:]

def power_spectrum_density(data, dt, ndivide=1, window=scipy.signal.windows.hann, overwrap_half=False):
    """Calculate power spectrum density of data.
    引数
    ----
    data         : input data
    dt           : time between each data
    ndivide      : do averaging(split data into ndivide, get psd of each, and average them)
    overwrap_half: a boolean, split data to half-overwrapped regions

    戻り値
    ------
    (frequencies, psd)
    """
    if overwrap_half:
        step = int(len(data)/(ndivide + 1))
        size = step*2
    else:
        step = int(len(data)/ndivide)
        size = step
    if bin(size).count("1") != 1:
        if overwrap_half:
            logger.warning('((length of data)/(ndivide+1))*2 is not power of 2: %d' % size)
        else:
            logger.warning('(length of data)/ndivide is not power of 2: %d' % size)
    psd = np.zeros(size)
    T = (size - 1)*dt
    vs = 1/dt
    vk_ = scipy.fft.fftfreq(size, dt)
    vk = vk_[np.where(vk_ >= 0)]
    for i in range(ndivide):
        d = data[i*step:i*step+size]
        if window is None:
            w    = np.ones(size)
            corr = 1.0
        else:
            w    = window(size)
            corr = np.mean(w**2)
        psd  = psd + 2*(np.abs(scipy.fft.fft(d*w)))**2/size*dt/corr
    return vk, psd[:len(vk)]/ndivide

class SweepData():
    """mkid_data.data.SweepDataを単純化したクラス"""
    def __init__(self, i=None, q=None, freq=None):
        self.i = i
        self.q = q
        self.freq = freq
        return

    @property
    def x(self):
        GHz = 1e9
        return self.freq/GHz

    @property
    def amplitude(self):
        return np.sqrt(self.i**2 + self.q**2)

    @property
    def iq(self):
        return self.i + 1j*self.q

    @property
    def deg(self):
        return np.arctan2(self.q, self.i)*180/np.pi

    def __len__(self):
        return len(self.freq)

class FixedData():
    """mkid_data.data.FixedDataを単純化したクラス"""
    def __init__(self):
        self.timestamp = None
        self.freq = None
        self._frequency = None
        self._i = None
        self._q = None
        self.fftgain = None
        self.info = None
        self.rate = 1
        return

    @property
    def frequency(self):
        GHz = 1e9
        return self._frequency/GHz

    @property
    def iq(self):
        return self.i + 1j*self.q

    @property
    def i(self):
        return self._i/self.fftgain

    @property
    def q(self):
        return self._q/self.fftgain

    @property
    def t(self):
        return self.timestamp/self.rate

class FixedFitData():
    """mkid_data.data.FixedFitDataを単純化したクラス"""
    def __init__(self, freq=None, timestamp=None, ampl=None, phase=None, info=None):
        self._freq = freq
        self._timestamp = timestamp
        self._ampl = ampl
        self._phase = phase
        self._info = info
        return

    @property
    def frequency(self):
        GHz = 1e9
        return self._freq/GHz

    @property
    def amplitude(self):
        return self._ampl

    @property
    def phase(self):
        return self._phase

    @property
    def t(self):
        return self._timestamp

    @property
    def info(self):
        return self._info

    def unpack(self):
        return (self._timestamp, self._ampl, self._phase)

class TODs(Mapping):
    """TODデータを保持するクラス
    注意
    ----
    このクラスのインスタンスを作成すると自動的にfitsファイルが開かれる。
    TODに関する処理が終わったら必ずclose()メソッドを明示的に呼び出すこと。

    コピー元
    --------
    bbsweeplib.kids.py fits_tods_fits
    """
    chunksize = int(1e6)
    kSampleRate = 2e9
    def __init__(self, infile):
        """コンストラクタ
        引数
        ----
        string TODが保存されたfitsファイルの名前

        戻り値
        ------
        なし
        """
        self.infile  = infile
        self.open()

    def open(self):
        self.hud = fits.open(self.infile)
        bintable = self.hud[1]
        self.fftgain  = bintable.header['fftgain']
        self.framert  = bintable.header['framert']
        self.nbins    = bintable.header['nbins']
        self.npoints  = bintable.header['npoints']
        self.lofreq   = bintable.header['lofreq']
        self.bins     = [_to_nbit_signed(bintable.header['BIN%d' % i], self.npoints)
                         for i in range(self.nbins)]
        self.if_freq      = self.kSampleRate * np.array(self.bins) / 2**self.npoints
        self.carrier_freq = self.if_freq + self.lofreq

        bindata = bintable.data
        timestamp_ = bindata.field('timestamp')
        framenr_ = bindata.field('framenr')
        self.timestamp = timestamp_
        self.framenr = framenr_
        self._read_bins = {}

        self.offset = -100 # remove last part of data (sometimes they behave bad)

    def close(self):
        self.hud.close()

    def __len__(self):
        return len(self.bins)

    def __contains__(self, key):
        return key in self.bins

    def __iter__(self):
        for key in sorted(self.bins):
            yield key

    def __getitem__(self, key):
        if key in self._read_bins:
            return self._read_bins[key]
        else:
            ind = self.bins.index(key)
            chunksize=int(1e4)
            name = 'kid_%04d' % ind
            rawdata = self.hud[1].data['data']
            read_I_ = rawdata[:self.offset, 2*ind    ]
            read_Q_ = rawdata[:self.offset, 2*ind + 1]

            info = dict()
            info['bins'] = self.bins
            info['freqs'] = self.carrier_freq[ind]
            info['header'] = self.hud[1].header
            d = FixedData()
            d.timestamp = self.timestamp[:self.offset]
            d._i = read_I_
            d._q = read_Q_
            d.fftgain = self.fftgain
            d._frequency = self.carrier_freq[ind]
            d.info = info
            self._read_bins[key] = d
            return d

    def __getstate__(self):
        if hasattr(self, 'hud'):
            dic = self.__dict__.copy()
            dic['_read_bins'] = {}
            dic['hud'] = None
            return dic
        else:
            return self.__dict__

    def __setstate__(self, dic):
        self.__dict__ = dic.copy()
        if 'hud' in dic:
            self.open()

def _to_nbit_signed(x, n):
    if x > 2**(n-1):
        return -((~x & (2**n-1))+1)
    else:
        return x

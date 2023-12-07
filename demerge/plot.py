"""1つのKIDについてSweepのフィットの結果をプロットする。

Filename: plot.py
(C) 2021 内藤システムズ

使い方
------
第一引数 フィットの結果が格納された一時ファイル(*.pkl)の名前
第二引数 グラフを保存するディレクトリ名

例)
  $ python plot.py cache/kid00001.pkl graph
"""
# standard library
import sys
import pickle
from logging import basicConfig, getLogger


# dependencies
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from . import demerge


# module logger
logger = getLogger(__name__)


def plot(kid, oldkid=None):
    if kid['enabled'] == False:
        return None, None
    if kid['resonance_params'] is None:
        return None, None
    arga = kid['resonance_params']['arga']
    absa = kid['resonance_params']['absa']
    tau = kid['resonance_params']['tau']
    fr = kid['resonance_params']['fr'].value # GHz
    Qr = kid['resonance_params']['Qr'].value
    Qc = kid['resonance_params']['Qc'].value
    phi0 = kid['resonance_params']['phi0']
    c = kid['resonance_params']['c']

    fitrange = kid['fitrange']
    rewound_localsweep = demerge.gaolinbg_rewind(kid['localsweep'].x, kid['localsweep'].iq, arga, absa, tau, fr, Qr, Qc, phi0, c)
    fitted_localsweep = demerge.gaolinbg_function()(kid['localsweep'].x, arga, absa, tau, fr, Qr, Qc, phi0, c)
    limited_fitted_localsweep = demerge.gaolinbg_function()(kid['localsweep'].x[fitrange], arga, absa, tau, fr, Qr, Qc, phi0, c)
    rewound_fitted_localsweep = demerge.gaolinbg_rewind(kid['localsweep'].x, fitted_localsweep, arga, absa, tau, fr, Qr, Qc, phi0, c)
    fitted_background = demerge.gaolinbg_bg_function()(kid['localsweep'].x, arga, absa, tau, fr, Qr, Qc, phi0, c)
    limited_fitted_background = demerge.gaolinbg_bg_function()(kid['localsweep'].x[fitrange], arga, absa, tau, fr, Qr, Qc, phi0, c)
    rewound_tod = demerge.gaolinbg_rewind(kid['tod'].frequency, kid['tod'].iq[::100], arga, absa, tau, fr, Qr, Qc, phi0, c)
    fitted_peak = demerge.gaolinbg_function()(fr, arga, absa, tau, fr, Qr, Qc, phi0, c)

    plt.ioff()
    fig1 = plt.figure(figsize=(16,10))
    ax1 = plt.subplot(231)
    ax2 = plt.subplot(232)
    ax3 = plt.subplot(233)
    ax4 = plt.subplot(212)

    # Sweepのプロット(KidFitResult::plot()から抜粋)
    ax1.plot(kid['localsweep'].x, abs(kid['localsweep'].iq),'.g', label='Ampl.')
    ax1.plot(kid['localsweep'].x, kid['localsweep'].i,'.b', label='I')
    ax1.plot(kid['localsweep'].x, kid['localsweep'].q,'.r', label='Q')
    ax1.plot(kid['localsweep'].x, abs(fitted_background), '-', color='gray', label='bg')
    ax1.plot(kid['localsweep'].x[fitrange], abs(limited_fitted_localsweep),'-y')
    ax1.plot(kid['localsweep'].x[fitrange], np.real(limited_fitted_localsweep),'-c')
    ax1.plot(kid['localsweep'].x[fitrange], np.imag(limited_fitted_localsweep),'-m')
    ax1.plot(fr, np.abs(fitted_peak),'r*')
    ax1.set_title('Amplitude vs Freq')
    ax1.set_xlabel('Frequency [GHz]')
    ax1.set_ylabel('Amplitude')
    ax1.axhline(color='k')
    ax1.legend(loc='best')
    ax1.grid()

    #KidFitResult::plot()から抜粋
    ax2.plot(kid['localsweep'].i, kid['localsweep'].q, '.b', label='data')
    ax2.plot(np.real(limited_fitted_localsweep), np.imag(limited_fitted_localsweep), '-c', label='fit')
    ax2.plot(np.real(fitted_peak), np.imag(fitted_peak), 'r*')
    ax2.plot(np.real(fitted_background), np.imag(fitted_background), '-', color='gray', label='bg')
    if not oldkid is None:
        ax2.plot(oldkid['I'], oldkid['Q'], '.b', label='data')
        ax2.plot(np.real(oldkid['fity']), np.imag(oldkid['fity']), '-c', label='fit')
        ax2.plot(np.real(fitted_peak), np.imag(fitted_peak), 'r*')
        ax2.plot(np.real(oldkid['bg']), np.imag(oldkid['bg']), '-', color='gray', label='bg')

    ax2.set_title('KID[{}], I vs Q'.format(kid['kidid']))
    ax2.set_xlabel('I')
    ax2.set_ylabel('Q')
    ax2.axhline(color='k')
    ax2.axvline(color='k')
    ax2.legend(loc='best')
    ax2.axis('equal')
    ax2.grid()

    #FitSweep()から抜粋
    ax3.set_title('rewind')
    ax3.plot(np.real(rewound_localsweep), np.imag(rewound_localsweep), label='Sweep')
    ax3.plot(np.real(rewound_tod), np.imag(rewound_tod), '.', alpha=0.2, label='TOD')
    ax3.plot(np.real(rewound_fitted_localsweep), np.imag(rewound_fitted_localsweep), 'y', lw=3, alpha=0.5, label='fit')
    ax3.set_xlabel('Re')
    ax3.set_ylabel('Im')
    ax3.axhline(color='k')
    ax3.axvline(color='k')
    ax3.legend(loc='best')
    ax3.axis('equal')
    ax3.grid()

    # TODのプロット
    phase_tod = -np.angle(-rewound_tod)
    phase_localsweep = -np.angle(-rewound_localsweep)
    phase_fitted_localsweep = -np.angle(-rewound_fitted_localsweep)
    ax4.plot(kid['localsweep'].x, phase_localsweep, label='Sweep')
    ax4.plot(kid['localsweep'].x, phase_fitted_localsweep, 'y', lw=3, alpha=0.5, label='fit')

    ## spline interpolation
    try:
        tck = interpolate.splrep(phase_fitted_localsweep, kid['localsweep'].x, s=0)
        x = interpolate.splev(phase_tod, tck, der=0)
        ax4.plot(x, phase_tod, '.', alpha=0.2, label='TOD')
    except ValueError:
        logger.warning('KID[{}] TOD interpolation failed..'.format(kid['kidid']))

    tck = interpolate.splrep(kid['localsweep'].x, phase_localsweep, s=0)
    y = interpolate.splev(kid['tod'].frequency, tck, der=0)
    ax4.plot(kid['tod'].frequency, y, 'go', label='carrier f')
    y = interpolate.splev(fr, tck, der=0)
    ax4.plot(fr, y, 'r*', label='resonance f')
    ax4.set_title('KID[{}], Phase vs Freq'.format(kid['kidid']))
    ax4.set_xlabel('Frequency [GHz]')
    ax4.set_ylabel('Phase [rad]')
    ax4.set_xlim(min(kid['localsweep'].x), max(kid['localsweep'].x))
    ax4.legend(loc='best')
    ax4.grid()

    #SaveFits.py::Plot()から抜粋
    fc = kid['tod'].frequency #carrier frequency
    tck = interpolate.splrep(kid['localsweep'].x, phase_localsweep, s=0)
    yfc = interpolate.splev(fc, tck, der=0) # phase of carrier f
    yfr = interpolate.splev(fr, tck, der=0) # phase of resonance f

    tck = interpolate.splrep(phase_fitted_localsweep, kid['localsweep'].x, s=0)
    f = interpolate.splev(yfc, tck, der=0)
    linyfc = 4.0*Qr*(f - fr)/fr

    rewound_tod_amplitude = kid['rewound_tod'].amplitude
    rewound_tod_phase = kid['rewound_tod'].phase

    f = interpolate.splev(rewound_tod_phase, tck, der=0)
    linphase = 4.0*Qr*(f - fr)/fr

    ts = kid['tod_timestamp'][:-100]#データの最後の部分にはよくないデータが含まれる可能性があるので、最後の100個ぐらいを除外する(make_reduced_fitsで処理されている)。
    ts = ts - ts[0]

    fig2 = plt.figure(figsize=(8,10))
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312)
    ax3 = plt.subplot(313)

    ax1.plot(ts, rewound_tod_amplitude, 'r', label='ampl. (from cache)')
    ax1.plot(ts, linphase, 'b', label='linearized phase (from cache)')
    ax1.plot(ts, rewound_tod_phase, 'b', lw=3, alpha=0.3, label='phase (from cache)')
    ax1.plot(ts, np.ones(len(ts))*linyfc, 'c--', lw=1, label='carrier linphase (from cache)')
    ax1.plot(ts, np.ones(len(ts))*yfc, 'c--', lw=3, alpha=0.3, label='carrier phase (from cache)')
    #ax1.set_ylim(-2, 2)
    ax1.set_title('KID[{}]'.format(kid['kidid']))
    ax1.set_ylabel('Normalized Response')
    ax1.legend(loc='best')
    ax1.grid()

    ampl = kid['amp']
    phase = kid['ph']
    linphase = kid['linph']
    yfc = kid['yfc']
    linyfc = kid['linyfc']

    ax2.plot(ts, ampl, 'r', label='ampl. (deglitch from fits)')
    ax2.plot(ts, linphase, 'b', label='linearized phase (deglitch from fits)')
    ax2.plot(ts, phase, 'b', lw=3, alpha=0.3, label='phase (deglitch from fits)')
    ax2.plot(ts, np.ones(len(ts))*linyfc, 'c--', lw=1, label='carrier linphase (from fits)')
    ax2.plot(ts, np.ones(len(ts))*yfc, 'c--', lw=3, alpha=0.3, label='carrier phase (from fits)')
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel('Normalized Response')
    ax2.legend(loc='best')
    ax2.grid()

    # rebin = hdu['READOUT'].header['DSAMPLE']
    # dt = 1/hdu['READOUT'].header['FRAMERT'] * rebin
    rebin = 1
    dt = 1/kid['framert']*rebin
    plot_tod_ratio = 0.25
    size = 2**int(np.floor(np.log2(len(ts)*plot_tod_ratio)))
    f_, ampl_ = demerge.power_spectrum_density(ampl[:size], dt, 7, window=None, overwrap_half=True)
    f_, phase_ = demerge.power_spectrum_density(phase[:size], dt, 7, window=None, overwrap_half=True)
    f_, linphase_ = demerge.power_spectrum_density(linphase[:size], dt, 7, window=None, overwrap_half=True)
    ax3.semilogx(f_, np.log10(ampl_)*10.0, 'r', label='ampl. (from fits)')
    ax3.semilogx(f_, np.log10(linphase_)*10.0, 'b', label='linearized phase (from fits)')
    ax3.semilogx(f_, np.log10(phase_)*10.0, 'g', lw=3, alpha=0.3, label='phase (from fits)')
    ax3.set_xlabel('Frequency [Hz]')
    ax3.set_ylabel('PSD [dBk/Hz]')
    ax3.set_ylim(-110.,-50.)
    ax3.legend(loc='best')
    ax3.grid()

    return fig1, fig2


def main() -> None:
    """
    コマンドライン引数
    ------------------
    args[1] string フィットの結果が格納された一時ファイル(*.pkl)の名前
    args[2] string グラフを保存するディレクトリ名
    """
    basicConfig(
        datefmt='%Y-%m-%d %H:%M:%S',
        format='[%(asctime)s %(name)s %(levelname)s] %(message)s',
    )

    args = sys.argv
    filename = args[1]
    graph_dir = args[2]
    with open(filename, 'rb') as f:
        kid = pickle.load(f)
    fig1, fig2 = plot(kid)
    if not fig1 is None:
        plt.figure(fig1.number)
        plt.savefig('{}/sweep_kid{:05}.png'.format(graph_dir, kid['kidid']), transparent=False)
        plt.close()
    if not fig2 is None:
        plt.figure(fig2.number)
        plt.savefig('{}/tod_kid{:05}.png'.format(graph_dir, kid['kidid']), transparent=False)
        plt.close()


if __name__=='__main__':
    main()

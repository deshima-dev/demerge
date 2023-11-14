"""フィットの結果からphaseなどを計算してfitsに格納する

Filename: make_reduced_fits.py
(C) 2021 内藤システムズ

使い方
------
このスクリプトは第一引数与えられたディレクトリの中のkidXXXXX.pklファイルをすべて読み込む。
各pklファイルにはフィットの結果が格納されていることが正常に動作するための前提条件となる。
解析結果をまとめて第二引数に与えられるファイル名で保存する。

例) カレントディレクトリにあるcacheディレクトリの中身をhoge.fitsへまとめる。
  $ python make_reduced_fits.py ./cache hoge.fits
"""
import os
import sys
import glob
import shutil
import pickle
import numpy as np
from astropy.io import fits
from scipy import interpolate
from . import demerge

def make_reduced_fits(kidfiles, output_filename):
    fitinfo  = [] # [[kidid, pread, fc, yfc, linyfc, fr, dfr, Qr, dQr, Qc, dQc, Qi, dQi], ...]
    readout  = demerge.readout_dict()
    kidsinfo = demerge.kids_dict()
    readout['hdr_val_lis'][1]  = 'hoge'
    kidsinfo['hdr_val_lis'][1] = 'localsweep.sweep'
    downsampling_rate          = 1

    framert = None
    npoints = None
    for kidfile in kidfiles:
        with open(kidfile, 'rb') as f:
            kid = pickle.load(f)

        timestamp = kid['tod_timestamp'][:-100] #データの最後の部分にはよくないデータが含まれる可能性があるので、最後の100個ぐらいを除外する。
        ts        = demerge.rebin_array(timestamp, downsampling_rate) #1はdown samplingしないことを意味する
        fc        = kid['tod'].frequency

        if kid.get('resonance_params') is None:
            fr          = np.nan
            dfr         = np.nan
            Qr          = np.nan
            dQr         = np.nan
            Qc          = np.nan
            dQc         = np.nan
            Qi          = np.nan
            dQi         = np.nan
            yfc         = np.nan
            linyfc      = np.nan
            ampl        = np.array([np.nan for j in range(len(ts))])
            phase       = np.array([np.nan for j in range(len(ts))])
            linphase    = np.array([np.nan for j in range(len(ts))])
            rewound_tod = None
        else:
            fr  = kid['resonance_params']['fr'].value # GHz
            dfr = kid['resonance_params']['fr'].stderr # GHz
            Qr  = kid['resonance_params']['Qr'].value
            dQr = kid['resonance_params']['Qr'].stderr
            Qc  = kid['resonance_params']['Qc'].value
            dQc = kid['resonance_params']['Qc'].stderr
            Qi  = kid['resonance_params']['Qi'].value
            dQi = kid['resonance_params']['Qi'].stderr

            arga = kid['resonance_params']['arga']
            absa = kid['resonance_params']['absa']
            tau  = kid['resonance_params']['tau']
            phi0 = kid['resonance_params']['phi0']
            c = kid['resonance_params']['c']

            rewound_localsweep        = demerge.gaolinbg_rewind(kid['localsweep'].x, kid['localsweep'].iq, arga, absa, tau, fr, Qr, Qc, phi0, c)
            fitted_localsweep         = demerge.gaolinbg_function()(kid['localsweep'].x, arga, absa, tau, fr, Qr, Qc, phi0, c)
            rewound_fitted_localsweep = demerge.gaolinbg_rewind(kid['localsweep'].x, fitted_localsweep, arga, absa, tau, fr, Qr, Qc, phi0, c)

            phase_localsweep        = -np.angle(-rewound_localsweep)
            phase_fitted_localsweep = -np.angle(-rewound_fitted_localsweep)

            ## spline interpolation
            tck = interpolate.splrep(kid['localsweep'].x, phase_localsweep, s=0)
            yfc = interpolate.splev(fc, tck, der=0) # phase of carrier f
            #yfr = interpolate.splev(fr, tck, der=0) # phase of resonance f

            ## linphaseの計算(linyfc = k.convert_to_fshift(yfc, opt='linphase'))
            tck    = interpolate.splrep(phase_fitted_localsweep, kid['localsweep'].x, s=0)
            f      = interpolate.splev(yfc, tck, der=0)
            linyfc = 4.0*Qr*(f - fr)/fr

            # 抜粋:mkid_data.data.FixedData::calibrate_with_blind_tones()メソッド
            GHz = 1e9
            calibrated = demerge.calibrate_with_blind_tones(kid['tod'], kid['blind_tone_left'], kid['blind_tone_right'])
            #calibrated_tod = FixedData(kid['tod'].frequency*GHz, kid['tod'].t, np.real(calibrated), np.imag(calibrated), info=kid['tod'].info)
            # calibrated_tod = demerge.FixedData()
            # calibrated_tod.timestamp = kid['tod'].t
            # calibrated_tod.i = np.real(calibrated)
            # calibrated_tod.q = np.imag(calibrated)
            # calibrated_tod._frequency = kid['tod'].frequency*GHz
            # calibrated_tod.info = kid['tod'].info

            # 抜粋:mkid_data.kidfit.KidFitResult::rewound_ampl_phase()メソッド
            rw = demerge.gaolinbg_rewind(kid['tod'].frequency, calibrated, arga, absa, tau, fr, Qr, Qc, phi0, c)
            #rw = demerge.gaolinbg_rewind(calibrated_tod.frequency, calibrated_tod.iq, arga, absa, tau, fr, Qr, Qc, phi0, c)
            ampl      = 2*np.abs(rw)
            phase     = -np.angle(-rw)
            index_ar0 = np.where(abs(phase)>3.1)
            if len(index_ar0[0])>0:
                index0  = index_ar0[0][0]
                initial = phase[index0]
                if initial>0.:
                    index_ar1 = np.where(phase<0.)
                    factor    = 2.*np.pi
                else:
                    index_ar1 = np.where(phase>0.)
                    factor    = -2.*np.pi
                for index in index_ar1[0]:
                    if index>index0: phase[index] += factor

            #rewound_tod = FixedFitData(calibrated_tod.frequency, calibrated_tod.t, ampl, phase, info=calibrated_tod.info)
            #rewound_tod = FixedFitData(kid['tod'].frequency, kid['tod'].t, ampl, phase, info=kid['tod'].info)
            rewound_tod = demerge.FixedFitData(kid['tod'].frequency, kid['tod'].t, ampl, phase, info=kid['tod'].info)

            #ts, ampl, phase = k.get_cache('deglitch').unpack()
            baseline_thresh   = 6.0
            glitch_thresh     = 5.0
            clusterize_thresh = 2
            interp_offset     = 0
            bad = demerge.find_glitch([rewound_tod.phase], baseline_thresh, glitch_thresh, clusterize_thresh, interp_offset)

            # 抜粋:bbsweeplib.kids.KID::deglitch()
            ampl  = demerge.interpolate_bad(rewound_tod.amplitude, bad)
            phase = demerge.interpolate_bad(rewound_tod.phase, bad)
            d     = demerge.FixedFitData(rewound_tod.frequency*GHz, rewound_tod.t, ampl, phase)
            ts, ampl, phase = d.unpack()

            #tck = interpolate.splrep(phase_fitted_localsweep, kid['localsweep'].x, s=0)
            f        = interpolate.splev(phase, tck, der=0)
            linphase = 4.0*Qr*(f - fr)/fr
            ts       = demerge.rebin_array(ts, downsampling_rate)
            ampl     = demerge.rebin_array(ampl, downsampling_rate)
            phase    = demerge.rebin_array(phase, downsampling_rate)

        fitinfo.append([kid['kidid'], kid['readpower'], fc, yfc, linyfc, fr, dfr, Qr, dQr, Qc, dQc, Qi, dQi])

        if len(readout['cols_data_lis']) == 0:
            readout['cols_data_lis'].append(timestamp)
            readout['cols_data_lis'].append(np.ones(len(timestamp))*kid['pixelid'])

            framert = kid['framert']
            npoints = kid['npoints']

        readout['cols_key_lis'].append('Amp, Ph, linPh {}'.format(kid['kidid']))
        readout['cols_data_lis'].append(np.array([ampl, phase, linphase]).T)
        readout['tform'].append('3E')
        readout['tunit'].append(None)

        # plot.pyのtodのプロットで使うために保存
        kid['amp']         = ampl
        kid['ph']          = phase
        kid['linph']       = linphase
        kid['yfc']         = yfc
        kid['linyfc']      = linyfc
        kid['rewound_tod'] = rewound_tod
        with open(kidfile, 'wb') as f:
            pickle.dump(kid, f)

        readout['hdr_com_lis'].append('label for field {}'.format(kid['kidid'] + 2))
        readout['hdr_com_lis'].append('data format of field {}'.format(kid['kidid'] + 2))
        continue

    pixelid = 0
    #nkid = len(glob.glob('cache/kid*.pkl'))
    nkid = len(kidfiles)
    readout['hdr_val_lis'][2] = framert
    readout['hdr_val_lis'][3] = npoints
    readout['hdr_val_lis'][4] = downsampling_rate
    readout['hdr_key_lis'].append('NKID{}'.format(pixelid))
    readout['hdr_val_lis'].append(nkid) #全KID数
    readout['hdr_com_lis'].insert(5, 'number of KIDs (pixel {})'.format(pixelid))

    kidsinfo['hdr_key_lis'].append('NKID{}'.format(pixelid))
    kidsinfo['hdr_val_lis'].append(nkid)
    kidsinfo['hdr_com_lis'].insert(2, 'number of KIDs (pixel {})'.format(pixelid))

    kidid_, pread_, fc_, yfc_, linyfc_, fr_, dfr_, Qr_, dQr_, Qc_, dQc_, Qi_, dQi_ = [np.array(fitinfo)[:,i] for i in range(13)]
    kidsinfo['cols_data_lis'].append(np.ones(nkid)*pixelid)
    kidsinfo['cols_data_lis'].append(kidid_)
    kidsinfo['cols_data_lis'].append(pread_)
    kidsinfo['cols_data_lis'].append(fc_)
    kidsinfo['cols_data_lis'].append(np.array([yfc_, linyfc_]).T)
    kidsinfo['cols_data_lis'].append(np.array([fr_, dfr_]).T)
    kidsinfo['cols_data_lis'].append(np.array([Qr_, dQr_]).T)
    kidsinfo['cols_data_lis'].append(np.array([Qc_, dQc_]).T)
    kidsinfo['cols_data_lis'].append(np.array([Qi_, dQi_]).T)

    hdus = fits.HDUList()
    hdus.append(fits.PrimaryHDU())
    hdus.append(demerge.createBinTableHDU(kidsinfo))
    hdus.writeto(output_filename)
    hdus.close()

    readout_hdu = demerge.createBinTableHDU(readout)
    fits.append(output_filename, readout_hdu.data, readout_hdu.header, memmap=True)
    return

def main() -> None:
    """
    コマンドライン引数
    ------------------
    args[1] string pklファイルが格納されているcacheディレクトリへのパス
    args[2] string 出力ファイル名
    """
    args            = sys.argv
    cache_dir       = args[1]
    output_filename = args[2]

    kidfiles = sorted(glob.glob(cache_dir + '/kid*.pkl')) #ファイル名で昇順で並べ替える

    make_reduced_fits(kidfiles, output_filename)

if __name__=='__main__':
    main()

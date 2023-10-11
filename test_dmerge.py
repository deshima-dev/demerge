"""DESHIMA解析プログラムのテスト

File name: test_dmerge.py
Python 3.7
(C) 2021 内藤システムズ
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

script_dir = os.path.dirname(os.path.abspath(__file__))
libpath = os.path.join(os.path.dirname(script_dir), "libs")
sys.path.append(libpath)

import dmerge
import dmerge_deprecated
import mkid_data
from bbsweeplib.files import read_kidslist, read_localsweep, read_fits
from bbsweeplib.kids import fits_tods_fits
from mkid_data.data import FixedData, FixedFitData


class DmergeTestDrive(unittest.TestCase):
    """dmerge.pyモジュールの単体テスト"""

    def setUp(self):
        return

    def load_data(self):
        """SweepデータやTODを読み込む"""
        kidslist = dmerge.load_kidslist(
            "../raw_dataset/obs/cosmos_20171103184836/kids.list"
        )
        localsweeps = dmerge.load_localsweep(
            "../raw_dataset/obs/cosmos_20171103184836/localsweep.sweep",
            framelen=kidslist[0]["framelen"],
        )
        tods = dmerge.TODs(
            "../raw_dataset/obs/cosmos_20171103184836/20171103184836.fits"
        )
        return kidslist, localsweeps, tods

    def test_SweepData(self):
        """SweepDataクラスのテスト"""
        sweep = dmerge_deprecated._read_sweep(
            "../raw_dataset/obs/cosmos_20171103184836/localsweep.sweep"
        )
        lofreqs, bins, data = sweep
        framelen = 19
        dfreq = 2e9 / (2**framelen)
        GHz = 1e9
        for b, d in zip(bins, data.T):
            freq = lofreqs + b * dfreq
            sd = dmerge.SweepData(np.real(d), np.imag(d), freq)
            amplitude = np.sqrt(np.real(d) ** 2 + np.imag(d) ** 2)
            iq = d
            deg = np.arctan2(np.imag(d), np.real(d)) * 180 / np.pi
            self.assertTrue((sd.x == freq / GHz).all())
            self.assertTrue((sd.amplitude == amplitude).all())
            self.assertTrue((sd.iq == iq).all())
            self.assertTrue((sd.deg == deg).all())
            self.assertEqual(len(sd), len(freq))
            break
        return

    def test_FixedData(self):
        """FixedDataクラスのテスト"""
        if os.path.exists("cache"):
            shutil.rmtree("cache")
        cmd = " ".join(
            [
                "python make_divided_data.py",
                "../raw_dataset/obs/cosmos_20171103184836/kids.list",
                "../raw_dataset/obs/cosmos_20171103184836/localsweep.sweep",
                "../raw_dataset/obs/cosmos_20171103184836/20171103184836.fits",
                "cache",
            ]
        )
        os.system(cmd)
        cmd = " ".join(["python", "calc_resonance_params.py", "cache/kid00000.pkl"])
        os.system(cmd)
        with open("cache/kid00000.pkl", "rb") as f:
            kid = pickle.load(f)
        GHz = 1e9
        calibrated = dmerge.calibrate_with_blind_tones(
            kid["tod"], kid["blind_tone_left"], kid["blind_tone_right"]
        )
        tods = dmerge.TODs(
            "../raw_dataset/obs/cosmos_20171103184836/20171103184836.fits"
        )
        index = tods.bins.index(-201203)
        chunksize = int(1e4)
        rawdata = tods.hud[1].data["data"]
        read_I_ = rawdata[: tods.offset, 2 * index]
        read_Q_ = rawdata[: tods.offset, 2 * index + 1]
        info = dict()
        info["bins"] = tods.bins
        info["freqs"] = tods.carrier_freq[index]
        info["header"] = tods.hud[1].header
        expected = FixedData(
            "I-Q-Gain",
            (tods.timestamp[: tods.offset],),
            (read_I_, read_Q_, tods.fftgain),
            tods.carrier_freq[index],
            info=info,
        )
        # expected = FixedData(kid['tod'].frequency*GHz, kid['tod'].t, np.real(calibrated), np.imag(calibrated), info=kid['tod'].info)

        result = dmerge.FixedData()
        result.timestamp = kid["tod"].t
        result.fftgain = tods.fftgain
        result._i = read_I_
        result._q = read_Q_
        result._frequency = tods.carrier_freq[index]
        result.info = kid["tod"].info

        self.assertTrue((result.i == expected.i).all())
        self.assertTrue((result.q == expected.q).all())
        self.assertTrue((result.t == expected.t).all())
        self.assertTrue((result.iq == expected.iq).all())
        self.assertTrue((result.frequency == expected.frequency).all())
        return

    def test_FixedFitData(self):
        """単純にしたFixedFitDataのテスト"""
        if os.path.exists("cache"):
            shutil.rmtree("cache")
        cmd = " ".join(
            [
                "python make_divided_data.py",
                "../raw_dataset/obs/cosmos_20171103184836/kids.list",
                "../raw_dataset/obs/cosmos_20171103184836/localsweep.sweep",
                "../raw_dataset/obs/cosmos_20171103184836/20171103184836.fits",
                "cache",
            ]
        )
        os.system(cmd)
        cmd = " ".join(["python", "calc_resonance_params.py", "cache/kid00000.pkl"])
        os.system(cmd)
        with open("cache/kid00000.pkl", "rb") as f:
            kid = pickle.load(f)
        GHz = 1e9
        fr = kid["resonance_params"]["fr"].value  # GHz
        Qr = kid["resonance_params"]["Qr"].value
        Qc = kid["resonance_params"]["Qc"].value
        arga = kid["resonance_params"]["arga"]
        absa = kid["resonance_params"]["absa"]
        tau = kid["resonance_params"]["tau"]
        phi0 = kid["resonance_params"]["phi0"]
        c = kid["resonance_params"]["c"]
        calibrated = dmerge.calibrate_with_blind_tones(
            kid["tod"], kid["blind_tone_left"], kid["blind_tone_right"]
        )
        rw = dmerge.gaolinbg_rewind(
            kid["tod"].frequency, calibrated, arga, absa, tau, fr, Qr, Qc, phi0, c
        )
        ampl = 2 * np.abs(rw)
        phase = -np.angle(-rw)
        index_ar0 = np.where(abs(phase) > 3.1)
        if len(index_ar0[0]) > 0:
            index0 = index_ar0[0][0]
            initial = phase[index0]
            if initial > 0.0:
                index_ar1 = np.where(phase < 0.0)
                factor = 2.0 * np.pi
            else:
                index_ar1 = np.where(phase > 0.0)
                factor = -2.0 * np.pi
            for index in index_ar1[0]:
                if index > index0:
                    phase[index] += factor
        expected = FixedFitData(
            kid["tod"].frequency, kid["tod"].t, ampl, phase, info=kid["tod"].info
        )
        result = dmerge.FixedFitData(
            kid["tod"].frequency, kid["tod"].t, ampl, phase, kid["tod"].info
        )
        self.assertTrue((result.amplitude == expected.amplitude).all())
        self.assertTrue((result.phase == expected.phase).all())
        self.assertTrue((result.t == expected.t).all())
        self.assertEqual(result.info, expected.info)
        self.assertTrue((result.frequency == expected.frequency).all())
        result_timestamp, result_ampl, result_phase = result.unpack()
        expected_timestamp, expected_ampl, expected_phase = expected.unpack()
        self.assertTrue((result_timestamp == expected_timestamp).all())
        self.assertTrue((result_ampl == expected_ampl).all())
        self.assertTrue((result_phase == expected_phase).all())
        return

    def test_Expr_with_args(self):
        """Expr_with_argsの仕組みを調べる"""
        from collections import namedtuple

        class Expr_with_args(
            namedtuple("Expr_with_args", "expr, arg_symbols, arg_names")
        ):
            def __call__(self, x, *args, **kws):
                return expr, arg_symbols, arg_names

        a, b, c = Expr_with_args("expr", "symbols", "names")
        self.assertEqual(a, "expr")
        self.assertEqual(b, "symbols")
        self.assertEqual(c, "names")
        return

    def test_divide_data(self):
        """KID毎にSweepデータやTODを分割して辞書に格納する機能の試験
        機能
        ----
        divide_data()関数はsweep_data, kidlist, todを受け取り、
        それらをKID毎に分割して辞書に格納する。
        その辞書を適当に分割してファイルに保存すれば、
        並列処理に利用できる。
        """
        kidslist, localsweeps, tods = self.load_data()
        result = dmerge.divide_data(kidslist, localsweeps, tods)
        tods.close()
        self.assertEqual(len(result), 63)
        nbin = -196483
        self.assertEqual(result[0]["bin"], nbin)
        self.assertEqual(type(result[0]["localsweep"]), dmerge.SweepData)
        self.assertEqual(type(result[0]["localsweep"].x), np.ndarray)
        self.assertEqual(result[0]["readpower"], -75)
        # self.assertEqual(result[0]['blind_tone_left'], -201203)
        # self.assertEqual(result[0]['blind_tone_right'], -192949)
        self.assertEqual(type(result[0]["blind_tone_left"]), dmerge.FixedData)
        self.assertEqual(type(result[0]["blind_tone_right"]), dmerge.FixedData)
        self.assertEqual(type(result[0]["tod"]), dmerge.FixedData)
        return

    def test_make_divided_data(self):
        """make_divided_data.pyスクリプトの試験"""
        if os.path.exists("cache"):
            shutil.rmtree("cache")
        cmd = " ".join(
            [
                "python make_divided_data.py",
                "../raw_dataset/obs/cosmos_20171103184836/kids.list",
                "../raw_dataset/obs/cosmos_20171103184836/localsweep.sweep",
                "../raw_dataset/obs/cosmos_20171103184836/20171103184836.fits",
                "cache",
            ]
        )
        os.system(cmd)
        self.assertEqual(len(glob.glob("cache/kid*.pkl")), 63)
        return

    def test_kids_with_both_blinds(self):
        """kids_with_both_blinds()関数のテストと単機能の関数のテスト"""
        kidslist = dmerge.load_kidslist(
            "../raw_dataset/obs/cosmos_20171103184836/kids.list"
        )
        expected = dmerge_deprecated.kids_with_both_blinds(
            kidslist, allow_without_blind=True
        )
        result = dmerge.find_nearest_blind_tones(kidslist)
        self.assertTrue((result == expected).all())
        return

    def test_read_kidslist(self):
        """新しいread_kidslist()関数をテストする"""
        expected = dmerge_deprecated.read_kidslist(
            "../raw_dataset/obs/cosmos_20171103184836/kids.list"
        )
        info, kids, blinds, powers = expected
        result = dmerge.load_kidslist(
            "../raw_dataset/obs/cosmos_20171103184836/kids.list"
        )
        r_info, r_kids, r_blinds, r_powers = result
        self.assertTrue(info == r_info)
        self.assertTrue(powers == r_powers)
        self.assertTrue((kids == r_kids).all())
        self.assertTrue((blinds == r_blinds).all())
        return

    def test_load_localsweep(self):
        """ローカルsweepデータをロードするテスト(ロード結果が古い関数と同じであることを確認する)"""
        info, kids, blinds, powers = dmerge.load_kidslist(
            "../raw_dataset/obs/cosmos_20171103184836/kids.list"
        )
        filename = "../raw_dataset/obs/cosmos_20171103184836/localsweep.sweep"
        framelen = info["framelen"]
        expected = dmerge_deprecated.read_localsweep(filename, framelen=framelen)
        result = dmerge.load_localsweep(filename, framelen)
        self.assertEqual(result.keys(), expected.keys())
        for b in result.keys():
            with self.subTest(i=b):
                self.assertTrue((result[b].x == expected[b].x).all())
                self.assertTrue((result[b].amplitude == expected[b].amplitude).all())
                self.assertTrue((result[b].iq == expected[b].iq).all())
                self.assertTrue((result[b].deg == expected[b].deg).all())
                self.assertEqual(len(result[b]), len(expected[b]))
        return

    def test_fit_sweep(self):
        """1つのKIDのSweepデータをフィットする試験"""
        kidslist, localsweeps, tods = self.load_data()
        kids = dmerge.divide_data(kidslist, localsweeps, tods)
        tods.close()
        # kidid, fr, dfr, Qr, dQr, Qc, dQc, Qi, dQi
        expected = np.load(
            "workspace/analysis/ana_data/cosmos_20171103184836_01011_ana/FitSweep_fit.npy"
        )
        offset = 0  # フィットに失敗した場合そのKID IDは欠番になる。その番号を調整するためのオフセット。
        for kidid in range(len(kids)):
            with self.subTest(i=kidid):
                fc = kids[kidid]["tod"].frequency
                Q_search = 100
                peaks = dmerge.find_peaks(
                    kids[kidid]["localsweep"].x,
                    kids[kidid]["localsweep"].amplitude,
                    fc,
                    minq=Q_search,
                    maxratio=1,
                )
                try:
                    result, fitrange = dmerge.fit_onepeak(
                        kids[kidid]["localsweep"], peaks, nfwhm=5
                    )
                except Exception:
                    offset += 1
                    continue
                fr = result.params["fr"].value  # GHz
                dfr = result.params["fr"].stderr  # GHz
                Qr = result.params["Qr"].value
                dQr = result.params["Qr"].stderr
                Qc = result.params["Qc"].value
                dQc = result.params["Qc"].stderr
                Qi = result.params["Qi"].value
                dQi = result.params["Qi"].stderr

                if dfr / fr < 0.0 or dQr / Qr < 0.0 or dQc / Qc < 0.0 or dQi / Qi < 0.0:
                    offset += 1
                    continue
                elif fr != fr or dfr != dfr or Qr != Qr or dQr != dQr:  # reject Nan
                    offset += 1
                    continue
                elif abs(Qr) == float("inf") or abs(dQr) == float(
                    "inf"
                ):  # reject +/-inf
                    offset += 1
                    continue

                kidid = kidid - offset
                self.assertEqual(result.params["fr"].value, expected[1][kidid])
                self.assertEqual(result.params["fr"].stderr, expected[2][kidid])
                self.assertEqual(result.params["Qr"].value, expected[3][kidid])
                self.assertEqual(result.params["Qr"].stderr, expected[4][kidid])
                self.assertEqual(result.params["Qc"].value, expected[5][kidid])
                self.assertEqual(result.params["Qc"].stderr, expected[6][kidid])
                self.assertEqual(result.params["Qi"].value, expected[7][kidid])
        return

    def test_calc_resonance_params(self):
        """calc_resonance_params.pyスクリプトの試験"""
        cmd = " ".join(
            [
                "python make_divided_data.py",
                "../raw_dataset/obs/cosmos_20171103184836/kids.list",
                "../raw_dataset/obs/cosmos_20171103184836/localsweep.sweep",
                "../raw_dataset/obs/cosmos_20171103184836/20171103184836.fits",
                "cache",
            ]
        )
        os.system(cmd)
        cmd = " ".join(["python", "calc_resonance_params.py", "cache/kid00001.pkl"])
        os.system(cmd)
        with open("cache/kid00001.pkl", "rb") as f:
            result = pickle.load(f)
        self.assertEqual(type(result["resonance_params"]), lmfit.Parameters)
        self.assertTrue("fitrange" in result)
        return

    def test_parallel_calc_resonance_params(self):
        """calc_resonance_params.pyスクリプトを用いた並列処理の試験"""
        cmd = " ".join(
            [
                "python make_divided_data.py",
                "../raw_dataset/obs/cosmos_20171103184836/kids.list",
                "../raw_dataset/obs/cosmos_20171103184836/localsweep.sweep",
                "../raw_dataset/obs/cosmos_20171103184836/20171103184836.fits",
                "cache",
            ]
        )
        os.system(cmd)
        cmd = " ".join(
            ["ls cache/*.pkl | xargs -P7 -n1", "python", "calc_resonance_params.py"]
        )
        os.system(cmd)
        for filename in os.listdir("cache"):
            with open("cache/" + filename, "rb") as f:
                result = pickle.load(f)
                with self.subTest(i=filename):
                    if result["enabled"]:
                        self.assertEqual(
                            type(result["resonance_params"]), lmfit.Parameters
                        )
                    else:
                        self.assertEqual(result["resonance_params"], None)
        return

    def test_rewind(self):
        """rewind()関数の試験
        方法
        ----
        rewind()関数を用いてキャリア周波数と共振周波数の位相を計算する。
        reduced_*.fitsに保存されている位相と比べる。
        """
        if not os.path.exists("cache"):
            cmd = " ".join(
                [
                    "python make_divided_data.py",
                    "../raw_dataset/obs/cosmos_20171103184836/kids.list",
                    "../raw_dataset/obs/cosmos_20171103184836/localsweep.sweep",
                    "../raw_dataset/obs/cosmos_20171103184836/20171103184836.fits",
                    "cache",
                ]
            )
            os.system(cmd)
        for kidid in range(len(glob.glob("cache/kid*.pkl"))):
            with open("cache/kid{:05}.pkl".format(kidid), "rb") as f:
                kid = pickle.load(f)
            fc = kid["tod"].frequency
            peaks = dmerge.find_peaks(
                kid["localsweep"].x,
                kid["localsweep"].amplitude,
                fc,
                minq=100,
                maxratio=1,
            )
            kidfitresult, fitrange = dmerge.fit_onepeak(
                kid["localsweep"], peaks, nfwhm=5
            )

            fr = kidfitresult.params["fr"].value  # GHz
            dfr = kidfitresult.params["fr"].stderr  # GHz
            Qr = kidfitresult.params["Qr"].value
            dQr = kidfitresult.params["Qr"].stderr
            Qc = kidfitresult.params["Qc"].value
            dQc = kidfitresult.params["Qc"].stderr
            Qi = kidfitresult.params["Qi"].value
            dQi = kidfitresult.params["Qi"].stderr

            if dfr / fr < 0.0 or dQr / Qr < 0.0 or dQc / Qc < 0.0 or dQi / Qi < 0.0:
                continue
            elif fr != fr or dfr != dfr or Qr != Qr or dQr != dQr:  # reject Nan
                continue
            elif abs(Qr) == float("inf") or abs(dQr) == float("inf"):  # reject +/-inf
                continue

            arga = kidfitresult.params["arga"]
            absa = kidfitresult.params["absa"]
            tau = kidfitresult.params["tau"]
            fr = kidfitresult.params["fr"]
            Qr = kidfitresult.params["Qr"]
            Qc = kidfitresult.params["Qc"]
            phi0 = kidfitresult.params["phi0"]
            c = kidfitresult.params["c"]

            rewound_localsweep = dmerge.gaolinbg_rewind(
                kid["localsweep"].x,
                kid["localsweep"].iq,
                arga,
                absa,
                tau,
                fr,
                Qr,
                Qc,
                phi0,
                c,
            )
            fitted_localsweep = dmerge.gaolinbg_function()(
                kid["localsweep"].x, arga, absa, tau, fr, Qr, Qc, phi0, c
            )
            rewound_fitted_localsweep = dmerge.gaolinbg_rewind(
                kid["localsweep"].x,
                fitted_localsweep,
                arga,
                absa,
                tau,
                fr,
                Qr,
                Qc,
                phi0,
                c,
            )

            phase_localsweep = -np.angle(-rewound_localsweep)
            phase_fitted_localsweep = -np.angle(-rewound_fitted_localsweep)

            ## spline interpolation
            tck = interpolate.splrep(kid["localsweep"].x, phase_localsweep, s=0)
            yfc = interpolate.splev(fc, tck, der=0)  # phase of carrier f
            # yfr = interpolate.splev(fr, tck, der=0) # phase of resonance f

            ## linphaseの計算(linyfc = k.convert_to_fshift(yfc, opt='linphase'))
            tck = interpolate.splrep(phase_fitted_localsweep, kid["localsweep"].x, s=0)
            f = interpolate.splev(yfc, tck, der=0)
            linyfc = 4.0 * Qr * (f - fr) / fr

            filename = (
                "../raw_dataset/obs/cosmos_20171103184836/reduced_20171103184836.fits"
            )
            error_limit = 0.0001  # 0.01%
            with fits.open(filename, memmap=False) as hdul:
                with self.subTest(i=kidid):
                    expected = hdul["KIDSINFO"].data["yfc, linyfc"][kidid][0]
                    if np.isnan(expected):
                        if kid["enabled"]:
                            self.assertTrue(np.isnan(yfc))
                        else:
                            self.assertTrue(True)
                    else:
                        error = np.abs(yfc - expected) / expected
                        self.assertLess(error, error_limit, "yfc")

                    expected = hdul["KIDSINFO"].data["yfc, linyfc"][kidid][1]
                    if np.isnan(expected):
                        self.assertTrue(np.isnan(linyfc))
                    else:
                        error = np.abs(linyfc - expected) / expected
                        self.assertLess(error, error_limit, "linyfc")
        return

    def test_make_reduced_fits(self):
        """make_reduced_fits.pyスクリプトの試験
        方法
        ----
        make_reduced_fits.pyで生成されるreduced.fitsファイルと旧プログラムで生成されるものと中身を比較する。
        """
        kidfile = "cache/kid00000.pkl"
        if os.path.exists(kidfile):
            with open(kidfile, "rb") as f:
                kid = pickle.load(f)
            if not "resonance_params" in kid:
                shutil.rmtree("cache")
        if not os.path.exists("cache/kid00000.pkl"):
            cmd = " ".join(
                [
                    "python make_divided_data.py",
                    "../raw_dataset/obs/cosmos_20171103184836/kids.list",
                    "../raw_dataset/obs/cosmos_20171103184836/localsweep.sweep",
                    "../raw_dataset/obs/cosmos_20171103184836/20171103184836.fits",
                    "cache",
                ]
            )
            os.system(cmd)
            cmd = " ".join(
                ["ls cache/*.pkl | xargs -P7 -n1", "python", "calc_resonance_params.py"]
            )
            os.system(cmd)
        if os.path.exists("cache/reduced.fits"):
            os.remove("cache/reduced.fits")

        cmd = " ".join(["python make_reduced_fits.py", "cache", "cache/reduced.fits"])
        os.system(cmd)
        self.assertTrue(os.path.exists("cache/reduced.fits"))
        hdul_expected = fits.open(
            "../raw_dataset/obs/cosmos_20171103184836/reduced_20171103184836.fits",
            memmap=False,
        )
        hdul_result = fits.open("cache/reduced.fits", memmap=False)

        extname = "KIDSINFO"
        fields = [
            "pixelid",
            "kidid",
            "Pread",
            "fc",
            "yfc, linyfc",
            "fr, dfr (300K)",
            "Qr, dQr (300K)",
            "Qc, dQc (300K)",
            "Qi, dQi (300K)",
        ]
        for field in fields:
            j = 0
            for expected in hdul_expected[extname].data[field]:
                with self.subTest(i=field, j=j):
                    if type(expected) == np.ndarray:
                        if np.isnan(expected[0]):
                            self.assertTrue(
                                np.isnan(hdul_result[extname].data[field][j][0])
                            )
                        else:
                            self.assertEqual(
                                hdul_result[extname].data[field][j][0], expected[0]
                            )

                        if np.isnan(expected[1]):
                            self.assertTrue(
                                np.isnan(hdul_result[extname].data[field][j][1])
                            )
                        else:
                            self.assertEqual(
                                hdul_result[extname].data[field][j][1], expected[1]
                            )

                    elif np.isnan(expected):
                        self.assertTrue(np.isnan(hdul_result[extname].data[field][j]))
                    else:
                        self.assertEqual(hdul_result[extname].data[field][j], expected)
                j += 1

        extname = "READOUT"  # 注意!!!最後に空白が1文字入っている!!!
        fields = ["timestamp", "pixelid"]
        for field in fields:
            j = 0
            for expected in hdul_expected[extname].data[field]:
                with self.subTest(i=field, j=j):
                    if type(expected) == np.ndarray:
                        if np.isnan(expected[0]):
                            self.assertTrue(
                                np.isnan(hdul_result[extname].data[field][j][0])
                            )
                        else:
                            self.assertEqual(
                                hdul_result[extname].data[field][j][0], expected[0]
                            )

                        if np.isnan(expected[1]):
                            self.assertTrue(
                                np.isnan(hdul_result[extname].data[field][j][1])
                            )
                        else:
                            self.assertEqual(
                                hdul_result[extname].data[field][j][1], expected[1]
                            )

                    elif np.isnan(expected):
                        self.assertTrue(np.isnan(hdul_result[extname].data[field][j]))
                    else:
                        self.assertEqual(hdul_result[extname].data[field][j], expected)
                j += 1

        extname = "READOUT"  # 注意!!!最後に空白が1文字入っている!!!
        nmax = 63
        fields = ["Amp, Ph, linPh {}".format(i) for i in range(nmax)]
        for field in fields:
            j = 0
            for expected in hdul_expected[extname].data[field]:
                with self.subTest(i=field, j=j):
                    if np.isnan(expected[0]):
                        self.assertTrue(
                            np.isnan(hdul_result[extname].data[field][j][0]),
                            "Amp is nan",
                        )
                        self.assertTrue(
                            np.isnan(hdul_result[extname].data[field][j][1]),
                            "Ph is nan",
                        )
                        self.assertTrue(
                            np.isnan(hdul_result[extname].data[field][j][2]),
                            "linPh is nan",
                        )
                    else:
                        self.assertEqual(
                            hdul_result[extname].data[field][j][0], expected[0], "Amp"
                        )
                        self.assertEqual(
                            hdul_result[extname].data[field][j][1], expected[1], "Ph"
                        )
                        self.assertEqual(
                            hdul_result[extname].data[field][j][2], expected[2], "linPh"
                        )
                j += 1

        hdul_result.close()
        hdul_expected.close()
        return

    def test_mergefits(self):
        """MergeFitsのテスト
        TODO:
        (v) とりあえず、reduced.fitsをそのままマージできるかどうか、試してみる。
        (v) マージの結果が旧プログラムと同じことを確認する。
        """
        # reducedfitsfile = 'cache/reduced_20171103184836.fits'
        reducedfitsfile = "cache/reduced.fits"
        if not os.path.exists(reducedfitsfile):
            if os.path.exists("cache"):
                shutil.rmtree("cache")
            cmd = " ".join(
                [
                    "python make_divided_data.py",
                    "../raw_dataset/obs/cosmos_20171103184836/kids.list",
                    "../raw_dataset/obs/cosmos_20171103184836/localsweep.sweep",
                    "../raw_dataset/obs/cosmos_20171103184836/20171103184836.fits",
                    "cache",
                ]
            )
            os.system(cmd)
            cmd = " ".join(
                ["ls cache/*.pkl | xargs -P7 -n1", "python", "calc_resonance_params.py"]
            )
            os.system(cmd)
            cmd = " ".join(
                ["python make_reduced_fits.py", "cache", "cache/reduced.fits"]
            )
            os.system(cmd)
        dfitsfile = "cache/dfits_20171103184836.fits.gz"
        if os.path.exists(dfitsfile):
            os.remove(dfitsfile)

        from merge_to_dfits import MergeToDfits

        ddb_fits = "../caldb/DDB_20180619.fits.gz"
        dfitsdict = "dfits_dict.yaml"
        obsinst = "../raw_dataset/obs/cosmos_20171103184836/20171103184836.obs"
        antennalog = "../raw_dataset/obs/cosmos_20171103184836/20171103184836.ant"
        weatherlog = "../raw_dataset/obs/cosmos_20171103184836/20171103184836.wea"
        cabinlog = "cabin.db"

        mtd = MergeToDfits(
            ddbfits=ddb_fits,
            dfitsdict=dfitsdict,
            obsinst=obsinst,
            antennalog=antennalog,
            rout_data=reducedfitsfile,
            weatherlog=weatherlog,
            cabinlog=cabinlog,
        )
        dfits_hdus = mtd.dfits
        dfits_hdus.writeto(dfitsfile)
        mtd.kidsinfo_hdus.close()
        self.assertTrue(os.path.exists(dfitsfile))
        return

    def test_plot_sweep(self):
        """plot_sweep()関数を用いて1つのKIDについてプロットを作成するテスト"""
        os.system("rm -rf tmp_graph")
        os.system("mkdir tmp_graph")
        os.system("python plot_sweep.py cache/kid00001.pkl tmp_graph")
        self.assertTrue(os.path.exists("graphs"))
        self.assertTrue(os.path.exists("graphs/sweep_kid00001.png"))
        self.assertTrue(os.path.exists("graphs/tod_kid00001.png"))
        return

    def test_parallel_plot_sweep(self):
        """plot_sweep.pyスクリプトを用いた並列処理の試験"""
        os.system("rm -rf cache tmp_graph")
        os.system("mkdir cache tmp_graph")
        cmd = " ".join(
            [
                "python make_divided_data.py",
                "../raw_dataset/obs/cosmos_20171103184836/kids.list",
                "../raw_dataset/obs/cosmos_20171103184836/localsweep.sweep",
                "../raw_dataset/obs/cosmos_20171103184836/20171103184836.fits",
                "cache",
            ]
        )
        os.system(cmd)
        cmd = " ".join(
            ["ls cache/*.pkl | xargs -P7 -n1", "python", "calc_resonance_params.py"]
        )
        os.system(cmd)
        cmd = " ".join(["python make_reduced_fits.py", "cache", "cache/reduced.fits"])
        os.system(cmd)
        s = ""
        for filename in glob.glob("cache/kid*.pkl"):
            s += filename + " tmp_graph "
        cmd = " ".join(['echo "', s, '"', "| xargs -P7 -n2", "python", "plot_sweep.py"])
        os.system(cmd)
        self.assertEqual(
            len(os.listdir("graphs")), 2 * (63 - 1)
        )  # 63個のKIDの内46番だけフィットできていないのでプロットされない。
        return


if __name__ == "__main__":
    unittest.main()

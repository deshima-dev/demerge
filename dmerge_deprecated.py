"""dmerge.pyの内、新しい関数やクラスに置き換えられ、テストのために消さずにおきたいものをこのファイルに記す。

(C) 2021 内藤システムズ
"""
import numpy as np


def kids_with_both_blinds(kidslist, tod=None, allow_without_blind=False):
    """
    ** bbsweeplib.TODcalibからコピー **
    search kids its nearest blind tones.

    kidslist :: a tuple of (info, kids, blinds, powers),
                that read_kidslist returns.
    tod      :: None or a OrderedDict.
                if not None, only search KID in `tod`.

    Return a list of list of form
      [[signal0, left0, right0],
       [signal1, left1, right1],
         :
         :
       [signalp, leftp, rightp]],
      where
        signalN : carrier bin number for KID,
        leftN   : carrier bin number for nearest left blind tone,
        rightN  : carrier bin number for nearest right blind tone,
      respectively.
      this list is sorted by bin number in lowest-first order.
    """
    info, kids, blinds, powers = kidslist
    if tod:
        todkeys = tod.keys()
        kids = np.intersect1d(kids, todkeys)

    allkeys = sorted(list(kids) + list(blinds))
    results = []
    for k in sorted(kids):
        pair = []
        pos = allkeys.index(k) - 1
        while pos >= 0:
            key = allkeys[pos]
            if ((key in blinds) and not tod) or (
                (key in blinds) and tod and (key in todkeys)
            ):
                pair.append(allkeys[pos])
                break
            pos = pos - 1
        else:
            if allow_without_blind:
                pair.append(None)
            else:
                continue
        pos = allkeys.index(k) + 1
        while pos < len(allkeys):
            key = allkeys[pos]
            if ((key in blinds) and not tod) or (
                (key in blinds) and tod and (key in todkeys)
            ):
                pair.append(key)
                break
            pos = pos + 1
        else:
            if allow_without_blind:
                pair.append(None)
            else:
                continue
        results.append([k] + pair)
    return np.array(results)


def read_kidslist(fname):
    """Read kidslist file. ** bbsweeplib.files.py **

    Return (info, kids, blinds, powers)
      where info   is a dict with keys 'LO' (LO freq in Hz) and
                      framelen (frame count);
            kids   is a 1D list of bin indices for KIDs;
            blinds is a 1D list of bin indices for blind tones;
            powers is a OrderedDict, whose key is bin indices
                      and value is power for that bin.
    """
    mode = ["header"]
    kids = []
    blinds = []
    powers = dict()

    def checkmode(line):
        if line[:5] == "#KIDs":
            mode[0] = "kids"
            return True
        elif line[:7] == "#blinds":
            mode[0] = "blinds"
            return True
        else:
            return False

    with open(fname) as f:
        info = dict()
        for l in f:
            l = l.strip()
            if not l:
                continue
            elif checkmode(l):
                continue
            if mode[0] == "header":
                key, val = l[1:].split(":")
                if key == "LO":
                    info["LO"] = float(val) * 1e6
                elif key == "framelen":
                    info["framelen"] = int(val)
                else:
                    raise RuntimeError("error")
            else:
                values = l.split()
                nbin = int(values[0])
                power = float(values[1])
                if mode[0] == "kids":
                    kids.append(nbin)
                    powers[nbin] = power
                else:
                    blinds.append(nbin)
                    powers[nbin] = power
    return info, np.array(kids), np.array(blinds), powers


def read_localsweep(sweepfname, framelen):
    """Read local sweep file. ** bbsweeplib.files.py **
    Returns a OrderedDict of md.SweepData's, keyed by bin number.
    It is sorted by bin number, in lowest-first order.
    """
    sweep = _read_sweep(sweepfname)
    from collections import OrderedDict

    lofreqs, bins, data = sweep
    sweeps = OrderedDict()
    dfreq = 2e9 / (2**framelen)
    from dmerge import SweepData

    for b, d in zip(bins, data.T):
        swpdata = SweepData(np.real(d), np.imag(d), lofreqs + b * dfreq)
        sweeps[b] = swpdata

    result = OrderedDict(sorted(sweeps.items(), key=lambda x: x[0]))  # sort by bin
    return result


def _read_sweep(fname):
    """Read local sweep file. ** bbsweeplib.files.py **
    Return (lofreqs, bins, data)
      where lofreqs is an 1D array of frequency in Hz;
            bins is an 1D array of bin indices;
            data is an 2D array of complex demodulated amplitudes,
            one row per LO frequency, one column per bin.
    """
    rawdata = np.loadtxt(fname)
    nrow, ncol = rawdata.shape
    bins = map(int, rawdata[0, 1::3])
    data = 1.0 * rawdata[:, 2::3] + 1.0j * rawdata[:, 3::3]
    lofreqs = rawdata[:, 0] * 1e6
    framelen = 2**16
    return lofreqs, bins, data

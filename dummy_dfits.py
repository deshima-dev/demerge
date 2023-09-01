"""テスト用のダミーデータからなるdfits.gzファイルを生成する"""

from astropy.io import fits
from datetime   import datetime, timedelta

import math
import yaml
import numpy          as np
import merge_function as mf

if __name__ == '__main__':
    with open('dfits_dict.yaml', 'r') as f:
        dfits_dict = yaml.load(f, Loader=yaml.Loader)


    od = dfits_dict['obsinfo_dict']
    
    od['hdr_vals']['FITSTYPE'] = 'DESHIMAv2'
    od['hdr_vals']['TELESCOP'] = 'dummy_telescope'
    od['hdr_vals']['SITELON']  = 1.0
    od['hdr_vals']['SITELAT']  = 2.0
    od['hdr_vals']['DATE-OBS'] = '2023-08-31T14:14:14'
    od['hdr_vals']['OBSERVER'] = 'dummy_observer'
    od['hdr_vals']['OBJECT']   = 'dummy_object'
    od['hdr_vals']['RA']       = 3.0
    od['hdr_vals']['DEC']      = 4.0
    od['hdr_vals']['EQUINOX']  = 1000
    od['hdr_vals']['DDBID']    = '20230831'

    od['col_vals']['interval']  = np.array([0.1])
    od['col_vals']['integtime'] = np.array([0.2])
    od['col_vals']['beamsize']  = np.array([0.005]) # 18 arcsec
    od['col_vals']['pixelid']   = np.array([99])
    od['col_vals']['offsetaz']  = np.array([1.1, 1.2, 1.3])
    od['col_vals']['offsetel']  = np.array([2.1, 2.2, 2.3])
    od['col_vals']['gain']      = np.array([3.1, 3.2, 3.3])
    od['col_vals']['masterids'] = np.array([[1, 2, 3]])
    od['col_vals']['kidids']    = np.array([[10, 20, 30]])
    od['col_vals']['kidtypes']  = np.array([[0, 1, 2]])
    od['col_vals']['kidfreqs']  = np.array([[300.1, 300.2, 300.3]])
    od['col_vals']['kidQs']     = np.array([[100.1, 200.1, 300.1]])

    obsinfo = mf.create_bintablehdu(od)

    # 時刻を生成する
    now = datetime.now()

    n = 10
    timestamps = [(now + timedelta(minutes=i)).isoformat() for i in range(n)] # 10分間
    
    ad = dfits_dict['antenna_dict']

    ad['hdr_vals']['FILENAME']  = 'antenna.log'
    ad['col_vals']['time']      = np.array(timestamps)
    ad['col_vals']['scantype']  = np.array(['GRAD']*n)
    ad['col_vals']['az']        = np.array([0.1*i for i in range(n)])
    ad['col_vals']['el']        = np.array([0.2*i for i in range(n)])
    ad['col_vals']['ra']        = np.array([0.3*i for i in range(n)])
    ad['col_vals']['dec']       = np.array([0.4*i for i in range(n)])
    ad['col_vals']['az_center'] = np.array([0.5*i for i in range(n)])
    ad['col_vals']['el_center'] = np.array([0.6*i for i in range(n)])

    antenna = mf.create_bintablehdu(ad)

    rd = dfits_dict['readout_dict']

    n = math.floor(10*60/0.005) # 10分間の点
    starttimes = [(now + timedelta(milliseconds=(i*5))).isoformat() for i in range(n)] # 10分間

    rd['hdr_vals']['FILENAME']  = 'readoutfile'
    rd['col_vals']['starttime'] = np.array(starttimes).astype(np.datetime64)
    rd['col_vals']['pixelid']   = np.array([99])
    rd['col_vals']['lin_phase'] = np.array([[1.1, 1.2, 1.3]]*n).astype(np.float64)
    rd['col_vals']['Tsignal']   = np.array([[2.1, 2.2, 2.3]]*n).astype(np.float64)

    readout = mf.create_bintablehdu(rd)

    n = 10*6 # 10秒間隔で10分間
    timestamps = [datetime.strftime(now + timedelta(seconds=i*10), '%Y%m%d%H%M%S') for i in range(n)] # 10分間

    wd = dfits_dict['weather_dict']
    wd['hdr_vals']['FILENAME']       = 'weather.log'
    wd['col_vals']['time']           = np.array(timestamps)
    wd['col_vals']['temperature']    = np.array([6.5]*n).astype(np.float64)
    wd['col_vals']['pressure']       = np.array([570.5]*n).astype(np.float64)
    wd['col_vals']['vapor-pressure'] = np.array([0.57]*n).astype(np.float64)
    wd['col_vals']['windspd']        = np.array([19.5]*n).astype(np.float64)
    wd['col_vals']['winddir']        = np.array([250.0]*n).astype(np.float64)

    weather = mf.create_bintablehdu(wd)


    n = 10 # 1分間隔で10分間
    timestamps = [(now + timedelta(minutes=i)).isoformat() for i in range(n)] # 10分間
    
    cd = dfits_dict['cabin_t_dict']
    cd['hdr_vals']['FILENAME']    = 'cabin.log'
    cd['col_vals']['time']        = np.array(timestamps)
    cd['col_vals']['upper_cabin'] = np.array([6.5]*n).astype(np.float64)
    cd['col_vals']['main_cabin']  = np.array([16.5]*n).astype(np.float64)
    cabin_t = mf.create_bintablehdu(cd)

    n = math.floor(10*60/0.001) # 1ミリ秒間隔で10分間の点
    datetimes = [(now + timedelta(milliseconds=(i))).timestamp() for i in range(n)] # 10分間

    # 半分の時間がたったら状態が変化するようにstatesを生成する
    states = []
    for i in range(n):
        if i < (n/2):
            states.append(0)
        else:
            states.append(1)

    sd = dfits_dict['skychop_dict']

    sd['hdr_vals']['FILENAME'] = 'skychop.log'
    sd['col_vals']['time']     = np.array(datetimes).astype(np.float64)
    sd['col_vals']['state']    = np.array(states).astype(np.int8)

    skychop = mf.create_bintablehdu(sd)
    
    # KIDSINFO(dfits_dictにKIDSINFOの定義はないので最初から作る)
    header = fits.Header()
    header['EXTNAME']  = 'KIDSINFO', 'name of binary data'
    header['FILENAME'] = 'localsweep.sweep', 'localsweep filename'

    columns = [
        fits.Column(name='pixelid',        format='I',  array=[99, 99, 99]),
        fits.Column(name='kidid',          format='I',  array=[1, 2, 3]),
        fits.Column(name='Pread',          format='E',  array=[1.1, 2.2, 3.3], unit='dBm'),
        fits.Column(name='fc',             format='E',  array=[1.1, 2.2, 3.3], unit='GHz'),
        fits.Column(name='yfc, linyfc',    format='2E', array=[[1.0, 1.1], [2.0, 2.2], [3.0, 3.3]]),
        fits.Column(name='fr, dfr (300K)', format='2E', array=[[1.0, 1.1], [2.0, 2.2], [3.0, 3.3]]),
        fits.Column(name='Qr, dQr (300K)', format='2E', array=[[1.0, 1.1], [2.0, 2.2], [3.0, 3.3]]),
        fits.Column(name='Qc, dQc (300K)', format='2E', array=[[1.0, 1.1], [2.0, 2.2], [3.0, 3.3]]),
        fits.Column(name='Qi, dQi (300K)', format='2E', array=[[1.0, 1.1], [2.0, 2.2], [3.0, 3.3]]),
    ]
    kidsinfo = fits.BinTableHDU.from_columns(columns, header)

    hdul = fits.HDUList()
    hdul.append(fits.PrimaryHDU())
    hdul.append(obsinfo)
    hdul.append(antenna)
    hdul.append(readout)
    hdul.append(weather)
    hdul.append(skychop)
    hdul.append(cabin_t)
    hdul.append(kidsinfo)
    hdul.writeto('dfits_dummy.fits.gz', overwrite=True)

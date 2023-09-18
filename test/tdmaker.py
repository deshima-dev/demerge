"""単純なテストデータを作成するクラス

python 3.9

(C) 内藤システムズ
"""
from datetime      import datetime, timedelta, timezone
from astropy.io    import fits, ascii
from astropy.table import Table

import numpy as np
import math
import sys

class TestDataMaker():
    def __init__(self, time):
        """コンストラクタ

        引数
        ====
        integer 測定時間(分)

        戻り値
        ======
        なし
        """
        self.time       = time
        self.over       = 1 # READOUTより少し多めの観測時間(分)
        measure_time    = timedelta(minutes=self.time + self.over) # READOUTの時間より環境測定の時間の方が長くなるように設定
        readout_time    = timedelta(minutes=self.time)
        self.begin_time = datetime.now(tz=timezone.utc)
        self.end_time   = self.begin_time + measure_time

        # 搭載されているMKIDの数
        self.n_kid = 63

        # データ取得周期(秒)
        self.T_antenna = 0.1
        self.T_readout = 0.00625
        self.T_skychop = 0.001
        self.T_weather = 10
        self.T_misti   = 0.1
        self.T_cabin   = 60

        # データ点数
        measure_time = measure_time.total_seconds()
        readout_time = readout_time.total_seconds()
        self.n_readout = math.floor(readout_time/self.T_readout)
        self.n_antenna = math.floor(measure_time/self.T_antenna)
        self.n_skychop = math.floor(measure_time/self.T_skychop)
        self.n_weather = math.floor(measure_time/self.T_weather)
        self.n_misti   = math.floor(measure_time/self.T_misti)
        self.n_cabin   = math.floor(measure_time/self.T_cabin)
        print('n_readout: {}'.format(self.n_readout))
        print('n_skychop: {}'.format(self.n_skychop))
        return

    def generate_all(self):
        prefix = 'testdata'
        self.antenna.write('{}.ant'.format(prefix),     format='ascii.commented_header', overwrite=True)
        self.skychop.write('{}.skychop'.format(prefix), format='ascii.commented_header', overwrite=True)
        self.weather.write('{}.wea'.format(prefix),     format='ascii.commented_header', overwrite=True)
        self.cabin.write('{}.cabin'.format(prefix),     format='ascii.commented_header', overwrite=True)
        self.misti.write('{}.misti'.format(prefix),     format='ascii.no_header',        overwrite=True)
        self.ddb.writeto('{}_DDB.fits.gz'.format(prefix),              overwrite=True)
        self.readout.writeto('{}_reduced_readout.fits'.format(prefix), overwrite=True)
        self.dfits.writeto('{}_dfits.fits.gz'.format(prefix),          overwrite=True)
        return

    @property
    def antenna(self):
        antenna_table = Table()
        # ANTENNA時刻は秒が少数第一桁まで。そのため下2桁から6桁までを[:-5]を用いて文字列として削除している。
        antenna_table['time'] = [(self.begin_time + timedelta(milliseconds=self.T_antenna*1e3*i)).strftime('%Y%m%d%H%M%S.%f')[:-5] for i in range(self.n_antenna)]
        bias = 2.1
        dummy = np.array([1.1 for i in range(self.n_antenna)])
        antenna_table['ra-prg']          = dummy
        antenna_table['dec-prg']         = dummy
        antenna_table['az-prg']          = dummy + bias
        antenna_table['el-prg']          = dummy + bias
        antenna_table['az-real']         = dummy
        antenna_table['el-real']         = dummy
        antenna_table['x']               = dummy
        antenna_table['y']               = dummy
        antenna_table['z']               = dummy
        antenna_table['xt']              = dummy
        antenna_table['yt']              = dummy
        antenna_table['zt']              = dummy
        antenna_table['lst']             = dummy
        antenna_table['az-prg(no-col)']  = dummy
        antenna_table['el-prog(no-col)'] = dummy
        antenna_table['az-prog(center)'] = dummy
        antenna_table['el-prog(center)'] = dummy
        antenna_table['type']            = ['GRAD']*self.n_antenna

        # READOUTの時間で丁度半分の時刻でONに切り替えるためにantennaの打刻数を比を使って補正
        #
        #                     self.time                             self.over
        # |------------------------------------------------|-------------------------|
        # start                                            end of readout            end of antenna
        #
        # 全打点数をNとするreadoutが終わるまでの時刻の打点数nは
        #           n = N * (time / (time + over))
        # で表せる。
        #
        n = math.floor(self.n_antenna*(self.time/(self.time + self.over)))
        antenna_table['type'][math.floor(n/2):] = 'ON'
        return antenna_table

    @property
    def skychop(self):
        skychop_table = Table()

        skychop_table['ts']    = [(self.begin_time + timedelta(seconds=self.T_skychop*i)).timestamp() for i in range(self.n_skychop)]
        skychop_table['state'] = [1]*self.n_skychop

        # READOUTの時間で丁度半分の時刻で0に切り替えるためにskychopの打刻数を比を使って補正(理屈はantennaの補正と同じ)
        n = math.floor(self.n_skychop*(self.time/(self.time + self.over)))
        skychop_table['state'][math.floor(n/2):] = 0
        return skychop_table

    @property
    def weather(self):
        weather_table = Table()

        dummy = [15.0 for i in range(self.n_weather)]

        weather_table['time'] = [(self.begin_time + timedelta(seconds=self.T_weather*i)).strftime('%Y%m%d%H%M%S') for i in range(self.n_weather)]
        weather_table['tmperature']     = dummy # 実データのtypoをそのまま再現する
        weather_table['presure']        = dummy
        weather_table['vapor-pressure'] = dummy
        weather_table['aux1']           = dummy
        weather_table['aux2']           = dummy
        weather_table['aux3']           = dummy
        return weather_table

    @property
    def misti(self):
        misti_table = Table()

        misti_table['date']    = [(self.begin_time + timedelta(seconds=self.T_misti*i)).strftime('%Y/%m/%d') for i in range(self.n_misti)]
        misti_table['time']    = [(self.begin_time + timedelta(seconds=self.T_misti*i)).strftime('%H:%M:%S.%f') for i in range(self.n_misti)]
        misti_table['az']      = [180.0]*self.n_misti
        misti_table['el']      = [90.0]*self.n_misti
        misti_table['pwv']     = [610.0]*self.n_misti
        misti_table['Tground'] = [250.0]*self.n_misti
        return misti_table

    @property
    def cabin(self):
        cabin_table = Table()

        dummy = np.array([15.0 for i in range(self.n_cabin)])
        bias  = 2.0

        cabin_table['date']  = [(self.begin_time + timedelta(seconds=self.T_cabin*i)).strftime('%Y/%m/%d') for i in range(self.n_cabin)]
        cabin_table['time']  = [(self.begin_time + timedelta(seconds=self.T_cabin*i)).strftime('%H:%M') for i in range(self.n_cabin)]
        cabin_table['col3']  = dummy + bias
        cabin_table['col4']  = dummy
        cabin_table['col5']  = dummy
        cabin_table['col6']  = dummy
        cabin_table['col7']  = dummy
        cabin_table['col8']  = dummy
        cabin_table['col9']  = dummy
        cabin_table['col10'] = dummy
        cabin_table['col11'] = dummy
        cabin_table['col12'] = dummy
        cabin_table['col13'] = dummy
        cabin_table['col14'] = dummy
        cabin_table['col15'] = dummy
        cabin_table['col16'] = dummy
        cabin_table['col17'] = dummy
        cabin_table['col18'] = dummy
        cabin_table['col19'] = dummy
        cabin_table['col20'] = dummy
        cabin_table['col21'] = dummy
        cabin_table['col22'] = dummy
        return cabin_table

    @property
    def ddb(self):
        n_data = 31944
        n_kid  = 66

        header = fits.Header()
        header['EXTNAME']  = 'KIDDES', 'name of binary data'
        header['FILENAME'] = 'LT119_FB2.2G_49ch.csv', 'input filename'
        header['PIXEL0']   = 'LT119_FB2.2G_49ch', 'name of pixel 0'
        header['NKID0']    = n_kid, 'number of KIDs (pixel 0)'
        columns = [
            fits.Column(name='pixelid',   format='I',  array=[0]*n_kid),
            fits.Column(name='kidid',     format='I',  array=[i for i in range(n_kid)]),
            fits.Column(name='masterid',  format='I',  array=[1 for i in range(n_kid)]),
            fits.Column(name='attribute', format='8A', array=['filter']*n_kid),
        ]
        kiddes = fits.BinTableHDU.from_columns(columns, header)

        header = fits.Header()
        header['EXTNAME']  = 'KIDFILT', 'name of binary data'
        header['FILENAME'] = 'filter_table_DDBXXX.npy', 'localsweep filename'
        header['JSONNAME'] = 'kid_corresp_XXX.json', 'localsweep filename'
        header['NKID0']    = self.n_kid, 'number of KIDs (pixel 0)'
        dummy = (1.5, 0.5)
        columns = [
            fits.Column(name='pixelid',  format='I',  array=[0]*self.n_kid),
            fits.Column(name='kidid',    format='I',  array=[i for i in range(self.n_kid)]),
            fits.Column(name='masterid', format='I',  array=[1 for i in range(self.n_kid)]),
            fits.Column(name='F_filter, dF_filter', format='2E', array=[dummy]*self.n_kid),
            fits.Column(name='Q_filter, dQ_filter', format='2E', array=[dummy]*self.n_kid),
        ]
        kidfilt = fits.BinTableHDU.from_columns(columns, header)

        header = fits.Header()
        header['EXTNAME']  = 'KIDRESP', 'name of binary data'
        header['FILENAME'] = 'responsibity_table_DDBXXX.npy', 'localsweep filename'
        header['JSONNAME'] = 'kid_corresp_XXX.json', 'localsweep filename'
        dummy = (1.0, 0.5, 1.0) # p0, etaf, T0
        columns = [
            fits.Column(name='pixelid',    format='I',  array=[0]*self.n_kid),
            fits.Column(name='kidid',      format='I',  array=[i for i in range(self.n_kid)]),
            fits.Column(name='cal params', format='3E', array=[dummy]*self.n_kid),
        ]
        kidresp = fits.BinTableHDU.from_columns(columns, header)

        hdul = fits.HDUList()
        hdul.append(fits.PrimaryHDU())
        hdul.append(kiddes)
        hdul.append(kidfilt)
        hdul.append(kidresp)
        return hdul

    @property
    def readout(self):
        header = fits.Header()
        header['EXTNAME']  = 'KIDSINFO', 'name of binary data'
        header['FILENAME'] = 'localsweep.sweep', 'localsweep filename'
        header['NKID0']    = self.n_kid, 'number of KIDs (pixel 0)'
        dummy = (1.1, 0.2)
        columns = [
            fits.Column(name='kidid',          format='I',  array=[i for i in range(self.n_kid)]),
            fits.Column(name='pixelid',        format='I',  array=[0]*self.n_kid),
            fits.Column(name='yfc, linyfc',    format='2E', array=[dummy]*self.n_kid),
            fits.Column(name='fr, dfr (300K)', format='2E', array=[dummy]*self.n_kid),
            fits.Column(name='Qr, dQr (300K)', format='2E', array=[dummy]*self.n_kid),
            fits.Column(name='Qc, dQc (300K)', format='2E', array=[dummy]*self.n_kid),
            fits.Column(name='Qi, dQi (300K)', format='2E', array=[dummy]*self.n_kid),
        ]
        kidsinfo = fits.BinTableHDU.from_columns(columns, header)

        header = fits.Header()
        header['EXTNAME']  = 'READOUT', 'name of binary data'
        header['FILENAME'] = 'fuga', 'localsweep filename'
        header['NKID0']    = self.n_kid, 'number of KIDs (pixel 0)'

        columns = [
            fits.Column(name='timestamp', format='D', array=[(self.begin_time + timedelta(microseconds=self.T_readout*1e6*i)).timestamp() for i in range(self.n_readout)]),
            fits.Column(name='pixelid',   format='I', array=[0]*self.n_readout),
        ]

        dummy = (1.0, 1.0, 1.0)
        for i in range(self.n_kid):
            name   = 'Amp, Ph, linPh {}'.format(i)
            column = fits.Column(name=name, format='3E', array=[dummy]*self.n_readout)
            columns.append(column)

        readout = fits.BinTableHDU.from_columns(columns, header)

        hdul = fits.HDUList()
        hdul.append(fits.PrimaryHDU())
        hdul.append(kidsinfo)
        hdul.append(readout)
        return hdul

    @property
    def dfits(self):
        header = fits.Header()
        header['EXTNAME'] = 'OBSINFO', 'name of binary data'
        header['RA']      = 1.2, 'right ascension of the object in units of dec'
        header['DEC']     = 2.3, 'declination of the object in units of deg'
        dummy = 1.1
        columns = [
            fits.Column(name='masterids', format='63K', array=[i for i in range(self.n_kid)]),
            fits.Column(name='kidids',    format='63K', array=[i for i in range(self.n_kid)]),
            fits.Column(name='kidtypes',  format='63K', array=[1]*self.n_kid),
            fits.Column(name='kidfreqs',  format='63D', array=[dummy]*self.n_kid),
        ]
        obsinfo = fits.BinTableHDU.from_columns(columns, header)

        header = fits.Header()
        header['EXTNAME'] = 'READOUT', 'name of binary data'
        rmin = 100
        rmax = 300
        dr   = (rmax - rmin)/self.n_readout
        dummy = (1.1)*self.n_kid
        columns = [
            fits.Column(name='starttime', format='26A', array=[(self.begin_time + timedelta(microseconds=self.T_readout*1e6*i)).isoformat() for i in range(self.n_readout)]),
            fits.Column(name='Tsignal',   format='63D', array=[rmin + dr*i for i in range(self.n_readout)]),
        ]
        readout = fits.BinTableHDU.from_columns(columns, header)

        header = fits.Header()
        header['EXTNAME'] = 'ANTENNA', 'name of binary data'
        dummy = 1.1
        columns = [
            fits.Column(name='time',      format='26A', array=[(self.begin_time + timedelta(milliseconds=self.T_antenna*1e3*i)).isoformat() for i in range(self.n_antenna)]),
            fits.Column(name='az',        format='D',   array=[1.5 ]*self.n_antenna),
            fits.Column(name='el',        format='D',   array=[1.25]*self.n_antenna),
            fits.Column(name='az_center', format='D',   array=[0.5 ]*self.n_antenna),
            fits.Column(name='el_center', format='D',   array=[0.25]*self.n_antenna),
            fits.Column(name='ra',        format='D',   array=[2.5 ]*self.n_antenna),
            fits.Column(name='dec',       format='D',   array=[2.25]*self.n_antenna),
            fits.Column(name='scantype',  format='4A',  array=['GRAD']*self.n_antenna),
        ]
        antenna = fits.BinTableHDU.from_columns(columns, header)

        header = fits.Header()
        header['EXTNAME'] = 'WEATHER', 'name of binary data'
        dummy = 1.1
        columns = [
            fits.Column(name='time',           format='26A', array=[(self.begin_time + timedelta(seconds=self.T_weather*i)).isoformat() for i in range(self.n_weather)]),
            fits.Column(name='temperature',    format='D',   array=[1.5 ]*self.n_weather),
            fits.Column(name='pressure',       format='D',   array=[1.25]*self.n_weather),
            fits.Column(name='vapor-pressure', format='D',   array=[0.5 ]*self.n_weather),
            fits.Column(name='windspd',        format='D',   array=[0.25]*self.n_weather),
            fits.Column(name='winddir',        format='D',   array=[2.5 ]*self.n_weather),
        ]
        weather = fits.BinTableHDU.from_columns(columns, header)

        hdul = fits.HDUList()
        hdul.append(fits.PrimaryHDU())
        hdul.append(obsinfo)
        hdul.append(readout)
        hdul.append(antenna)
        hdul.append(weather)
        return hdul
    
if __name__ == '__main__':
    """
    説明
    ====
    並列処理を行ってデータを生成する時はコマンドラインの第一引数にデータ名を指定する。
    """
    tdm = TestDataMaker(3) # 3分間の測定データを生成する

    args = sys.argv
    if (len(args) == 1):
        tdm.generate_all()
        sys.exit(0)

    if (len(args) > 1):
        data_name = args[1]
        prefix    = 'testdata'
        if (data_name == 'antenna'):
            tdm.antenna.write('{}.ant'.format(prefix), format='ascii.commented_header', overwrite=True)
            sys.exit(0)
        if (data_name == 'skychop'):
            tdm.skychop.write('{}.skychop'.format(prefix), format='ascii.commented_header', overwrite=True)
            sys.exit(0)
        if (data_name == 'weather'):
            tdm.weather.write('{}.wea'.format(prefix), format='ascii.commented_header', overwrite=True)
            sys.exit(0)
        if (data_name == 'misti'):
            tdm.misti.write('{}.misti'.format(prefix), format='ascii.no_header', overwrite=True)
            sys.exit(0)
        if (data_name == 'cabin'):
            tdm.cabin.write('{}.cabin'.format(prefix), format='ascii.commented_header', overwrite=True)
            sys.exit(0)
        if (data_name == 'ddb'):
            tdm.ddb.writeto('{}_DDB.fits.gz'.format(prefix), overwrite=True)
            sys.exit(0)
        if (data_name == 'readout'):
            tdm.readout.writeto('{}_reduced_readout.fits'.format(prefix), overwrite=True)
            sys.exit(0)
        if (data_name == 'dfits'):
            tdm.dfits.writeto('{}_dfits.fits.gz'.format(prefix), overwrite=True)
            sys.exit(0)
        


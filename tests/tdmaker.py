"""単純なテストデータを作成するクラス

python 3.9

(C) 内藤システムズ
"""
from datetime      import datetime, timedelta, timezone
from astropy.io    import fits, ascii
from astropy.table import Table

import numpy as np
import lzma
import math
import sys
import argparse

class TestDataMaker():
    def __init__(self, time, **kwargs):
        """データ取得周期や打点数の計算

        引数
        ====
        integer 測定時間(分)

        戻り値
        ======
        なし
        """
        # 引数の処理
        self.time             = time
        self.p0               = kwargs.pop('p0',               1.0)
        self.etaf             = kwargs.pop('etaf',             0.5)
        self.T0               = kwargs.pop('T0',               1.0)
        self.linyfc           = kwargs.pop('linyfc',           0.25)
        self.Qr               = kwargs.pop('Qr',               1.1)
        self.lower_cabin_temp = kwargs.pop('lower_cabin_temp', 15)
        self.linear_readout   = kwargs.pop('linear_readout',   '') # inc, dec 増加か減少を選ぶ
        self.all_grad         = kwargs.pop('all_grad',         False)
        self.linear_antenna   = kwargs.pop('linear_antenna',   False)
        self.measure_time     = kwargs.pop('measure_time',     time) # 測定時間を意図的に変更する場合に使う

        readout_time = timedelta(minutes=self.time)
        measure_time = timedelta(minutes=self.measure_time)

        self.begin_time = datetime(year=2023, month=9, day=22, hour=9, minute=00, tzinfo=timezone.utc)
        self.end_time   = self.begin_time + measure_time

        # 見つかったMKIDの数
        self.n_kid = 63

        # データ取得周期(秒)
        self.T_antenna = 0.1
        self.T_readout = 0.00625
        self.T_skychop = 0.001
        self.T_weather = 10
        self.T_misti   = 0.1
        self.T_cabin   = 60

        # データ点数
        #
        # !!注意!!
        # measure_timeとreadout_timeはここで秒に置き換えられる!
        #
        # 打点数は割り算の結果を小数点以下切り捨てで計算されるため、
        # end_timeにぴったりの打点数以下になる。
        # そのためREADOUTの時刻に併せて補間すると毎回時刻の終わりに外挿が発生しNaNが格納される。
        #
        self.measure_time = measure_time.total_seconds()
        self.readout_time = readout_time.total_seconds()
        self.n_readout = math.floor(self.readout_time/self.T_readout)
        self.n_antenna = math.floor(self.measure_time/self.T_antenna)
        self.n_skychop = math.floor(self.measure_time/self.T_skychop)
        self.n_weather = math.floor(self.measure_time/self.T_weather)
        self.n_misti   = math.floor(self.measure_time/self.T_misti)
        self.n_cabin   = math.floor(self.measure_time/self.T_cabin)

        self.over = (self.measure_time - self.readout_time)/60 # 秒から分へ変換
        if self.over < 0:
            self.over = 0
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
        bias = 2.2
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

        # lonがs線形に0から180まで変化するデータを作る
        if self.linear_antenna:
            lon_max = 180
            lon = [lon_max/self.measure_time*(i*self.T_antenna) for i in range(self.n_antenna)]

            antenna_table['az-prg']          = [0.0]*self.n_antenna
            antenna_table['el-prg']          = [0.0]*self.n_antenna
            antenna_table['az-real']         = [0.0]*self.n_antenna
            antenna_table['el-real']         = [0.0]*self.n_antenna
            antenna_table['az-prg(no-col)']  = lon
            antenna_table['el-prog(no-col)'] = [0.0]*self.n_antenna

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
        if self.all_grad == False:
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
        cabin_table['col4']  = self.lower_cabin_temp
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
        header['DDB_ID'] = 'YYYYmmdd'
        primary = fits.PrimaryHDU(header=header)

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
        dummy = (self.p0, self.etaf, self.T0)
        columns = [
            fits.Column(name='pixelid',    format='I',  array=[0]*self.n_kid),
            fits.Column(name='kidid',      format='I',  array=[i for i in range(self.n_kid)]),
            fits.Column(name='cal params', format='3E', array=[dummy]*self.n_kid),
        ]
        kidresp = fits.BinTableHDU.from_columns(columns, header)

        hdul = fits.HDUList()
        hdul.append(primary)
        hdul.append(kiddes)
        hdul.append(kidfilt)
        hdul.append(kidresp)
        return hdul

    @property
    def readout(self, **kwargs):
        header = fits.Header()
        header['EXTNAME']  = 'KIDSINFO', 'name of binary data'
        header['FILENAME'] = 'localsweep.sweep', 'localsweep filename'
        header['NKID0']    = self.n_kid, 'number of KIDs (pixel 0)'
        dummy = (1.5, 0.25)
        columns = [
            fits.Column(name='kidid',          format='I',  array=[i for i in range(self.n_kid)]),
            fits.Column(name='pixelid',        format='I',  array=[0]*self.n_kid),
            fits.Column(name='yfc, linyfc',    format='2E', array=[(1.0, self.linyfc)]*self.n_kid),
            fits.Column(name='fr, dfr (300K)', format='2E', array=[dummy]*self.n_kid),
            fits.Column(name='Qr, dQr (300K)', format='2E', array=[(self.Qr, 0.25)]*self.n_kid),
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

        dummy = None
        linPh_max = 300
        if self.linear_readout == 'inc':
            linPh = np.sqrt(linPh_max/self.readout_time)
            dummy = [(1.0, 1.0, linPh*np.sqrt(i*self.T_readout)) for i in range(self.n_readout)]
        elif self.linear_readout == 'dec':
            dummy = [(1.0, 1.0, np.sqrt(-linPh_max/self.readout_time*i*self.T_readout + linPh_max)) for i in range(self.n_readout)]
        else:
            dummy = [(1.0, 1.0, 1.0)]*self.n_readout

        for i in range(self.n_kid):
            name   = 'Amp, Ph, linPh {}'.format(i)
            column = fits.Column(name=name, format='3E', array=dummy)
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
        rmin = 0
        rmax = 300
        dr   = (rmax - rmin)/self.n_readout
        dummy = (1.1)*self.n_kid

        if self.linear_readout == 'inc':
            Tsignal = [rmin + dr*i for i in range(self.n_readout)]
        elif self.linear_readout == 'dec':
            Tsignal = [rmax - dr*i for i in range(self.n_readout)]
        else:
            Tsignal = [rmin + dr*i for i in range(self.n_readout)]


        columns = [
            fits.Column(name='starttime', format='26A', array=[(self.begin_time + timedelta(microseconds=self.T_readout*1e6*i)).isoformat() for i in range(self.n_readout)]),
            fits.Column(name='Tsignal',   format='63D', array=Tsignal),
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
    引数
    ====

    説明
    ====
    並列処理を行ってデータを生成する時はコマンドラインの第一引数にデータ名を指定する。
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('data_name',          type=str,   default='', nargs='?')
    parser.add_argument('--time',             type=int,   default=3,          help='測定時間(分)を整数で指定して下さい')
    parser.add_argument('--p0',               type=float, default=1.0,        help='p0をfloatで指定して下さい')
    parser.add_argument('--etaf',             type=float, default=0.5,        help='etafをfloatで指定して下さい')
    parser.add_argument('--T0',               type=float, default=1.0,        help='T0をfloatで指定して下さい')
    parser.add_argument('--Qr',               type=float, default=1.1,        help='Qrをfloatで指定して下さい')
    parser.add_argument('--linyfc',           type=float, default=0.25,       help='linyfcをfloatで指定して下さい')
    parser.add_argument('--linear_readout',   type=str,   default='',         help='readoutの値を線形に変化させる場合はinc/decのいずれかを指定して下さい')
    parser.add_argument('--linear_antenna',   type=bool,  default=False,      help='antennaのlonを線形に変化させる場合はTrueを指定して下さい')
    parser.add_argument('--all_grad',         type=bool,  default=False,      help='すべてのSCAN状態をGRADにする場合はTrueを指定して下さい')
    parser.add_argument('--lower_cabin_temp', type=float, default=15,         help='MainCabinの温度(degC)をfloatで指定して下さい')
    parser.add_argument('--prefix',           type=str,   default='testdata', help='生成されるファイル名のprefixを指定して下さい')
    parser.add_argument('--measure_time',     type=int,   default=None,       help='環境測定時間(分)を整数で指定して下さい')
    parser.add_argument('--xz',               type=bool,  default=False,      help='Trueを指定するとskychopファイルをxzで圧縮する')
    parser.add_argument('--gz',               type=bool,  default=False,      help='Trueを指定するとreadoutファイルをgzで圧縮する')
    a = parser.parse_args()

    if a.measure_time == None:
        a.measure_time = a.time

    tdm = TestDataMaker(time            =a.time,
                        p0              =a.p0,
                        etaf            =a.etaf,
                        T0              =a.T0,
                        Qr              =a.Qr,
                        linyfc          =a.linyfc,
                        lower_cabin_temp=a.lower_cabin_temp,
                        linear_readout  =a.linear_readout,
                        linear_antenna  =a.linear_antenna,
                        all_grad        =a.all_grad,
                        measure_time    =a.measure_time,
                        )

    if a.data_name == '':
        tdm.generate_all()
        sys.exit(0)

    if (a.data_name == 'antenna'):
        tdm.antenna.write('{}.ant'.format(a.prefix), format='ascii.commented_header', overwrite=True)
        sys.exit(0)
    if (a.data_name == 'skychop'):
        if (a.xz):
            with lzma.open('{}.skychop.dat.xz'.format(a.prefix), 'wt') as f:
                tdm.skychop.write(f, format='ascii.commented_header', overwrite=True)
        else:
            tdm.skychop.write('{}.skychop'.format(a.prefix), format='ascii.commented_header', overwrite=True)
        sys.exit(0)
    if (a.data_name == 'weather'):
        tdm.weather.write('{}.wea'.format(a.prefix), format='ascii.commented_header', overwrite=True)
        sys.exit(0)
    if (a.data_name == 'misti'):
        tdm.misti.write('{}.misti'.format(a.prefix), format='ascii.no_header', overwrite=True)
        sys.exit(0)
    if (a.data_name == 'cabin'):
        tdm.cabin.write('{}.cabin'.format(a.prefix), format='ascii.commented_header', overwrite=True)
        sys.exit(0)
    if (a.data_name == 'ddb'):
        tdm.ddb.writeto('{}_DDB.fits.gz'.format(a.prefix), overwrite=True)
        sys.exit(0)
    if (a.data_name == 'readout'):
        if (a.gz):
            tdm.readout.writeto('{}_reduced_readout.fits.gz'.format(a.prefix), overwrite=True)
        else:
            tdm.readout.writeto('{}_reduced_readout.fits'.format(a.prefix), overwrite=True)
        sys.exit(0)
    if (a.data_name == 'dfits'):
        tdm.dfits.writeto('{}_dfits.fits.gz'.format(a.prefix), overwrite=True)
        sys.exit(0)

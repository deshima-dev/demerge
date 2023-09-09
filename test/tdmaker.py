"""テストデータを作成するクラス

(C) 内藤システムズ
"""
from datetime      import datetime, timedelta
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
        measure_time = timedelta(minutes=time)
        self.begin_time   = datetime.now()
        self.end_time     = self.begin_time + measure_time

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
        self.n_antenna = math.floor(measure_time/self.T_antenna)
        self.n_readout = math.floor(measure_time/self.T_readout)
        self.n_skychop = math.floor(measure_time/self.T_skychop)
        self.n_weather = math.floor(measure_time/self.T_weather)
        self.n_misti   = math.floor(measure_time/self.T_misti)
        self.n_cabin   = math.floor(measure_time/self.T_cabin)
        return

    def generate_all(self):
        prefix = 'testdata'
        self.antenna.write('{}.ant'.format(prefix),     format='ascii.commented_header', overwrite=True)
        self.skychop.write('{}.skychop'.format(prefix), format='ascii.commented_header', overwrite=True)
        self.weather.write('{}.wea'.format(prefix),     format='ascii.commented_header', overwrite=True)
        self.cabin.write('{}.cabin'.format(prefix),     format='ascii.commented_header', overwrite=True)
        self.misti.write('{}.misti'.format(prefix),     format='ascii.no_header', overwrite=True)
        self.ddb.writeto('{}_DDB.fits.gz'.format(prefix),              overwrite=True)
        self.readout.writeto('{}_reduced_readout.fits'.format(prefix), overwrite=True)
        return

    @property
    def antenna(self):
        antenna_table = Table()
        # ANTENNA時刻は秒が少数第一桁まで。そのため下2桁から6桁までを[:-5]を用いて文字列として削除している。
        antenna_table['time'] = [(self.begin_time + timedelta(seconds=self.T_antenna*i)).strftime('%Y%m%d%H%M%S.%f')[:-5] for i in range(self.n_antenna)]
        tmp = [0.0 for i in range(self.n_antenna)]
        antenna_table['ra-prg']          = tmp
        antenna_table['dec-prg']         = tmp
        antenna_table['az-prg']          = tmp
        antenna_table['el-prg']          = tmp
        antenna_table['az-real']         = tmp
        antenna_table['el-real']         = tmp
        antenna_table['x']               = tmp
        antenna_table['y']               = tmp
        antenna_table['z']               = tmp
        antenna_table['xt']              = tmp
        antenna_table['yt']              = tmp
        antenna_table['zt']              = tmp
        antenna_table['lst']             = tmp
        antenna_table['az-prg(no-col)']  = tmp
        antenna_table['el-prog(no-col)'] = tmp
        antenna_table['az-prog(center)'] = tmp
        antenna_table['el-prog(center)'] = tmp
        antenna_table['type']            = ['GRAD']*self.n_antenna
        return antenna_table

    @property
    def skychop(self):
        skychop_table = Table()

        skychop_table['ts']    = [(self.begin_time + timedelta(seconds=self.T_skychop*i)).timestamp() for i in range(self.n_skychop)]
        skychop_table['state'] = [1]*self.n_skychop
        return skychop_table

    @property
    def weather(self):
        weather_table = Table()

        tmp = [0.0 for i in range(self.n_weather)]

        weather_table['time'] = [(self.begin_time + timedelta(seconds=self.T_weather*i)).strftime('%Y%m%d%H%M%S') for i in range(self.n_weather)]
        weather_table['temperature']    = tmp
        weather_table['presure']        = tmp
        weather_table['vapor-pressure'] = tmp
        weather_table['aux1']           = tmp
        weather_table['aux2']           = tmp
        weather_table['aux3']           = tmp
        return weather_table

    @property
    def misti(self):
        misti_table = Table()

        tmp = [0.0 for i in range(self.n_misti)]

        misti_table['YYYY']  = [(self.begin_time + timedelta(seconds=self.T_misti*i)).strftime('%Y') for i in range(self.n_misti)]
        misti_table['mm']    = [(self.begin_time + timedelta(seconds=self.T_misti*i)).strftime('%m') for i in range(self.n_misti)]
        misti_table['dd']    = [(self.begin_time + timedelta(seconds=self.T_misti*i)).strftime('%d') for i in range(self.n_misti)]
        misti_table['HH']    = [(self.begin_time + timedelta(seconds=self.T_misti*i)).strftime('%H') for i in range(self.n_misti)]
        misti_table['MM']    = [(self.begin_time + timedelta(seconds=self.T_misti*i)).strftime('%M') for i in range(self.n_misti)]
        misti_table['SS.SS'] = [(self.begin_time + timedelta(seconds=self.T_misti*i)).strftime('%S.%f')[:-4] for i in range(self.n_misti)]
        misti_table['az']                       = tmp
        misti_table['el']                       = tmp
        misti_table['power']                    = tmp
        misti_table['hot_load_temp']            = tmp
        misti_table['receiver_room_temp']       = tmp
        misti_table['primary_mirror_room_temp'] = tmp
        misti_table['chopper_mirror_status']    = [1]*self.n_misti
        return misti_table

    @property
    def cabin(self):
        cabin_table = Table()

        tmp = [0.0 for i in range(self.n_cabin)]

        cabin_table['date']  = [(self.begin_time + timedelta(seconds=self.T_cabin*i)).strftime('%Y/%m/%d') for i in range(self.n_cabin)]
        cabin_table['time']  = [(self.begin_time + timedelta(seconds=self.T_cabin*i)).strftime('%H:%M') for i in range(self.n_cabin)]
        cabin_table['col3']  = tmp
        cabin_table['col4']  = tmp
        cabin_table['col5']  = tmp
        cabin_table['col6']  = tmp
        cabin_table['col7']  = tmp
        cabin_table['col8']  = tmp
        cabin_table['col9']  = tmp
        cabin_table['col10'] = tmp
        cabin_table['col11'] = tmp
        cabin_table['col12'] = tmp
        cabin_table['col13'] = tmp
        cabin_table['col14'] = tmp
        cabin_table['col15'] = tmp
        cabin_table['col16'] = tmp
        cabin_table['col17'] = tmp
        cabin_table['col18'] = tmp
        cabin_table['col19'] = tmp
        cabin_table['col20'] = tmp
        cabin_table['col21'] = tmp
        cabin_table['col22'] = tmp
        return cabin_table

    @property
    def ddb(self):
        n_data = 31944

        header = fits.Header()
        header['EXTNAME']  = 'KIDFILT', 'name of binary data'
        header['FILENAME'] = 'filter_table_DDBXXX.npy', 'localsweep filename'
        header['JSONNAME'] = 'kid_corresp_XXX.json', 'localsweep filename'
        columns = [
            fits.Column(name='pixelid',  format='I',  array=[0]*self.n_kid),
            fits.Column(name='kidid',    format='I',  array=[i for i in range(self.n_kid)]),
            fits.Column(name='masterid', format='I',  array=[1 for i in range(self.n_kid)]),
        ]
        kidfilt = fits.BinTableHDU.from_columns(columns, header)

        header = fits.Header()
        header['EXTNAME']  = 'KIDRESP', 'name of binary data'
        header['FILENAME'] = 'responsibity_table_DDBXXX.npy', 'localsweep filename'
        header['JSONNAME'] = 'kid_corresp_XXX.json', 'localsweep filename'
        columns = [
            fits.Column(name='pixelid',    format='I',  array=[0]*self.n_kid),
            fits.Column(name='kidid',      format='I',  array=[i for i in range(self.n_kid)]),
            fits.Column(name='cal params', format='3E', array=[(1.0, 1.0, 1.0)]*self.n_kid),
        ]
        kidresp = fits.BinTableHDU.from_columns(columns, header)

        hdul = fits.HDUList()
        hdul.append(fits.PrimaryHDU())
        hdul.append(kidfilt)
        hdul.append(kidresp)
        return hdul

    @property
    def readout(self):
        header = fits.Header()
        header['EXTNAME']  = 'READOUT', 'name of binary data'
        header['FILENAME'] = 'fuga', 'localsweep filename'
        columns = [
            fits.Column(name='timestamp', format='D',  array=[(self.begin_time + timedelta(seconds=self.T_readout*i)).timestamp() for i in range(self.n_readout)]),
            fits.Column(name='pixelid',   format='I',  array=[0]*self.n_readout),
        ]

        for i in range(self.n_kid):
            name   = 'Amp, Ph, linPh {}'.format(i)
            column = fits.Column(name=name, format='3E', array=[(1.0, 1.0, 1.0)]*self.n_readout)
            columns.append(column)

        readout = fits.BinTableHDU.from_columns(columns, header)

        hdul = fits.HDUList()
        hdul.append(fits.PrimaryHDU())
        hdul.append(readout)
        return hdul

if __name__ == '__main__':
    """
    説明
    ====
    並列処理を行ってデータを生成する時はコマンドラインの第一引数にデータ名を指定する。
    """
    tdm = TestDataMaker(10) # 10分間の測定データを生成する

    args = sys.argv
    if (len(args) == 1):
        tdm.generate_all()
        sys.exit(0)

    if (len(args) > 1):
        data_name = args[1]
        if (data_name == 'antenna'):
            tdm.antenna.write('testdata.ant', format='ascii.commented_header', overwrite=True)
            sys.exit(0)
        if (data_name == 'skychop'):
            tdm.skychop.write('testdata.skychop', format='ascii.commented_header', overwrite=True)
            sys.exit(0)
        if (data_name == 'weather'):
            tdm.weather.write('testdata.wea', format='ascii.commented_header', overwrite=True)
            sys.exit(0)
        if (data_name == 'misti'):
            tdm.misti.write('testdata.misti', format='ascii.no_header', overwrite=True)
            sys.exit(0)
        if (data_name == 'cabin'):
            tdm.cabin.write('testdata.cabin', format='ascii.commented_header', overwrite=True)
            sys.exit(0)
        if (data_name == 'ddb'):
            tdm.ddb.writeto('testdata_DDB.fits.gz', overwrite=True)
            sys.exit(0)
        if (data_name == 'readout'):
            tdm.readout.writeto('testdata_reduced.fits', overwrite=True)
            sys.exit(0)
        


"""テストデータを作成する

(C) 内藤システムズ
"""
from datetime      import datetime, timedelta
from astropy.io    import fits, ascii
from astropy.table import Table

import numpy as np
import math

measure_time = timedelta(minutes=10)
begin_time   = datetime.now()
end_time     = begin_time + measure_time

# データ取得周期(秒)
T_antenna = 0.1
T_readout = 0.00625
T_skychop = 0.001
T_weather = 10
T_misti   = 0.1
T_cabin   = 60

# データ点数
measure_time = measure_time.total_seconds()
n_antenna = math.floor(measure_time/T_antenna)
n_readout = math.floor(measure_time/T_readout)
n_skychop = math.floor(measure_time/T_skychop)
n_weather = math.floor(measure_time/T_weather)
n_misti   = math.floor(measure_time/T_misti)
n_cabin   = math.floor(measure_time/T_cabin)

# antenna

antenna_table = Table()

# ANTENNA時刻は秒が少数第一桁まで。そのため下2桁から6桁までを[:-5]を用いて文字列として削除している。
antenna_table['time'] = [(begin_time + timedelta(seconds=T_antenna*i)).strftime('%Y%m%d%H%M%S.%f')[:-5] for i in range(n_antenna)]

tmp = [0.0 for i in range(n_antenna)]

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
antenna_table['type']            = ['GRAD']*n_antenna
antenna_table.write('testdata.ant', format='ascii.commented_header', overwrite=True)

# skychop
skychop_table = Table()

skychop_table['ts']    = [(begin_time + timedelta(seconds=T_skychop*i)).timestamp() for i in range(n_skychop)]
skychop_table['state'] = [1]*n_skychop
skychop_table.write('testdata.skychop', format='ascii.commented_header', overwrite=True)

# weather
weather_table = Table()

tmp = [0.0 for i in range(n_weather)]

weather_table['time'] = [(begin_time + timedelta(seconds=T_weather*i)).strftime('%Y%m%d%H%M%S') for i in range(n_weather)]
weather_table['temperature']    = tmp
weather_table['presure']        = tmp
weather_table['vapor-pressure'] = tmp
weather_table['aux1']           = tmp
weather_table['aux2']           = tmp
weather_table['aux3']           = tmp
weather_table.write('testdata.wea', format='ascii.commented_header', overwrite=True)

# misti
misti_table = Table()

tmp = [0.0 for i in range(n_misti)]

misti_table['YYYY']  = [(begin_time + timedelta(seconds=T_misti*i)).strftime('%Y') for i in range(n_misti)]
misti_table['mm']    = [(begin_time + timedelta(seconds=T_misti*i)).strftime('%m') for i in range(n_misti)]
misti_table['dd']    = [(begin_time + timedelta(seconds=T_misti*i)).strftime('%d') for i in range(n_misti)]
misti_table['HH']    = [(begin_time + timedelta(seconds=T_misti*i)).strftime('%H') for i in range(n_misti)]
misti_table['MM']    = [(begin_time + timedelta(seconds=T_misti*i)).strftime('%M') for i in range(n_misti)]
misti_table['SS.SS'] = [(begin_time + timedelta(seconds=T_misti*i)).strftime('%S.%f')[:-4] for i in range(n_misti)]
misti_table['az']                       = tmp
misti_table['el']                       = tmp
misti_table['power']                    = tmp
misti_table['hot_load_temp']            = tmp
misti_table['receiver_room_temp']       = tmp
misti_table['primary_mirror_room_temp'] = tmp
misti_table['chopper_mirror_status']    = [1]*n_misti
misti_table.write('testdata.misti', format='ascii.no_header', overwrite=True)

# cabin
cabin_table = Table()

tmp = [0.0 for i in range(n_cabin)]

cabin_table['date']  = [(begin_time + timedelta(seconds=T_cabin*i)).strftime('%Y/%m/%d') for i in range(n_cabin)]
cabin_table['time']  = [(begin_time + timedelta(seconds=T_cabin*i)).strftime('%H:%M') for i in range(n_cabin)]
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
cabin_table.write('testdata.cabin', format='ascii.commented_header', overwrite=True)

# ddb
n_data = 31944
n_kid  = 63

header = fits.Header()
header['EXTNAME']  = 'KIDFILT', 'name of binary data'
header['FILENAME'] = 'filter_table_DDBXXX.npy', 'localsweep filename'
header['JSONNAME'] = 'kid_corresp_XXX.json', 'localsweep filename'
columns = [
    fits.Column(name='pixelid',  format='I',  array=[0]*n_kid),
    fits.Column(name='kidid',    format='I',  array=[i for i in range(n_kid)]),
    fits.Column(name='masterid', format='I',  array=[1 for i in range(n_kid)]),
]
kidfilt = fits.BinTableHDU.from_columns(columns, header)

header = fits.Header()
header['EXTNAME']  = 'KIDRESP', 'name of binary data'
header['FILENAME'] = 'responsibity_table_DDBXXX.npy', 'localsweep filename'
header['JSONNAME'] = 'kid_corresp_XXX.json', 'localsweep filename'
columns = [
    fits.Column(name='pixelid',    format='I',  array=[0]*n_kid),
    fits.Column(name='kidid',      format='I',  array=[i for i in range(n_kid)]),
    fits.Column(name='cal params', format='3E', array=[(1.0, 1.0, 1.0)]*n_kid),
]
kidresp = fits.BinTableHDU.from_columns(columns, header)

hdul = fits.HDUList()
hdul.append(fits.PrimaryHDU())
hdul.append(kidfilt)
hdul.append(kidresp)
hdul.writeto('testdata_DDB.fits.gz', overwrite=True)

# readout

header = fits.Header()
header['EXTNAME']  = 'READOUT', 'name of binary data'
header['FILENAME'] = 'fuga', 'localsweep filename'
columns = [
    fits.Column(name='timestamp', format='D',  array=[(begin_time + timedelta(seconds=T_readout*i)).timestamp() for i in range(n_readout)]),
    fits.Column(name='pixelid',   format='I',  array=[0]*n_readout),
]

for i in range(n_kid):
    name   = 'Amp, Ph, linPh {}'.format(i)
    column = fits.Column(name=name, format='3E', array=[(1.0, 1.0, 1.0)]*n_readout)
    columns.append(column)

readout = fits.BinTableHDU.from_columns(columns, header)

hdul = fits.HDUList()
hdul.append(fits.PrimaryHDU())
hdul.append(readout)
hdul.writeto('testdata_reduced_readout.fits', overwrite=True)

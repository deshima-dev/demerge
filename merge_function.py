"""merge_to_dfits.py: Read logging data and merge them into a FITS object

 Author : Tetsutaro Ueda, Junya Suzuki, Kenichi Karatsu, Tatsuya Takekoshi
 Created: 2017/11/02
 Revision History:
     2018/02/08 - KK - rewrite using class.
     2018/06/08 - TT - apply new calibration method.
     2021         NAITO systems modfied.
"""
from __future__ import print_function

__all__ = [
    'FORM_FITSTIME',
    'FORM_FITSTIME_P',
    'DEFAULT_ROOM_T',
    'CabinTempDB',
    'create_bintablehdu',
    'load_obsinst',
    'get_maskid_corresp'
    'Tlos_model',
    'calibrate_to_power',
    'get_Troom',
    'convert_asciitime',
    'convert_timestamp'
]

from datetime import datetime
from datetime import timedelta
from calendar import timegm
import numpy as np
import scipy.interpolate
from astropy.io import fits
import sqlite3
import sys
import dfits2dems

#-------------------------------- CONSTANTS
FORM_FITSTIME   = '%Y-%m-%dT%H:%M:%S'                          # YYYY-mm-ddTHH:MM:SS
FORM_FITSTIME_P = '%Y-%m-%dT%H:%M:%S.%f'                       # YYYY-mm-ddTHH:MM:SS.ss

CABIN_Q_MARGIN  = 5*60 # seconds. Margin for cabin data query.
DEFAULT_ROOM_T  = 17. + 273. # Kelvin
DEFAULT_AMB_T  = 0. + 273. # Kelvin

class CabinTempDB:
    def __init__(self, db_path):
        self._conn = sqlite3.connect(db_path)
        self._cur = self._conn.cursor()
        self.db_path = db_path

    def get_between(self, start_str, end_str):
        sql = '''
        SELECT datetime, upper_cabin, main_cabin FROM cabin_temp
        WHERE
            datetime >= '{0}'
            AND
            datetime <= '{1}';
        '''.format(start_str, end_str)        
        self._cur.execute(sql)
        return self._cur.fetchall()

    def __del__(self):
        self._conn.close()

def create_bintablehdu(hd):
    """Create Binary Table HDU from 'hdu_dict'"""
    header = fits.Header()
    for (i, j) in zip(hd['hdr_vals'].items(), hd['hdr_coms'].items()):
        header[i[0]] = i[1], j[1]
    columns = [
        fits.Column(name=i[0], format=j[1], array=i[1], unit=k[1])
        for (i, j, k) in zip(
            hd['col_vals'].items(),
            hd['col_form'].items(),
            hd['col_unit'].items()
        )
    ]
    hdu = fits.BinTableHDU.from_columns(columns, header)
    for i in hd['hdr_coms'].items():
        hdu.header.comments[i[0]] = i[1]
    return hdu

def load_obsinst(obsinst):
    """Get data for 'OBSINFO'"""
    if not '.obs' in obsinst:
        raise ValueError('The input file must be an observational instruction!!')

    with open(obsinst, 'r') as f:
        equinox = 2000  # Default parameter
        for line in f:
            if 'SET ANTENNA_G TRK_TYPE' in line:
                trktype = line.split()[-1].strip('\'')
            elif 'SET ANTENNA_G SRC_NAME' in line:
                obs_object = line.split()[-1].strip('\'')
            elif 'SET ANTENNA_G SRC_POS' in line:
                srcpos = [float(c) for c in line.split()[-1].strip('()').split(',')]
            elif 'SET ANTENNA_G EPOCH' in line:
                equinox = line.split()[-1].strip('\'JB')
            elif 'SET DES OBS_USER' in line:
                observer = line.split()[-1].strip('\'')
    if trktype == 'RADEC':
        ra  = srcpos[0]
        dec = srcpos[1]
    else:
        ra  = 0
        dec = 0
    return {'observer': observer, 'obs_object': obs_object,  'ra': ra, 'dec': dec, 'equinox': equinox}

def get_maskid_corresp(pixelid, ddb):
    """Get Correspondance of 'master' and 'kid'"""
    nkid = ddb['KIDFILT'].header['NKID%d' %pixelid]
    kiddict, kidfilt = {}, {}
    for (i, j, k, l) in zip(ddb['KIDFILT'].data['kidid'],
                            ddb['KIDFILT'].data['masterid'],
                            ddb['KIDFILT'].data['F_filter, dF_filter'],
                            ddb['KIDFILT'].data['Q_filter, dQ_filter']):
        kiddict[i] = j
        kidfilt[i] = (k[0], l[0])
    kidname = {}
    for (i, j) in zip(ddb['KIDDES'].data['masterid'], ddb['KIDDES'].data['attribute']):
        kidname[i] = j

    masterids, kidids, kidtypes, kidfreqs, kidQs = [], [], [], [], []
    for i in range(nkid):
        masterid = kiddict[i]
        if masterid < 0:
            kind = 'unknown'
        else:
            kind = kidname[masterid]

        if kind=='wideband': attr = 0
        elif kind=='filter': attr = 1
        elif kind=='blind':  attr = 2
        elif kind=='Al':     attr = 3
        elif kind=='NbTiN':  attr = 4 
        else:                attr = -1

        masterids.append( masterid )
        kidids.append( i )        
        kidtypes.append( attr )
        kidfreqs.append( kidfilt[i][0] )
        kidQs.append( kidfilt[i][1] )
    return masterids, kidids, kidtypes, kidfreqs, kidQs

def Tlos_model(dx, p0, etaf, T0, Troom, Tamb):
    """Calibrate 'amplitude' and 'phase' to 'power'"""
    return (dx + p0*np.sqrt(Troom+T0))**2 / (p0**2 * etaf) - T0/etaf - (1-etaf)/etaf*Tamb

def calibrate_to_power(pixelid, Troom, Tamb, rhdus, ddb):
    nkid = rhdus['READOUT'].header['NKID%d' %pixelid]
    kiddict = {}
    for (i, j) in zip(ddb['KIDFILT'].data['kidid'], ddb['KIDFILT'].data['masterid']):
        kiddict[i] = j

    linphase = np.transpose([rhdus['READOUT'].data['Amp, Ph, linPh %d' %i].T[2] for i in range(nkid)])
    linyfc   = rhdus['KIDSINFO'].data['yfc, linyfc'].T[1]
    Qr       = rhdus['KIDSINFO'].data['Qr, dQr (300K)'].T[0]
    fshift = np.array((linphase - linyfc)/ (4.*Qr)).T
    fshift_err = np.zeros( len(fshift) )
    #---- Responsivity curve
    (p0, etaf, T0) = ddb['KIDRESP'].data['cal params'].T
    Tsignal = []
    for i in range(nkid):
        masterid = kiddict[i]
        if masterid<0:
            Tsignal.append( [np.nan for j in range( len(fshift[i]) )] )
            continue
        #---- Convert to power
        Tsignal.append(Tlos_model(fshift[i], p0[i], etaf[i], T0[i], Troom, Tamb))
    return np.array(Tsignal).T

# def get_Troom(rhdus, cabin_db):
#     """Get Troom from cabin temperature database"""
#     starttime = convert_timestamp(rhdus['READOUT'].data['timestamp'])
#     margin = timedelta(seconds=CABIN_Q_MARGIN)
#     db_start = datetime.strptime(starttime[0], FORM_FITSTIME_P)
#     db_start_mod = db_start - timedelta(microseconds=db_start.microsecond) - margin
#     db_end = datetime.strptime(starttime[-1], FORM_FITSTIME_P)
#     db_end_mod = db_end - timedelta(microseconds=db_end.microsecond) + margin        
#     db_ret = cabin_db.get_between(db_start_mod.isoformat(), db_end_mod.isoformat())
#     if len(db_ret)==0:
#         print('Cannot find cabin temperature data in', cabin_db.db_path, file=sys.stderr)
#         print('Instead, Troom', DEFAULT_ROOM_T, 'K will be used', file=sys.stderr)
#         Troom = DEFAULT_ROOM_T
#     else:
#         # Check whether the database covers the measurement 
#         db_dts = [datetime.strptime(entry[0], FORM_FITSTIME) for entry in db_ret]
#         main_t = [entry[2] for entry in db_ret]
#         if (db_start < db_dts[0]) or (db_dts[-1] < db_end):
#             print('DB data does not cover the measurement span', file=sys.stderr)
#             print('Temperature at the edge(s) will be used to cover the range', file=sys.stderr)
#             # Expand
#             if (db_start < db_dts[0]):
#                 db_dts = [db_start] + db_dts
#                 main_t = [main_t[0]] + main_t
#             if (db_dts[-1] < db_end):
#                 db_dts = db_dts + [db_end]
#                 main_t = main_t + [main_t[-1]]

#         starttime_ut = rhdus['READOUT'].data['timestamp']
#         db_uts = [timegm(dt.timetuple()) for dt in db_dts]
#         ## linear interpolate 
#         int_func = scipy.interpolate.interp1d(db_uts, main_t)
#         Troom = int_func(starttime_ut) + 273.15
#     return Troom, db_ret

def convert_asciitime(asciitime, form_fitstime):
    """Ascii time"""
    asciitime = [datetime.strptime('%14.6f' %t, '%Y%m%d%H%M%S.%f') for t in asciitime]
    asciitime = [datetime.strftime(t, form_fitstime) for t in asciitime]
    return np.array(asciitime)

def convert_timestamp(timestamp):
    """Timestamp"""
    timestamp = [datetime.utcfromtimestamp(t) for t in timestamp]
    timestamp = [datetime.strftime(t, FORM_FITSTIME_P) for t in timestamp]
    return np.array(timestamp)

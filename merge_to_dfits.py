"""merge_to_dfits.py: Read logging data and merge them into a FITS object

 Author : Tetsutaro Ueda, Junya Suzuki, Kenichi Karatsu, Tatsuya Takekoshi
 Created: 2017/11/02
 Revision History:
     2018/02/08 - KK - rewrite using class.
     2018/06/08 - TT - apply new calibration method.
     2021         NAITO systems modfied.
"""
from __future__ import print_function
import os
import sys
import yaml
import numpy as np
from astropy.io import fits
from astropy.io import ascii
from astropy import coordinates
from astropy import units
import merge_function as fc
import dfits2dems

# -------------------------------- CONSTANTS
TELESCOP = "ASTE"
D_ASTE = (10.0 * units.m).value  # Diameter  of the ASTE
LON_ASTE = coordinates.Angle("-67d42m11.89525s").deg  # Longitude of the ASTE
LAT_ASTE = coordinates.Angle("-22d58m17.69447s").deg  # Latitude  of the ASTE
FORM_FITSTIME = fc.FORM_FITSTIME  # YYYY-mm-ddTHH:MM:SS
FORM_FITSTIME_P = fc.FORM_FITSTIME_P  # YYYY-mm-ddTHH:MM:SS.ss
DEFAULT_ROOM_T = fc.DEFAULT_ROOM_T  # Kelvin
DEFAULT_AMB_T = fc.DEFAULT_AMB_T  # Kelvin


class MergeToDfits:
    """Read logging data of ASTE and merge them into a FITS object

    Args.:
        ddbfits    (str): File name of the DESHIMA Database FITS.
        dfitsdict  (str): File name of the DFITS Dictionay file (yaml).
        obsinst    (str): File name of the Observation Instruction.
        antennalog (str): File name of the Antennalogging.
        rout_data  (str): File name of the Readout Data (FITS).
        skychoplog (str): File name of the Skychoplogging.
        mistilog   (str): File name of the Mistilogging.
        weatherlog (str): File name of the Weatherlogging.
                          Default is None.
        cabinlog   (str): File name of the cabin temperature log database.
                          Default is None.


    Example:
        "ddbfits"   : DDB_20180208.fits
        "dfitsdict" : dfits_dict.yaml
        "obsinst"   : 20171024012916.obs
        "antennalog": 20171024012916.ant
        "rout_data" : reduced_20171024012916.fits
        "skychoplog": 20171024012916.sky
        "mistilog"  : 20171024012916.misti
        "weatherlog": 20171024012916.wea
        "cabinlog"  : 20171024012916.cabin

        Get HDUList of DFITS:
        >>> mtd = MergeToDfits(
                ddbfits=ddbfits,
                dfitsdict=dfitsdict,
                obsinst=obsinst,
                antennalog=antennalog,
                rout_data=rout_data,
                skychoplog=skychoplog
                mistilog=mistilog
            )
        >>> dfits = mtd.dfits
        >>> dfits.indo()
        Filename: (No file associated with this HDUList)
        No.    Name         Type      Cards   Dimensions   Format
        0  PRIMARY     PrimaryHDU       4   ()
        1  OBSINFO     BinTableHDU     56   1R x 12C   ['K', 'D', 'D', 'D', 'D', 'D', 'D', '63K', '63K', '63K', '63D', '63D']
        2  ANTENNA     BinTableHDU     32   3854R x 8C   ['26A', '4A', 'D', 'D', 'D', 'D', 'D', 'D']
        3  KIDSINFO    BinTableHDU     29   63R x 8C   ['I', 'I', 'E', '2E', '2E', '2E', '2E', '2E']
        4  READOUT     BinTableHDU     23   29234R x 5C   ['26A', 'K', '63D', '63D', '63D']

        Get each HDU:
        >>> mtd = MergeToDfits(
                ddbfits=ddbfits,
                dfitsdict=dfitsdict,
                obsinst=obsinst,
                antennalog=antennalog,
                rout_data=rout_data,
                skychoplog=skychoplog,
                weatherlog=weatherlog,
                cabinlog=cabindb
            )
        >>> obsinfo = mtd.obsinfo
        >>> antenna = mtd.antenna
        >>> readout = mtd.readout
        >>> weather = mtd.weather
    """

    def __init__(
        self,
        ddbfits,
        dfitsdict,
        obsinst,
        antennalog,
        rout_data,
        skychoplog,
        mistilog,
        weatherlog=None,
        cabinlog=None,
    ):
        # -------- Path
        self.ddbfits = os.path.expanduser(ddbfits)
        self.dfitsdict = os.path.expanduser(dfitsdict)
        self.obsinst = os.path.expanduser(obsinst)
        self.antennalog = os.path.expanduser(antennalog)
        self.skychoplog = os.path.expanduser(skychoplog)
        self.mistilog = os.path.expanduser(mistilog)
        self.rout_data = os.path.expanduser(rout_data)
        if weatherlog is None:
            self.weatherlog = None
        else:
            self.weatherlog = os.path.expanduser(weatherlog)
        if cabinlog is None:
            self.cabinlog = None
        else:
            # -------- Define Troom
            self.cabinlog = os.path.expanduser(cabinlog)
            (
                self.cabin_datetimes,
                upper_cabin_temps,
                lower_cabin_temps,
            ) = dfits2dems.retrieve_cabin_temps(self.cabinlog)
            self.upper_cabin_temps = upper_cabin_temps + 273.15
            self.lower_cabin_temps = lower_cabin_temps + 273.15
        self.pixelid = 0
        self.db_ret = 0
        # -------- DFITS Dictionary
        with open(self.dfitsdict, "r") as f:
            self.dfits_dict = yaml.load(f, Loader=yaml.Loader)

        # -------- Read 'antennalog'
        antlog_data = ascii.read(self.antennalog)[:-1]
        self.ant_time = fc.convert_asciitime(antlog_data["time"], FORM_FITSTIME_P)

    @property
    def dfits(self):
        """HDU list of DFITS(Merge to DFITS)"""
        hdus = fits.HDUList()
        hdus.append(fits.PrimaryHDU())  # PRIMARY
        hdus.append(self.obsinfo)  # OBSINFO
        hdus.append(self.antenna)  # ANTENNA
        hdus.append(self.kidsinfo)  # KIDSINFO
        hdus.append(self.readout)  # READOUT: must be before CABIN_T
        hdus.append(self.skychop)
        hdus.append(self.misti)
        if not self.weatherlog is None:
            hdus.append(self.weather)  # WEATHER
        if not self.cabinlog is None:
            hdus.append(self.cabin_t)  # CABIN_T
        return hdus

    @property
    def obsinfo(self):
        """HDU of 'OBSINFO'"""
        # -------- Get the Dicitinary of 'OBSINFO': 'obsinfo_dict'
        od = self.dfits_dict["obsinfo_dict"]
        # -------- Get Header Values
        obsinst_dict = fc.load_obsinst(self.obsinst)
        od["hdr_vals"]["FITSTYPE"] = "DESHIMAv1"
        od["hdr_vals"]["TELESCOP"] = TELESCOP
        od["hdr_vals"]["SITELON"] = LON_ASTE
        od["hdr_vals"]["SITELAT"] = LAT_ASTE
        od["hdr_vals"]["DATE-OBS"] = self.ant_time[0][:19]
        od["hdr_vals"]["OBSERVER"] = obsinst_dict["observer"]
        od["hdr_vals"]["OBJECT"] = obsinst_dict["obs_object"]
        od["hdr_vals"]["RA"] = obsinst_dict["ra"]
        od["hdr_vals"]["DEC"] = obsinst_dict["dec"]
        od["hdr_vals"]["EQUINOX"] = obsinst_dict["equinox"]
        # /*--------------------------- Not confirmed ----------------------------*/
        od["col_vals"]["interval"] = np.array([1.0 / 196])
        od["col_vals"]["integtime"] = np.array([1.0 / 196])
        od["col_vals"]["beamsize"] = np.array([0.005])  # 18 arcsec
        # /*----------------------------------------------------------------------*/
        # -------- Get DDBID and Values for Columns
        ddb = fits.open(self.ddbfits)
        od["hdr_vals"]["DDBID"] = ddb["PRIMARY"].header["DDB_ID"]
        od["col_vals"]["pixelid"] = ddb["KIDRCP"].data["pixelid"]
        od["col_vals"]["offsetaz"] = ddb["KIDRCP"].data["offsetaz"]
        od["col_vals"]["offsetel"] = ddb["KIDRCP"].data["offsetel"]
        od["col_vals"]["gain"] = ddb["KIDRCP"].data["gain"]
        mas_kid_corresp = fc.get_maskid_corresp(self.pixelid, ddb)
        od["col_vals"]["masterids"] = np.array([mas_kid_corresp[0]])
        od["col_vals"]["kidids"] = np.array([mas_kid_corresp[1]])
        od["col_vals"]["kidtypes"] = np.array([mas_kid_corresp[2]])
        od["col_vals"]["kidfreqs"] = np.array([mas_kid_corresp[3]])
        od["col_vals"]["kidQs"] = np.array([mas_kid_corresp[4]])
        ddb.close()
        return fc.create_bintablehdu(od)

    @property
    def antenna(self):
        """HDU of 'ANTENNA'"""
        # -------- Get the Dicitinary of 'ANTENNA': 'antenna_dict'
        ad = self.dfits_dict["antenna_dict"]
        # -------- Get Header Values
        ad["hdr_vals"]["FILENAME"] = os.path.basename(self.antennalog)
        # -------- Read 'antennalog'
        antlog_data = ascii.read(self.antennalog)[:-1]
        # -------- Get Values for Columns
        ad["col_vals"]["time"] = self.ant_time
        ad["col_vals"]["scantype"] = antlog_data["type"]
        try:
            ad["col_vals"]["az"] = (
                antlog_data["az-prg(no-cor)"]
                + antlog_data["az-real"]
                - antlog_data["az-prg"]
            )
            ad["col_vals"]["el"] = (
                antlog_data["el-prog(no-cor)"]
                + antlog_data["el-real"]
                - antlog_data["el-prg"]
            )
        except:
            ad["col_vals"]["az"] = (
                antlog_data["az-prg(no-col)"]
                + antlog_data["az-real"]
                - antlog_data["az-prg"]
            )
            ad["col_vals"]["el"] = (
                antlog_data["el-prog(no-col)"]
                + antlog_data["el-real"]
                - antlog_data["el-prg"]
            )
        ad["col_vals"]["ra"] = antlog_data["ra-prg"]
        ad["col_vals"]["dec"] = antlog_data["dec-prg"]
        ad["col_vals"]["az_center"] = antlog_data["az-prog(center)"]
        ad["col_vals"]["el_center"] = antlog_data["el-prog(center)"]
        return fc.create_bintablehdu(ad)

    @property
    def readout(self):
        """HDU of 'READOUT'"""
        # -------- Get the Dicitinary of 'READOUT': 'readout_dict'
        rd = self.dfits_dict["readout_dict"]

        # -------- Get Header Values
        rd["hdr_vals"]["FILENAME"] = os.path.basename(self.rout_data)

        # -------- Open 'DDB' and 'rout_data'
        ddb = fits.open(self.ddbfits)
        rhdus = fits.open(self.rout_data)

        # -------- Define Tamb
        if self.weatherlog:
            wlog_data = ascii.read(self.weatherlog)
            Tamb = np.nanmean(wlog_data["tmperature"]) + 273.15
        else:
            print("Weatherlog DB is not specified.", file=sys.stderr)
            print("Tamb", DEFAULT_AMB_T, "K will be used", file=sys.stderr)
            Tamb = DEFAULT_AMB_T

        # -------- Get Values for Columns
        nkid = rhdus["READOUT"].header["NKID%d" % self.pixelid]
        reduce_data = np.transpose(
            [rhdus["READOUT"].data["Amp, Ph, linPh %d" % i].T for i in range(nkid)]
        )

        rd["col_vals"]["starttime"] = fc.convert_timestamp(
            rhdus["READOUT"].data["timestamp"]
        )
        rd["col_vals"]["pixelid"] = rhdus["READOUT"].data["pixelid"]
        rd["col_vals"]["lin_phase"] = reduce_data[:, 2]

        # (Comment for dfits 20180703: Troom is fixed to Troom[0], because we should not change Troom from base temperature (but Tamb should be in future...)
        Tsignal = fc.calibrate_to_power(
            self.pixelid, self.lower_cabin_temps[0], Tamb, rhdus, ddb
        )
        rd["col_vals"]["Tsignal"] = Tsignal
        ddb.close()
        rhdus.close()
        return fc.create_bintablehdu(rd)

    @property
    def kidsinfo(self):
        """HDU of 'KIDSINFO'"""
        self.kidsinfo_hdus = fits.open(self.rout_data)
        return self.kidsinfo_hdus["KIDSINFO"]

    @property
    def weather(self):
        """HDU of 'WEATHER'"""
        # -------- Error Handling: Case of 'weatherlog' is None
        if self.weatherlog is None:
            raise ValueError('No "weatherlog" is inputed!!')
        # -------- Get the Dicitinary of 'WEATHER': 'weather_dict'
        wd = self.dfits_dict["weather_dict"]
        # -------- Get Header Values
        wd["hdr_vals"]["FILENAME"] = os.path.basename(self.weatherlog)
        # -------- Read 'weatherlog'
        wlog_data = ascii.read(self.weatherlog)
        # -------- Get Values for Columns
        wd["col_vals"]["time"] = fc.convert_asciitime(wlog_data["time"], FORM_FITSTIME)
        wd["col_vals"]["temperature"] = wlog_data["tmperature"]
        wd["col_vals"]["pressure"] = wlog_data["presure"]
        wd["col_vals"]["vapor-pressure"] = wlog_data["vapor-pressure"]
        wd["col_vals"]["windspd"] = wlog_data["aux1"]
        wd["col_vals"]["winddir"] = wlog_data["aux2"]
        return fc.create_bintablehdu(wd)

    @property
    def cabin_t(self):
        """HDU of 'CABIN_T'"""
        if self.cabinlog is None:
            raise ValueError('No "cabinlog" is inputed!!')
        # -------- Get the Dicitinary of 'CABIN_T': 'cabin_t_dict'
        cabin_t_dict = self.dfits_dict["cabin_t_dict"]
        # -------- Set Data to the Dictinary 'cabin_t_dict'
        cabin_t_dict["hdr_vals"]["FILENAME"] = os.path.basename(self.cabinlog)
        cabin_t_dict["col_vals"]["time"] = self.cabin_datetimes
        cabin_t_dict["col_vals"]["upper_cabin"] = self.upper_cabin_temps
        cabin_t_dict["col_vals"]["main_cabin"] = self.lower_cabin_temps
        return fc.create_bintablehdu(cabin_t_dict)

    @property
    def skychop(self):
        skychop_dict = self.dfits_dict["skychop_dict"]
        datetimes, states = dfits2dems.retrieve_skychop_states(self.skychoplog)
        # -------- Set Data to the Dictinary 'cabin_t_dict'
        skychop_dict["hdr_vals"]["FILENAME"] = os.path.basename(self.skychoplog)
        skychop_dict["col_vals"]["time"] = datetimes
        skychop_dict["col_vals"]["state"] = states
        return fc.create_bintablehdu(skychop_dict)

    @property
    def misti(self):
        md = self.dfits_dict["misti_dict"]
        datetimes, az, el = dfits2dems.retrieve_misti_log(self.mistilog)
        # -------- Set Data to the Dictinary 'cabin_t_dict'
        md["hdr_vals"]["FILENAME"] = os.path.basename(self.mistilog)
        md["col_vals"]["time"] = datetimes
        md["col_vals"]["az"] = az
        md["col_vals"]["el"] = el
        return fc.create_bintablehdu(md)

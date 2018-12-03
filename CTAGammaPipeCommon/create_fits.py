# ==========================================================================
# Import fits file in database
#
# Copyright (C) 2018 Nicolo' Parmiggiani, Andrea Bulgarelli
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==========================================================================


import sys
import mysql.connector as mysql
#from conf import get_pipedb_conf
#from conf import get_evtdb_conf
from GammaPipeCommon.utility import *
import numpy as np
from numpy import rec
#from pyfits import Column
import os
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.io import fits

# create FITS with event in time window from observation
def write_fits(tstart_tt,tstop_tt,observationid,datarepositoryid,path_base_fits,tref_mjd,obs_ra,obs_dec,emin,emax,fov,instrumentname):

    tstart = float(tstart_tt)
    tstop = float(tstop_tt)
    emin = str(emin)
    emax = str(emax)
    obs_ra = str(obs_ra)
    obs_dec = str(obs_dec)
    fov = str(fov)

    #connect to database
    conf_dictionary = get_pipedb_conf()

    pipedb_hostname = conf_dictionary['host']
    pipedb_username = conf_dictionary['username']
    pipedb_password = conf_dictionary['password']
    pipedb_port = conf_dictionary['port']
    pipedb_database = conf_dictionary['database']

    #connect to database
    conf_dictionary = get_evtdb_conf()

    evtdb_hostname = conf_dictionary['host']
    evtdb_username = conf_dictionary['username']
    evtdb_password = conf_dictionary['password']
    evtdb_port = conf_dictionary['port']
    evtdb_database = conf_dictionary['database']

    # get events list
    conn = mysql.connect(host=evtdb_hostname, user=evtdb_username, passwd=evtdb_password, db=evtdb_database)
    cursor = conn.cursor(dictionary=True)

    query = "SELECT * FROM streaming_evt WHERE timerealtt > "+str(tstart)+" AND timerealtt < "+str(tstop)+" AND observationid = "+str(observationid)+" AND datarepositoryid = "+str(datarepositoryid)
    print(query)
    cursor.execute(query)

    events = cursor.fetchall()
    print(len(events))

    cursor.close()
    conn.close()


    hdulist = fits.open(path_base_fits)

    primary_hdu = hdulist[0]

    events_hdu = hdulist[1]

    tref_mjd = float(tref_mjd)

    for x in events:
        time_real_mjd = Utility.convert_tt_to_mjd(x['timerealtt'])
        #print(time_real_mjd)
        #print(tref_mjd)
        time_real_seconds = str((float(time_real_mjd)-float(tref_mjd))*86400)
        #print(time_real_seconds)
        x['timerealtt'] = time_real_seconds

    # CREATE EVENTS data table HDU

    c_e_1 = fits.Column(name = 'EVENT_ID', format = '1J', bscale = 1, bzero = 2147483648, array=np.array([x['evtid'] for x in events]))
    c_e_2 = fits.Column(name = 'TIME',format = '1D', unit = 'day', array=np.array([x['timerealtt'] for x in events]))
    c_e_3 = fits.Column(name = 'RA',format = '1E', unit = 'deg', array=np.array([x['ra_deg'] for x in events]))
    c_e_4 = fits.Column(name = 'DEC', format = '1E', unit = 'deg', array=np.array([x['dec_deg'] for x in events]))
    c_e_5 = fits.Column(name = 'ENERGY', format = '1E', unit = 'TeV', array=np.array([x['energy'] for x in events]))
    c_e_6 = fits.Column(name = 'DETX', format = '1E', unit = 'deg', array=np.array([x['detx'] for x in events]))
    c_e_7 = fits.Column( name = 'DETY', format = '1E', unit = 'deg', array=np.array([x['dety'] for x in events]))
    c_e_8 = fits.Column(name = 'MC_ID', format = '1J', array=np.array([x['mcid'] for x in events]))

    coldefs = fits.ColDefs([c_e_1, c_e_2,c_e_3,c_e_4,c_e_5,c_e_6,c_e_7,c_e_8])

    data_tbhdu = fits.BinTableHDU.from_columns(coldefs)
    data_tbhdu.header = hdulist[1].header



    #change header content
    data_tbhdu.header['NAXIS2'] = len(events)
    data_tbhdu.header['DSVAL2'] = emin+":"+emax
    data_tbhdu.header['DSVAL3'] = "CIRCLE("+obs_ra+","+obs_dec+","+fov+")"
    data_tbhdu.header['MMN00001'] = "None"
    data_tbhdu.header['MMN00002'] = "None"
    data_tbhdu.header['TELESCOP'] = str(instrumentname)
    data_tbhdu.header['OBS_ID'] = str(observationid)

    # # time_second/86400 + 51544.5 = time_mjd
    # tstart_mjd = str(float(tstart)/86400 + 51544.5)
    # tstop_mjd = str(float(tstop)/86400 + 51544.5)
    # print tstart_mjd
    # print tstop_mjd
    # # UTC timescale
    # tstart_astropy = Time(tstart_mjd, format='mjd')
    # tstop_astropy = Time(tstop_mjd, format='mjd')
    #
    # #convert to ISO
    # tstart_astropy = tstart_astropy.iso
    # tstop_astropy = tstop_astropy.iso

    #date time
    data_tbhdu.header['DATE_OBS'] = ""
    data_tbhdu.header['TIME_OBS'] = ""
    data_tbhdu.header['DATE_END'] = ""
    data_tbhdu.header['TIME_END'] = ""
    data_tbhdu.header['TSTART'] = str(tstart)
    data_tbhdu.header['TSTOP'] = str(tstop)

    int_tref = int(tref_mjd)
    refint,refdecimal = int(tref_mjd),float(tref_mjd)-int(tref_mjd)

    data_tbhdu.header['MJDREFI'] = str(refint)
    data_tbhdu.header['MJDREFF'] = str(refdecimal)

    data_tbhdu.header['TELAPSE'] = str(tstop-tstart)
    data_tbhdu.header['ONTIME'] = str(tstop-tstart)
    data_tbhdu.header['LIVETIME'] = str(tstop-tstart)
    data_tbhdu.header['DEADC'] = "1"
    data_tbhdu.header['TIMEDEL'] = "1"

    data_tbhdu.header['RA_PNT'] = obs_ra
    data_tbhdu.header['DEC_PNT'] = obs_dec

    # CREATE GTI data table HDU

    time_start_mjd = Utility.convert_tt_to_mjd(tstart)
    time_stop_mjd = Utility.convert_tt_to_mjd(tstop)

    time_start_second = str((float(time_start_mjd)-float(tref_mjd))*86400)
    time_stop_second = str((float(time_stop_mjd)-float(tref_mjd))*86400)

    gti_tstart = np.array([time_start_second])
    gti_tstop = np.array([time_stop_second])
    c1 = fits.Column(name='START', format='1D', array=gti_tstart)
    c2 = fits.Column(name='STOP', format='1D', array=gti_tstop)
    coldefs = fits.ColDefs([c1, c2])

    gti_tbhdu = fits.BinTableHDU.from_columns(coldefs)
    gti_tbhdu.header = hdulist[2].header

    thdulist = fits.HDUList([primary_hdu,data_tbhdu,gti_tbhdu])

    filename = "obs_"+str(observationid)+"_"+str(tstart)+"_"+str(tstop)+".fits"
    if os.path.exists(filename):
        os.unlink(filename)
    thdulist.writeto(filename)

    hdulist_new = fits.open(filename)

    #print hdulist_new[2].header

    hdulist_new.close()

    return filename


if __name__ == '__main__':

    # crete the XML for the specific observation
    tstart_tt = sys.argv[1]
    tstop_tt = sys.argv[2]
    observationid = sys.argv[3]
    tref_mjd = sys.argv[4]
    obs_ra = sys.argv[5]
    obs_dec = sys.argv[6]
    emin = sys.argv[7]
    emax = sys.argv[8]
    fov = sys.argv[9]
    instrumentname = sys.argv[10]
    path_base_fits = "templates/base_empty.fits"

    filename = write_fits(tstart_tt,tstop_tt,observationid,path_base_fits,tref_mjd,obs_ra,obs_dec,emin,emax,fov,instrumentname)

    print(filename)

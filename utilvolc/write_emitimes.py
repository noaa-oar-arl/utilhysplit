# write_emitimes.py
# Writes a HYSPLIT EMITIMES.txt file for a cylindrical volcanic source
# and for inserting volcat data into hysplit
from utilhysplit import cylsource
from utilhysplit import emitimes
from datetime import datetime
import xarray as xr
import numpy as np
import numpy.ma as ma
from utilvolc import volcat
from os import path

"""
This script contains functions (class) that writes emitimes files
for cylinder source and volcat data insertion hysplit runs.
--------------
Functions:
--------------
write_cyl_file: writes emitimes file to working directory
Class: InsertVolcat

--------------
"""


def write_cyl_file(wdir, date_time, lat, lon, volcname, radius, dr, duration, pollnum, pollpercents, altitude, umbrella):
    """ Write EMITIMES.txt file for cylinder source hysplit runs
    Inputs:
    wdir: working directory (string)
    date_time: date and time of emission start (datetime object)
    lat: latitude for center of cylinder (float)
    lon: longitude for center of cylinder (float)
    volcname: name of volcano (string)
    radius: radius of cylinder around volcano vent in meters (integer)
    dr: will determine number of concentric circles in column,
    must be evenly divisible into radius (integer)
    duration: hour and minute (HHMM) of emission (string)
    pollnum: number of particle size bins (integer)
    pollpercents: percentage of particles for each size bin,
    represented as values from [0] to [1] (list)
    altitude: height of column in meters [bottom, top] (list)
    umbrella: 1 - uniform column (integer)

    Outputs:
    filename: location of newly written emitimes file (string)
    """

    dt = date_time
    fname = 'EMIT_'+volcname+'_cyl_'+dt.strftime("%Y%m%d%H%M")+'_'+duration+'hrs.txt'

    filename = wdir+fname

    # Creating emitimes files
    efile = cylsource.EmitCyl(filename=filename)
    latlist, lonlist = efile.calc_cyl(lat, lon, radius, dr)
    nrecords = efile.calc_nrecs(latlist, pollnum=pollnum, umbrella=umbrella)
    efile.write_data(dt, latlist, lonlist, nrecords, pollnum, duration=duration,
                     pollpercents=pollpercents, height=altitude)

    print('File '+filename+' created')

    return filename


class InsertVolcat:

    def __init__(self, wdir, 
                 vdir, date_time, duration,
                 pollpercents=1,
                 pollnum=1,  
                 vname=None):
        """
        Class of tools for inserting volcat data into hysplit
        -------------
        Inputs:
        wdir: working directory - where emitimes file will be located (string)
        vdir: volcat directory - where volcat data files are located (string)
        date_time: date and time of volcat file to use (datetime object)
        pollnum: number of particle size bins (integer)
        pollpercents: percentage of particles for each size bin,
        represented as values from [0] to [1] (list) 
        duration: hour and minute (HHMM) of emission (string)
        --------------
        """

        if wdir[-1] != "/":
            wdir += "/"
        self.wdir = wdir
        if vdir[-1] != "/":
            vdir += "/"
        self.vdir = vdir
        self.date_time = date_time
        self.pollnum = pollnum
        self.pollpercents = pollpercents
        self.duration = duration
        self.vname = vname

    def add_vname(self, vname):
        self.vname = vname

    def find_fname(self):
        """ Determines filename for volcat file"""
        volcfile = 'SCOPE_NWC_ASH-LS-ASH_PRODUCTS-HIMAWARI8_NOAA-RAIKOKE-' + \
            self.date_time.strftime('%Y%m%d-%H%M%S')+'-fv2.nc'
        fname = self.vdir+volcfile
        return fname

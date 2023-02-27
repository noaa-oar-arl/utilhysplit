"""
This script contains functions (class) that writes emitimes files
for cylinder source and volcat data insertion hysplit runs.
--------------
Functions:
--------------
find_emit
write_cyl_file: writes emitimes file to working directory

Class: InsertVolcat
     add_vname: adds volcano name to instance variables
     find_fnames: finds list of volcat files based on datetime object 
     (can also use volcano ID)
     make_match: makes a string containing event datetime, image datetime, and 
     volcano id which is used for generating emitimes filenames and area filenames
     get_area: calculates area of lat/lon grid
     make_1D: makes 1D arrays of lat, lon, ash height, ash mass, and area
     write_emit: writes emitimes files
Class: EmitName 

--------------
Change log
"""

#import os
#from datetime import datetime

#import numpy as np
#import numpy.ma as ma
#import pandas as pd
#import xarray as xr
#from utilhysplit import emitimes

#from utilvolc import volcat
#from utilvolc.volcat import VolcatName
from utilhysplit import cylsource


def write_cyl_file(
    wdir,
    date_time,
    lat,
    lon,
    volcname,
    radius,
    dr,
    duration,
    pollnum,
    pollpercents,
    altitude,
    umbrella,
):
    """Write EMITIMES.txt file for cylinder source hysplit runs
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
    fname = (
        "EMIT_"
        + volcname
        + "_cyl_"
        + dt.strftime("%Y%m%d%H%M")
        + "_"
        + duration
        + "hrs.txt"
    )

    filename = wdir + fname

    # Creating emitimes files
    efile = cylsource.EmitCyl(filename=filename)
    latlist, lonlist = efile.calc_cyl(lat, lon, radius, dr)
    nrecords = efile.calc_nrecs(latlist, pollnum=pollnum, umbrella=umbrella)
    efile.write_data(
        dt,
        latlist,
        lonlist,
        nrecords,
        pollnum,
        duration=duration,
        pollpercents=pollpercents,
        height=altitude,
    )

    print("File " + filename + " created")

    return filename

# write_cylsource_emitimes.py
# Writes a HYSPLIT EMITIMES.txt file for a cylindrical volcanic source
from utilhysplit import cylsource
from utilhysplit import emitimes
from datetime import datetime

"""
This script contains a function that writes emitimes files
for cylinder source hysplit runs.
--------------
Functions:
--------------
write_cyl_file: writes emitimes file to working directory
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
    efile = cylsource.EmitTimes(filename=filename)
    latlist, lonlist = efile.calc_cyl(lat, lon, radius, dr)
    nrecords = efile.calc_nrecs(latlist, pollnum=pollnum, umbrella=umbrella)
    efile.write_data(dt, latlist, lonlist, nrecords, pollnum, duration=duration,
                     pollpercents=pollpercents, height=altitude)

    print('File '+filename+' created')

    return filename

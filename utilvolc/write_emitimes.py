# write_emitimes.py
# Writes a HYSPLIT EMITIMES.txt file for a cylindrical volcanic source
# and for inserting volcat data into hysplit
from utilhysplit import cylsource
from utilhysplit import emitimes
from datetime import datetime
from glob import glob
import xarray as xr
import numpy as np
import numpy.ma as ma
from utilvolc import volcat
from os import path
from math import pi, cos

"""
This script contains functions (class) that writes emitimes files
for cylinder source and volcat data insertion hysplit runs.
--------------
Functions:
--------------
write_cyl_file: writes emitimes file to working directory
Class: InsertVolcat
     add_vname: adds volcano name to instance variables
     find_match: finds the indentifying string of date, time, volcano ID
     find_fname: finds volcat file based on directory and datetime object
     (can also use volcano ID)
     get_area: calculates area of lat/lon grid
     make_1D: makes 1D arrays of lat, lon, ash height, ash mass, and area
     write_emit: writes emitimes files
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

    def __init__(self, wdir, vdir, date_time,
                 duration='0010',
                 pollpercents=[1],
                 pollnum=1,
                 vname=None,
                 vid=None):
        """
        Class of tools for inserting volcat data into hysplit
        -------------
        Inputs:
        wdir: working directory - where emitimes file will be located (string)
        vdir: volcat directory - where volcat data files are located (string)
        date_time: date and time of volcat file to use (datetime object)
        duration: hour and minute (HHMM) of emission - default is '0010' (string)
        pollnum: number of particle size bins - default is 1 (integer)
        pollpercents: percentage of particles for each size bin,
        represented as values from [0] to [1] - default is [1] (list)
        vname: volcano name - default is None (string)
        vid: volcano id - default is None (string)
        Outputs:
        --------------
        Functions:
        add_vname: adds volcano name
        find_match: finds string with date_time and vid to match
        find_fname: finds volcat filename
        get_area: calculates the domain area, gridded
        make_1D: makes 1D array of lat, lon, ash height, ash mass, area
        write_emit: writes emitimes files
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
        self.vid = vid
        self.find_fname()

    def add_vname(self, vname):
        """ Adds volcano name"""
        self.vname = vname

    def find_match(self):
        """Determines matching string to identify file"""
        if self.vid != None:
            match = self.date_time.strftime('%Y%j_%H%M%S_v')+self.vid
        else:
            match = self.date_time.strftime('%Y%j_%H%M%S')
        return match

    def find_fname(self):
        """ Determines filename for volcat file based on vdir, date_time and vid"""
        vfiles = '*.nc'
        volclist = glob(self.vdir+vfiles)
        match = self.find_match()
        fname = [f for f in volclist if match in f]
        self.fname = fname[0]
        return self.fname

    def get_area(self, write=False, correct_parallax=True, clip=False):
        """Calculates the area (km^2) of each volcat grid cell
        Converts degress to meters using a radius of 6378.137km.
        Input:
        write: boolean (default: False) Write area to file
        correct_parallax: boolean (default: True) Use parallax correct lat/lon values
        clip: boolean (default: False) Use clipped array around data, reduces domain
        output:
        area: xarray containing gridded area values
        """
        d2r = pi / 180.0  # convert degress to radians
        d2km = 6378.137 * d2r  # convert degree latitude to kilometers

        if self.fname:
            dset = volcat.open_dataset(self.fname, correct_parallax=correct_parallax)
        else:
            print('ERROR: Need volcat filename!')

        # Extracts ash mass array (two methods - one is smaller domain around feature)
        # Removes time dimension
        if clip == True:
            mass = volcat.get_mass(dset)[0, :, :]
        else:
            mass = dset.ash_mass_loading[0, :, :]

        lat = mass.latitude
        lon = mass.longitude
        latrad = lat * d2r  # Creating latitude array in radians
        coslat = np.cos(latrad) * d2km * d2km
        shape = np.shape(mass)

        # Make shifted lat and shifted lon arrays to use for calculations
        lat_shift = lat[1:, :].values
        lon_shift = lon[:, 1:].values
        # Adding row/column of nans to shifted lat/lon arrays
        to_add_lon = np.empty([shape[0]]) * np.nan
        to_add_lat = np.empty([shape[1]]) * np.nan
        # Back to xarray for calculations
        lat2 = xr.DataArray(np.vstack((lat_shift, to_add_lat)), dims=['y', 'x'])
        lon2 = xr.DataArray(np.column_stack((lon_shift, to_add_lon)), dims=['y', 'x'])

        # area calculation
        area = abs(lat-lat2) * abs(abs(lon)-abs(lon2)) * coslat
        area.name = 'area'
        area.attrs['long_name'] = 'area of each lat/lon grid box'
        area.attrs['units'] = 'km^2'
        # Reformatting array attributes
        if write == True:
            directory = self.vdir+'Area/'
            match = self.find_match()
            if correct_parallax == True:
                areafname = 'area_'+match+'_pc.nc'
            else:
                areafname = 'area_'+match+'.nc'
            print(directory+areafname)
            area.to_netcdf(directory+areafname)
        return area

    def make_1D(self, correct_parallax=True):
        """ Makes compressed 1D arrays of latitude, longitude, ash height,
        mass emission rate, and area
        For use in writing emitimes files. See self.write_emit()
        Input:
        area: xarray dataset of area between lat/lon grid (from netcdf files)
        correct_parallax: use parallax corrected lat/lon (default = True)
        Output:
        lat: 1D array of latitude
        lon: 1D array of longitude
        hgt: 1D array of ash top height
        mass: 1D array of ash mass
        area: 1D array of area
        """
        import numpy.ma as ma

        if self.fname:
            dset = volcat.open_dataset(self.fname, correct_parallax=correct_parallax)
        else:
            print('ERROR: Need volcat filename!')

        # Extracts ash mass array - Removes time dimension
        mass0 = dset.ash_mass_loading[0, :, :]
        height0 = dset.ash_cloud_height[0, :, :]
        # Making 0. in arrays nans
        mass = mass0.where(mass0 > 0.)
        height = height0.where(height0 > 0.)

        # Finds and area files. If no file detected, produces warning
        areadir = self.vdir+'Area/'
        match = self.find_match()
        if correct_parallax == True:
            arealist = glob(areadir+'*_pc.nc')
        else:
            arealist = glob(areadir + '*0.nc')
        areafile = [f for f in arealist if match in f]
        if areafile:
            areafile = areafile[0]
            areaf = xr.open_dataset(areafile).area
        else:
            print('No area file detected! Please use self.get_area()')

        # Calculating mass - rate is (mass/hr) in HYSPLIT
        # area needs to be in m^2 not km^2!
        ashmass = mass * areaf * 1000. * 1000.
        latitude = mass.latitude
        longitude = mass.longitude

        # Creating compressed 1-d arrays (removing nans) to prepare for file writing
        hgt_nan = ma.compressed(height)  # arra with nan values to remove
        hgt = hgt_nan[~np.isnan(hgt_nan)]  # removing nan values
        mass = ma.compressed(ashmass)[~np.isnan(hgt_nan)]  # removing nan values
        lat = ma.compressed(latitude)[~np.isnan(hgt_nan)]  # removing nan values
        lon = ma.compressed(longitude)[~np.isnan(hgt_nan)]  # removing nan values
        area = ma.compressed(areaf)[~np.isnan(hgt_nan)]  # removing nan values
        area = area * 1000. * 1000  # Moving area from (km^2) to (m^2)
        hgt = hgt * 1000.  # Moving hgt from (km) to (m)

        return lat, lon, hgt, mass, area

    def write_emit(self, correct_parallax=True, heat='0.00e+00'):
        """ Writes emitimes file from volcat data.
        Inputs are created using self.make_1D()
        Uses instance variables: date_time, duration, par,
        Inputs:
        lat: 1D array of latitude
        lon: 1D array of longitude
        hgt: 1D array of ash top height
        mass: 1D array of ash mass
        area: 1D array of area
        correct_parallax: (boolean) use parallax corrected lat lon values
        heat: (string) default=0.00e+00
        Output:
        emitimes file located in wdir
        """
        # Call make_1D() array to get lat, lon, height, mass, and area arrays
        lat, lon, hgt, mass, area = self.make_1D(correct_parallax=correct_parallax)
        # Need to write some code to handle different durations and
        # resulting constants for mass rate (g/hr) determination
        if self.duration == '0010':
            const = 6
        match = self.find_match()
        filename = 'VOLCAT_'+match+'_par'+str(self.pollnum)
        f = open(self.wdir + filename, 'w')
        f.write('YYYY MM DD HH        DURATION(HHMM) #RECORDS \n')
        f.write('YYYY MM DD HH MM DURATION(HHMM) LAT LON HGT(m) RATE(g/hr) AREA(m2) HEAT(w) \n')
        f.write('{:%Y %m %d %H  } {} {}\n'.format(self.date_time,
                                                  self.duration, np.shape(hgt)[0] * self.pollnum))
        h = 0
        while h < len(hgt):
            i = 0
            while i < len(self.pollpercents):
                f.write('{:%Y %m %d %H %M} {} {:9.6f} {:10.6f} {:8.2f} {:.2E} {:.2E} {} \n'.format(
                    self.date_time, self.duration, lat[h], lon[h], hgt[h], mass[h]*const*self.pollpercents[i], area[h], heat))
                i += 1
            h += 1
        f.close()
        return('Emitimes file written: '+self.wdir+filename)

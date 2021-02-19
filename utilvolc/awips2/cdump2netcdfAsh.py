import datetime
import os
import sys
import xarray as xr
import numpy as np
import logging
import pandas as pd
from netCDF4 import Dataset
import hysplit
# 01/28/2020 AMC cdump2awips created to make a netcdf file appropriate for input into AWIPS
# hysplit.py was modified in the forked version of MONET to make this work.

# combine_cdump creates a 6 dimensional xarray dataarray object from cdump files.
#
# Notes 
# global attribute time_origin=2020-02-04 16:20:47
# ensid variable is an integer.
# "ensemble_tag" variable is a string
# dimension order is time, ensid, latitude, longitude
# one time period per file.

# AWIPS2 expects netcdf files to have only one time period per file.
# It can read a zipped file consisting of multiple netcdf files, one for each time period.

logger = logging.getLogger(__name__)

def mass_loading(xrash):
    # returns mass loading in each level.
    # return hysplit.calc_aml(xrash)
    return hysplit.hysp_massload(xrash)

def meters2FL(meters):
    flight_level = meters * 3.28084 / 100.0
    return int(np.ceil(flight_level / 10.0) * 10)

class Cdump2Awips:
    def __init__(
        self,
        xrash,
        outname,
        munit="unit",
        fileformat="NETCDF4",
        jobid="unknown",
        globalhash={},
    ):
        """
        xrash1 : xarray data-array output by combine_datset of hysplit.py in monetio
                 module.
        outname : str
        munit : str
        fileformat : str
        source_file_description : str
        source_location: tuple (float,float) (latitude,longitude)

        globalhash : dictionary used in add_global_attributes method.

        OUTPUT
        netcdf files for AWIPS2.
        coordinates are latitude, longitude, time, ensid.
        Although time is a variable only one time per file.
        Each level containing concentrations is a variable.
     
        """
        # self.dt = dt
        self.add_probs = False 
          
        self.dfmt = "%Y-%m-%d %H:%M:%S"
        self.outname = outname
        self.munit = munit  # mass unit for concentration
        self.massunit = munit  # mass unit for column mass loading.
        self.fileformat = fileformat
        self.jobid = jobid
        # self.coordlist = ('time', 'ensid', 'longitude', 'latitude')
        self.coordlist = ("time", "ensid", "latitude", "longitude")

        self.globalhash = globalhash

        if "sample time hours" in xrash.attrs:
            dt = xrash.attrs["sample time hours"]
        else:
            dt = 1

        # mass loading should be in g/m2 to compare to satellite.
        # concentration should be in mg/m3 to compare to threshold levels.
        self.sample_time = np.timedelta64(int(dt), "h")
        # stack the ensemble and source dimensions so it is one dimension
        self.xrash = xrash.stack(ensemble=("ens", "source"))
        # put dimensionsin correct order.
        self.xrash = self.xrash.transpose(
            "time", "ensemble", "x", "y", "z", transpose_coords=False
        )
        # array with mass loading rather than concentration.
        mass = mass_loading(xrash)
        # if concentration  unit is mg then make mass loading unit in g.
        if self.munit == "mg":
            mass = mass / 1000.0
            self.massunit = "g"
        self.mass = mass.stack(ensemble=("ens", "source"))
        self.mass = self.mass.transpose(
            "time", "ensemble", "x", "y", transpose_coords=False
        )
        self.levelra = xrash.z.values
        self.nra = xrash.values

    def create_all_files(self):
        """
        returns list of filenames.
        """
        flist = [
            self.make_dataset(x) for x in np.arange(0, len(self.xrash.time.values))
        ]
        return flist
        # for iii, tm in enumerate(self.xrash.time.values):
        #    self.filenamelist.append(self.make_dataset(iii))

    def add_global_attributes(self, fid):
        # GLOBAL ATTRIBUTES
        #logger.debug("Adding global attributes")
        fid.jobid = self.jobid
        fid.time_origin = datetime.datetime.now().strftime(self.dfmt)
        fid.concentration_mass_unit = self.munit
        fid.massloading_mass_unit = self.massunit

        if "MER" in self.globalhash.keys():
            fid.MER = self.globalhash["MER"]
        if "MER_unit" in self.globalhash.keys():
            fid.MER_unit = self.globalhash["MER_unit"]
        if "source_latitude" in self.globalhash.keys():
            fid.soure_latitude = self.globalhash["source_latitude"]
        if "source_longitude" in self.globalhash.keys():
            fid.soure_longitude = self.globalhash["source_longitude"]
        if "source_name" in self.globalhash.keys():
            fid.source_name = self.globalhash["source_name"]
        if "emission_start" in self.globalhash.keys():
            fid.emission_start = self.globalhash["emission_start"].strftime(self.dfmt)
        if "emission_duration_hours" in self.globalhash.keys():
            fid.emission_duration_hours = self.globalhash["emission_duration_hours"]

        return fid

    def make_conc_level(self, fid, variable_name, min_level, max_level):
        coordlist = self.coordlist
        concid = fid.createVariable(variable_name, "f4", coordlist)
        concid.units = self.munit + "/m3"
        concid.long_name = "Concentration Array"
        concid.bottomlevel = min_level
        concid.toplevel = max_level
        return concid

    def make_dataset(self, iii):
        """
        Inputs:
           iii : integer
        Returns :
           outname : str

        uses attributes: 
            self.xrash
            self.mass
            self.sample_time
            self.munit
            self.massunit
            self.coordlist
            self.add_probs
 
        uses methods:
             self.add_global_attributes         
             self.make_conc_level

        uses function:
             makeconc
        """
        xrash = self.xrash.copy()
        munit = self.munit
        sample_time = self.sample_time
        date1 = xrash.time[iii].values
        datestr = pd.to_datetime(date1).strftime("%Y%m%d_%H%MUTC")
        # outname = '{}_{}_t{:02d}.nc'.format(self.outname, datestr,iii)
        outname = "{}_{}.nc".format(self.outname, iii)
        fid = Dataset(outname, "w", format="NETCDF4")
        fid = self.add_global_attributes(fid)

        # DEFINE DIMENSIONS
        lat_shape = xrash.shape[3]
        lon_shape = xrash.shape[2]
        lat = fid.createDimension("latitude", lat_shape)
        lon = fid.createDimension("longitude", lon_shape)
        # level = fid.createDimension('levels',len(levelra))

        clevs = [0.02,0.2, 2, 5, 10, 100]
        clevels = fid.createDimension("contour_levels", len(clevs))
        ens_shape = xrash.coords["ensemble"].shape[0]
        # add ensemble mean, ensemble standard deviation,
        # probability of exceedances for 0.02, 0.2, 2, 5, 10 mg/m3.
        if self.add_probs: 
           ens_shape += 7
        ensemble = fid.createDimension("ensid", ens_shape)

        time = fid.createDimension("time", 1)  # one time per file
        bnds = fid.createDimension("bnds", 2)  # two bounds per time.

        # origin = fid.createDimension('origins',hxr.attrs['Number Start Locations'])
        origin = fid.createDimension("origins", 1)
        # Scalar variables

        # latra, lonra = hf.getlatlon(hxr)
        latra = xrash.latitude[:, 0]
        lonra = xrash.longitude[0]

        # DEFINE A DIFFERENT VARIABLE for EACH LEVEL.
        # DIMENSIONS are ensemble tag, latitude, longitude
        levs = xrash.z.values
        # lev1, lev2, lev3 = handle_levels(levs)

        # concid_list = [self.make_conc_level(x) for x in xrash.z.values]
        concid_list = []
        for jjj, level in enumerate(xrash.z.values):
            maxlev = "FL{}".format(meters2FL(level))
            if jjj > 0:
                minlev = "FL{}".format(meters2FL(xrash.z.values[jjj - 1]))
            else:
                minlev = "SFC"
            varname = maxlev
            concid_list.append(self.make_conc_level(fid, varname, minlev, maxlev))

        massid = fid.createVariable("MassLoading", "f4", self.coordlist)
        massid.units = self.massunit + "/m2"
        massid.long_name = "Mass Loading from SFC to " + maxlev

        # Standard Contour levels for concentration in mg/m3
        clevelid = fid.createVariable("Contour_levels", "f4", ("contour_levels"))
        clevelid[:] = clevs

        # Dimension with different ensemble members.
        ensembleid = fid.createVariable("ensemble", "str", ("ensid"))
        ensid = fid.createVariable("ensid", "i4", ("ensid"))
        sourceid = fid.createVariable("source", "str", ("ensid"))

        latid = fid.createVariable("latitude", "f4", ("latitude"))
        latid.long_name = "latitude degrees north from the equator"
        latid.units = "degrees_north"
        latid.point_spacing = "even"
        lonid = fid.createVariable("longitude", "f4", ("longitude"))
        lonid.long_name = "longitude degrees east from the greenwhich meridian"
        lonid.units = "degrees_east"
        lonid.point_spacing = "even"

        timeid = fid.createVariable("time", "f4", ("time"))
        # attributes for time grid.
        timeid.units = "days since 1970-01-01 00:00:00"
        timeid.standard_name = "time"
        timeid.bounds = "time_bnds"
        timeid.calendar = "gregorian"

        time_bnds = fid.createVariable("time_bnds", "f4", ("time", "bnds"))

        # Put data into variables
        # only one time per file.

        epoch = np.datetime64("1970-01-01T00:00:00Z")
        # date1 = xrash.time[iii].values
        t1 = (xrash.time[iii].values - epoch) / np.timedelta64(1, "s")
        # change seconds to days
        t1 = t1 / (24.0 * 60 * 60)
        t2 = ((xrash.time[iii].values + sample_time) - epoch) / np.timedelta64(1, "s")

        t2 = t2 / (24.0 * 60 * 60)

        temp = xrash.loc[dict(time=date1)]

        mult = 1
        for jjj, concid in enumerate(concid_list):
            # logger.debug('adding concentration info')
            lev = self.xrash.z.values[jjj]
            concid[:] = makeconc(
                self.xrash.copy(), date1, lev, dotranspose=True, mult=mult
            )
        # logger.debug('adding massloading info')
        massid[:] = makeconc(self.mass, date1, dotranspose=True, level=None)

        latid[:] = latra
        lonid[:] = lonra
        # levelid[:] = levelra
        timeid[:] = t1
        time_bnds[:] = [[t1, t2]]
        # these may be duplicated since ensemble and source
        # dimensions are stacked.
        ensid_list = self.xrash.coords["ens"].values
        if self.add_probs:
           ensid_list.extend(['ensemble mean', 'ensemble standard deviation'])
           ensid_list.extend(['Probability of exceedance 0.02 mg/m3'])
           ensid_list.extend(['Probability of exceedance 0.2 mg/m3'])
           ensid_list.extend(['Probability of exceedance 2 mg/m3'])
           ensid_list.extend(['Probability of exceedance 5 mg/m3'])
           ensid_list.extend(['Probability of exceedance 10 mg/m3'])
        ensembleid[:] =  ensid_list
        sourceid[:] = self.xrash.coords["source"].values
        # this does not need to be changed with add_prob.
        ensid[:] = np.arange(1, ens_shape + 1)
        fid.close()
        logger.debug("made file {}".format(outname))
        return outname


def makeconc(xrash, date1, level, mult=1, dotranspose=False, verbose=False):
    """
    INPUTS
    xrash : xarray data-array
    date1 : datetime.datetime object
    level : list of level names
    RETURNS 
    c1 : data array with concentration from multiple levels combined.
    """
    # this is for mass loading.
    if not level:
        c1 = mult * xrash.sel(time=date1)
    else:
        c1 = mult * xrash.sel(time=date1, z=level)
    # this line is for netcdf awips output
    if dotranspose:
        c1 = c1.transpose("time", "ensemble", "y", "x", transpose_coords=True)
    if verbose:
        print("C1", c1)
    if verbose:
        print(c1.shape)
    return c1

def maketestblist(dname='./'):
    # Need list of tuples. (filename, sourcetag, mettag)
    blist = []
    dname = dname
    fname = "cdump.Aegec00"
    blist.append((os.path.join(dname, fname), "S1", "gec00"))
    fname = "cdump.Aegep01"
    blist.append((os.path.join(dname, fname), "S1", "gep01"))
    return blist

def maketestncfile():
    blist = maketestblist()
    # base name of the netcdf file.
    oname = "out.nc"
    # xarray dataset produced by hysplit.combine_dataset.
    xrash = maketestra()
    # 
    Cdump2Awips(xrash, oname)
    fnames = c2n.create_all_files()

def maketestra():
    blist = maketestblist()
    # xrash is an xarray dataset which can be input into
    # Cdump2Awips class initialization.
    xrash = hysplit.combine_dataset(blist)
    return xrash

# volcat.py
# A reader for VOLCAT ash data using xarray
import sys
import os
from os import walk
import datetime
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import logging
import numpy as np
import numpy.ma as ma
import pandas as pd
import monet
from  utilvolc.get_area import get_area
from utilhysplit import hysplit_gridutil
from monetio.models import hysplit

logger = logging.getLogger(__name__)

# from pyresample.bucket import BucketResampler

# change log
# 2022 Nov 17 AMC updated correct_pc with better regrid support.
# 2022 Nov 17 AMC updated open_dataset
# 2022 Nov 22 AMC correct_pc need to use np.nanmin to get min values

"""
This script contains routines that open/read VOLCAT data in xarray format,
manipulate arrays as necessary, and plots desirable variables.
-------------
Functions:
-------------

For opening files
    open_dataset: opens single VOLCAT file
    get_volcat_list: returns list of data-arrays with volcat data
    find_volcat: finds volcat files in designated directory
    get_volcat_name_df: puts parts of volcat file name in pandas dataframe


Regridding and parallax correction
    correct_pc: corrects parallax
    write_parallax_corrected_files: writes new files with parallax corrected lat/lon

Helper functions
    __pc_loop
    __matchvals:
    _get_time: set time dimension for VOLCAT data
    _get_time2: set time dimension for VOLCAT data
    _get_latlon: rename lat/lon, set coordinates of VOLCAT data

Utilities
    combine_regridded : utilized in ash_inverse.py
    test_volcat: tests to see if parallax corrected lat/lon values exist
    bbox: finds bounding box around data - used for trimming arrays


Extracting variables

    get_data: extracts desired data from large volcat file
    check_names:
    get_pc_latitude: uses parallax corrected latitude
    get_pc_longitude: uses parallax corrected longitude
    get_height: returns array of ash top height from VOLCAT
    get_radius: returns array of ash effective radius from VOLCAT
    get_total_mass: returns total mass in volcat file
    get_mass:  returns array of ash mass loading from VOLCAT
    get_ashdet:
    mass_sum:
    get_time:
    get_atherr: returns array of ash top height error from VOLCAT

Classes: 
VolcatName
     compare: compares volcat file names - shows difference
     parse: parses name into components
     create_name: in progress
------------

Functions moved to volcat_legacy.py

open_dataset2: opens single volcat file in Reventador format
open_mfdataset: opens multiple VOLCAT files
open_hdf: opens single NRT HDF VOLCAT file
create_netcdf: creates netcdf of import variables from NRT HDF
regrid_volcat
regrid_volcat_xesmf
regrid_volcat2
write_regridded_files: writes regridded files
_make2d_latlon
average_volcat :
average_volcat_new :

"""


def open_dataset(
    fname,
    gridspace=None,
    correct_parallax=False,
    mask_and_scale=True,
    decode_times=False,
):
    """
    Opens single VOLCAT file
    gridspace: only necessary if doing parallax correction
    mask_and_scale : needs to be set to True for Bezymianny data.
    decode_times   : needs to be True for parallax corrected files and some of hdf data.

    """
    # 03/07/2021 The Bezy height data has a fill value of -1,
    # scale_factor of 0.01 and offset of 0.
    # The scale factor needs to be applied to get output in km.

    # ash_mass_loading has no scale_factor of offset and fill value is -999.
    if "pc.nc" in fname or "rg.nc" in fname:
        dset = xr.open_dataset(fname, mask_and_scale=mask_and_scale, decode_times=True)
        return dset
    else:
        dset = xr.open_dataset(fname, mask_and_scale=mask_and_scale, decode_times=decode_times)
        # not needed for new Bezy data.
    if 'Dim1' in dset.dims.keys() and 'Dim2' in dset.dims.keys(): 
        dset = dset.rename({"Dim1": "y", "Dim0": "x"})
    #if "some_vars.nc" in fname:
    #    pass
    #elif "Karymsky" in fname:
    #    pass
    #else:
        # use parallax corrected if available and flag is set.
    dset = _get_latlon(dset, "latitude", "longitude")
    dset = _get_time(dset)
    if "pc_latitude" in dset.data_vars and correct_parallax:
        if not gridspace:
            dset = correct_pc(dset)
        else:
            dset = correct_pc(dset, gridspace=gridspace)
        dset.attrs.update({"parallax corrected coordinates": "True"})
    elif "pc_latitude" not in dset.data_vars and correct_parallax:
        print("WARNING: cannot correct parallax. Data not found in file")
        dset.attrs.update({"parallax corrected coordinates": "False"})
    else:
        dset.attrs.update({"parallax corrected coordinates": "False"})
    return dset

def get_volcat_name_df(tdir, daterange=None, vid=None, fid=None, include_last=False):
    """
    Returns dataframe with columns being the information in the vhash
    dictionary of the VolcatName class. This is all the information collected from the filename.
    """
    tlist = find_volcat(tdir, vid=None, daterange=None, return_val=2)
    # vlist is a list of dictionaries with information from the name.
    vlist = [x.vhash for x in tlist]
    if not vlist:
        return pd.DataFrame()
    temp = pd.DataFrame(vlist)
    if isinstance(daterange, (list, np.ndarray)):
        temp = temp[temp["idate"] >= daterange[0]]
        if include_last:
            temp = temp[temp["idate"] <= daterange[1]]
        else:
            temp = temp[temp["idate"] < daterange[1]]
    if vid:
        temp = temp[temp["volcano id"] == vid]
    if fid:
        temp = temp[temp["fid"] == fid]
    if 'volcano id' in temp.columns and 'edate' in temp.columns:
        if "fid" in temp.columns:
            temp = temp.sort_values(["volcano id", "fid", "edate"], axis=0)
        else:
            temp = temp.sort_values(["volcano id", "edate"], axis=0)
    return temp


def get_volcat_list(tdir, daterange=None, vid=None, fid=None,
                    fdate=None, flist=None, return_val=2,
                    correct_parallax=True, mask_and_scale=True,
                    decode_times=True, verbose=False,
                    include_last=True,):
    """
    returns list of data-arrays with volcat data.
    Inputs:
    tdir: string - directory of volcat files
    daterange: datetime object -  [datetime0, datetime1] or none
    vid: string - volcano ID
    fid: feature ID
    fdate: Feature datetime - edate in tframe (datetime object or datetime64)
    flist: list of filenames
    return_val: integer (1,2,3) - see find_volcat() for explanation
    correct_parallax: boolean
    mask_and_scale: boolean
    decode_times: boolean
    verbose: boolean
    include_last: boolean
    Outputs:
    das: list of datasets
    """
    if isinstance(flist,(list,np.ndarray)):
        filenames = flist
    else:
        tframe = get_volcat_name_df(
            tdir, vid=vid, fid=fid, daterange=daterange, include_last=include_last
        )
        if fdate:
            eventd = pd.to_datetime(fdate).to_datetime64()
            filenames = tframe.loc[tframe['edate'] == eventd, 'filename'].tolist()
        else:
            filenames = tframe.filename.values
    das = []
    for nnn, iii in enumerate(filenames):
        if not os.path.isfile(os.path.join(tdir,iii)):
           logger.warning('file not found, skipping file {}'.format(iii))
           print('file not found, skipping file {}'.format(iii))
           continue
        if verbose:
           logger.info('working on {} {} out of {}'.format(nnn, iii, len(filenames)))
           print('working on {} {} out of {}'.format(nnn, iii, len(filenames)))
        # opens volcat files using volcat.open_dataset
        if not "_pc" in iii:
            if verbose: print('opening', iii, correct_parallax, mask_and_scale, decode_times)
            das.append(
                open_dataset(
                    os.path.join(tdir, iii),
                    correct_parallax=correct_parallax,
                    mask_and_scale=mask_and_scale,
                    decode_times=decode_times,
                )
            )
        else:
            das.append(xr.open_dataset(os.path.join(tdir, iii)))
    return das


def write_parallax_corrected_files(
    tdir,
    wdir,
    vid=None,
    daterange=None,
    verbose=False,
    flist=None,
    gridspace=None,
    tag="pc",
):
    """
    ***If flist is not specified, this does not work. There are folders in the tdir and they
    cause a problem with the function. Flist must not include directories, just file names***
    tdir : str : location of volcat files.
    wdir : str : location to write new files
    vid : volcano id : if None will find all
    daterange : [datetime, datetime] : if None will find all.
    verbose: boolean
    flist: list of files? ***NEED TO SPECIFY FILE LIST***
    gridspace: float : grid size of pc array
    tag: used to create filename of new file.

    creates netcdf files with parallax corrected values.
    files have same name with _{tag}.nc added to the end.
    Current convention is to use tag=pc.
    These will be needed for input into MET.

    Currently no overwrite option exists in this function. If the file
    already exists, then this function returns a message to that effect and
    does not overwrite the file.
    """
    logger.info('volcat write_parallax_corrected_files function')
    anum = 0
    newnamelist = []
    if isinstance(flist,(list,np.ndarray)):
        if verbose: print('Using filenamse from list')
        vlist = flist
    else:
        if verbose: print('Finding filenames')
        vlist = find_volcat(tdir, vid, daterange, verbose=verbose, return_val=2)

    for iii, val in enumerate(vlist):
        print(type(val))
        if verbose: print('working on {}'.format(val))
        if isinstance(val, str):
            fname = val
        else:
            fname = val.fname
        new_fname = fname.replace(".nc", "_{}.nc".format(tag))
        newname = os.path.join(wdir,new_fname)
        #print('wdir {}'.format(newname))
        if os.path.isfile(os.path.join(wdir, new_fname)):
            if verbose: print(
                "Netcdf file exists {} in directory {} cannot write ".format(
                    new_fname, wdir
                )
            )
            anum+=1
        else:
            if verbose:
                print("writing {} to {}".format(new_fname, wdir))
            dset = open_dataset(
                os.path.join(tdir, fname), gridspace=gridspace, correct_parallax=True
            )
            newname = os.path.join(wdir,new_fname)
            #print('wdir {}'.format(newname))
            dset.to_netcdf(newname)
            if not os.path.isfile(os.path.join(wdir, new_fname)):
                logger.warning('Warning: file did not write {} {}'.format(wdir,new_fname))
            else:  
                logger.info('file did write {} {}'.format(wdir,new_fname))
    #print('Number of files which were already written {}'.format(anum))

def find_volcat(
    tdir, vid=None, daterange=None, return_val=2, verbose=False, include_last=False
):
    ##NOT WORKING FOR NISHINOSHIMA DATA##
    """
    Locates files in tdir which follow the volcat naming
    convention as defined in VolcatName class.
    If a daterange is defined will return only files

    Inputs:
    tdir : string - volcat files directory
    vid: string - volcano id
    daterange : [datetime, datetime] or None
    include_last : boolean
               True - includes volcat data with date = daterange[1]
               False - only include data with date < daterange[1]
    return_val : integer
               1 - returns dictionary
               2-  returns list of VolcatName objects.
               3 - returns list of filenames
    Returns:
               1 - returns dictionary. key is date. values is VolcatName object.
               2 - returns list of VolcatName objects.
               3 - returns list of filenames
    """
    import sys

    vhash = {}  # dictionary
    nflist = []  # list of filenames
    vnlist = []  # list of filenames
    if not os.path.isdir(tdir):
        print("directory not valid {}".format(tdir))
    for fln in os.listdir(tdir):
        try:
            vn = VolcatName(fln)
        except:
            if verbose:
                print("Not VOLCAT filename {}".format(fln))
            continue
        if daterange and include_last:
            if vn.image_date < daterange[0] or vn.image_date > daterange[1]:
                if verbose:
                    print("date not in range", vn.image_date, daterange[0], daterange[1])
                continue
        elif daterange and not include_last:
            if vn.image_date < daterange[0] or vn.image_date >= daterange[1]:
                if verbose:
                    print("date not in range", vn.image_date, daterange[0], daterange[1])
                continue
        if vid and vn.vhash["volcano id"] != vid:
            continue
        if return_val == 1:
            if vn.image_date not in vhash.keys():
                vhash[vn.image_date] = vn
            else:
                print("two files with same date")
                print(vhash[vn.image_date].compare(vn))
        elif return_val == 2:
            vnlist.append(vn)
        elif return_val == 3:
            nflist.append(fln)
    if return_val == 1:
        return vhash
    elif return_val == 2:
        return vnlist
    elif return_val == 3:
        return nflist


def test_volcat(tdir, daterange=None, verbose=True):
    """
    checks the pc_latitude field for values greater than 0.
    """
    vnlist = find_volcat(tdir, daterange, verbose)
    for key in vnlist.keys():
        vname = vnlist[key].fname
        dset = open_dataset(os.path.join(tdir, vname), pc_correct=False)
        if np.max(dset.pc_latitude) > 0:
            print("passed")
        else:
            print("failed")


class VolcatName:
    """
    12/18/2020 works with 'new' data format.
    parse the volcat name to get information.
    attributes:
    self.fname name of file
    self.date date associated with file
    self.vhash is a dictionary which contains info
    gleaned from the naming convention.

    methods:
    compare: returns what is different between two file names.
    """

    def __init__(self, fname):
        # if full directory path is input then just get the filename
        if "/" in fname:
            temp = fname.split("/")
            self.fname = temp[-1]
        else:
            self.fname = fname
        self.vhash = {}
        self.date = None
        self.image_date = None
        self.event_date = None
        self.image_dtfmt = "s%Y%j_%H%M%S"
        self.event_dtfmt = "b%Y%j_%H%M%S"

        self.make_keylist()
        self.make_datekeys()

        self.pc_corrected = False
        # parse only if a string is given.
        if isinstance(fname, str):
            self.parse(self.fname)
        self.vhash["filename"] = fname

    def make_datekeys(self):
        self.datekeys = [3, 4, 10, 11]


    def make_keylist(self):
        self.keylist = ["algorithm name"]
        self.keylist.append("satellite platform")
        self.keylist.append("event scanning strategy")
        self.keylist.append("image date")  # should be image date (check)
        self.keylist.append("image time")
        self.keylist.append("fid")
        self.keylist.append("volcano id")
        self.keylist.append("description")
        self.keylist.append("WMO satellite id")
        self.keylist.append("image scanning strategy")
        self.keylist.append("event date")  # should be event date (check)
        self.keylist.append("event time")
        self.keylist.append("feature id")

    def __lt__(self, other):
        """
        sort by
        volcano id first.
        event date
        image date
        feature id if it exists.
        """
        if self.vhash["volcano id"] < other.vhash["volcano id"]:
            return True
        if "fid" in self.vhash.keys() and "fid" in other.vhash.keys():
            if self.vhash["fid"] < other.vhash["fid"]:
                return True
        if self.event_date < other.event_date:
            return True
        if self.image_date < other.image_date:
            return True
        sortlist = [
            "feature id",
            "image scanning strategy",
            "WMO satellite id",
            "description",
            "event scanning strategy",
            "satellite platform",
            "algorithm name",
        ]
        for key in sortlist:
            if key in other.vhash.keys() and key in self.vhash.keys():
                if self.vhash[key] < other.vhash[key]:
                    return True

    def compare(self, other):
        """
        other is another VolcatName object.
        Returns
        dictionary of information which is different.
        values is a  tuple of (other value, self value).
        """
        diffhash = {}
        for key in self.keylist:
            if key in other.vhash.keys() and key in self.vhash.keys():
                if other.vhash[key] != self.vhash[key]:
                    diffhash[key] = (other.vhash[key], self.vhash[key])
        return diffhash

    def __str__(self):
        # 2023 14 Jan (amc) make sure keys are in the dictionary.
        keys = self.vhash.keys()
        keylist = [x for x in self.keylist if x in keys]
        val = [self.vhash[x] for x in keylist]
        return str.join("_", val)

    def parse(self, fname):
        # if full_disk in filename replace with fulldisk because _ is used as separator
        fname = fname.replace('Full_Disk', 'FullDisk')
        fname = fname.replace('FULL_DISK', 'FullDisk')
        temp = fname.split("_")
        if "pc" in temp[-1]:
            self.pc_corrected = True
        jjj = 0
        for iii, key in enumerate(self.keylist):
            val = temp[jjj]
            # nishinoshima files have a g00? code before the volcano id.
            if key == "fid":
                if val[0] == "g":
                    self.vhash[key] = val
                else:
                    continue
            self.vhash[key] = val
            jjj += 1
        # Image date marks date of the data collection
        dk = self.datekeys
        dstr = "{}_{}".format(self.vhash[self.keylist[dk[0]]], self.vhash[self.keylist[dk[1]]])
        self.image_date = datetime.datetime.strptime(dstr, self.image_dtfmt)
        # Event date is start of event
        dstr = "{}_{}".format(
            self.vhash[self.keylist[dk[2]]], self.vhash[self.keylist[dk[3]]]
        )
        self.event_date = datetime.datetime.strptime(dstr, self.event_dtfmt)
        self.vhash[self.keylist[dk[3]]] = self.vhash[self.keylist[dk[3]]].replace(".nc", "")

        self.vhash["idate"] = self.image_date
        self.vhash["edate"] = self.event_date
        self.date = self.image_date
        return self.vhash

    def create_name(self):
        """
        To do: returns filename given some inputs.
        """
        return -1



def bbox(darray, fillvalue):
    """Returns bounding box around data
    Input: Must be dataarray
    Outupt: Lower left corner, upper right corner of bounding box
    around data.
    if fillvalue is None then assume Nan's.
    """
    try:
        arr = darray[0, :, :].values
    except:
        print('arr failed', darray)
    if fillvalue:
        a = np.where(arr != fillvalue)
    else:
        a = np.where(~np.isnan(arr))
    if np.min(a[0]) != 0.0 and np.min(a[1]) != 0.0:
        bbox = (
            [np.min(a[0] - 3), np.min(a[1]) - 3],
            [np.max(a[0] + 3), np.max(a[1]) + 3],
        )
    else:
        bbox = ([np.min(a[0]), np.min(a[1])], [np.max(a[0]), np.max(a[1])])
    return bbox


def _get_latlon(dset, name1="latitude", name2="longitude"):
    dset = dset.set_coords([name1, name2])
    return dset


def _get_time(dset):
    import pandas as pd

    temp = dset.attrs["time_coverage_start"]
    time = pd.to_datetime(temp)
    dset["time"] = time
    dset = dset.expand_dims(dim="time")
    dset = dset.set_coords(["time"])
    return dset


def _get_time2(dset):
    import pandas as pd

    date = "20" + str(dset.attrs["Image_Date"])[1:]
    time1 = str(dset.attrs["Image_Time"])
    if len(time1) == 5:
        time1 = "0" + str(dset.attrs["Image_Time"])
    time = pd.to_datetime(date + time1, format="%Y%j%H%M%S", errors="ignore")
    dset["time"] = time
    dset = dset.expand_dims(dim="time")
    dset = dset.set_coords(["time"])
    return dset


# Extracting variables


def get_data(dset, vname, clip=True):
    gen = dset.data_vars[vname]
    atvals = gen.attrs
    fillvalue = None
    if "_FillValue" in gen.attrs:
        fillvalue = gen._FillValue
        gen = gen.where(gen != fillvalue)
        fillvalue = None
    if clip:
        status = True
        try:
            box = bbox(gen, fillvalue)
        except:
            print('volcat get_data bbox for clipping failed')
            status = False
        if status:
            gen = gen[:, box[0][0]: box[1][0], box[0][1]: box[1][1]]
            if "_FillValue" in gen.attrs:
                gen = gen.where(gen != fillvalue)
            else:
                gen = gen.where(gen)
    # applies scale_factor and offset if they are in the attributes.
    if "scale_factor" in gen.attrs:
        gen = gen * gen.attrs["scale_factor"]
    if "offset" in gen.attrs:
        gen = gen + gen.attrs["offset"]
    if "add_offset" in gen.attrs:
        gen = gen + gen.attrs["add_offset"]
    # keep relevant attributes.
    new_attr = {}
    for key in atvals.keys():
        if key not in ["_FillValue", "add_offset", "offset", "scale_factor"]:
            new_attr[key] = atvals[key]
    gen.attrs = new_attr
    return gen


def check_names(dset, vname, checklist, clip=True):
    if vname:
        return get_data(dset, vname, clip=clip)
    for val in checklist:
        if val in dset.data_vars:
            return get_data(dset, val, clip=clip)
    return xr.DataArray()



def get_pc_latitude(dset, vname=None, clip=True):
    """Returns array with retrieved height of the highest layer of ash."""
    """Default units are km above sea-level"""
    checklist = ["pc_latitude"]
    return check_names(dset, vname, checklist, clip=clip)


def get_pc_longitude(dset, vname=None, clip=True):
    """Returns array with retrieved height of the highest layer of ash."""
    """Default units are km above sea-level"""
    checklist = ["pc_longitude"]
    return check_names(dset, vname, checklist, clip=clip)


def get_height(dset, vname=None, clip=True):
    """Returns array with retrieved height of the highest layer of ash.
    Default units are km above sea-level"""
    checklist = ["ash_cth", "ash_cloud_height"]
    return check_names(dset, vname, checklist, clip=clip)


def get_radius(dset, vname=None, clip=True):
    """Returns 2d array of ash effective radius
    Default units are micrometer"""
    checklist = ["ash_r_eff", "effective_radius_of_ash"]
    return check_names(dset, vname, checklist, clip=clip)


def check_total_mass(dset):
    mass = get_mass(dset)
    area = get_area(mass)
    # area returned in km2.
    # mass in g/m2
    # 1e6 to convert are to m2
    # 1e-9 to convert g to Tg
    masstot = mass * area * 1e-6
    masstot = masstot.sum().values  
    # return unit is in Tg
    return  masstot
   

def get_total_mass(dset):
    # unit is in Tg.
    """Units are in Tg"""
    return dset.ash_mass_loading_total_mass.values[0]


def get_mass(dset, vname=None, clip=True):
    """Returns 2d array of ash mass loading
    Default units are grams / meter^2"""
    checklist = ["ash_mass", "ash_mass_loading"]
    return check_names(dset, vname, checklist, clip=clip)


def get_ashdet(dset, vname=None, clip=True):
    """Returns 2d array of detected ash
    Values > 0.0 = detected ash
    Values < 0.0 = no detected ash
    Can be used to determine if ash was detected, but ash mass or ash height was not"""
    checklist = ["ash_spectral_signature_strength"]
    return check_names(dset, vname, checklist, clip=clip)


def mass_sum(dset):
    mass = get_mass(dset)
    mass2 = mass.where(mass > 0.0, 0.0).values
    mass_sum = np.sum(mass2)
    return mass_sum


def get_time(dset):
    time = dset.time_coverage_start
    return time


def get_atherr(dset):
    """Returns array with uncertainty in ash top height from VOLCAT."""
    """Default units are km above sea-level"""
    height_err = dset.ash_cth_uncertainty
    height_err = height_err.where(height_err != height_err._FillValue, drop=True)
    return height_err


def __matchvals(pclat, pclon, massra, height,area):
    # pclat : xarray DataArray
    # pclon : xarray DataArray
    # mass  : xarray DataArray
    # height : xarray DataArray
    # used in correct_pc
    # returns 1D list of tuples of values in the 4 DataArrays
    pclon = pclon.values.flatten()
    pclat = pclat.values.flatten()
    mass = massra.values.flatten()
    height = height.values.flatten()
    area = area.values.flatten()
    tlist = list(zip(pclat, pclon, mass, height,area))
    # only return tuples in which mass has a valid value
    if "_FillValue" in massra.attrs:
        fill = massra.attrs["_FillValue"]
        tlist = [x for x in tlist if x[2] != fill]
    else:
        # get rid of Nans.
        tlist = [x for x in tlist if ~np.isnan(x[2])]
    return tlist


def determine_pc_grid_space():
    return 1

# Notes
# to improve speed of parallax correction
# could look at getting the difference between parallax corrected coordinates
# and actual coordinates. Hopefully difference would be the same across all points.
# then just shift the whole grid by that much.
# then additional regridding would have to be performed.


def correct_pc(dset, gridspace=None):
    """
    moves mass and height values into the coordinate values closest
    to the parallax corrected values. Results in dataset with mass and height shifted
    to parallax corrected positions.

    gridspace : float
                if gridspace is not None then also regrids onto grid with this spacing.
    """
    # 06/02/2021 amc commented out use of the ashdet field.
    # AMR: Added ability to grid parallax corrected data to regular grid
    # 2023 Jan 14 (amc) make sure grids stay in range -180,180 and -90,90.
    # 2023 Jan 14 (amc) replace np.arange with np.linspace 

    mass = get_mass(dset, clip=False)
    height = get_height(dset, clip=False)
    effrad = get_radius(dset, clip=False)
    # ashdet = get_ashdet(dset, clip=False)

    if not gridspace:
        newmass = xr.zeros_like(mass.isel(time=0))
        newnum = xr.zeros_like(mass.isel(time=0))
        newhgt = xr.zeros_like(mass.isel(time=0))
        newrad = xr.zeros_like(mass.isel(time=0))
        lon = dset.ash_mass_loading.longitude.values
        lat = dset.ash_mass_loading.latitude.values
    else:
        # AMR: Making regular grid for parallax correction
        # Changes pushed 8/26/2021
        # Mass, height, and effrad should all have the same lat/lon
        # if this is not the case, need to adjust this portion of the function
        # maybe loop through the variables to calculate max,min?
        # NEED to make this more general, for various grid spaces
        #latmin = round(mass.latitude.values.min())
        #latmax = round(mass.latitude.values.max()) + 1.0
        #lonmin = round(mass.longitude.values.min())
        #lonmax = round(mass.longitude.values.max()) + 1.0
        # some files had nan's at the edges for lat/lon values
        # so change to using this.
        latmin = np.nanmin(mass.latitude.values)
        latmax = np.nanmax(mass.latitude.values) + 1.0
        lonmin = np.nanmin(mass.longitude.values)
        lonmax = np.nanmax(mass.longitude.values) + 1.0
  
        latmin = int(np.floor(latmin))
        latmax =  int(np.ceil(latmax))
        lonmin = int(np.floor(lonmin))
        lonmax = int(np.ceil(lonmax))
        if latmax > 90.0: latmax = 90
        if latmin < -90.0: latmin = -90
        if lonmax > 180.0: lonmax = 180
        if lonmin < -180.0: lonmin = -180

        # --------------------------------------------
        def set_array(lmin,lmax,ddd,tp): 
            # would like to decrease res but linspace
            res = 0.001
            if tp=='latitude':
               tmax = 90.0
               tmin = -90
            else:
               tmax =  180.0
               tmin = -180.0
            numd = int(np.ceil((lmax-lmin)/ddd))
            
            # linspace was giving step sizes of 0.101.
            outra,step = np.linspace(lmin,lmax,num=numd,retstep=True)

            if step - ddd > res:
                print('WARNING  spacing not correct')
                print(outra[1]-outra[0], ddd)
                print(lmin,lmax,numd)
                print(outra)
            
            if outra[-1] > tmax+res: 
                print('WARNING, value greater than max')
                print(outra[-10:])
            return outra
        # --------------------------------------------

        lats = set_array(latmin,latmax,gridspace,'latitude')
        lons = set_array(lonmin,lonmax,gridspace,'longitude')

        #lats = np.arange(latmin, latmax, gridspace)
        #lats = np.arange(latmax, latmin, gridspace * -1)
        #lons = np.arange(lonmin, lonmax, gridspace)
        #lons = [x for x in lons if x<=180.001]
        lon, lat = np.meshgrid(lons, lats)
        #lon = longitude
        #lat = latitude
        tmp = np.zeros_like(lat)
        # Making zeros like arrays
        das = xr.DataArray(
            data=tmp,
            dims=["y", "x"],
            coords=dict(
                latitude=(["y", "x"], lat), longitude=(["y", "x"], lon)
            ),
        )
        newmass = das.copy()
        newhgt = das.copy()
        newrad = das.copy()
        # END of Additional code - AMR
    newmass.attrs = mass.attrs
    newhgt.attrs = height.attrs
    newrad.attrs = effrad.attrs
    area = get_area(mass)
    newarea = get_area(newmass)
    time = mass.time
    pclat = get_pc_latitude(dset, clip=False)
    pclon = get_pc_longitude(dset, clip=False)
    tlist = np.array(__matchvals(pclon, pclat, mass, height,area))

    indexlist = []
    prev_point = 0


    def check_compat(t1, t2):
        if not np.array_equal(t1.longitude.values,t2.longitude.values,equal_nan=True):
           return False
        else:
           return True

    newmass,newhgt,newrad = __pc_loop(tlist,newmass,newhgt,newrad)

    # TODO currently the fill value is 0.
    # possibly change to nan or something else?
    newmass = newmass.assign_attrs({"_FillValue": 0})
    newhgt = newhgt.assign_attrs({"_FillValue": 0})
    newrad = newrad.assign_attrs({"_FillValue": 0})

    newmass = newmass.expand_dims("time")
    newhgt = newhgt.expand_dims("time")
    newrad = newrad.expand_dims("time")

    newmass = newmass.transpose("time", "y", "x", transpose_coords=True)
    newhgt = newhgt.transpose("time", "y", "x", transpose_coords=True)
    newrad = newrad.transpose("time", "y", "x", transpose_coords=True)
    # keep original names for mass and height.
 
    # ---------------------------------------
    # this block resets the longitude and latitude values
    # becase they are floats they change in the pc_loop function
    # slightly and the the Dataset creation throws an error. 
    newmass = newmass.drop("longitude")
    newhgt = newhgt.drop("longitude")
    newrad = newrad.drop("longitude") 
    newmass = newmass.drop("latitude")
    newhgt = newhgt.drop("latitude")
    newrad = newrad.drop("latitude") 

    newmass = newmass.assign_coords(longitude=(("y","x"),lon))
    newmass = newmass.assign_coords(latitude=(("y","x"),lat))
    newhgt = newhgt.assign_coords(longitude=(("y","x"),lon))
    newhgt = newhgt.assign_coords(latitude=(("y","x"),lat))
    newrad = newrad.assign_coords(longitude=(("y","x"),lon))
    newrad = newrad.assign_coords(latitude=(("y","x"),lat))

    # adjust the mass to mass loading by dividing by new area.
    newmass = xr.where(newmass>0, newmass/newarea, newmass)
    newmass = xr.where(newmass==0, np.nan, newmass)

    if not check_compat(newmass, newhgt):
       print('warning: coordinate values are not the same')
    # ---------------------------------------

    dnew = xr.Dataset(
        {
            "ash_mass_loading": newmass,
            "ash_cloud_height": newhgt,
            "effective_radius_of_ash": newrad,
            "ash_mass_loading_total_mass": dset.ash_mass_loading_total_mass,
            "feature_area": dset.feature_area,
            "feature_age": dset.feature_age,
            "feature_id": dset.feature_id,
        }
    )

    dnew.ash_mass_loading.attrs.update(dset.ash_mass_loading.attrs)
    dnew.ash_cloud_height.attrs.update(dset.ash_cloud_height.attrs)
    dnew.effective_radius_of_ash.attrs.update(dset.effective_radius_of_ash.attrs)
    dnew.time.attrs.update({"standard_name": "time"})
    #dnew.latitude.attrs.update({"standard_name": "latitude"})
    #dnew.longitude.attrs.update({"standard_name": "longitude"})
    #dnew = dnew.assign_attrs(dset.attrs)

    checklon = dnew.isel(y=0).longitude.values
    if np.any(checklon > 180.001):
       print('warning longitude values above 180 found')
    if np.any(checklon < -180.001):
       print('warning longitude values below -180 found')


    return dnew

def __pc_loop(tlist,newmass,newhgt,newrad):
    testmass = newmass.copy()
    indexhash = {}
    hthash = {}
    #aratio = 1
    for point in tlist:
        lon = point[0]
        if lon>180.0001: 
           print('changing lon', lon, lon-360.0)
           lon=lon-360.0
        iii = newmass.monet.nearest_ij(lat=point[1], lon=point[0])
        area = point[4]
        hgt = point[3]
        totmass = point[2] * area
        if iii in indexhash.keys():
            indexhash[iii] += totmass
            hthash[iii] = np.max((hthash[iii],hgt))
        else:
            indexhash[iii] = totmass
            hthash[iii] = hgt

        newmass = xr.where(
            (newmass.coords["x"] == iii[0]) & (newmass.coords["y"] == iii[1]),
            indexhash[iii],
            newmass,
        )

        newhgt = xr.where(
            (newhgt.coords["x"] == iii[0]) & (newhgt.coords["y"] == iii[1]),
            hthash[iii],
            newhgt,
        )
        # AMR: Need to adjust this for effective radius - not currently in tlist, and therefore not in iii
        newrad = xr.where(
            (newrad.coords["x"] == iii[0]) & (newrad.coords["y"] == iii[1]),
            point[3],
            newrad,
        )
    return newmass,newhgt, newrad


def combine_regridded(vlist):
    """
    vlist : list of xarray DataArray or DataSet. 
            volcat data which has been regridded onto
            regular lat-lon grid.
            Each object needs to have the same grid spacing.
            If objects have a time dimension then length of time dimension should be 1.

    Returns :
            xarray DataArray or Dataset which is combination of
            those in the input list.
            If input objects have time dimension or coordinate then output will
            be concatenated along the time dimension. If they do not then an
            ens dimension will be created.

    Possible Issues :
           grid corners are set at (0,0) or (-180,0). If input arrays are not
           compatible with those grid corners (e.g. spacing of 0.5 and original grid
           corner of 20.25) then this will fail or produce unexpected results.
    """
    emptyvalue = xr.DataArray()

    # check range of longitude values.
    # from 0 to 360 or -180 to 180?
    # also find gid spacing and make sure grid spacing
    # same for all grids.
    minlon = 9000
    maxlon = -9000
    dlat = []
    dlon = []
    for iii, das in enumerate(vlist):
        lonra = das.isel(y=0).longitude.values
        small = np.min(lonra)
        large = np.max(lonra)
        minlon = np.min([minlon,small])
        maxlon = np.max([maxlon,large])

        temp = das.copy()
        xval = temp.x.values+1
        yval = temp.y.values+1
        temp = temp.assign_coords(x=("x",xval))
        temp = temp.assign_coords(y=("y",yval))
        atemp = hysplit_gridutil.find_grid_specs(temp)
        dlat.append(atemp['Latitude Spacing'])
        dlon.append(atemp['Longitude Spacing'])

    # ----------------------------------
    # set corner longitude accordingly
    # should be -180 to -180
    # crnrlon should always be -180.
    # TODO - this block can be streamlined. 
    if maxlon <= 180.0:
       crnrlon = -180
    elif minlon >=0 and maxlon <=360:
       print('warning lon points > 180 detected')
       crnrlon = -180
    else:
       print('grids not using same longitude range')
       print('longitudes from {} to {}'.format(minlon, maxlon))
       #print('warning: cannot combine')
       crnrlon = -180
       #return emptyvalue
    # ----------------------------------

    # check that all grids have same spacing.
    tol = 0.01
    dlat.sort()
    check1 = dlat[-1]-dlat[0]<=tol 
    if not check1:
       print('Grids do not have the same spacing')
       print(str.join(' ',[str(x) for x in dlat]))
       print('warning: cannot combine')
       return emptyvalue
    dlon.sort()
    check2 = dlon[-1]-dlon[0]<=tol 
    if not check2:
       print('Grids do not have the same spacing')
       print(str.join(' ',[str(x) for x in dlon]))
       print('warning: cannot combine')
       return emptyvalue

    dlat = dlat[0]
    dlon = dlon[0] 
    logger.info('grid spacing is {} {}'.format(dlat,dlon))

    # check if time is a dimension

    rlist = []
    dimvalue = []
    for iii, das in enumerate(vlist):

        # check if time is a dimension
        if 'time' in das.dims:
            das = das.isel(time=0)
            dimvalue.append(das.time.values) 
            newdim = 'time' 
        elif 'time' in das.coords:
            dimvalue.append(das.time.values)
            newdim = 'time' 
        else:
            newdim = 'ens'
            dimvalue.append(iii)

        # put in new x and y indices for global grid.
        #latra = das.sel(x=1).latitude.values
        #lonra = das.sel(y=1).longitude.values
        #ynew, xnew =  hysplit_gridutil.get_new_indices(latra,lonra,attrs)
        #das = das.assign_coords(y=ynew)
        #das = das.assign_coords(x=xnew)

        # create global grid and use that. 
        # this changes the x and y values that correspond to the lat lon values.
        # important because aligning with other grids is done with x and y since
        # they are integers. It is difficult to align with floats.
        attrs = hysplit_gridutil.create_global_grid(dlat,dlon,crnrlon)
        try:
            das = hysplit_gridutil.change_grid(das,attrs)
        except:
            print('volcat not compatible with global grid: ', attrs)
        rlist.append(das)

        # xnew will encompass areas of all the data-arrays.
        if iii==0:
           dnew = das.copy()
        else:
           aaa, dnew = xr.align(das,dnew,join="outer")
           dnew = dnew.fillna(0) 
   
    # now expand each array to size of xnew.   
    blist = []
    for iii, das in enumerate(rlist):
        bbb, junk = xr.align(das,dnew,join="outer")
        bbb.expand_dims(newdim)
        bbb[newdim] = dimvalue[iii]
        blist.append(bbb) 
    # concatenate the expanded arrays which all have same size.
    newra = xr.concat(blist,newdim)
    # add grid information to dataset
    newra.attrs.update(attrs)
    # if don't do this then some of the lat lon coords will be nan.
    newra = hysplit.reset_latlon_coords(newra)
    if 'time' in newra.dims:
        newra = newra.sortby('time',ascending=True)
    hysplit.reset_latlon_coords(newra)
    return newra

 


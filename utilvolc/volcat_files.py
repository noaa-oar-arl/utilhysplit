# volcat.py
# A reader for VOLCAT data using xarray
# For use with MONET
import sys
import os
from os import walk
import datetime
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import numpy as np
import numpy.ma as ma
import pandas as pd
# from pyresample.bucket import BucketResampler


"""
This script contains routines that open/read VOLCAT data in xarray format,
manipulate arrays as necessary, and plots desirable variables.
-------------
Functions:
-------------
open_dataset: opens single VOLCAT file
open_hdf: opens single NRT HDF VOLCAT file
create_netcdf: creates netcdf of import variables from NRT HDF
open_mfdataset: opens multiple VOLCAT files
regrid_volcat: regrids with monet.remap_nearest
regrid_volcat2: regrids with monet.remap_xesmf
average_volcat_new: regrids using regrid_volcat() and adds avg mass load, max ash height variables to dataset
average_volcat: averages VOLCAT files over designated time period
get_volcat_name_df: puts parts of volcat file name in pandas dataframe
get_volcat_list: returns list of data-arrays with volcat data
write_regridded_files: writes regridded files
write_parallax_corrected_files: writes new files with parallax corrected lat/lon
find_volcat: finds volcat files in designated directory
test_volcat: tests to see if parallax corrected lat/lon values exist
open_dataset2: opens single volcat file in Reventador format
bbox: finds bounding box around data - used for trimming arrays
_get_time: set time dimension for VOLCAT data
_get_latlon: rename lat/lon, set coordinates of VOLCAT data
get_data: extracts desired data from large volcat file
check_names:
create_pc_plot: plots parallax corrected vs uncorrected values
compare_pc: compares corrected vs uncorrected values
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
plot_height: plots ash top height from VOLCAT
plot_radius: plots ash effective radius from VOLCAT
plot_mass: plots ash mass loading from VOLCAT
plot_gen: generates quick plot, not saved
matchvals:
matchvals2:
find_iii:
correct_pc: corrects parallax
Class: VolcatName
     compare: compares volcat file names - shows difference
     parse: parses name into components
     create_name: in progress
------------
"""


def open_dataset(fname,
                 correct_parallax=False,
                 mask_and_scale=True, decode_times=False):
    """
    Opens single VOLCAT file
    mask_and_scale : needs to be set to True for Bezymianny data.
    decode_times   : needs to be True for some of the hdf data.

    """
    # 03/07/2021 The Bezy height data has a fill value of -1,
    # scale_factor of 0.01 and offset of 0.
    # The scale factor needs to be applied to get output in km.

    # ash_mass_loading has no scale_factor of offset and fill value is -999.
    dset = xr.open_dataset(fname, mask_and_scale=mask_and_scale,
                           decode_times=decode_times)
    # not needed for new Bezy data.
    try:
        dset = dset.rename({"Dim1": 'y', "Dim0": 'x'})
    except:
        pass
    if 'some_vars.nc' in fname:
        pass
    elif 'pc.nc' in fname or 'rg.nc' in fname:
        return dset
    else:
        # use parallax corrected if available and flag is set.
        dset = _get_latlon(dset, 'latitude', 'longitude')
        dset = _get_time(dset)
    if 'pc_latitude' in dset.data_vars and correct_parallax:
        print('correcting pc')
        dset = correct_pc(dset)
        dset.attrs.update({'parallax corrected coordinates': 'True'})
    elif 'pc_latitude' not in dset.data_vars and correct_parallax:
        print('WARNING: cannot correct parallax. Data not found in file')
        dset.attrs.update({'parallax corrected coordinates': 'False'})
    else:
        dset.attrs.update({'parallax corrected coordinates': 'False'})
    return dset


def open_hdf(fname):
    """ Opens single HDF NRT VOLCAT file"""
    # 4/12/2021 Creating reader for HDF files - NRT, not processed
    dset = xr.open_dataset(fname, mask_and_scale=False, decode_times=False)
    dset = dset.rename({"lines": 'y', "elements": 'x'})
    # create 2d lat lon arrays
    lon2d, lat2d = _make2d_latlon(dset)
    # Rename, add variable and assing to coordinates
    lon = xr.DataArray(lon2d, name='longitude', dims=['y', 'x'])
    lat = xr.DataArray(lat2d, name='latitude', dims=['y', 'x'])
    attrs = dset.attrs
    dset = xr.merge([dset, lat, lon])
    dset = dset.set_coords(['latitude', 'longitude'])
    dset.attrs = attrs
    return dset


def create_netcdf(fname1, fname2):
    """ Creates netcdf of important variables from L1 and L2 VOLCAT hdf files
    Writes to same directory as fname2 files"""
    dset1 = xr.open_dataset(fname1, mask_and_scale=False, decode_times=False)
    lat = dset1.pixel_latitude.rename({'lines': 'y', 'elements': 'x'}).rename('latitude')
    lon = dset1.pixel_longitude.rename({'lines': 'y', 'elements': 'x'}).rename('longitude')
    # Ash Top Height, Ash Mass, Ash Effective Radius
    dset2 = xr.open_dataset(fname2, mask_and_scale=False, decode_times=False)
    attrs = dset2.attrs
    namestr = dset2.attrs['Default_Name_ash_ret']
    mass = dset2[namestr +
                 '_ash_mass_loading'].rename({'lines': 'y', 'elements': 'x'}).rename('ash_mass_loading')
    height = dset2[namestr +
                   '_ash_top_height'].rename({'lines': 'y', 'elements': 'x'}).rename('ash_cloud_height')
    radius = dset2[namestr +
                   '_ash_effective_radius'].rename({'lines': 'y', 'elements': 'x'}).rename('ash_effective_radius')

    # Creating netcdf of important variables
    dset = xr.merge([mass, height, radius, lat, lon])
    dset = dset.set_coords(['latitude', 'longitude'])
    dset.attrs = attrs
    dset = _get_time2(dset)
    dset.to_netcdf(fname2[:-11]+'_'+fname2[-10:-3]+'some_vars.nc')
    return(print(fname2[:-11]+'_'+fname2[-10:-3]+'some_vars.nc created!'))


def open_mfdataset(fname):
    # 12/1/2d020 Not modified for new files (Bezy)
    # TO DO - modify for new files.
    """Opens multiple VOLCAT files"""
    print(fname)
    # dset = xr.open_mfdataset(fname, concat_dim='time', decode_times=False, mask_and_scale=False)
    from glob import glob
    from numpy import sort
    files = sort(glob(fname))
    das = []
    for i in files:
        das.append(open_dataset(i))
    dset = xr.concat(das, dim='time')
    dset = _get_latlon(dset)
    dset = dset.rename({"lines": 'y', "elements": 'x'})
    return dset


def regrid_volcat(das, cdump):
    """
    Regridding with monet.remap_nearest
    das : list of volcat xarrays
    cdump : xarray DataArray with latitude/longitude to regrid to.
    returns xarray with volcat data with dimension of time and regridded to match cdump.
    """
    # In progress.
    # das is list of volcat datasets.
    # cdump is dataset with appropriate grid.
    # This function maps to new grid.

    # remap_nearest may not be what we want to use. Seems that some small areas with low
    # mass 'disappear' using this regridding scheme. May want to look into using pyresample.bucket or other.
    rai = 1e5
    mlist = []
    hlist = []
    total_mass = []
    feature_area = []
    feature_id = []

    for iii, dset in enumerate(das):
        print('time in loop', dset.time.values)
        near_mass = cdump.monet.remap_nearest(
            dset.ash_mass_loading.isel(time=0), radius_of_influence=rai)
        near_height = cdump.monet.remap_nearest(
            dset.ash_cloud_height.isel(time=0), radius_of_influence=rai)
        near_mass = near_mass.compute()
        near_height = near_height.compute()
        mlist.append(near_mass)
        hlist.append(near_height)
        total_mass.append(dset.ash_mass_loading_total_mass)
        feature_area.append(dset.feature_area)
    newmass = xr.concat(mlist, dim='time')
    newhgt = xr.concat(hlist, dim='time')
    totmass = xr.concat(total_mass, dim='time')
    farea = xr.concat(feature_area, dim='time')
    dnew = xr.Dataset({'ash_mass_loading': newmass,
                       'ash_cloud_height': newhgt,
                       # 'effective_radius_of_ash': newrad,
                       'ash_mass_loading_total_mass': totmass,
                       'feature_area': farea})
    # 'feature_area': dset.feature_area,
    # 'feature_age': dset.feature_age,
    # 'feature_id': dset.feature_id})
    # add global attributes.
    dnew = dnew.assign_attrs(dset.attrs)
    dnew.ash_mass_loading.attrs.update(dset.ash_mass_loading.attrs)
    dnew.ash_cloud_height.attrs.update(dset.ash_cloud_height.attrs)
    dnew.time.attrs.update({'standard_name': 'time'})
    # propogate attributes on latitude and longitude
    dnew.latitude.attrs.update(dset.latitude.attrs)
    dnew.longitude.attrs.update(dset.longitude.attrs)
    dnew.attrs.update({'Regrid Method': 'remap_nearest'})
    return dnew

def regrid_volcat2(das, cdump):
    """
    Regridding with monet.remap_xesmf
    das : list of volcat xarrays
    returns xarray with volcat data with dimension of time and regridded to match cdump.
    """
    # In progress.
    # das is list of volcat datasets.
    # cdump is dataset with appropriate grid.
    # This function maps to new grid.
    # mass 'disappear' using this regridding scheme. May want to look into using pyresample.bucket or other
    mlist = []
    hlist = []
    total_mass = []
    for iii, dset in enumerate(das):
        near_mass = cdump.p006.squeeze().monet.remap_xesmf(
            dset.ash_mass_loading.squeeze(), method='bilinear')
        near_height = cdump.p006.monet.squeeze().remap_xesmf(
            dset.ash_cloud_height.squeeze(), method='bilinear')
        near_mass = near_mass.compute()
        near_height = near_height.compute()
        mlist.append(near_mass)
        hlist.append(near_height)
        total_mass.append(dset.ash_mass_loading_total_mass)
    newmass = xr.concat(mlist, dim='time')
    newhgt = xr.concat(hlist, dim='time')
    totmass = xr.concat(total_mass, dim='time')
    dnew = xr.Dataset({'ash_mass_loading': newmass,
                       'ash_cloud_height': newhgt,
                       # 'effective_radius_of_ash': newrad,
                       'ash_mass_loading_total_mass': totmass})
    # 'feature_area': dset.feature_area,
    # 'feature_age': dset.feature_age,
    # 'feature_id': dset.feature_id})
    dnew.time.attrs.update({'standard_name': 'time'})
    dnew.latitude.attrs.update({'standard_name': 'latitude'})
    dnew.longitude.attrs.update({'standard_name': 'longitude'})
    dnew.attrs.update({'Regrid Method': 'remap_xesmf bilinear'})
    return dnew


def regrid_volcat_xesmf(das, cdump, method):
    """
    Regridding with xesmf - first gridding irregular volcat grid to
    regular grid (same resolution) using nearest neighbor
    Then gridding to regular cdump grid using conservative method
    das: list of volcat xarrays
    cdump: hysplit cdump xarray
    method: string of possible xesmf regridding techniques
    returns xarray with volcat data with dimension of time and regridded to match cdump.
    """
    import xesmf as xe
    import numpy as np
    # In progress.

    def regrid(ds_source, ds_target, da_source, method):
        # ds_source: dataset of data to be regridded
        # ds_target: dataset of target array
        # da_source: dataarray of data to be regridded
        # NOTE: latitude and longitude must be named lat lon for this to work - use rename
        regridder = xe.Regridder(ds_source, ds_target, method, periodic=True)
        da_target = regridder(da_source)
        regridder.clean_weight_file()
        return da_target

    def rename(xra):
        # Xarray to rename latitude and longitude
        newx = xra.rename({'latitude': 'lat', 'longitude': 'lon'})
        return newx

    def make_grid(xra, d_lon, d_lat):
        # xra: xarray
        # d_lon: delta lon
        # d_lat: delta lat
        xra = rename(xra)
        grid = xe.util.grid_2d(np.min(xra.lon)-1., np.max(xra.lon)+1., d_lon,
                               np.min(xra.lat)-1., np.max(xra.lat)+1., d_lat)
        return grid

    mlist = []
    hlist = []
    total_mass = []
    vgrid = make_grid(das[-1], 0.05, 0.05)
    # vgrid = make_grid(cdump, 0.1, 0.1)
    for iii, dset in enumerate(das):
        dset2 = rename(dset)
        ashmass = regrid(dset2, vgrid, dset.ash_mass_loading, 'nearest_s2d')
        height = regrid(dset2, vgrid, dset.ash_cloud_height, 'nearest_s2d')
        mlist.append(ashmass)
        hlist.append(height)
        total_mass.append(dset.ash_mass_loading_total_mass)

    mass = xr.concat(mlist, dim='time')
    hgt = xr.concat(hlist, dim='time')
    totmass = xr.concat(total_mass, dim='time')

    # Regridding to cdump array - conservative method needs box bounds
    cgrid = make_grid(cdump, 0.1, 0.1)
    newmass = regrid(vgrid, cgrid, mass, method)
    newhgt = regrid(vgrid, cgrid, hgt, method)
    newmass = mass
    newhgt = hgt

    dnew = xr.Dataset({'ash_mass_loading': newmass,
                       'ash_cloud_height': newhgt,
                       # 'effective_radius_of_ash': newrad,
                       'ash_mass_loading_total_mass': totmass})
    # 'feature_area': dset.feature_area,
    # 'feature_age': dset.feature_age,
    # 'feature_id': dset.feature_id})
    dnew.time.attrs.update({'standard_name': 'time'})
    return dnew


def average_volcat_new(das, cdump, skipna=False, convert_nans=False):
    # STILL IN PROGRESS
    """
    Very similar to average_volcat() except it regrids using regrid_volcat()
    Output contains full array from regrid_volcat() with average/maximum
    arrays added

    Inputs:
    das: list of volcat datasets
    cdump: xarray of target grid (hysplit cdump usually)
    skipna: boolean - skip nans when taking mean/max
    convert_nans: boolean - convert nans to 0. in all xarrays BEFORE regridding
    outputs:
    dsetnew: dataset with added ash mass mean/ash height max
    """
    fill = 'nan'
    if convert_nans:
        fill = 'No fill value'
        dset = []
        i = 0
        while i < len(das):
            dset.append(das[i].fillna(0.))
            i += 1
        hxr = cdump.fillna(0)
    else:
        dset = das
        hxr = cdump
    dnew = regrid_volcat(dset, hxr)

    avgmass = dnew.ash_mass_loading.mean(dim='time', skipna=skipna, keep_attrs=True)
    maxhgt = dnew.ash_cloud_height.max(dim='time', skipna=skipna, keep_attrs=True)
    # renaming variable
    avgmass = avgmass.load().rename('ash_mass_avg')
    maxhgt = maxhgt.load().rename('ash_height_max')
    # Adding time dimension, changing long_name
    avgmass = avgmass.assign_coords(time=dnew.time[-1]).expand_dims('time')
    maxhgt = maxhgt.assign_coords(time=dnew.time[-1]).expand_dims('time')
    avgmass.attrs['long_name'] = 'Average total column loading of ash in the highest continuous ash layer for the previous hour'
    avgmass.attrs['fill_value'] = fill
    maxhgt.attrs['long_name'] = 'Maximum cloud top height of the highest continuous ash layer for the previous hour'
    maxhgt.attrs['fill_value'] = fill

    # Merging datasets
    dsetnew = xr.merge([dnew, avgmass, maxhgt], combine_attrs='drop_conflicts')

    return dsetnew


def average_volcat(das, cdump, skipna=False, convert_nans=False):
    # In progress.
    """
    Function first regrids to new grid, then calculates then average
    mass loading and maximum ash height.

    Inputs:
    das: list of volcat datasets
    cdump: is dataset with appropriate grid
    skipna: boolean - flag to skip nan values in grid when calculating mean or max
    convert_nans: noolean - flag to convert nans in datasets to 0
    Output:
    avgmass: mean of volcat ash mass loading, from das
    maxhgt: maximum of volcat ash cloud height

    Notes:
    remap_nearest may not be what we want to use.
    Seems that some small areas with low
    mass 'disappear' using this regridding scheme.
    May want to look into using pyresample.bucket or other.
    """
    rai = 1e5
    mlist = []
    hlist = []
    for iii, dset in enumerate(das):
        near_mass = cdump.monet.remap_nearest(
            dset.ash_mass_loading.isel(time=0), radius_of_influence=rai).load()
        near_height = cdump.monet.remap_nearest(
            dset.ash_cloud_height.isel(time=0), radius_of_influence=rai).load()
        mlist.append(near_mass)
        hlist.append(near_height)
    newmass = xr.concat(mlist, dim='time')
    newhgt = xr.concat(hlist, dim='time')
    # when averaging the mass need to convert nan's to zero?
    if convert_nan:
        newmass = newmass.fillna(0.)
        newhgt = newhgt.fillna(0.)
    # option to skip nans
    avgmass = newmass.mean(dim='time', skipna=skipna)
    # note that averaging the height is not correct, better to take maximum along time
    maxhgt = newhgt.max(dim='time', skipna=skipna)
    return avgmass, maxhgt



def get_volcat_name_df(tdir, 
                       daterange=None, 
                       vid=None,
                       fid=None,
                       include_last=False):
    """
    Returns dataframe with columns being the information in the vhash
    dictionary of the VolcatName class. This is all the information collected from the filename.
    """
    tlist = find_volcat(tdir, vid=None, daterange=None, return_val=2)
    vlist = [x.vhash for x in tlist]
    temp = pd.DataFrame(vlist)
    if isinstance(daterange, (list,np.ndarray)):
        temp = temp[temp['edate'] >= daterange[0]]  
        if include_last:
            temp = temp[temp['edate'] <= daterange[1]]  
        else:
            temp = temp[temp['edate'] < daterange[1]]  
    if vid:
        temp = temp[temp['volcano id']==vid]
    if fid:
        temp = temp[temp['fid']==fid]

    if 'fid' in temp.columns:
        temp = temp.sort_values(['volcano id','fid','edate'],axis=0)
    else:
        temp = temp.sort_values(['volcano id','edate'],axis=0)
    return temp

# two json files.
# the first one 
def read_event_summary():
    # list of all active eruptions in the satellite image.
    # Return dataframe made from dictionary in the event summary.
    return df

def summarize_files(volcat_event_df):
    """
    what volcano id's are available.
    what time periods.
    """
    # return volcano id's.
    return -1

def choose_files(volcat_event_df, vid, frequency=10):
    """
    volcat_event_df with columns specified by VolcatName dictionary.
    frequency: how far apart should files to be spaced (minutes)
    """
    return -1




def get_volcat_list(tdir, 
                    daterange=None, 
                    vid=None, 
                    fid=None,
                    flist=None, 
                    return_val=2, 
                    correct_parallax=True, 
                    mask_and_scale=True, 
                    decode_times=True, 
                    verbose=False, 
                    include_last=True):
    """
    returns list of data-arrays with volcat data.
    Inputs:
    tdir: string - directory of volcat files
    daterange: datetime object -  [datetime0, datetime1] or none
    vid: string - volcano ID
    return_val: integer (1,2,3) - see find_volcat() for explanation
    correct_parallax: boolean
    mask_and_scale: boolean
    decode_times: boolean
    verbose: boolean
    include_last: boolean
    Outputs:
    das: list of datasets
    """
    if flist:
       filnames = flist
    else:
       tframe = get_volcat_name_df(tdir,vid=vid,fid=fid,daterange=daterange,include_last=include_last)
       filenames = tframe.filename.values 
    das = []
    for iii in filenames:
        # opens volcat files using volcat.open_dataset
        if not '_pc' in iii:
            das.append(open_dataset(os.path.join(tdir, iii),
                                    correct_parallax=correct_parallax,
                                    mask_and_scale=mask_and_scale,
                                    decode_times=decode_times))
        else:
            das.append(xr.open_dataset(os.path.join(tdir, iii)))
    return das

def find_volcat(tdir, vid=None, daterange=None,
                return_val=2, verbose=False, include_last=False):
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
        print('directory not valid {}'.format(tdir))
    for fln in os.listdir(tdir):
        try:
            vn = VolcatName(fln)
        except:
            if verbose:
                print('Not VOLCAT filename {}'.format(fln))
            continue
        if daterange and include_last:
            if vn.date < daterange[0] or vn.date > daterange[1]:
                if verbose:
                    print('date not in range', vn.date, daterange[0], daterange[1])
                continue
        elif daterange and not include_last:
            if vn.date < daterange[0] or vn.date >= daterange[1]:
                if verbose:
                    print('date not in range', vn.date, daterange[0], daterange[1])
                continue
        if vid and vn.vhash['volcano id'] != vid:
            continue
        if return_val == 1:
            if vn.date not in vhash.keys():
                vhash[vn.date] = vn
            else:
                print('two files with same date')
                print(vhash[vn.date].compare(vn))
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
        if '/' in fname:
            temp = fname.split('/')
            self.fname = temp[-1]
        else:
            self.fname = fname
        self.vhash = {}
        self.date = None
        self.dtfmt = "s%Y%j_%H%M%S"
        self.image_dtfmt = "b%Y%j_%H%M%S"

        self.keylist = ['algorithm name']
        self.keylist.append('satellite platform')
        self.keylist.append('event scanning strategy')
        self.keylist.append('event date')
        self.keylist.append('event time')
        self.keylist.append('fid')
        self.keylist.append('volcano id')
        self.keylist.append('description')
        self.keylist.append('WMO satellite id')
        self.keylist.append('image scanning strategy')
        self.keylist.append('image date')
        self.keylist.append('image time')
        self.keylist.append('feature id')

        self.pc_corrected = False
        self.parse(self.fname)
        self.vhash['filename'] = fname

    def __lt__(self, other):
        """
        sort by 
        volcano id first.
        date
        feature id if it exists. 
        """
        if self.vhash['volcano id'] < other.vhash['volcano id']:
            return True
        if 'fid' in self.vhash.keys() and 'fid' in other.vhash.keys():
            if self.vhash['fid'] < other.vhash['fid']:
                return True
        if self.date < other.date:
            return True
        if self.image_date < other.image_date:
            return True
        sortlist = ['feature id', 'image scanning strategy',
                    'WMO satellite id', 'description', 'event scanning strategy',
                    'satellite platform', 'algorithm name']
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
        val = [self.vhash[x] for x in self.keylist]
        return str.join('_', val)

    def parse(self, fname):
        temp = fname.split('_')
        if 'pc' in temp[-1]:
            self.pc_corrected = True
        jjj = 0
        for iii, key in enumerate(self.keylist):
            val = temp[jjj]
            # nishinoshima files have a g00? code before the volcano id.
            if key == 'fid':
                if val[0] == 'g':
                    self.vhash[key] = val
                else:
                    continue
            self.vhash[key] = val
            jjj += 1
        # Event date marks date of the data collection
        dstr = '{}_{}'.format(self.vhash[self.keylist[3]],
                              self.vhash[self.keylist[4]])
        self.date = datetime.datetime.strptime(dstr, self.dtfmt)
        # Image date is ?
        dstr = '{}_{}'.format(self.vhash[self.keylist[10]],
                              self.vhash[self.keylist[11]])
        self.image_date = datetime.datetime.strptime(dstr, self.image_dtfmt)
        self.vhash[self.keylist[11]] = self.vhash[self.keylist[11]].replace('.nc', '')

        self.vhash['idate'] = self.image_date
        self.vhash['edate'] = self.date

        return self.vhash

    def create_name(self):
        """
        To do: returns filename given some inputs.
        """
        return -1



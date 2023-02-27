# volcat_legacy.py
# Some functions from volcat.py that are no longer maintained.
import sys
import os
#from os import walk
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
from utilvolc import find_volcat, bbox, _get_time, _get_latlon

logger = logging.getLogger(__name__)

# from pyresample.bucket import BucketResampler

# change log
# 2022 Nov 17 AMC updated correct_pc with better regrid support.
# 2022 Nov 17 AMC updated open_dataset
# 2022 Nov 22 AMC correct_pc need to use np.nanmin to get min values

"""
This script contains legacy routines that open/read VOLCAT data in xarray format,
manipulate arrays as necessary, and plots desirable variables.

Functions

open_dataset2: opens single volcat file in Reventador format
open_hdf: opens single NRT HDF VOLCAT file
create_netcdf: creates netcdf of import variables from NRT HDF
open_mfdataset: opens multiple VOLCAT files
regrid_volcat :
regrid_volcat2 :
regrid_volcat_xesmf :
average_volcat_new :
average_volcat :
write_regridded_files :
_make2d_latlon :
"""


#def open_hdf(fname):
#    """Opens single HDF NRT VOLCAT file"""
#    # 4/12/2021 Creating reader for HDF files - NRT, not processed
#    dset = xr.open_dataset(fname, mask_and_scale=False, decode_times=False)
#    dset = dset.rename({"lines": "y", "elements": "x"})
#    # create 2d lat lon arrays
#    lon2d, lat2d = _make2d_latlon(dset)
#    # Rename, add variable and assing to coordinates
#    lon = xr.DataArray(lon2d, name="longitude", dims=["y", "x"])
#    lat = xr.DataArray(lat2d, name="latitude", dims=["y", "x"])
#    attrs = dset.attrs
#    dset = xr.merge([dset, lat, lon])
#    dset = dset.set_coords(["latitude", "longitude"])
#    dset.attrs = attrs
#    return dset


#def create_netcdf(fname1, fname2):
#    """Creates netcdf of important variables from L1 and L2 VOLCAT hdf files
#    Writes to same directory as fname2 files"""
#    dset1 = xr.open_dataset(fname1, mask_and_scale=False, decode_times=False)
#    lat = dset1.pixel_latitude.rename({"lines": "y", "elements": "x"}).rename(
#        "latitude"
#    )
#    lon = dset1.pixel_longitude.rename({"lines": "y", "elements": "x"}).rename(
#        "longitude"
#    )
    # Ash Top Height, Ash Mass, Ash Effective Radius
#    dset2 = xr.open_dataset(fname2, mask_and_scale=False, decode_times=False)
#    attrs = dset2.attrs
#    namestr = dset2.attrs["Default_Name_ash_ret"]
#    mass = (
#        dset2[namestr + "_ash_mass_loading"]
#        .rename({"lines": "y", "elements": "x"})
#        .rename("ash_mass_loading")
#    )
#    height = (
#        dset2[namestr + "_ash_top_height"]
#        .rename({"lines": "y", "elements": "x"})
#        .rename("ash_cloud_height")
#    )
#    radius = (
#        dset2[namestr + "_ash_effective_radius"]
#        .rename({"lines": "y", "elements": "x"})
#        .rename("ash_effective_radius")
#    )
#
    # Creating netcdf of important variables
#    dset = xr.merge([mass, height, radius, lat, lon])
#    dset = dset.set_coords(["latitude", "longitude"])
#    dset.attrs = attrs
#    dset = _get_time2(dset)
#    dset.to_netcdf(fname2[:-11] + "_" + fname2[-10:-3] + "some_vars.nc")
#    return print(fname2[:-11] + "_" + fname2[-10:-3] + "some_vars.nc created!")


#def open_mfdataset(fname):
    # 12/1/2d020 Not modified for new files (Bezy)
    # TO DO - modify for new files.
#    """Opens multiple VOLCAT files"""
#    print(fname)
    # dset = xr.open_mfdataset(fname, concat_dim='time', decode_times=False, mask_and_scale=False)
#    from glob import glob
#    from numpy import sort

#    files = sort(glob(fname))
#    das = []
#    for i in files:
#        das.append(open_dataset(i))
#    dset = xr.concat(das, dim="time")
#    dset = _get_latlon(dset)
#    dset = dset.rename({"lines": "y", "elements": "x"})
#    return dset


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
        print("time in loop", dset.time.values)
        near_mass = cdump.monet.remap_nearest(
            dset.ash_mass_loading.isel(time=0), radius_of_influence=rai
        )
        near_height = cdump.monet.remap_nearest(
            dset.ash_cloud_height.isel(time=0), radius_of_influence=rai
        )
        near_mass = near_mass.compute()
        near_height = near_height.compute()
        mlist.append(near_mass)
        hlist.append(near_height)
        total_mass.append(dset.ash_mass_loading_total_mass)
        feature_area.append(dset.feature_area)
    newmass = xr.concat(mlist, dim="time")
    newhgt = xr.concat(hlist, dim="time")
    totmass = xr.concat(total_mass, dim="time")
    farea = xr.concat(feature_area, dim="time")
    dnew = xr.Dataset(
        {
            "ash_mass_loading": newmass,
            "ash_cloud_height": newhgt,
            'effective_radius_of_ash': newrad,
            "ash_mass_loading_total_mass": totmass,
            "feature_area": farea,
        }
    )
    # 'feature_area': dset.feature_area,
    # 'feature_age': dset.feature_age,
    # 'feature_id': dset.feature_id})
    # add global attributes.
    dnew = dnew.assign_attrs(dset.attrs)
    dnew.ash_mass_loading.attrs.update(dset.ash_mass_loading.attrs)
    dnew.ash_cloud_height.attrs.update(dset.ash_cloud_height.attrs)
    dnew.time.attrs.update({"standard_name": "time"})
    # propogate attributes on latitude and longitude
    dnew.latitude.attrs.update(dset.latitude.attrs)
    dnew.longitude.attrs.update(dset.longitude.attrs)
    dnew.attrs.update({"Regrid Method": "remap_nearest"})
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
            dset.ash_mass_loading.squeeze(), method="bilinear"
        )
        near_height = cdump.p006.monet.squeeze().remap_xesmf(
            dset.ash_cloud_height.squeeze(), method="bilinear"
        )
        near_mass = near_mass.compute()
        near_height = near_height.compute()
        mlist.append(near_mass)
        hlist.append(near_height)
        total_mass.append(dset.ash_mass_loading_total_mass)
    newmass = xr.concat(mlist, dim="time")
    newhgt = xr.concat(hlist, dim="time")
    totmass = xr.concat(total_mass, dim="time")
    dnew = xr.Dataset(
        {
            "ash_mass_loading": newmass,
            "ash_cloud_height": newhgt,
            # 'effective_radius_of_ash': newrad,
            "ash_mass_loading_total_mass": totmass,
        }
    )
    # 'feature_area': dset.feature_area,
    # 'feature_age': dset.feature_age,
    # 'feature_id': dset.feature_id})
    dnew.time.attrs.update({"standard_name": "time"})
    dnew.latitude.attrs.update({"standard_name": "latitude"})
    dnew.longitude.attrs.update({"standard_name": "longitude"})
    dnew.attrs.update({"Regrid Method": "remap_xesmf bilinear"})
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
        newx = xra.rename({"latitude": "lat", "longitude": "lon"})
        return newx

    def make_grid(xra, d_lon, d_lat):
        # xra: xarray
        # d_lon: delta lon
        # d_lat: delta lat
        xra = rename(xra)
        grid = xe.util.grid_2d(
            np.min(xra.lon) - 1.0,
            np.max(xra.lon) + 1.0,
            d_lon,
            np.min(xra.lat) - 1.0,
            np.max(xra.lat) + 1.0,
            d_lat,
        )
        return grid

    mlist = []
    hlist = []
    total_mass = []
    vgrid = make_grid(das[-1], 0.05, 0.05)
    # vgrid = make_grid(cdump, 0.1, 0.1)
    for iii, dset in enumerate(das):
        dset2 = rename(dset)
        ashmass = regrid(dset2, vgrid, dset.ash_mass_loading, "nearest_s2d")
        height = regrid(dset2, vgrid, dset.ash_cloud_height, "nearest_s2d")
        mlist.append(ashmass)
        hlist.append(height)
        total_mass.append(dset.ash_mass_loading_total_mass)

    mass = xr.concat(mlist, dim="time")
    hgt = xr.concat(hlist, dim="time")
    totmass = xr.concat(total_mass, dim="time")

    # Regridding to cdump array - conservative method needs box bounds
    cgrid = make_grid(cdump, 0.1, 0.1)
    newmass = regrid(vgrid, cgrid, mass, method)
    newhgt = regrid(vgrid, cgrid, hgt, method)
    newmass = mass
    newhgt = hgt

    dnew = xr.Dataset(
        {
            "ash_mass_loading": newmass,
            "ash_cloud_height": newhgt,
            # 'effective_radius_of_ash': newrad,
            "ash_mass_loading_total_mass": totmass,
        }
    )
    # 'feature_area': dset.feature_area,
    # 'feature_age': dset.feature_age,
    # 'feature_id': dset.feature_id})
    dnew.time.attrs.update({"standard_name": "time"})
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
    fill = "nan"
    if convert_nans:
        fill = "No fill value"
        dset = []
        i = 0
        while i < len(das):
            dset.append(das[i].fillna(0.0))
            i += 1
        hxr = cdump.fillna(0)
    else:
        dset = das
        hxr = cdump
    dnew = regrid_volcat(dset, hxr)

    avgmass = dnew.ash_mass_loading.mean(dim="time", skipna=skipna, keep_attrs=True)
    maxhgt = dnew.ash_cloud_height.max(dim="time", skipna=skipna, keep_attrs=True)
    # renaming variable
    avgmass = avgmass.load().rename("ash_mass_avg")
    maxhgt = maxhgt.load().rename("ash_height_max")
    # Adding time dimension, changing long_name
    avgmass = avgmass.assign_coords(time=dnew.time[-1]).expand_dims("time")
    maxhgt = maxhgt.assign_coords(time=dnew.time[-1]).expand_dims("time")
    avgmass.attrs[
        "long_name"
    ] = "Average total column loading of ash in the highest continuous ash layer for the previous hour"
    avgmass.attrs["fill_value"] = fill
    maxhgt.attrs[
        "long_name"
    ] = "Maximum cloud top height of the highest continuous ash layer for the previous hour"
    maxhgt.attrs["fill_value"] = fill

    # Merging datasets
    dsetnew = xr.merge([dnew, avgmass, maxhgt], combine_attrs="drop_conflicts")

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
            dset.ash_mass_loading.isel(time=0), radius_of_influence=rai
        ).load()
        near_height = cdump.monet.remap_nearest(
            dset.ash_cloud_height.isel(time=0), radius_of_influence=rai
        ).load()
        mlist.append(near_mass)
        hlist.append(near_height)
    newmass = xr.concat(mlist, dim="time")
    newhgt = xr.concat(hlist, dim="time")
    # when averaging the mass need to convert nan's to zero?
    if convert_nan:
        newmass = newmass.fillna(0.0)
        newhgt = newhgt.fillna(0.0)
    # option to skip nans
    avgmass = newmass.mean(dim="time", skipna=skipna)
    # note that averaging the height is not correct, better to take maximum along time
    maxhgt = newhgt.max(dim="time", skipna=skipna)
    return avgmass, maxhgt



def write_regridded_files(
    cdump, tdir, wdir, tag="rg", vid=None, daterange=None, verbose=False
):
    """
    cdump : xarray DataArray with latitude and longitude grid for regridding.
    tdir : str : location of volcat files.
    wdir : str : location to write new files
    vid : volcano id : if None will find all
    daterange : [datetime, datetime] : if None will find all.
    verbose: boolean
    tag: used to create filename of new file.

    creates netcdf files regridded values.
    files have same name with _{tag}.nc added to the end.
    Current convention is to use tag=pc.
    These will be needed for input into MET.

    Currently no overwrite option exists in this function. If the file
    already exists, then this function returns a message to that effect and
    does not overwrite the file.
    """
    vlist = find_volcat(tdir, vid, daterange, verbose=verbose, return_val=2)
    for iii, val in enumerate(vlist):
        fname = val.fname
        new_fname = fname.replace(".nc", "_{}.nc".format(tag))
        if os.path.isfile(os.path.join(wdir, new_fname)):
            print(
                "Netcdf file exists {} in directory {} cannot write ".format(
                    new_fname, wdir
                )
            )
        else:
            if verbose:
                print("writing {} to {}".format(new_fname, wdir))
            dset = open_dataset(
                os.path.join(tdir, fname), correct_parallax=False, decode_times=True
            )
            dnew = regrid_volcat([dset], cdump)
            dnew.to_netcdf(os.path.join(wdir, new_fname))




def open_dataset2(fname):
    """Opens single VOLCAT file in reventador format"""
    print(fname)
    dset = xr.open_dataset(fname, mask_and_scale=False, decode_times=False)
    # dset = dset.rename({"Dim1":'y',"Dim0":'x'})
    # dset = _get_latlon(dset)
    # dset = _get_time(dset)
    return dset




def _make2d_latlon(dset):
    lon = np.linspace(
        dset.attrs["Longitude_Range"][0],
        dset.attrs["Longitude_Range"][1],
        dset.attrs["Last_Element_Processed"],
    )
    lat = np.linspace(
        dset.attrs["Latitude_Range"][1],
        dset.attrs["Latitude_Range"][0],
        dset.attrs["Line_Segment_Size"],
    )
    lon2d, lat2d = np.meshgrid(lon, lat)
    return lon2d, lat2d











 


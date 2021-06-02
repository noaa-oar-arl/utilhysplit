import os
import xarray as xr
import numpy as np
import datetime
import monetio.models.hysplit as hysplit
from utilhysplit.evaluation.statmain import cdf
from utilhysplit.evaluation.statmain import pixel_matched_cdf

"""
FUNCTIONS
The functions in this file are for
manipulating concentration xarray objects which have dimensions of lat,lon,z,time,ens,source
manipulating massloading   xarray objects which have dimensions of lat,lon,time,ens,source

ATL       : applied threshold level
APLra     : applied percentile level
APL       : applied percentile level
ens_cdf   : cumulative distribution functions.


# Needs more work:
make_ATL_hysp:
topheight : returns height of level

"""

# _______________________________________________________________________

# 2021 Jun 2 amc some functions moved to ensemble_tools_plotting.py
# 2021 Jun 2 amc replaced choose_and_stack function with preprocess.


def APL(indra, problev=50, enslist=None, sourcelist=None):
    """
    Applied Percentile Level
    """
    newra = APLra(indra, enslist, sourcelist)
    vpi = np.where(newra.percent_level.values <= problev)
    # if pick value lower than one available then just provide lowest level.
    try:
        pindex = vpi[0][-1]
    except:
        pindex=0
    return newra.sel(index=pindex)


def APLra(indra, enslist=None, sourcelist=None):
    """
    Applied Percentile Level

    indra: xr dataarray produced by combine_dataset or by hysp_massload
            it must have 'ens' dimension, 'source' dimension or both.
            time and z dimensions are optional.
 
    enslist : list of values to use for 'ens' coordinate
    sourcelist : list of values to use for 'source' coordinate

    Returns:
           xarray data array sorted along the 'ensemble' dimension.
          "ens" and/or 'source'  dimension is replaced by "percent_level" dimension
           and 'index' coordinate.
    """
    # combine 'source' and 'ens' dimensions if applicable.
    dra2, dim = preprocess(indra,enslist,sourcelist)
    coords = dra2.coords
    coords2 = dict(coords)
    # remove 'ens' from the coordinate dictionary.
    coords2.pop('ens')
    dims = list(dra2.dims)
    # find which dimension is the 'ens' dimension
    dii = dims.index(dim)
    # sort along the ensemble axis.
    # sort is an inplace operation. returns an empty array.
    # this type of sort only available on numpy array.
    # does not work in xarray dataarray.
    dvalues = dra2.values.copy()
    dvalues.sort(axis=dii)
    # instead of an 'ensemble' dimension, now have an 'index' dimension.
    dims[dii] = 'index'
    newra = xr.DataArray(dvalues, dims=dims, coords=coords2)
    percents = 100*(newra.index.values+1) / len(newra.index.values)
    newra = newra.assign_coords(percent_level=('index', percents))
    return newra

def volcATL(indra):
    newra = indra.MER * indra
    return ATL(indra)


def preprocess(indra,enslist=None,sourcelist=None):
    """
    indra : xarray dataArray
    enslist : list or np.ndarray
    sourcelist : list or np.ndarray

    Returns:
    dra : xarray datatArray: same dimensions as indra except will stack 'source' and 'ens'
          dimensions if both present.
    dim : str : indicates whether 'ens' or 'source' dimension is to be used.
    """

    dra = indra.copy()

    # handle whether using 'ens' dimension, 'source' dimension or both.
    dim = 'ens'

    # if lists are an np.ndarray then testing them with not
    # returns an error.
    if isinstance(sourcelist, np.ndarray):
        pass
    elif "source" in dra.dims and not sourcelist:
        sourcelist = dra.source.values 

    if isinstance(enslist, np.ndarray):
        pass
    elif "ens" in dra.dims and not enslist:
        enslist = dra.ens.values 


    if "source" in dra.dims and 'ens' in dra.dims:
        dra = dra.sel(ens=enslist)
        dra = dra.sel(source=sourcelist)
        dra = dra.stack(ens=("ens", "source"))

    elif "source" in dra.dims:
        dra = dra.sel(source=sourcelist)
        dim = 'source'
    elif "ens" in dra.dims:
        dra = dra.sel(ens=enslist)
    else:
        print('Warning: could not find source or ens dimension')
        dim = None
    return dra, dim



def ATL(indra, enslist=None, sourcelist=None,
        thresh=0,  norm=False, weights=None):
    """
     Applied Threshold Level (also ensemble frequency of exceedance).

     indra: xr dataarray produced by combine_dataset or by hysp_massload
            it must have 'ens' dimension, 'source' dimension or both.
            time and z dimensions are optional.
 
     enslist : list of values to use for 'ens' coordinate
     sourcelist : list of values to use for 'source' coordinate

     thresh : int or float. If 0 then use > for test. otherwise use >=.

     weights : numpy array of same length as enslist + sourcelist containing weight for
               each member.

     Returns array with number of ensemble members above
     given threshold at each location.
     dimensions will be same as input except no 'ens' or 'source' dimension.

     norm : boolean
            if True return percentage of members.
            if False return number of members.
            using norm is same as applying uniform weighting.
            should not be used with weights.

     """
    # handle whether using 'ens' dimension, 'source' dimension or both.
    dra, dim = preprocess(indra,enslist,sourcelist)
    
    if thresh == 0: 
        dra2 = xr.where(dra>thresh,1,0)
    else: 
        dra2 = xr.where(dra>=thresh,1,0)

    if isinstance(weights,np.ndarray) or isinstance(weights,list):
       if len(weights) != len(dra2[dim].values):
          print('WARNING : weights do not coorespond to values')
          print('weights :', weights)
          print('values : ', dra2[dim].values)
       else:
          wra = xr.DataArray(weights,dims=dim)
          dra2 = wra * dra2 
    dra2 = dra2.sum(dim=[dim])
    if norm:
        nmembers = len(enslist)
        dra2 = dra2/ nmembers
    return dra2



def make_ATL_hysp(xra, variable='p006', threshold=0., MER=None):
    # amr version. maybe not needed now that ATL has been modified.
    """
    Uses threshold value to make binary field of ensemble members.
    For use in statistics calculations. Based on mass loading, no vertical component
    Inputs:
    xra: xarray dataset - netcdf files from hysplit.combine_dataset
    variable: variable name (string) from xra
    threshold: ash threshold - default=0. (float/integer)
    MER: mass eruption rate - default is Mastin equation MER from xra attributes, units are  (float/integer)
    Output:
    xrabin: binary xarray dataarray (source, x, y)
    """

    xra2 = xra[variable] * 1000.  # Each vertical layer is 1000m - MAY WANT TO ALLOW INPUT
    xra2 = xra2.sum(dim='z')  # Summing along z makes the units g/m^2
    # May want to add loops for time and ensemble dimension in the future
    # MER must used for two ensemble members
    if MER == None:
        MER = xra.attrs['Fine Ash MER']  # Allowing for MER input - default is Mastin equation MER

    xra3 = xra2
    a = 0
    while a < len(xra2.source):
        if xra2[a, :, :].source in ('CylinderSource', 'LineSource'):
            xra3[a, :, :] = xra2[a, :, :] * MER
        else:
            xra3[a, :, :] = xra2[a, :, :]
        a += 1

    xrabin = xr.where(xra3 >= threshold, 1., 0.)
    return xrabin


def topheight(draash, time, ens, level, thresh=0.01, tp=0):
    # more work needs to be done on this or there may be
    # a similar function elsewhere.

    source = 0
    if isinstance(level, int):
        level = [level]
    iii = 0
    for lev in level:
        lev_value = draash.z.values[lev]
        r2 = draash.sel(time=time)
        r2 = r2.isel(ens=ens)
        r2 = r2.isel(source=source, z=lev)
        # place zeros where it is below threshold
        r2 = r2.where(r2 >= thresh)
        r2 = r2.fillna(0)
        # place ones where it is above threshold
        r2 = r2.where(r2 < thresh)
        r2 = r2.fillna(lev_value)
        if iii == 0:
            rtot = r2
            rtot.expand_dims('z')
        else:
            r2.expand_dims('z')
            rtot = xr.concat([rtot, r2], 'z')
        rht = rtot.max(dim='z')
    return rht




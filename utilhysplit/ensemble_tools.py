""" 
manipulating xarray objects with dimensions of lat,lon,time,ens,source,z(optional)
FUNCTIONS
The functions in this file are for manipulating xarray objects which represent an ensemble of HYSPLIT runs
with concentration or mass loading information.

ATL       : applied threshold level
ATLra     : applied threshold level

APL       : applied percentile level
APLra     : applied percentile level

preprocess : stacks ens and source dimensions

listvals : returns 1d list of values in data-array

# Needs checking
    ATLra output for the problev 50 and 99 levels.
"""
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr

# _______________________________________________________________________

# 2021 Jun 2  amc some functions moved to ensemble_tools_plotting.py
# 2021 Jun 2  amc replaced choose_and_stack function with preprocess.
# 2022 Jan 10 amc moved functions for fss and afss to ensemble_stat.py to remove dependency on plume_stat.py
# 2023 Feb 28 commented out volcATL and make_ATL_hysp
# 2024 Mar 15 moved some functions to anothe file.

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
        pindex = 0
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
    dra2, dim = preprocess(indra, enslist, sourcelist)
    coords = dra2.coords
    coords2 = dict(coords)
    # remove 'ens' from the coordinate dictionary.
    coords2.pop("ens")
    coords2.pop("source")
    dims = list(dra2.dims)
    # find which dimension is the 'ens' dimension
    dii = dims.index(dim)
    # sort along the ensemble axis from smallest to greatest.
    # sort is an inplace operation. returns an empty array.
    # this type of sort only available on numpy array.
    # does not work in xarray dataarray.
    dvalues = dra2.values.copy()
    dvalues.sort(axis=dii)
    # instead of an 'ensemble' dimension, now have an 'index' dimension.
    dims[dii] = "index"
    newra = xr.DataArray(dvalues, dims=dims, coords=coords2)
    # gives percentage of concentrations which are smaller than or equal to the one given.
    percents = 100 * (newra.index.values + 1) / len(newra.index.values)
    #percents = 100-percents
    newra = newra.assign_coords(percent_level=("index", percents))
    return newra


# def volcATL(indra):
#    # 2021 Jun 2 amc was returning ATL(indra) instead of ATL(newra)
#    newra = indra.MER * indra
#    return ATL(newra)


def preprocess(indra, enslist=None, sourcelist=None):
    """
    indra : xarray dataArray
    enslist : list or np.ndarray
    sourcelist : list or np.ndarray

    Returns:
    dra : xarray datatArray: same dimensions as indra except will stack 'source' and 'ens'
          dimensions if both present.
    """

    dra = indra.copy()

    # dim : str : indicates whether 'ens' or 'source' dimension is to be used.
    # handle whether using 'ens' dimension, 'source' dimension or both.
    dim = "ens"

    # if lists are an np.ndarray then testing them with not
    # returns an error.
    if isinstance(sourcelist, np.ndarray):
       sourcelist = list(sourcelist)

    # if no sourcelist is passed and source dimension is only length 1,
    # then 'squeeze' it rather than stack it.
    elif "source" in dra.dims and len(dra.source.values) == 1:
        dra = dra.isel(source=0)
    elif "source" in dra.dims and not sourcelist:
        sourcelist = dra.source.values

    if isinstance(enslist, np.ndarray):
        pass
    # if no enslist is passed and source dimension is only length 1,
    # then 'squeeze' it rather than stack it.
    elif "ens" in dra.dims and len(dra.ens.values) == 1:
        dra = dra.isel(ens=0)
    elif "ens" in dra.dims and not enslist:
        enslist = dra.ens.values

    if "source" in dra.dims and "ens" in dra.dims:
        dra = dra.rename({"ens": "metens"})
        dra = dra.sel(metens=enslist)
        dra = dra.sel(source=sourcelist)
        dra = dra.stack(ens=("metens", "source"))
    elif "source" in dra.dims:
        dra = dra.sel(source=sourcelist)
        dim = "source"
    elif "ens" in dra.dims:
        dra = dra.sel(ens=enslist)
    else:
        # print("Warning: could not find source or ens dimension")
        dim = None
    return dra, dim


def ATLra(
    indra, enslist, sourcelist, threshlist, norm=True, weights=None, include_zero=False,
    **kwargs):

    # ------------------------------------------------
    # making frequency of exceedances
    kwargkeys = list(kwargs.keys())
    kwargkeys = [x.lower() for x in kwargkeys]

    atl_list = []
    for thresh in threshlist:
        temp = ATL(indra, enslist, sourcelist, thresh, norm, weights, include_zero)
        temp = temp.assign_coords(thresh=thresh)
        temp = temp.expand_dims("thresh")
        atl_list.append(temp)
    new = xr.concat(atl_list, dim="thresh")

    atthash = {}
    atthash["Description"] = "Ensemble relative frequency of exceedance"
    atthash["thresh unit"] = "mg/m3"
    atthash["unit"] = "Fraction of ensemble members above threshold"
    new.attrs.update(atthash)

    # ------------------------------------------------
    # making gridded concentration
    problev = 50
    slen = len(sourcelist)
    #problev = np.floor(slen/2)/slen
    #print('USING problev ', problev)
    apl = APL(indra, problev=problev,sourcelist=sourcelist,enslist=enslist)
    descrip = "Concentrations at {} percentile level".format(problev)
    atthash = {}
    atthash["Description"] = descrip
    atthash["units"] = "mg/m3"
    apl.attrs.update(atthash)

    problev = 99
    apl90 = APL(indra, problev=problev,sourcelist=sourcelist,enslist=enslist)
    descrip = "Concentrations at {} percentile level".format(problev)
    atthash = {}
    atthash["Description"] = descrip
    atthash["units"] = "mg/m3"
    apl90.attrs.update(atthash)

    # ------------------------------------------------

    # ------------------------------------------------
    # making dataset
    dhash = {}
    dhash["FrequencyOfExceedance"] = new
    dhash["Concentration"] = apl.drop('percent_level',dim=None)
    dhash["Concentration99"] = apl90.drop('percent_level',dim=None)
    dset = xr.Dataset(data_vars=dhash)

    #nhash = {}
    #for key in ['vaac','volcano name','meteorological model', 'model']:
    #    if key in kwargkeys:
    #       nhash[key] = kwargs[key]
    #    else:
    #       nhash[key] = 'unknown'


    dset.attrs.update(kwargs)

    return dset


def ATL(
    indra,
    enslist=None,
    sourcelist=None,
    thresh=0,
    norm=True,
    weights=None,
    include_zero=False,
):
    """
    Applied Threshold Level (also ensemble frequency of exceedance).

    indra: xr dataarray produced by combine_dataset or by hysp_massload
           it must have 'ens' dimension, 'source' dimension or both.
           time and z dimensions are optional.

    enslist : list of values to use for 'ens' coordinate
    sourcelist : list of values to use for 'source' coordinate

    thresh : int or float. If 0 then use > for test. otherwise use >=.
    thresh : list or npndarray of 2 floats. Number of ensemble members between the two values.

    include_zero : boolean. If true then use >= when threshold is zero.
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

    if isinstance(weights, (list, np.ndarray)):
        norm = False

    dra, dim = preprocess(indra, enslist, sourcelist)

    # allow for multiple category forecasts.
    # probability that value is between two values.
    if isinstance(thresh, list) or isinstance(thresh, np.ndarray):
        threshmax = float(thresh[1])
        thresh = float(thresh[0])
    else:
        thresh = float(thresh)
        threshmax = None

    # place ones where above threshold. zeros elsewhere.
    if thresh == 0 and not include_zero:
        dra2 = xr.where(dra >= thresh, 1, 0)
    else:
        dra2 = xr.where(dra >= thresh, 1, 0)

    # if this threshold is not None.
    if threshmax:
        # place ones where below max threshold.
        dra3 = xr.where(dra <= threshmax, 1, 0)
        # multiply with dra2 to get where above and below threshold.
        dra2 = dra3 * dra2

    if isinstance(weights, (np.ndarray, list)):
        if len(weights) != len(dra2[dim].values):
            print("WARNING : weights do not coorespond to values")
            print("weights :", weights)
            print("values : ", dra2[dim].values)
        else:
            wra = xr.DataArray(weights, dims=dim)
            dra2 = wra * dra2
    dra2 = dra2.sum(dim=[dim])
    if norm:
        nmembers = len(dra[dim].values)
        dra2 = dra2 / nmembers
    return dra2

def listvals(dra):
    """
    returns 1d list of values in the data-array
    used in ens_cdf function.
    """
    vals = dra.values
    vshape = np.array(vals.shape)
    return vals.reshape(np.prod(vshape))

def topheight(inash, time, level=None, dlev=0,enslist=None,sourcelist=None, thresh=0.01):
    """
    inash - xarray
    dlev - level thickness. To be added on to get the top height values.

    Returns
    rht - xarray with highest level that contains concentration above threshold
    rbt - xarray with lowest level that contains concentration above threshold
    """
    if 'time' in inash.dims:
        inash = inash.sel(time=time)
    revash,dim = preprocess(inash,enslist,sourcelist)
    if isinstance(level, int):
        level = [level]
    if not isinstance(level,list):
        level = np.arange(0,len(revash.z.values))

    for iii, lev in enumerate(level):
        lev_value = revash.z.values[lev]
        rr2 = revash.isel(z=lev)
        #if dim in rr2.dims:
        #   rr2 = rr2.max(dim=dim)
        # place zeros where it is below threshold
        # place level value + level thickness in meters where above threshold.
        # assume that level height indicates the bottom of the level. 
        rr2 = xr.where(rr2>thresh,lev_value+dlev,0)

        if iii == 0:
            rtot = rr2
            rtot = rtot.expand_dims('z')
        else:
            rr2 = rr2.expand_dims('z')
            rtot = xr.concat([rtot, rr2], 'z')
    rht = rtot.max(dim='z')

    rht2 = xr.where(rtot==0,1e6,rtot)
    rbt = rht2.min(dim='z')
    rbt = xr.where(rbt<1e5,rbt,0)
    return rht, rbt



def get_hull(z,thresh1=0.1,thresh2=1000):
    lat = z.longitude.values.flatten()
    lon = z.latitude.values.flatten()
    zzz = z.values.flatten()
    tlist = list(zip(lat,lon,zzz))
    print('MAX tlist', np.max(tlist))
    #thresh=0.1
    tlist = [x for x in tlist if ~np.isnan(x[2])]
    tlist = [x for x in tlist if x[2]>=thresh1]
    tlist = [x for x in tlist if x[2]<=thresh2]
    print(tlist[0:5])
    lon = [x[1] for x in tlist]
    lat = [x[0] for x in tlist]
    mpts = geotools.make_multi(lon,lat)
    ch, ep = geotools.concave_hull(mpts,alpha=1)
    return ch, ep




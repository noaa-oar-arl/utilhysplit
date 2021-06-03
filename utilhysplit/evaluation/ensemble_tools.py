import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import pandas as pd
from utilhysplit.evaluation.statmain import cdf
from utilhysplit.evaluation.statmain import pixel_matched_cdf
from utilhysplit.evaluation import plume_stat

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
    dims[dii] = "index"
    newra = xr.DataArray(dvalues, dims=dims, coords=coords2)
    percents = 100 * (newra.index.values + 1) / len(newra.index.values)
    newra = newra.assign_coords(percent_level=("index", percents))
    return newra


def volcATL(indra):
    # 2021 Jun 2 amc was returning ATL(indra) instead of ATL(newra)
    newra = indra.MER * indra
    return ATL(newra)


def preprocess(indra, enslist=None, sourcelist=None):
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
    dim = "ens"

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

    if "source" in dra.dims and "ens" in dra.dims:
        dra = dra.sel(ens=enslist)
        dra = dra.sel(source=sourcelist)
        dra = dra.stack(ens=("ens", "source"))

    elif "source" in dra.dims:
        dra = dra.sel(source=sourcelist)
        dim = "source"
    elif "ens" in dra.dims:
        dra = dra.sel(ens=enslist)
    else:
        print("Warning: could not find source or ens dimension")
        dim = None
    return dra, dim


def ATL(indra, enslist=None, sourcelist=None, thresh=0, norm=False, weights=None):
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
    dra, dim = preprocess(indra, enslist, sourcelist)

    if thresh == 0:
        dra2 = xr.where(dra > thresh, 1, 0)
    else:
        dra2 = xr.where(dra >= thresh, 1, 0)

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
        nmembers = len(enslist)
        dra2 = dra2 / nmembers
    return dra2


def listvals(dra):
    """
    returns 1d list of values in the data-array
    """
    vals = dra.values
    vshape = np.array(vals.shape)
    return vals.reshape(np.prod(vshape))

def ens_time_fss(
    indralist,
    obsralist,
    enslist=None,
    sourcelist=None,
    #timelist=None,
    neighborhoods = [1,3,5,7],
    threshold=0,
    plot=True,
    #pixel_match=None,
):
    """
    RETURNS
    pandas dataframe with columns
    Nlen, FBS, FBS_ref, FSS, ens, time
    """
    dflist = []
    for pairs in zip(indralist, obsralist):
        df = ens_fss(pairs[0],pairs[1],enslist,sourcelist,neighborhoods,
                     threshold,plot)
        dflist.append(df)
    return pd.concat(dflist)

def ens_fss(
    indra,
    obsra,
    enslist=None,
    sourcelist=None,
    #timelist=None,
    neighborhoods = [1,3,5,7],
    threshold=0,
    plot=True,
    #pixel_match=None,
):
    """
    indra and obsra need to be same time period and grid.

    RETURNS
    pandas dataframe with columns
    Nlen, FBS, FBS_ref, FSS, ens, time

    """
    dra, dim = preprocess(indra, enslist, sourcelist)
    dflist = []
    for ens in dra[dim].values:
        if dim == "ens":
            subdra = dra.sel(ens=ens)
        elif dim == "source":
            subdra = dra.sel(source=ens)
        scores = plume_stat.CalcScores(obsra, subdra,threshold=threshold)
        df1 = scores.calc_fss(makeplots=False,szra=neighborhoods)
        df1['ens'] = ens
        dflist.append(df1)

    # calculate fss for ensemble mean
    meanra = dra.mean(dim=dim)
    mean_scores = plume_stat.CalcScores(obsra, meanra,threshold=threshold)
    df1 = scores.calc_fss(makeplots=False,szra=neighborhoods)
    df1['ens'] = 'mean'
    # add time to dataframe.
    dflist.append(df1)
    dfall = pd.concat(dflist)
    if 'time' in indra.coords:
        dfall['time'] = pd.to_datetime(indra.coords['time'].values)
    return dfall

def plot_ens_fss_ts(ensdf, nval=5, sizemult=1, enslist=None):
    """
    Plot FSS on y and time on x.
    ensdf : pandas DataFrame output from  ens_time_fss
    nval : neighborhood size to plot.
    """
    if nval:
       tempdf = ensdf[ensdf['Nlen']==nval]
    fig, ax = plt.subplots(1,1)

    ensfss = tempdf.pivot(columns='ens',values='FSS',index='time')
    uniform = tempdf.pivot(columns='ens',values='uniform',index='time')
    ensfss.plot(ax=ax, legend=None,colormap='tab20')
    colA = uniform.columns[0] 
    uniform.plot(ax=ax, y=colA, legend=None,colormap='winter')
    if 'mean' in ensfss.columns:
        ensfss.plot(ax=ax, y='mean',LineWidth=5,colormap="winter")

def plot_ens_fss(ensdf, sizemult=1, timelist=None, enslist=None, nlist=None):
    """
    Plot FSS on y and Neighborhood length on x.

    ensdf : pandas DataFrame output from ens_fss function
    sizemult : set to value other than one to convert to degrees.
    """
    fig, ax = plt.subplots(1,1)
    if sizemult != 1:
       ensdf['length (degrees)'] = ensdf['Nlen']*sizemult
       ensfss = ensdf.pivot(columns='ens',values='FSS',index='length (degrees)')
    else:
       ensfss = ensdf.pivot(columns='ens',values='FSS',index='Nlen')
    ensfss.plot(ax=ax, legend=None,colormap='spring')
    random = ensdf['random'].unique()
    nmin = float(np.min(ensdf['Nlen']))* sizemult
    nmax = float(np.max(ensdf['Nlen']))* sizemult
    if 'mean' in ensfss.columns:
        ensfss.plot(ax=ax, y='mean',LineWidth=5,colormap="winter")
    # plot random forecast
    for randomval in random:
        plt.plot([nmin,nmax],[randomval, randomval], '--k')
    uniform = ensdf['uniform'].unique()
    # plot uniform forecast
    for uniformval in uniform:
        plt.plot([nmin,nmax],[uniformval, uniformval], '--r')
     


def ens_cdf(
    indra,
    enslist=None,
    sourcelist=None,
    timelist=None,
    threshold=0,
    plot=True,
    pixel_match=None,
):
    """
    produces plots of cumulative distribution functions.
    indra : xarray DataArray produced by combine_dataset function or hysp_massload function..
    timelist : list of times in the time coordinate to produce plots for
    threshold : float : produce CDF for values > threshold.
    pixel_match : Number of pixels in observation.
                  if not None then will use this instead of threshold.
                  Cut length of modeled data to same length as observed. 
    Returns:
    cdfhash : dictionary. key is the time. value is a tuple of the CDF (x,y)
    """
    # select sources of interest and stack
    if "source" in indra.dims:
        sourcekey = True
    else:
        sourcekey = False
    dra, dim = preprocess(indra, enslist, sourcelist)
    cdfhash = {}
    if plot:
        fig = plt.figure(1)
        ax = fig.add_subplot(1, 1, 1)
    
    # loop through ens/source members
    for ens in dra[dim].values:
        if dim == "ens":
            subdra = dra.sel(ens=ens)
        elif dim == "source":
            subdra = dra.sel(source=ens)

        # check if time is a dimension.
        if not isinstance(timelist, (list, np.ndarray)) and 'time' in subdra.dims:
            timelist = subdra.time.values
        elif 'time' not in subdra.dims:
            timelist = [0]

        # loop through times. If no time, then just go through once.
        for tm in timelist:
            # check if time is a dimension.
            if 'time' in subdra.dims:
                tvals = subdra.sel(time=tm)
            else:
                tvals = subdra 
            # create cdf from values above threshold
            if not pixel_match:
                sdata, y = cdf([x for x in listvals(tvals) if x > threshold])
            # create cdf from highest pixel_match values.
            else:
                sdata, y = pixel_matched_cdf(listvals(tvals), pixel_match)

            # create the key for the dictionary.
            if sourcekey:
                key = (tm, ens[0], ens[1])
            else:
                key = (tm, ens)
            cdfhash[key] = (sdata, y)
    if plot:
        plot_cdf(ax, cdfhash)
    return cdfhash


def plot_cdf(ax1, cdfhash):
    """
    Plots output from ens_cdf
    """
    clrs = ["r", "y", "g", "c", "b", "k"]
    for iii, key in enumerate(cdfhash.keys()):
        if iii > len(clrs) - 1:
            clrs.extend(clrs)
        ax1.step(cdfhash[key][0], cdfhash[key][1], "-" + clrs[iii])
    ax1.set_xscale("log")
    return ax1


def make_ATL_hysp(xra, variable="p006", threshold=0.0, MER=None):
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

    xra2 = (
        xra[variable] * 1000.0
    )  # Each vertical layer is 1000m - MAY WANT TO ALLOW INPUT
    xra2 = xra2.sum(dim="z")  # Summing along z makes the units g/m^2
    # May want to add loops for time and ensemble dimension in the future
    # MER must used for two ensemble members
    if not MER:
        MER = xra.attrs[
            "Fine Ash MER"
        ]  # Allowing for MER input - default is Mastin equation MER

    xra3 = xra2
    a = 0
    while a < len(xra2.source):
        if xra2[a, :, :].source in ("CylinderSource", "LineSource"):
            xra3[a, :, :] = xra2[a, :, :] * MER
        else:
            xra3[a, :, :] = xra2[a, :, :]
        a += 1

    xrabin = xr.where(xra3 >= threshold, 1.0, 0.0)
    return xrabin


def topheight(draash, time, ens, level, thresh=0.01):
    # more work needs to be done on this or there may be
    # a similar function elsewhere.

    source = 0
    if isinstance(level, int):
        level = [level]
    iii = 0
    for lev in level:
        lev_value = draash.z.values[lev]
        rr2 = draash.sel(time=time)
        rr2 = rr2.isel(ens=ens)
        rr2 = rr2.isel(source=source, z=lev)
        # place zeros where it is below threshold
        rr2 = rr2.where(rr2 >= thresh)
        rr2 = rr2.fillna(0)
        # place ones where it is above threshold
        rr2 = rr2.where(rr2 < thresh)
        rr2 = rr2.fillna(lev_value)
        if iii == 0:
            rtot = rr2
            rtot.expand_dims("z")
        else:
            rr2.expand_dims("z")
            rtot = xr.concat([rtot, rr2], "z")
        rht = rtot.max(dim="z")
    return rht

""" manipulating xarray objects with dimensions of lat,lon,time,ens,source,z(optional)
FUNCTIONS
The functions in this file are for manipulating xarray objects with concentration or mass loading information.


ens_cdf   : cumulative distribution functions.
plot_cdf  : plots outputs from ens_cdf

listvals : returns 1d list of values in data-array

get_pixel_match : finds threshold which would result in same number of pixels in two input arrays.

# Needs documentation 
    plot_ens_area
    plot_ens_accuracy
    ens_boxplot
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr
from utilhysplit.evaluation import hysplit_boxplots
from utilhysplit.evaluation.statmain import (
    cdf,
    get_pixel_matching_threshold,
    pixel_matched_cdf,
)

# _______________________________________________________________________

# 2021 Jun 2  amc some functions from ensemble_tools.py to here

def listvals(dra):
    """
    returns 1d list of values in the data-array
    used in ens_cdf function.
    """
    vals = dra.values
    vshape = np.array(vals.shape)
    return vals.reshape(np.prod(vshape))


def get_pixel_match(indra, obsra, thresh, return_binary=False):
    """
    Counts how many above threshold values in obsra.
    Sorts indra values form least to greatest.
    Finds threshold for indra ense members  which would result in same number of
    above threshold values as in obsra.

    Stores this threshold value in matchra and returns.
    applies thresholds to indra to return matchra which have same number of
    above threshold pixels as obsra.

    Inputs:
    Outputs:
    threshra : xarray dataArray with threshold for each ensemble value.
    matchra  : indra with
    """
    dra, dim = preprocess(indra)
    threshlist = []
    if dim:
        for ens in dra[dim].values:
            if dim == "ens":
                subdra = dra.sel(ens=ens)
            elif dim == "source":
                subdra = dra.sel(source=ens)
            pm_thresh = get_pixel_matching_threshold(obsra, subdra, thresh)
            threshlist.append(pm_thresh)
    else:
        subdra = dra
        pm_thresh = get_pixel_matching_threshold(obsra, subdra, thresh)
        threshlist.append(pm_thresh)
    threshra = xr.DataArray(threshlist, dims=dim)
    if return_binary:
        matchra = xr.where(indra >= threshra, 1, 0)
    else:
        matchra = xr.where(indra >= threshra, indra, 0)
    return threshra, matchra


# used in Bezy project
def plot_ens_area(
    ensdfin, ax=None, plotmean=False, legend=False, clrlist=None, enslist=None
):
    sns.set()
    sns.set_style("whitegrid")
    if not ax:
        fig, ax = plt.subplots(1, 1)
    sns.set_style("whitegrid")
    ensdf = ensdfin.copy()
    if "time" in ensdf.columns:
        val = ensdf.pivot(columns="ens", values="area_fc", index="time")
        if isinstance(enslist, list):
            val = val[enslist]
        obs = ensdf.pivot(columns="ens", values="area_obs", index="time")
        obs = obs[obs.columns[0]]
        ax.plot(
            obs.index.values,
            obs.values,
            linestyle="--",
            linewidth=10,
            alpha=0.5,
            label="obs",
        )
        if not clrlist:
            val.plot(ax=ax, legend=None, colormap="tab20")
        else:
            val.plot(ax=ax, legend=None, color=clrlist, alpha=0.5)
        # if 'mean' in val.columns and plotmean:
        #    val.plot(ax=ax, y='mean', legend=None, linewidth=5,colormap='winter')
    else:
        val = ensdf.pivot(columns="ens", values=cname, index="time")
        val = ensdf.pivot(columns="ens", values=cname, index="time")
    ax.set_ylabel("Area (number of pixels)")
    if legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels)
    return ax

# used in Bezy project
def plot_ens_accuracy(
    ensdfin, cname="MAE", plotmean=True, legend=False, clrlist=None, enslist=None
):
    if cname == "RMSE":
        rvalue = "RMSE"
        cname = "MSE"
    else:
        rvalue = cname
    ensdf = ensdfin.copy()
    sns.set()
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(1, 1)
    sns.set_style("whitegrid")
    if "time" in ensdf.columns:
        val = ensdf.pivot(columns="ens", values=cname, index="time")
        # this is for re-ordering the columns.
        if isinstance(enslist, list):
            val = val[enslist]
        if rvalue == "RMSE":
            val = val ** 0.5
        if clrlist:
            val.plot(ax=ax, legend=None, color=clrlist)
        else:
            val.plot(ax=ax, legend=None, colormap="tab20")
        # if 'mean' in val.columns and plotmean:
        #    val.plot(ax=ax, y='mean', legend=None, linewidth=5,colormap='winter')

    else:
        val = ensdf.pivot(columns="ens", values=cname, index="time")
        val = ensdf.pivot(columns="ens", values=cname, index="time")
    ax.set_ylabel(rvalue)
    if legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels)
    return ax


def ens_boxplot(
    indra,
    enslist=None,
    sourcelist=None,
    timelist=None,
    threshold=0,
    # plot=True,
    pixel_match=None,
    clist=None
    # xscale='log',
    # colors=None,
    # label='time'
):
    if "source" in indra.dims:
        sourcekey = True
    else:
        sourcekey = False
    dra, dim = preprocess(indra, enslist, sourcelist)

    for ens in dra[dim].values:
        vdata = []
        print(ens)
        if dim == "ens":
            subdra = dra.sel(ens=ens)
        elif dim == "source":
            subdra = dra.sel(source=ens)
        # check if time is a dimension.
        if not isinstance(timelist, (list, np.ndarray)) and "time" in subdra.dims:
            timelist = subdra.time.values
        elif "time" not in subdra.dims:
            timelist = [0]
        # loop through times. If no time, then just go through once.
        for tm in timelist:
            # check if time is a dimension.
            if "time" in subdra.dims:
                tvals = subdra.sel(time=tm)
            else:
                tvals = subdra
            # create list of above threshold values.
            vdata.append([x for x in listvals(tvals) if x > threshold])
            # create cdf from highest pixel_match values.
            # else:
            #    sdata, y = pixel_matched_cdf(listvals(tvals), pixel_match)
        print("boxplot", len(vdata), len(timelist))
        dj = hysplit_boxplots.prepare_boxplotdata(timelist, vdata)
        hysplit_boxplots.make_boxplot(dj, cols=clist)


def ens_cdf(
    indra,
    enslist=None,
    sourcelist=None,
    timelist=None,
    threshold=0,
    plot=True,
    pixel_match=None,
    xscale="log",
    colors=None,
    label="time",
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
        if not isinstance(timelist, (list, np.ndarray)) and "time" in subdra.dims:
            timelist = subdra.time.values
        elif "time" not in subdra.dims:
            timelist = [0]

        # loop through times. If no time, then just go through once.
        for tm in timelist:
            # check if time is a dimension.
            if "time" in subdra.dims:
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
        plot_cdf(ax, cdfhash, xscale, clrs=colors, label=label)
    return cdfhash


def plot_cdf(ax1, cdfhash, xscale="log", clrs=None, label="time"):
    """
    Plots output from ens_cdf
    """
    if isinstance(clrs, list):
        clrs = clrs
    else:
        clrs = ["r", "y", "g", "c", "b", "k"]
    for iii, key in enumerate(cdfhash.keys()):
        if iii > len(clrs) - 1:
            clrs.extend(clrs)
        if label == "time":
            lname = key[0]
        elif label == "ens":
            lname = key[1]
        ax1.step(cdfhash[key][0], cdfhash[key][1], ls="-", color=clrs[iii], label=lname)
    ax1.set_xscale(xscale)
    return ax1





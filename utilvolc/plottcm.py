import datetime
import logging
import os
from subprocess import call

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
import seaborn as sns
import xarray as xr

import utilhysplit.hysplit_gridutil as hgu
from monetio.models import hysplit
from utilhysplit import hcontrol
from utilhysplit.plotutils import colormaker
from utilvolc import volcat
from utilvolc.runhelper import Helper


logger = logging.getLogger(__name__)




def plot_outdat_profile_psize_function(
    dfdat, log=False, fignum=1, unit="kg/s", ax=None, cmap="viridis"
):
    # plots time series of MER. summed along column.
    # dfdat : pandas dataframe output by make_outdat_df function with
    #        part='basic'
    if not ax:
        sns.set()
        sns.set_style("whitegrid")
        fig = plt.figure(fignum, figsize=(10, 5))
        ax = fig.add_subplot(1, 1, 1)
    if isinstance(dfdat, str):
        df = pd.read_csv(dfdat, index_col=0, header=1)
        df = df.dropna()
    else:
        df = dfdat
    sns.set()
    plist = df.psize.unique()
    cm = colormaker.ColorMaker(cmap, len(plist), ctype="hex", transparency=None)
    cmlist = cm()
    for iii, psize in enumerate(plist):
        dfp = df[df.psize == psize]
        dfp = dfp.drop("psize", axis=1)
        # dfp = dfp.pivot(index='date',columns='ht')
        dfp = dfp.pivot(index="date", columns="ht")
        try:
            dfp = dfp.mass
        except:
            pass
        ts = dfp.sum()
        xval = ts.index.values * 1 / 1.0e3
        yvals = ts.values / 1.0e12
        ax.plot(yvals, xval, "#" + cmlist[iii], label=psize)
    ax.set_xlabel("Tg of mass emitted", fontsize=15)
    ax.set_ylabel("Height (km)", fontsize=15)
    # ax.set_ylabel('Mass {}'.format(unit),fontsize=15)
    return ax, ts


def plot_outdat_ts_psize_function(
    dfdat, log=False, fignum=1, unit="kg/s", ax=None, cmap="viridis"
):
    # plots time series of MER. summed along column.
    # dfdat : pandas dataframe output by make_outdat_df function with
    #        part='basic'

    if not ax:
        sns.set()
        sns.set_style("whitegrid")
        fig = plt.figure(fignum, figsize=(10, 5))
        ax = fig.add_subplot(1, 1, 1)
    if isinstance(dfdat, str):
        df = pd.read_csv(dfdat, index_col=0, header=1)
        df = df.dropna()
    else:
        df = dfdat
    sns.set()
    plist = df.psize.unique()
    cm = colormaker.ColorMaker(cmap, len(plist), ctype="hex", transparency=None)
    cmlist = cm()
    for iii, psize in enumerate(plist):
        dfp = df[df.psize == psize]
        dfp = dfp.drop("psize", axis=1)
        # dfp = dfp.pivot(index='date',columns='ht')
        dfp = dfp.pivot(index="ht", columns="date")
        try:
            dfp = dfp.mass
        except:
            pass
        ts = dfp.sum()
        if unit == "kg/s":
            yval = ts.values / 3.6e6
        elif unit == "g/h":
            yval = ts.values
        print(type(ts.index.values[0]), ts.index.values[0])
        xval = [pd.to_datetime(x) for x in ts.index.values]
        # ax.plot([x[0] for x in ts.index.values], yval, clr)
        print(cmlist[iii])
        ax.plot(xval, yval, "#" + cmlist[iii], label=psize)
        # fig.autofmt_xdate()
        ax.set_ylabel("MER {}".format(unit), fontsize=15)
    return ax, ts


def plot_emissions_timeseries(
    dfdat,
    log=False,
    fignum=1,
    unit="kg/s",
    ax=None,
    clr="k",
    label=None,
    alpha=1,
    lw=1,
    marker=None,
):
    # plots time series of MER. summed along column.
    # dfdat : pandas dataframe output by InverseOutDat class get_emis method.
    if not ax:
        sns.set()
        sns.set_style("whitegrid")
        fig = plt.figure(fignum, figsize=(10, 5))
        ax = fig.add_subplot(1, 1, 1)
    if isinstance(dfdat, str):
        df = pd.read_csv(dfdat, index_col=0, header=1)
        df = df.dropna()
    else:
        df = dfdat
    sns.set()
    if 'psize' in df.columns: df = df.drop("psize", axis=1)
    # dfp = dfp.pivot(index='date',columns='ht')
    df = df.pivot(index="date", columns="ht")
    try: 
       df = df.mass
    except: 
       pass
    ts = df.sum(axis=1)
    print('HERE', ts)
    if unit == "kg/s":
        yval = ts.values / 3.6e6
    elif unit == "g/h":
        yval = ts.values
        
    xval = [pd.to_datetime(x) for x in ts.index.values]
    #xval = [x for x in ts.index.values]
    print(xval)
    # ax.plot([x[0] for x in ts.index.values], yval, clr)
    ax.plot(
        xval, yval, label=label, color=clr, linestyle="-", marker=marker, alpha=alpha, linewidth=lw
    )
    # fig.autofmt_xdate()
    ax.set_ylabel("MER {}".format(unit), fontsize=15)
    return ax, df


def plot_outdat_profile_function(
    dfdat, fignum=1, unit="km", ax=None, clr="k", label=None, alpha=1, lw=1, marker="o"
):
    if not ax:
        sns.set()
        sns.set_style("whitegrid")
        fig = plt.figure(fignum, figsize=(10, 5))
        ax = fig.add_subplot(1, 1, 1)

    if isinstance(dfdat, str):
        df = pd.read_csv(dfdat, index_col=0, header=1)
        df = df.dropna()
    else:
        df = dfdat
 
    if 'psize' in df.columns: df = df.drop("psize", axis=1)
    # dfp = dfp.pivot(index='date',columns='ht')
    df2 = df.pivot(index="ht", columns="date")
    ts = df2.sum(axis=1)
    sns.set()
    xval = ts.values * 1 / 1e12
    try:
        yval = ts.index.values / 1000.0
    except:
        yval = list(map(float, list(ts.index.values)))
        yval = np.array(yval) / 1000.0
    if unit == "FL":
        yval = yval * 3280.84 / 100.0
    ax.plot(
        xval,
        yval,
        color=clr,
        marker=marker,
        label=label,
        alpha=alpha,
        linewidth=lw,
        linestyle="-",
    )
    ax.set_xlabel("Tg of mass emitted", fontsize=15)
    if unit == "FL":
        ax.set_ylabel("Height (FL)", fontsize=15)
    else:
        ax.set_ylabel("Height (km)", fontsize=15)
    totalmass = xval.sum()
    # print('total {} Tg'.format(totalmass))
    return ax, totalmass



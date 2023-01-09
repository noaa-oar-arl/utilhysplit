import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from utilhysplit.evaluation import web_ensemble_plots as wep
from utilhysplit.fixlondf import fixlondf
from utilhysplit.plotutils import map_util


def plottraj_map(tdump, clr="-c", ms=0.5, alpha=0.5, central_longitude=180):
    """
    central_longitude should be 180 if trajectories are crossing dateline.
    Otherwise use 0.
    """
    transform = wep.get_transform(central_longitude)
    fig, axra = plt.subplots(
        nrows=1,
        ncols=1,
        figsize=(10, 10),
        constrained_layout=False,
        subplot_kw={"projection": transform},
    )
    plottraj(tdump, axra, clr, alpha, ms)
    wep.format_plot(axra, transform)
    return axra


def plottraj(tdump2, ax, clr="-c", alpha=0.5, ms=0.5):
    sns.set()
    tnum = tdump2.traj_num.unique()
    tlist = tnum
    for tnum in tlist:
        temp = tdump2[tdump2["traj_num"] == tnum]
        temp = temp.sort_values(by="time")
        alt = temp.iloc[0].altitude
        xval = temp["longitude"].values
        xval = wep.shift_xvals(xval, 180)
        ax.plot(xval, temp["latitude"], clr, linewidth=0.1, markersize=ms, alpha=alpha)
    return ax


def frequency_plots_all(df, sdate, dtt, dres=1, vloc=None):
    """
    df : pandas aataframe
    sdate : datetime object
    dtt :
    """
    times = df.time.unique()
    times.sort()
    done = False
    if isinstance(dtt, int):
        dtt = datetime.timedelta(hours=dtt)
    while not done:
        onetime = df[df["time"] == sdate]
        onetime = fixlondf(onetime, colname="longitude", neg=False)
        frequency_plot(onetime, dres, vloc)
        plt.title(sdate)
        sdate += dtt
        if sdate > pd.to_datetime(times[-1]):
            print("done", sdate)
            print(times)
            done = True
        plt.show()


def frequency_plot(df, dres=1, vloc=None):
    """
    df : pandas dataframe with columns of longitude, latitude.
         weight column is optional. If it does not exist then all points are
         evenly weighted.
    dres : resolution of bins in degrees.
    """
    sns.set()
    if "weight" not in df.columns:
        df["weight"] = 1
    xmin = np.floor(np.min(df.longitude.values))
    xmax = np.ceil(np.max(df.longitude.values))
    ymin = np.floor(np.min(df.latitude.values))
    ymax = np.ceil(np.max(df.latitude.values))

    xlon = np.arange(int(xmin), int(xmax), dres)
    ylat = np.arange(int(ymin), int(ymax), dres)
    onetime = df.copy()
    onetime = fixlondf(onetime, colname="longitude", neg=False)
    hhh, xedges, yedges = np.histogram2d(
        onetime.longitude.values,
        onetime.latitude.values,
        weights=onetime.weight,
        bins=[xlon, ylat],
    )
    cb = plt.pcolormesh(xedges, yedges, hhh.T, cmap="viridis")
    plt.colorbar(cb)
    if isinstance(vloc, (list, tuple)):
        plt.plot(360 + vloc[1], vloc[0], "y^", markersize=10)
    return hhh, xedges, yedges

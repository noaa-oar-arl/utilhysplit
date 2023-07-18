""" creates ensemble relative frequency of exceedance and ensemble mean plots

Classes
    LabelData

Functions
    label_ax
    get_transform
    decide_central_longitude
    shift_sub
    shift_xvals
    set_levels
    massload_plot
    sub_massload_plot
    ATLtimeloop
    check_thresh
    set_ATL_text
    plot_ATL
    meterv2FL
    setup_plot
    reset_plots
    format_plot
    ATLsubplot
    massload_ensemble_mean

"""

import datetime
import logging

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.colors import BoundaryNorm

import cartopy
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from monetio.models import hysplit
from utilhysplit.evaluation.ensemble_tools import ATL
from utilhysplit.plotutils.colormaker import ColorMaker

logger = logging.getLogger(__name__)

def setup_logger_warning(level=logging.WARNING):
    MESSAGE_FORMAT_INFO = "%(asctime)s.%(msecs)03d %(levelname)s - %(message)s"
    MESSAGE_FORMAT_DEBUG = \
    "%(asctime)s.%(msecs)03d %(levelname)s %(name)s - %(message)s"
    import sys
    logging.basicConfig(
            stream=sys.stdout,
            level=level,
            format=MESSAGE_FORMAT_INFO,
            datefmt="%H:%M:%S")


class LabelData:
    def __init__(self, time, descrip, units, source="", tag=""):
        """
        time : str, numpy datetime64, or datetime.datetime
        descript : str
        """
        if type(time) == str:
            self.time = time
        else:
            timefmt = "%m/%d/%Y %H:%M UTC"
            self.model = "NOAA HYSPLIT"
            # somtetimes time is a numpy datetime64  object.
            # convert to datetime.
            if type(time) == np.datetime64:
                time = pd.to_datetime(str(time))
            self.time = time.strftime(timefmt)
            temp2 = time + datetime.timedelta(hours=1)
            self.time = "{} to {}".format(self.time, temp2.strftime(timefmt))

        self.size = 15
        self.descrip = descrip
        self.units = "units: {}".format(units)
        self.source = "source: {}".format(source)
        self.tag = tag


def label_ax(ax, label, transform):
    # add labels to plot meant for labels.
    gl = ax.gridlines(
        crs=transform, draw_labels=False, linewidth=0, color="gray", alpha=0
    )
    gl.top_labels = False
    gl.bottom_labels = False
    gl.right_labels = False
    gl.left_labels = False
    size = label.size
    va = "bottom"
    ha = "left"
    rot = "horizontal"
    xloc = 0.05
    yloc = 0.5
    yspc = 0.1
    ax.text(
        xloc,
        yloc,
        label.model,
        va=va,
        ha=ha,
        rotation=rot,
        rotation_mode="anchor",
        transform=ax.transAxes,
        size=size,
    )
    ax.text(
        xloc,
        yloc - yspc,
        label.time,
        va=va,
        ha=ha,
        rotation=rot,
        rotation_mode="anchor",
        transform=ax.transAxes,
        size=size,
    )
    ax.text(
        xloc,
        yloc - 2 * yspc,
        label.descrip,
        va=va,
        ha=ha,
        rotation=rot,
        rotation_mode="anchor",
        transform=ax.transAxes,
        size=size,
    )
    ax.text(
        xloc,
        yloc - 3 * yspc,
        label.units,
        va=va,
        ha=ha,
        rotation=rot,
        rotation_mode="anchor",
        transform=ax.transAxes,
        size=size,
    )
    ax.text(
        xloc,
        yloc - 4 * yspc,
        label.source,
        va=va,
        ha=ha,
        rotation=rot,
        rotation_mode="anchor",
        transform=ax.transAxes,
        size=size,
    )
    ax.text(
        xloc,
        yloc - 5 * yspc,
        label.tag,
        va=va,
        ha=ha,
        rotation=rot,
        rotation_mode="anchor",
        transform=ax.transAxes,
        size=size,
    )


def get_transform(central_longitude=0):
    transform = cartopy.crs.PlateCarree(central_longitude=central_longitude, globe=None)
    # transform = cartopy.crs.AzimuthalEquidistant(central_longitude=180)
    return transform


def decide_central_longitude(xorg):
    # currently works when lon values are 0-360.
    min_x = np.min(xorg)
    max_x = np.max(xorg)
    # see if values span 180.
    if min_x < 180 and max_x > 180:
        central_longitude = 180
    else:
        central_longitude = 0
    return central_longitude


def shift_sub(x):
    if x > 0:
        newx = x - 180
    else:
        newx = x + 180
    if newx > 360:
        return 360 - newx
    else:
        return newx


def shift_xvals(xorg, central_longitude):
    if central_longitude == 0:
        return xorg
    xnew = np.vectorize(shift_sub)(xorg)
    # xnew = np.vectorize(shift_sub)(xorg - central_longitude)
    # xnew = [x + central_longitude  for x in xorg]
    # xnew = [x-360 if x > 360 else x for x in xnew]
    return xnew


def set_levels(mass):
    levels = [0.2, 2, 5, 10]
    if np.max(mass) < 0.01:
        levels = [0.0005, 0.001, 0.005, 0.01]
    elif np.max(mass) < 0.1:
        levels = [0.005, 0.01, 0.05, 0.1]
    elif np.max(mass) < 1:
        levels = [0.01, 0.05, 0.1, 0.5, 1]
    elif np.max(mass) < 11:
        levels = levels
    elif np.max(mass) < 51:
        levels.extend([20, 30, 40, 50])
    elif np.max(mass) < 101:
        levels.extend([25, 50, 75, 100])
    elif np.max(mass) < 501:
        levels.extend([50, 100, 250, 500])
    elif np.max(mass) < 751:
        levels.extend([50, 100, 250, 500, 750])
    else:
        mlev = np.max([1000, int(np.max(mass))])
        levels.extend([50, 100, 500, 750, mlev])
    return levels


def massload_plot(revash, enslist, name="None", vlist=None):
    setup_plot()
    mass = massload_ensemble_mean(revash, enslist)
    transform = get_transform()
    iii = 0
    fignamelist = []
    for time in mass.time.values:
        source = "Longitude: {:0.2f} Latitude {:0.2f}".format(vlist[0], vlist[1])
        label = LabelData(
            time, "Column mass loading (ensemble mean)", "g/m$^2$", source
        )

        mass2 = mass.sel(time=time).isel(source=0)
        x = mass2.longitude
        y = mass2.latitude
        # z = mass2.where(mass2!=0)
        central_longitude = decide_central_longitude(x)
        # if central_longitude !=0:
        #    x = shift_xvals(x.values,central_longitude)
        #     if iii==0: vlist[0] = shift_xvals(vlist[0],central_longitude)
        transform = get_transform(central_longitude)

        figname = name.replace("zzz", "{}".format(iii))
        z = mass2 / 1000.0  # convert to g/m2
        levels = set_levels(z)
        try:
            sub_massload_plot(x, y, z, transform, levels, label, figname, vlist)
        except ValueError as ex:
            print("ERROR: cannot create massload plot at {}: {}".format(time, str(ex)))
        iii += 1
        fignamelist.append(figname)
    return fignamelist


def sub_massload_plot(x, y, z, transform, levels, labeldata, name="None", vlist=None):
    nrow = 2
    ncol = 1  # second row is for text information.
    fig, axra = plt.subplots(
        nrows=nrow,
        ncols=ncol,
        figsize=(10, 10),
        constrained_layout=False,
        subplot_kw={"projection": transform},
    )
    ax = axra[0]
    dax = axra[1]
    # cmap = plt.get_cmap('tab20b')
    cm = ColorMaker("plasma", len(levels), ctype="rgb")
    clrs = cm()
    # norm = BoundaryNorm(levels,ncolors=cmap.N,clip=False)
    cb2 = ax.contourf(x, y, z, levels=levels, colors=clrs, transform=transform)
    # cb2 = ax.pcolormesh(x,y,z,cmap=cmap,transform=transform,norm=norm)
    if isinstance(vlist,(np.ndarray,list,tuple)): 
        ax.plot(vlist[0], vlist[1], "r^")
    format_plot(ax, transform)
    label_ax(dax, labeldata, transform)
    cb = plt.colorbar(cb2)
    cb.set_label("g/m$^2$")
    if name != "None":
        plt.savefig(name)
        plt.close()


def ATLtimeloop(
    revash,
    enslist,
    thresh,
    vlist=None,
    name="None",
    norm=True,
    clevels=[1, 10, 20, 40, 60, 80, 95],
    title="HYSPLIT Applied Threshold Levels",
    adjust=0,
):
    """
    creates plot for each time period.
    revash : xarray DataArray
    vlist : list of volcano coordinates to plot [longitude,latitude]
    name : str
           name for saving figures. if include zzz in name will replace
           with number for each time period.

    adjust: if >0 then will adjust threshold downwards if max concentration
            below the threshold which would result in empty plots.
            New threshold will be maxval / adjust.

    norm : boolean
           if True plot percentage of members.
           if False plot number of members.
           
    """
    fignamelist = []
    iii = 0
    if adjust > 0:
        thresh = check_thresh(np.max(revash), thresh, adjust)

    for time in revash.time.values:
        revash2 = revash.sel(time=time)
        source=''
        if isinstance(vlist,(np.ndarray,list,tuple)): 
            source = "Longitude: {:0.2f} Latitude {:0.2f}".format(vlist[0], vlist[1])
        label = LabelData(time, "Probability of exceedance", "%", source)
        rtot = ATL(revash2, enslist=enslist, thresh=thresh,norm=norm)
        # figname = name.replace('zzz',"{:02d}".format(iii))
        figname = name.replace("zzz", "{}".format(iii))
        if norm:
            rtot = rtot * 100
        title2 = "{}\n{}".format(title, label.time)
        plotATL(rtot, vlist, name=figname, levels=clevels, thresh=thresh, title=title2)
        iii += 1
        fignamelist.append(figname)
    reset_plots()
    return fignamelist


def check_thresh(cmax, thresh, adjust):
    """
    set threshold at half the maximum value.
    """
    thresh2 = thresh
    if cmax < thresh:
        thresh2 = float(cmax / float(adjust))
    return thresh2


def set_ATL_text(nrow):
    if nrow >= 6:
        yplace = 0.8
        xplace = 0.95
        size = 10
    else:
        yplace = 0.90
        xplace = 0.95
        size = 20
    return xplace, yplace, size


def plotATL(
    rtot, vlist, name="None", levels=[1, 5, 10, 20], thresh=0.2, title="HYSPLIT"
):
    """
    vlist : [latitude, longitude] of volcano location
    plot ensemble relative frequency for concentration.
    creates plot for one time period.
    creates subplot for each vertical level.
    """

    setup_plot()
    x = rtot.longitude
    y = rtot.latitude

    if isinstance(vlist,(np.ndarray,list,tuple)): 
        vlist_copy = vlist.copy()
    else:
        vlist_copy = None

    temp = rtot.max(dim="z")
    zvals = rtot.z.values
    nz = len(zvals)
    if nz > 1:
        ncol = 2
        nrow = int(np.ceil(nz / ncol))
    else:
        nrow = 1
        ncol = 1
    print("NROW, NCOL", nrow, ncol)
    z = temp.where(temp != 0)
    central_longitude = decide_central_longitude(x)
    print("central longitude", central_longitude)
    # if central_longitude !=0:
    #   x = shift_xvals(x.values,central_longitude)
    #   vlist_copy[0] = shift_xvals(vlist[0],central_longitude)
    transform = get_transform(central_longitude)
    fig, axarr = plt.subplots(
        nrows=nrow,
        ncols=ncol,
        figsize=(10, 10),
        constrained_layout=False,
        subplot_kw={"projection": transform},
    )
    if nrow > 1:
        axlist = axarr.flatten()
    else:
        axlist = [axarr]
    iii = 0
    for ax in axlist:
        # if level doesn't exist then break.
        try:
            z = rtot.isel(z=iii)
        except IndexError:
            logger.warning("plotATL. no level with index {}".format(iii))
            break
        z = z.where(z != 0)
        label = meterev2FL(rtot.z.values[iii])
        cb = ATLsubplot(
            ax, x, y, z, transform, label, vlist_copy, levels=levels, nrow=nrow
        )
        iii += 1
    # if at least one subplot was created.
    if iii > 0:
        cb2 = plt.colorbar(cb)
        # cb2.set_label('mg/m$^3$')
        cb2.set_label("Percent", size=10)
        cbar_ax = fig.axes[-1]
        cbar_ax.tick_params(labelsize=10)
        fig.suptitle(title, size=10)
        print("NAME", name)
        if name != "None":
            plt.savefig(name)
            plt.close()
    else:
        logger.warning("plotATL. no plots created")
        print("no plots")


def meterev2FL(meters):
    return "FL{:2.0f}".format(meters / 30.48)


def setup_plot():
    mpl.rcParams.update(mpl.rcParamsDefault)
    mpl.rcParams["font.family"] = "sans-serif"
    # mpl.use('pdf')
    plt.style.use("seaborn-poster")
    # plt.style.use('fivethirtyeight')


def reset_plots():
    mpl.rcParams.update(mpl.rcParamsDefault)


def format_plot(ax, transform):
    setup_logger_warning(level=logging.WARNING)
    # ax.add_feature(cartopy.feature.LAND)
    #ax.add_feature(cartopy.feature.BORDERS)
    ax.coastlines("50m")
    # This allows latitude and longitude axis
    # to have different scales.
    ax.set_aspect("auto", adjustable=None)
    # this will increase data limit to keep aspect ratio 1
    # ax.set_aspect(1, adjustable='datalim')
    # this will adjust axxes to  keep aspect ratio 1
    # when this is used, text is often mis-placed.
    # ax.set_aspect(1, adjustable='box')

    # don't use a different central_longitude when
    # making the tick labels.
    #return 1
    gl = ax.gridlines(
        crs=transform,
        draw_labels=True,
        linewidth=1,
        color="gray",
        alpha=0.5,
        linestyle="--",
    )
    gl.top_labels = False
    gl.bottom_labels = True
    gl.right_labels = False
    gl.left_labels = True
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {"size": 10, "color": "gray"}
    gl.ylabel_style = {"size": 10, "color": "gray"}


def ATLsubplot(
    ax,
    x,
    y,
    z,
    transform,
    label="",
    vlist=None,
    levels=[1, 5, 10, 15, 20],
    name="None",
    nrow=6,
):
    """
    ax : axes
    x :
    y :
    z :
    transform :
    label : str
    vlist : [longiude, latitude]
    levels : list of floats/ints, levels for contouring.
    name : str save figure to this name.
    """
    xplace, yplace, size = set_ATL_text(nrow)
    cmap = plt.get_cmap("viridis")
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
    cb2 = ax.pcolormesh(x, y, z, cmap=cmap, transform=transform, norm=norm)
    # ax.pcolormesh(x, y, z, cmap=cmap, transform=transform, norm=norm)
    if isinstance(vlist,(np.ndarray,list,tuple)): 
        ax.plot(vlist[0], vlist[1], "r^")
    # except Exception as ex:
    #    print('exception {} {}'.format(type(ex).__name__, ex.args)
    format_plot(ax, transform)
    ax.text(
        xplace,
        yplace,
        label,
        va="bottom",
        ha="center",
        rotation="horizontal",
        rotation_mode="anchor",
        transform=ax.transAxes,
        size=size,
        backgroundcolor="white",
    )
    # ax.plot(vlist[0],vlist[1],'r^')
    if name != "None":
        plt.savefig(name)
        # plt.close()
    # plt.show()
    return cb2


def massload_ensemble_mean(revash, enslist):
    rev2 = revash.sel(ens=enslist)
    mass = hysplit.hysp_massload(rev2)
    massmean = mass.mean(dim="ens")
    return massmean

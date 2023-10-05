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

# 2023 Oct 3  added _bool_ method to LabelData class.
#             added get_empty_label function
#             added ax argument to hplot and sub_height_plot
#             added height plotting to the PlotVAA class       

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
from utilhysplit.evaluation.ensemble_tools import ATL, preprocess, topheight
from utilhysplit.plotutils.colormaker import ColorMaker
from utilhysplit.evaluation.ensemble_polygons import HeightPolygons

logger = logging.getLogger(__name__)


def setup_logger_warning(level=logging.WARNING):
    MESSAGE_FORMAT_INFO = "%(asctime)s.%(msecs)03d %(levelname)s - %(message)s"
    MESSAGE_FORMAT_DEBUG = (
        "%(asctime)s.%(msecs)03d %(levelname)s %(name)s - %(message)s"
    )
    import sys

    logging.basicConfig(
        stream=sys.stdout, level=level, format=MESSAGE_FORMAT_INFO, datefmt="%H:%M:%S"
    )

def get_empty_label():
    return LabelData('none',None,None,None,None)

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

    def __bool__(self):
        if not isinstance(self.units,str): return False
        if not isinstance(self.source,str): return False
        if not isinstance(self.time,str): return False
        if not isinstance(self.descrip,str): return False
        return True 


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
    # if min_x < 180 and max_x > 180:
    if min_x < -90 or max_x > 90:
        central_longitude = 180
    else:
        central_longitude = 0
    # return central_longitude
    return 0


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


def get_source_str(vlist):
    if vlist:
        source = "Longitude: {:0.2f} Latitude {:0.2f}".format(vlist[1], vlist[0])
    else:
        source = "Unknown"
    return source




def height_plot(revash, thresh=0.2, figname="None", vlist=None, unit="FL"):
    # TO DO add plotting of bottom height.
    setup_plot()
    source = get_source_str(vlist)
    figname = "None"
    if "source" in revash.dims:
        # plot for ensemble mean
        height_plot_time(
            revash.mean(dim="source"),
            thresh=thresh,
            figname=figname,
            vlist=vlist,
            label=str(source),
            unit=unit,
        )
        # plot for individual sources
        for source in revash.source.values:
            height_plot_time(
                revash.sel(source=source),
                thresh=thresh,
                figname=figname,
                vlist=vlist,
                label=str(source),
                unit=unit,
            )
        else:
            height_plot_time(
                revash,
                thresh=thresh,
                figname=figname,
                vlist=vlist,
                label="",
                unit=unit,
            )


def height_plot_time(
    revash, thresh=0.2, vlist=None, label="", figname="None", unit="FL"
):
    level = np.arange(0, len(revash.z.values))
    dlev = revash.z.values[-1] - revash.z.values[-2]  
 
    if 'time' in revash.dims:
        loopvals = revash.time.values
    else:
        loopvals = [999]
    for time in loopvals:
        print(time)
        top, bottom = topheight(revash, time, level=level, dlev=dlev, thresh=thresh)
        slabel = LabelData(
            time,
            "Top and bottom height of cloud using {} mg/m3 threshold".format(thresh),
            unit,
            label,
        )
        hplot(top,thresh=thresh,label=slabel,figname=figname,vlist=vlist,unit=unit)
    return -1

def hplot(
    top, thresh=0.2, vlist=None, label="", figname="None", unit="FL",ax=None
):
        x = top.longitude
        y = top.latitude
        z = top.values
        #central_longitude = decide_central_longitude(x)
        transform = get_transform(0)
        levels = set_height_levels(z)
        sub_height_plot(x, y, z, transform, levels, label, figname, vlist, unit=unit,ax=ax)


def set_height_levels(zvals):
    levs = list(set(list(zvals.flatten() - 1)))
    levs.sort()
    levs[0] = 0
    if len(levs) <= 1:
        levs.append(100)
    return levs


def massload_plot(revash, enslist=None, sourcelist=None, name="None", vlist=None):
    """
    revash: xarray data-array with dimensions source, ens, time, z, y, x.
            values are concentrations in mg/m3.
    vlist: tuple of [latitude, longitude] with volcano location
    enslist:
    """
    setup_plot()
    mass = massload_ensemble_mean(revash, enslist, sourcelist)
    fignamelist = []

    source = get_source_str(vlist)

    for iii, time in enumerate(mass.time.values):
        label = LabelData(
            time, "Column mass loading (ensemble mean)", "g/m$^2$", source
        )
        mass2 = mass.sel(time=time)
        x = mass2.longitude
        y = mass2.latitude
        # z = mass2.where(mass2!=0)
        central_longitude = decide_central_longitude(x)
        transform = get_transform(central_longitude)

        figname = name.replace("zzz", "{}".format(iii))
        z = mass2 / 1000.0  # convert to g/m2
        levels = set_levels(z)
        try:
            sub_massload_plot(x, y, z, transform, levels, label, figname, vlist)
        except ValueError as ex:
            print("ERROR: cannot create massload plot at {}: {}".format(time, str(ex)))
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
    data_transform = get_transform(central_longitude=0)
    cb2 = ax.contourf(x, y, z, levels=levels, colors=clrs, transform=data_transform)
    # cb2 = ax.pcolormesh(x,y,z,cmap=cmap,transform=transform,norm=norm)
    if isinstance(vlist, (np.ndarray, list, tuple)):
        ax.plot(vlist[1], vlist[0], "r^", transform=data_transform)
    format_plot(ax, data_transform)
    label_ax(dax, labeldata, data_transform)
    cb = plt.colorbar(cb2)
    cb.set_label("g/m$^2$")
    if name != "None":
        plt.savefig(name)
    plt.show()
    # plt.close()



class PlotVAA:

    def __init__(self):

        self._transform = get_transform()
        self.setup()
        self._timestep = 3

    def setup(self):
        nrow = 3
        ncol = 2 
        fig, axra = plt.subplots(
            nrows=nrow,
            ncols=ncol,
            figsize=(20, 20),
            constrained_layout=False,
            subplot_kw={"projection": self.transform},
        )
        self.fig = fig
        self.axra = axra.flatten() 

    @property
    def transform(self):
        return self._transform

    @property
    def timestep(self):
        return datetime.timedelta(hours=self._timestep)

    @property
    def model(self): 
        return self._model 

    @model.setter
    def model(self,dset):
        self._model = dset

    def plot_one(self,vloc,thresh,plotmass=True):
        timelist = self._model.time.values
        zlevs = np.arange(0,len(self._model.z.values))
        hhh = []
        lab = []
        for iii, time in enumerate(timelist):
            ax = self.axra[iii]
            size=15
            dset = self._model.sel(time=time) 
            top,bottom = topheight(dset,time=time,level=zlevs,thresh=thresh)
          
            cmap='Blues'
            top_poly = HeightPolygons(cmap=cmap)
            top_poly.process(top,alpha=0.1)
            lw = 5
            top_poly = top_poly.merge(key='high')

            bottom_poly = HeightPolygons(cmap=cmap) 
            bottom_poly.process(bottom,alpha=0.1)
            bottom_poly = bottom_poly.merge(key='low')
            tpoly = top_poly.merge_with(bottom_poly)
            handles,labels = tpoly.plot(ax=self.axra[iii],vloc=vloc,pbuffer=0.15,legend=False,linewidth=lw)
            format_plot(self.axra[iii], self.transform)
            handles,labels = self.sort_labels(handles,labels)
            self.axra[iii].legend(handles,labels,fontsize=20)
            if plotmass:
                self.plotmass_method(self.axra[iii],dset)
            time = pd.to_datetime(time) 
            time_label=time.strftime("%d %b %H:%M UTC")
            yplace = 1.
            xplace = 0.2
            size = 20
            self.axra[iii].text(
                xplace,
                yplace,
                time_label,
                va="bottom",
                ha="center",
                rotation="horizontal",
                rotation_mode="anchor",
                transform=ax.transAxes,
                size=size,
                backgroundcolor="white")
        plt.tight_layout()
        self.polygons=top_poly

    def plotmass_method(self,ax,dset):
        levels = [0.02,0.1,0.2,2,5,10,50]
        cm = 'Reds'
        mass_cmap = plt.get_cmap(cm)
        norm = BoundaryNorm(levels,ncolors=mass_cmap.N,clip=False,extend='both')
        cset = xr.where(dset==0,np.nan,dset)
        cset = cset.max(dim='z')
        x = cset.longitude.values
        y = cset.latitude.values
        z = cset.values
        cb2 = ax.pcolormesh(x,y,z,transform=self.transform,norm=norm,cmap=cm)
        #cset.plot.pcolormesh(ax=ax,x='longitude',y='latitude',transform=self.transform,cmap='Reds')     
        cb = self.fig.colorbar(cb2)
        cb.set_label('mg m$^{-3}$',fontsize=20)               
        cb.ax.tick_params(labelsize=20)
        latr, lonr = self.find_limit(cset)
        ax.set_xlim(lonr[0],lonr[1])
        ax.set_ylim(latr[0],latr[1])

    def plotheight_method(self,ax,top):
        x = top.longitude
        y = top.latitude
        z = top.values
        levels = set_height_levels(z)
        cmap = plt.get_cmap("binary")
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
        cb2 = ax.pcolormesh(x, y, z, cmap=cmap, transform=self.transform, norm=norm)
        cb = plt.colorbar(cb2, ticks=levels)
        tick_labels = [meterev2FL(x) for x in levels]
        cb.ax.set_yticklabels(tick_labels)
        cb.set_label("Height")
        latr, lonr = self.find_limit(top)
        ax.set_xlim(lonr[0],lonr[1])
        ax.set_ylim(latr[0],latr[1])


    def plot(self,vloc,thresh,plotmass=True,plotheight=False,ppp='top'):
        timelist = self._model.time.values
        zlevs = np.arange(0,len(self._model.z.values))
        hhh = []
        lab = []
        for iii, time in enumerate(timelist):
            ax = self.axra[iii]
            xplace, yplace, size = set_ATL_text(iii)
            size=15
            dset = self._model.sel(time=time) 
            top,bottom = topheight(dset,time=time,level=zlevs,thresh=thresh)
            if ppp=='top': cmap='cividis'
            else: cmap='Blues'
            top_poly = HeightPolygons(cmap=cmap)
            self.polygons=top_poly
            top_poly.process(top,alpha=0.1)
            lw = 5
            #top_poly = top_poly.merge(key='high')

            bottom_poly = HeightPolygons(cmap=cmap) 
            bottom_poly.process(bottom,alpha=0.1)
            #bottom_poly = bottom_poly.merge(key='low')
          
            if ppp=='top': tpoly = top_poly
            if ppp=='bottom': tpoly = bottom_poly


            handles,labels = tpoly.plot(ax=self.axra[iii],vloc=vloc,pbuffer=0.15,legend=False,linewidth=lw)
            format_plot(self.axra[iii], self.transform)
            hhh.extend(handles)
            lab.extend(labels) 
            time = pd.to_datetime(time) 
            time_label=time.strftime("%d %b %H:%M UTC")
            if not plotmass:
              self.axra[iii].text(
                xplace,
                yplace,
                time_label,
                va="bottom",
                ha="center",
                rotation="horizontal",
                rotation_mode="anchor",
                transform=ax.transAxes,
                size=size,
                backgroundcolor="white",
              )
            if plotmass:
                self.plotmass_method(ax,dset)
            elif plotheight:
                self.plotheight_method(ax,top)
        handles,labels = self.sort_labels(hhh,lab)
        self.axra[0].legend(handles,labels,fontsize=20)
        plt.tight_layout()
   
    def find_limit(self,cset,buf=1,cmin=0.0001):
        lat = cset.latitude.values
        lon = cset.longitude.values
        vals = cset.values
        a = zip(lat.flatten(),lon.flatten(),vals.flatten())
        b = [x for x in a if x[2]>cmin]
        lat,lon,val = zip(*b)
        minlat = np.nanmin(lat)
        maxlat = np.nanmax(lat)
        minlon = np.nanmin(lon)
        maxlon = np.nanmax(lon)
        return (minlat-buf,maxlat+buf), (minlon-buf,maxlon+buf) 

    def sort_labels(self,handles,labels):
        lset = list(set(labels))
        lset.sort()
        hset = []
        zzz = list(zip(labels,handles))
        for val in lset:
            for zval in zzz:
                if zval[0]==val:
                   hset.append(zval[1])
                   break
        return hset, lset


def sub_height_plot(
    x, y, z, transform, levels, labeldata, name="None", vlist=None, unit="FL",ax=None
):
    if not ax:
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
    cmap = plt.get_cmap("binary")
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
    data_transform = get_transform(central_longitude=0)
    #data_transform = transform
    # cb2 = ax.pcolormesh(x, y, z, levels=levels, cmap=cmap, transform=data_transform)
    cb2 = ax.pcolormesh(x, y, z, cmap=cmap, transform=transform, norm=norm)
    if isinstance(vlist, (np.ndarray, list, tuple)):
        ax.plot(vlist[1], vlist[0], "r^", transform=data_transform)

    format_plot(ax, data_transform)
    cb = plt.colorbar(cb2, ticks=levels)
    if labeldata:
        label_ax(dax, labeldata, data_transform)
    if unit.upper() == "FL":
        tick_labels = [meterev2FL(x) for x in levels]
        cb.ax.set_yticklabels(tick_labels)
        cb.set_label("Height")
    elif unit.lower() == "m":
        tick_labels = levels
        cb.ax.set_yticklabels(tick_labels)
        cb.set_label("Height (m)")
    elif unit.lower() == "km":
        tick_labels = [x / 1000.0 for x in levels]
        cb.ax.set_yticklabels(tick_labels)
        cb.set_label("Height (km)")

    if name != "None" and isinstance(name,str):
        plt.savefig(name)
    plt.show()
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
        source = ""
        if isinstance(vlist, (np.ndarray, list, tuple)):
            source = "Longitude: {:0.2f} Latitude {:0.2f}".format(vlist[0], vlist[1])
        label = LabelData(time, "Probability of exceedance", "%", source)
        rtot = ATL(revash2, enslist=enslist, thresh=thresh, norm=norm)
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

    if isinstance(vlist, (np.ndarray, list, tuple)):
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
    z = temp.where(temp != 0)
    central_longitude = decide_central_longitude(x)
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
            # plt.close()
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
    # ax.add_feature(cartopy.feature.BORDERS)
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
    # return 1
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
    gl.xlabel_style = {"size": 20, "color": "gray"}
    gl.ylabel_style = {"size": 20, "color": "gray"}


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
    if isinstance(vlist, (np.ndarray, list, tuple)):
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


def massload_ensemble_mean(revash, enslist=None, sourcelist=None):
    rev2, dim = preprocess(revash, enslist, sourcelist)
    mass = hysplit.hysp_massload(rev2)
    massmean = mass.mean(dim=dim)
    return massmean

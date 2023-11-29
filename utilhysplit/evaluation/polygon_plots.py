""" creates ensemble relative frequency of exceedance and ensemble mean plots

Classes
    PlotVAA

"""

# 2023 Nov 20 copied from web_ensemble_plots

import datetime
import logging

#import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import numpy as np
import pandas as pd
import xarray as xr

from utilhysplit.evaluation import web_ensemble_plots as wep
from utilhysplit.evaluation.ensemble_polygons import HeightPolygons
from utilhysplit.evaluation.ensemble_tools import topheight

logger = logging.getLogger(__name__)


class PlotVAA:
    def __init__(self):
        self._transform = wep.get_transform()
        self.setup()
        self._timestep = 3
        self._model = xr.DataArray()
        #self.polygons = HeightPolygons(cmap="viridis")
        self.polygons = None

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
    def model(self, dset):
        self._model = dset

    def plot_one(self, vloc, thresh, plotmass=True):
        timelist = self._model.time.values
        zlevs = np.arange(0, len(self._model.z.values))
        for iii, time in enumerate(timelist):
            ax = self.axra[iii]
            size = 15
            dset = self._model.sel(time=time)
            top, bottom = topheight(dset, time=time, level=zlevs, thresh=thresh)

            cmap = "Blues"
            top_poly = HeightPolygons(cmap=cmap)
            top_poly.process(top, alpha=0.1)
            lw = 5
            top_poly = top_poly.merge(key="high")

            bottom_poly = HeightPolygons(cmap=cmap)
            bottom_poly.process(bottom, alpha=0.1)
            bottom_poly = bottom_poly.merge(key="low")
            tpoly = top_poly.merge_with(bottom_poly)
            handles, labels = tpoly.plot(
                ax=self.axra[iii], vloc=vloc, pbuffer=0.15, legend=False, linewidth=lw
            )
            format_plot(self.axra[iii], self.transform)
            handles, labels = sort_labels(handles, labels)
            self.axra[iii].legend(handles, labels, fontsize=20)
            if plotmass:
                self.plotmass_method(self.axra[iii], dset)
            time = pd.to_datetime(time)
            time_label = time.strftime("%d %b %H:%M UTC")
            yplace = 1.0
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
                backgroundcolor="white",
            )
        plt.tight_layout()
        self.polygons = top_poly

    def plotmass_method(self, ax, dset):
        levels = [0.02, 0.1, 0.2, 2, 5, 10, 50]
        cmap = "Reds"
        mass_cmap = plt.get_cmap(cmap)
        norm = BoundaryNorm(levels, ncolors=mass_cmap.N, clip=False, extend="both")
        cset = xr.where(dset == 0, np.nan, dset)
        cset = cset.max(dim="z")
        x = cset.longitude.values
        y = cset.latitude.values
        z = cset.values
        cb2 = ax.pcolormesh(x, y, z, transform=self.transform, norm=norm, cmap=cm)
        # cset.plot.pcolormesh(ax=ax,x='longitude',y='latitude',transform=self.transform,cmap='Reds')
        cb = self.fig.colorbar(cb2)
        cb.set_label("mg m$^{-3}$", fontsize=20)
        cb.ax.tick_params(labelsize=20)
        latr, lonr = find_limit(cset)
        ax.set_xlim(lonr[0], lonr[1])
        ax.set_ylim(latr[0], latr[1])

    def plotheight_method(self, ax, top):
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
        latr, lonr = find_limit(top)
        ax.set_xlim(lonr[0], lonr[1])
        ax.set_ylim(latr[0], latr[1])

    def plot(self, vloc, thresh, plotmass=True, plotheight=False, ppp="top"):
        timelist = self._model.time.values
        zlevs = np.arange(0, len(self._model.z.values))
        hhh = []
        lab = []
        for iii, time in enumerate(timelist):
            ax = self.axra[iii]
            xplace, yplace, size = wep.set_ATL_text(iii)
            size = 15
            dset = self._model.sel(time=time)
            top, bottom = topheight(dset, time=time, level=zlevs, thresh=thresh)
            if ppp == "top":
                cmap = "cividis"
            else:
                cmap = "Blues"
            top_poly = HeightPolygons(cmap=cmap)
            self.polygons = top_poly
            top_poly.process(top, alpha=0.1)
            lw = 5
            # top_poly = top_poly.merge(key='high')

            bottom_poly = HeightPolygons(cmap=cmap)
            bottom_poly.process(bottom, alpha=0.1)
            # bottom_poly = bottom_poly.merge(key='low')

            if ppp == "top":
                tpoly = top_poly
            if ppp == "bottom":
                tpoly = bottom_poly

            handles, labels = tpoly.plot(
                ax=self.axra[iii], vloc=vloc, pbuffer=0.15, legend=False, linewidth=lw
            )
            wep.format_plot(self.axra[iii], self.transform)
            hhh.extend(handles)
            lab.extend(labels)
            time = pd.to_datetime(time)
            time_label = time.strftime("%d %b %H:%M UTC")
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
                self.plotmass_method(ax, dset)
            elif plotheight:
                self.plotheight_method(ax, top)
        handles, labels = sort_labels(hhh, lab)
        self.axra[0].legend(handles, labels, fontsize=20)
        plt.tight_layout()

def find_limit(cset, buf=1, cmin=0.0001):
    lat = cset.latitude.values
    lon = cset.longitude.values
    vals = cset.values
    a = zip(lat.flatten(), lon.flatten(), vals.flatten())
    b = [x for x in a if x[2] > cmin]
    lat, lon, _ = zip(*b)
    minlat = np.nanmin(lat)
    maxlat = np.nanmax(lat)
    minlon = np.nanmin(lon)
    maxlon = np.nanmax(lon)
    return (minlat - buf, maxlat + buf), (minlon - buf, maxlon + buf)

def sort_labels(handles, labels):
    lset = list(set(labels))
    lset.sort()
    hset = []
    zzz = list(zip(labels, handles))
    for val in lset:
        for zval in zzz:
            if zval[0] == val:
                hset.append(zval[1])
                break
    return hset, lset

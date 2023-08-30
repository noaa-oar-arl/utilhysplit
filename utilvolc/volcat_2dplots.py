# volcat_2dplots.py
# plot volcat data
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd

from utilvolc import volcat
from utilhysplit.plotutils import vtools

# from pyresample.bucket import BucketResampler


"""
This script contains routines that plot VOLCAT data in xarray format,
-------------
Functions:
-------------
create_pc_plot: plots parallax corrected vs uncorrected values
compare_pc: compares corrected vs uncorrected values
plot_height: plots ash top height from VOLCAT
plot_radius: plots ash effective radius from VOLCAT
plot_mass: plots ash mass loading from VOLCAT
plot_gen: generates quick plot, not saved
------------
"""

# 2023 19 July (amc) modified plot_gen to return axes and use the vtools formatting.


def create_pc_plot(dset):
    """
    creates plots of parallax corrected vs. uncorrected values.
    """

    def subfunc(ax, vals):
        ax.plot(vals[0], vals[1], "k.", markersize=1)
        # plot 1:1 line
        minval = np.min(vals[0])
        maxval = np.max(vals[0])
        ax.plot([minval, maxval], [minval, maxval], "--r.", markersize=1)

    latitude, longitude = compare_pc(dset)
    fig = plt.figure(1)
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)

    ax1.set_ylabel("uncorrected")
    ax2.set_ylabel("uncorrected")
    ax2.set_xlabel("corrected")

    subfunc(ax1, latitude)
    subfunc(ax2, longitude)
    return fig, ax1, ax2


def compare_pc(dset):
    """
    Returns:
    latitude : [list of parrallax corrected values, list of uncorrected values]
    longitude : [list of parrallax corrected values, list of uncorrected values]
    """

    def process(pc, val):
        # pair corrected and uncorrected values.
        pzip = list(zip(pc, val))
        # remove nans
        new = [x for x in pzip if not np.isnan(x[0])]
        return list(zip(*new))

    pc_lat = volcat.get_pc_latitude(dset)
    pc_lon = volcat.get_pc_longitude(dset)
    latvals = pc_lat.latitude.values.flatten()
    lonvals = pc_lon.longitude.values.flatten()
    pclat = pc_lat.values.flatten()
    pclon = pc_lon.values.flatten()

    latitude = process(pclat, latvals)
    longitude = process(pclon, lonvals)
    return latitude, longitude


def plot_height(dset):
    """Plots ash top height from VOLCAT
    Does not save figure - quick image creation"""
    fig = plt.figure("Ash_Top_Height")
    title = "Ash Top Height (km)"
    # ax = fig.add_subplot(1, 1, 1)
    ax = plot_gen(dset, val="height", time=None, plotmap=True, title=title)
    return ax


def plot_radius(dset):
    """Plots ash effective radius from VOLCAT
    Does not save figure - quick image creation"""
    # fig = plt.figure("Ash_Effective_Radius")
    title = "Ash effective radius ($\mu$m)"
    # ax = fig.add_subplot(1, 1, 1)
    ax = plot_gen(dset, val="radius", time=None, plotmap=True, title=title)
    return ax


def plot_mass(dset, central_longitude=0):
    # fig = plt.figure("Ash_Mass_Loading")
    # ax = fig.add_subplot(1, 1, 1)
    ax = plot_gen(
        dset,
        val="mass",
        time=None,
        plotmap=True,
        title="Ash_Mass_Loading",
        unit="g m$^{-2}$",
        central_longitude=central_longitude,
    )
    return ax


def plot_gen(
    dset,
    val="mass",
    time=None,
    plotmap=True,
    title=None,
    central_longitude=180,
    unit=None,
):
    """Plot ash mass loading from VOLCAT
    Does not save figure - quick image creation"""
    # lat=dset.latitude
    # lon=dset.longitude
    if val == "mass":
        mass = volcat.get_mass(dset)
    elif val == "radius":
        mass = volcat.get_radius(dset)
    elif val == "height":
        mass = volcat.get_height(dset)
    if time and "time" in mass.coords:
        mass = mass.sel(time=time)
    elif "time" in mass.coords:
        mass = mass.isel(time=0)
    lat = mass.latitude
    lon = mass.longitude
    if plotmap:
        transform = ccrs.PlateCarree(central_longitude=central_longitude)
        vtransform = ccrs.PlateCarree(central_longitude=0)
        fig, axx = plt.subplots(1, 1, subplot_kw={"projection": transform})
        # m.add_feature(cfeat.LAND)
        # m.add_feature(cfeat.COASTLINE)
        # m.add_feature(cfeat.BORDERS)
        cb = axx.pcolormesh(lon, lat, mass, transform=vtransform)
        vtools.format_plot(axx, transform)
    else:
        plt.pcolormesh(lon, lat, mass)
    cb2 = plt.colorbar(cb)
    if isinstance(unit, str):
        cb2.set_label(unit)
    dstr = str(mass.time.values)
    plt.title(title + " " + dstr[0:-8])
    return axx
    # plt.show()

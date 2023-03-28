# volcat_2dplots.py
# plot volcat data
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd

import volcat

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


def create_pc_plot(dset):
    """
    creates plots of parallax corrected vs. uncorrected values.
    """

    def subfunc(ax, vals):
        ax.plot(vals[0], vals[1], "k.", MarkerSize=1)
        # plot 1:1 line
        minval = np.min(vals[0])
        maxval = np.max(vals[0])
        ax.plot([minval, maxval], [minval, maxval], "--r.", MarkerSize=1)

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
    ax = fig.add_subplot(1, 1, 1)
    plot_gen(dset, ax, val="height", time=None, plotmap=True, title=title)


def plot_radius(dset):
    """Plots ash effective radius from VOLCAT
    Does not save figure - quick image creation"""
    fig = plt.figure("Ash_Effective_Radius")
    title = "Ash effective radius ($\mu$m)"
    ax = fig.add_subplot(1, 1, 1)
    plot_gen(dset, ax, val="radius", time=None, plotmap=True, title=title)


def plot_mass(dset):
    fig = plt.figure("Ash_Mass_Loading")
    ax = fig.add_subplot(1, 1, 1)
    plot_gen(dset, ax, val="mass", time=None, plotmap=True, title="Ash_Mass_Loading")


def plot_gen(dset, ax, val="mass", time=None, plotmap=True, title=None):
    """Plot ash mass loading from VOLCAT
    Does not save figure - quick image creation"""
    # lat=dset.latitude
    # lon=dset.longitude
    if val == "mass":
        mass = get_mass(dset)
    elif val == "radius":
        mass = get_radius(dset)
    elif val == "height":
        mass = get_height(dset)
    if time and "time" in mass.coords:
        mass = mass.sel(time=time)
    elif "time" in mass.coords:
        mass = mass.isel(time=0)
    lat = mass.latitude
    lon = mass.longitude
    if plotmap:
        m = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
        m.add_feature(cfeat.LAND)
        m.add_feature(cfeat.COASTLINE)
        m.add_feature(cfeat.BORDERS)
        plt.pcolormesh(lon, lat, mass, transform=ccrs.PlateCarree())
    else:
        plt.pcolormesh(lon, lat, mass)
    plt.colorbar()
    plt.title(title)
    plt.show()

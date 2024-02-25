import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
import matplotlib.ticker as mticker
import numpy as np

def set_ticks(xmin,xmax,nticks):
    if xmin<0 and xmax>0:
       xticks = set_ticks_dateline(xmax,xmin,nticks)
    elif xmin>0 and xmax<0:
       xticks = set_ticks_dateline(xmin,xmax,nticks)
    elif xmax < xmin:
       xticks = set_ticks_normal(xmax,xmin,nticks)
    else:
       xticks = set_ticks_normal(xmin,xmax,nticks)
    return xticks

def set_ticks_dateline(xmin,xmax,nticks=3):
    # transform to 0 to 360 coordinates
    xlimits = [x+360 if x<0 else x for x in [xmin,xmax]]
    xticks = set_ticks_normal(xlimits[0],xlimits[1])
    # transform back to -180 to 180 coordinates
    xticks = [x-360 if x>180 else x for x in xticks]
    return xticks


def set_ticks_normal(xmin,xmax,nticks=3):
    # if the span is less than one degree
    if np.abs(xmax-xmin) > 1:
       first = np.floor(xmin)
       last = np.ceil(xmax)
    # round to nearest tenth of degree
    else:
       first = np.floor(xmin*10)/10
       last = np.ceil(xmax*10)/10
    span = last-first
    dx = span/nticks
    print('DX', dx)
    # round to nearest 5 if greater than 10
    if np.abs(dx) > 10:
       print('here')
       dx = int(int(dx/5.0)*5)
    elif np.abs(dx) > 1:
       dx = int(np.floor(dx))
    else:
       dx = int(dx*10)/10.0
    print('DX', dx)
    rticks = np.arange(first,last+dx,dx)
    rticks = [int(x*100)/100 for x in rticks]
    return list(rticks)


def draw_map(fignum, ax = None, fs=20):
    proj = ccrs.PlateCarree(central_longitude=180)
    if not ax:
        fig = plt.figure(fignum, figsize=(16,9))
        ax = fig.add_subplot(1,1,1)
    ax = plt.axes(projection=proj)
    gl = ax.gridlines(draw_labels=True, linewidth = 0.2, color="gray")
    gl.labels_right=False
    gl.labels_top=False
    #ocean = cfeature.NaturalEarthFeature(category='physical',name='ocean',
    #        scale='50m', facecolor = cfeature.COLORE['water'])

    #ax.add_feature(states.edgecolor='gray')
    ax.add_feature(cfeature.LAND)
    return ax

def reset_plots():
    mpl.rcParams.update(mpl.rcParamsDefault)

def setup_plot():
    mpl.rcParams.update(mpl.rcParamsDefault)
    mpl.rcParams["font.family"] = "sans-serif"
    # mpl.use('pdf')
    # plt.style.use("seaborn-poster")
    # plt.style.use('fivethirtyeight')


def format_plot(ax, transform, 
                xticks=None,
                yticks=None, 
                ylabel=True,
                xlabel=True,
                fsz = 20,
                land=False, borders=False, coastlines=True):

    #setup_logger_warning(level=logging.WARNING)
    if land: ax.add_feature(cartopy.feature.LAND)
    if borders: ax.add_feature(cartopy.feature.BORDERS)
    if coastlines: ax.coastlines("50m")
    # This allows latitude and longitude axis
    # to have different scales.
    ax.set_aspect("auto", adjustable=None)
    # this will increase data limit to keep aspect ratio 1
    # ax.set_aspect(1, adjustable='datalim')
    # this will adjust axxes to  keep aspect ratio 1
    # when this is used, text is often mis-placed.
    # ax.set_aspect(1, adjustable='box')

    # in general the transform should have central_longitude=0
    # otherwise it will label 180 degrees as 0.

    gl = ax.gridlines(
        crs=transform,
        draw_labels=True,
        linewidth=1,
        color="gray",
        alpha=0.5,
        linestyle="--",
    )
    gl.top_labels = False
    gl.bottom_labels = xlabel
    gl.right_labels = False
    gl.left_labels = ylabel
    if isinstance(xticks,list):
       gl.xlocator = mticker.FixedLocator(xticks)
    if isinstance(yticks,list):
       gl.ylocator = mticker.FixedLocator(yticks)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {"size": fsz, "color": "gray"}
    gl.ylabel_style = {"size": fsz, "color": "gray"}

def get_transform(central_longitude=0):
    transform = cartopy.crs.PlateCarree(central_longitude=central_longitude, globe=None)
    # transform = cartopy.crs.AzimuthalEquidistant(central_longitude=180)
    return transform


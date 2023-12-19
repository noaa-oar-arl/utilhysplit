import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
import matplotlib.ticker as mticker


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
    plt.style.use("seaborn-poster")
    # plt.style.use('fivethirtyeight')


def format_plot(ax, transform, 
                xticks=None,
                yticks=None, 
                land=False, borders=False, coastlines=True):

    setup_logger_warning(level=logging.WARNING)
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
    gl.bottom_labels = True
    gl.right_labels = False
    gl.left_labels = True
    if isinstance(xticks,list):
       gl.xlocator = mticker.FixedLocator(xticks)
    if isinstance(yticks,list):
       gl.ylocator = mticker.FixedLocator(yticks)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {"size": 20, "color": "gray"}

def get_transform(central_longitude=0):
    transform = cartopy.crs.PlateCarree(central_longitude=central_longitude, globe=None)
    # transform = cartopy.crs.AzimuthalEquidistant(central_longitude=180)
    return transform


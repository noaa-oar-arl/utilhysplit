# plot_CrIS_so2.py
# Reads and plots CrIS SO2
# From Dave Hyman (dhyman2@wisc.edu)
# Modified by Allison Ring (allison.ring@noaa.gov)

import volcat_so2 as so2
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.colors as c
#import cartopy.crs as ccrs
#import cartopy.feature as cfeat
from scipy import integrate
import numpy as np
from glob import glob

# NEED THIS FUNCTION TO INTERPOLATE
# THE MEDIAN HEIGHT ON EACH CELL


def dist_percentile(X, CDF_X, percentile, axis=-1):
    F = percentile / 100.
    over = CDF_X > F
    x_idx_over = np.argmax(over, axis=axis)
    x_idx_under = np.argmax(over, axis=axis) - 1
    X1 = X[x_idx_over]
    X0 = X[x_idx_under]
    Y1 = 1. - np.max((1.-CDF_X)*over, axis=axis)
    Y0 = np.max(CDF_X*(1-over), axis=axis)
    dY = Y1-Y0
    dX = X1-X0
    return X0 + dX/dY * (F - Y0)


# SPECIFY GRID FILE TO USE
flist = glob('/pub/ECMWF/JPSS/VOLCAT/Raikoke/SO2/*.nc')
for filename in flist:
    dset = so2.open_dataset(filename)
    date = filename[-30:-3]
    # EXTRACT NETCDF VARIABLES
    #time_utc = file_dset['time_utc'][:]
    lat = dset.latitude.values
    lon = dset.longitude.values
    zscore_initial = dset.zscore_initial
    #height = file_dset['height'][:]
    #height_pdf = file_dset['height_pdf'][:]
    ml_mean = dset.mass_loading_mean
    #mass_loading_std = file_dset['mass_loading_std'][:]

    # CALCULATE AUXILLIARY VARIABLES
    # HEIGHT CDF
    # height_cdf = np.concatenate(
    #    (0.*height_pdf[:, :, 0:1], integrate.cumtrapz(height_pdf, height, axis=-1)), axis=-1)
    # POSITIVE DETECTIONS (TRUE / FALSE)
    detected = (zscore_initial >= 5.)
    # ORBITAL COVERAGE (TRUE / FALSE)
    orbital_coverage = np.isfinite(zscore_initial)
    # WHERE IS THE CLOUD (TRUE / FALSE)
    cloud = np.logical_and(detected, orbital_coverage)
    # MEDIAN HEIGHT
    #median_height = (1. * detected) * dist_percentile(height, height_cdf, 50., axis=- 1)
    # MEAN HEIGHT
    #mean_height = np.trapz(height_pdf * height, height, axis=-1)
    # MEAN MASS LOADING
    mass = ml_mean[:, :, -1]
    mass_mask = np.ma.masked_where(~cloud, ml_mean[:, :, -1])

    # SET UP BASEMAP PROJECTION
    proj = 'laea'  # LAMBERT EQUAL AREA
    wide = 13E6  # meters wide
    high = 12E6  # meters high
    # center_lat = 90  # projection center latitude
    center_lat = 65
    center_lon = 180.0
    # center_lon = 153.2540  # projection center longitude (Raikoke)
    mmap = Basemap(resolution='l', projection=proj, width=wide,
                   height=high, lat_0=center_lat, lon_0=center_lon)

    # TRANSFORM GRID ACCORDING TO PROJECTION
    X_grid, Y_grid = mmap(lon, lat)
    # MEAN GRID RESOLUTION
    dX = np.diff(X_grid, axis=-1).mean()
    dY = np.diff(Y_grid, axis=0).mean()

    # PLOTTING
    fig = plt.figure('SO2', figsize=(20, 10))
    #img_proj = ccrs.PlateCarree(central_longitude=180)
    #m = plt.axes(projection=img_proj)
    # m.add_feature(cfeat.LAND)
    # m.add_feature(cfeat.BORDERS)
    # m.coastlines('50m')
    # PLOT COASTLINES, LINES OF LONGITUDE,LATITUDE
    coastid = mmap.drawcoastlines()
    #parallels = np.arange(-90., 91, 5)
    #meridians = np.arange(-180, 181, 5)
    #parid = m.drawparallels(parallels, color='lightgrey', labels=[1, 1, 0, 0])
    #merid = m.drawmeridians(meridians, labels=[0, 0, 1, 1], fmt='%3d', color='lightgrey')

    # PLOT BACKGROUND
    # (IE: WHERE THERE IS ORBITAL COVERAGE IN THE INTERVAL)
    #bgid = plt.pcolormesh(X_grid-dX/2., Y_grid-dY/2., np.ma.masked_where(~orbital_coverage,0*X_grid), vmin=-2, vmax=1, cmap='gray')

    # SEVERAL DATA PLOTTING OPTIONS
    # ORBITAL COVERAGE (TRUE / FALSE)
    #pltid = plt.pcolormesh(X_grid-dX/2.,Y_grid-dY/2.,orbital_coverage,vmin = 0, vmax = 1, norm = None)
    # TOTAL COLUMN LOADING (DU)
    #plt.pcolormesh(lat, lon, mass, vmin=0.1, vmax=100., norm=c.LogNorm())
    plt.pcolormesh(X_grid-dX/2., Y_grid-dY/2., np.ma.masked_where(~cloud,
                                                                  ml_mean[:, :, 30]), vmin=0.1, vmax=100, norm=c.LogNorm())
    # INITIAL Z - SCORE
    #pltid = plt.pcolormesh(X_grid-dX/2.,Y_grid-dY/2.,zscore_initial,vmin = -5, vmax = 5, norm = None)
    # MEAN HEIGHT (km)
    #pltid = m.pcolormesh(X_grid-dX/2.,Y_grid-dY/2.,np.ma.masked_where(~cloud, mean_height), vmin = 0 ,vmax=20, norm = None)
    # MEDIAN HEIGHT (km)
    #pltid = m.pcolormesh(X_grid-dX/2., Y_grid-dY/2., np.ma.masked_where(~cloud,median_height), vmin=0, vmax=20, norm=None)
    cb = plt.colorbar()
    cb.set_label('DU', fontsize=12, rotation='270', labelpad=25)
    cb.ax.tick_params(labelsize=12)
    plt.scatter(-26.75, 48.292, s=100, c='black', marker='^')
    plt.title('Raikoke SO2 '+date, fontsize=16)
    plt.savefig('SO2_'+date+'.png')
    # plt.show()
    plt.close()

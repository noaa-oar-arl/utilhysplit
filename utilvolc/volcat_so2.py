# volcat_so2.py
# Reads and plots CrIS SO2
# From Dave Hyman (dhyman2@wisc.edu)
# Modified by Allison Ring (allison.ring@noaa.gov)

import numpy as np
import xarray as xr
from scipy import integrate

# Opens the dataset


def open_dataset(fname):
    """ Opens CrIS SO2 netcdf, sets height, lat, lon as coordinates
    Input: filename string
    Output: xarray dataset """
    print(fname)
    dset = xr.open_dataset(fname, mask_and_scale=False, decode_times=False)
    dset = _get_latlon(dset)
    dset = dset.set_coords(['height'])
    dset = dset.set_coords(['time_utc'])
    dset = dset.rename({"cols": 'x', "rows": 'y', "levels": 'z', "utc_tuple": 'time'})
    return dset


def _get_latlon(dset):
    """ Renames latitude and longitude, sets them as coordinates
         Input: xarray of data """
    dset = dset.rename({'lat': 'latitude'})
    dset = dset.rename({'lon': 'longitude'})
    dset = dset.set_coords(['latitude', 'longitude'])
    return dset

# Extracting various variables


def get_mass(dset):
    """ Extracts the mean mass loading of SO2
         Input: xarray of data
         Output: mass loading mean in g/m^2 (not DU) """
    mass_loading_mean = dset.mass_loading_mean
    # Return units of g/m^2 not DU
    conv = (2.6867E20) * (64.066) * (1/(6.022E23))
    mass_loading_mean = mass_load_mean * conv
    return mass_loading_mean


def get_meanheight(dset):
    """ Extracts the height variable
         Input: xarray DataSet
         Output: xarray DataArray of mean height in km """
    height = dset.height
    height_pdf = dset.height_pdf
    mean_height = np.trapz(height_pdf * height, height, axis=-1)
    mean_height = xr.DataArray(mean_height)
    mean_height = mean_height.rename({'dim_0': 'y', 'dim_1': 'x'})
    mean_height = mean_height.set_coords(['latitude', 'longitude'])

    zscore_initial = dset['zscore_initial'][:]

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
    dY = Y1 - Y0
    dX = X1 - X0
    tmp = X0 + dX/dY * (F - Y0)
    return tmp


# CALCULATE AUXILLIARY VARIABLES
# HEIGHT CDF
#height_cdf = np.concatenate((0.*height_pdf[:, :, 0:1], integrate.cumtrapz(height_pdf, height, axis=-1)), axis=-1)
# POSITIVE DETECTIONS (TRUE / FALSE)
#detected = (zscore_initial >= 5.)
# ORBITAL COVERAGE (TRUE / FALSE)
#orbital_coverage = np.isfinite(zscore_initial)
# WHERE IS THE CLOUD (TRUE / FALSE)
#cloud = np.logical_and(detected, orbital_coverage)
# MEDIAN HEIGHT
#median_height = (1. * detected) * dist_percentile(height, height_cdf, 50., axis=- 1)
# MEAN HEIGHT


# SET UP BASEMAP PROJECTION
# proj = 'laea'  # LAMBERT EQUAL AREA
# wide = 13E6  # meters wide
# high = 13E6  # meters high
# center_lat = 90  # projection center latitude
# center_lon = 153.2540  # projection center longitude (Raikoke)
#m = Basemap(resolution='l', projection=proj, width=wide,height=high, lat_0=center_lat, lon_0=center_lon)

# TRANSFORM GRID ACCORDING TO PROJECTION
#X_grid, Y_grid = m(lon, lat)
# MEAN GRID RESOLUTION
#dX = np.diff(X_grid, axis=-1).mean()
#dY = np.diff(Y_grid, axis=0).mean()

# PLOTTING
#fig = plt.figure(figsize=(12, 8))

# PLOT COASTLINES, LINES OF LONGITUDE,LATITUDE
#coastid = m.drawcoastlines()
#parallels = np.arange(-90., 91, 5)
#meridians = np.arange(-180, 181, 5)
#parid = m.drawparallels(parallels, color='lightgrey', labels=[1, 1, 0, 0])
#merid = m.drawmeridians(meridians, labels=[0, 0, 1, 1], fmt='%3d', color='lightgrey')

# PLOT BACKGROUND
# (IE: WHERE THERE IS ORBITAL COVERAGE IN THE INTERVAL)
#bgid = plt.pcolormesh(X_grid-dX/2., Y_grid-dY/2., np.ma.masked_where(~orbital_coverage,0*X_grid), vmin=-2, vmax=1, cmap='gray')

# SEVERAL DATA PLOTTING OPTIONS
# ORBITAL COVERAGE (TRUE / FALSE)
# pltid = plt.pcolormesh(X_grid-dX/2.,Y_grid-dY/2.,orbital_coverage,vmin = 0, vmax = 1, norm = None)
# TOTAL COLUMN LOADING (DU)
# pltid = plt.pcolormesh(X_grid-dX/2.,Y_grid-dY/2.,np.ma.masked_where(~cloud, mass_loading_mean[:,:,-1]),vmin = 0.1, vmax = 100., norm = c.LogNorm())
# INITIAL Z - SCORE
# pltid = plt.pcolormesh(X_grid-dX/2.,Y_grid-dY/2.,zscore_initial,vmin = -5, vmax = 5, norm = None)
# MEAN HEIGHT (km)
# pltid = m.pcolormesh(X_grid-dX/2.,Y_grid-dY/2.,np.ma.masked_where(~cloud, mean_height), vmin = 0 ,vmax=20, norm = None)
# MEDIAN HEIGHT (km)
#pltid = m.pcolormesh(X_grid-dX/2., Y_grid-dY/2., np.ma.masked_where(~cloud, median_height), vmin=0, vmax=20, norm=None)
#cbarid = plt.colorbar()
# plt.show()

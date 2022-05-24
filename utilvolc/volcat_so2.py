# volcat_so2.py
# Reads and plots CrIS SO2
# From Dave Hyman (dhyman2@wisc.edu)
# Modified by Allison Ring (allison.ring@noaa.gov)
import numpy as np
import xarray as xr
import pandas as pd
from scipy import integrate
import matplotlib.pyplot as plt
from utilhysplit.par2conc import fixlondf
# Opens the dataset

    

class volcatSO2L3:

    def __init__(self,fname):
        self.fname = fname
        # use to convert from DU to  g/m^2 
        self.conv = (2.6867E20) * (64.066) * (1/(6.022E23))

    def open_dataset(self):
        self.dset = xr.open_dataset(self.fname, decode_times=True)
        return self.dset 


    def plotscatter(self,interp=False):
        lon = self.pframe.lon.values
        lat = self.pframe.lat.values
        mass = self.pframe.mass.values
        plt.scatter(lon,lat,mass)
       

    def plotmass(self,interp=True):
        if interp:
           vals = self.mass_interp
        else:
           vals = self.mass
        cb = plt.pcolormesh(vals)
        plt.colorbar(cb)

    def get_points(self):
        self.mass = self.dset.so2_column_loading_fov
        self.height = self.dset.so2_height_fov
        self.mass_interp = self.dset.so2_column_loading_interp
        self.height_interp = self.dset.so2_height_interp
        self.time = self.dset.time_so2
        lon = self.dset.lon.values
        lat = self.dset.lat.values
        self.lon, self.lat = np.meshgrid(lon,lat)

    def points2frame(self):
        mass = self.mass.values.flatten()
        height = self.height.values.flatten()
        time = self.time.values.flatten()
        lon = self.lon.flatten()
        lat = self.lat.flatten()
        rlist = []
        for iii in np.arange(0,len(mass)):
            if not np.isnan(mass[iii]):
               rvalue = (self.conv*mass[iii], lon[iii], lat[iii], height[iii], time[iii])
               rlist.append(rvalue)
        pframe = pd.DataFrame(rlist)
        pframe.columns = ['mass', 'lon', 'lat','height','time']
        pframe = pframe.sort_values(['lon','lat','time'],ascending=True)
        #pframe = fixlondf(pframe,neg=False)  
        self.pframe = pframe
        return pframe


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

# Extracting  variables


def get_mass(dset):
    """ Extracts the mean mass loading of SO2
         Input: xarray of data
         Output: mass loading mean in g/m^2 (not DU) """
    mass_load_mean = dset.mass_loading_mean
    # Return units of g/m^2 not DU
    conv = (2.6867E20) * (64.066) * (1/(6.022E23))
    mass_loading_mean = mass_load_mean * conv
    mass_loading_mean.attrs['units'] = 'g/m^2'
    return mass_loading_mean


def get_meanhgt(dset):
    """ Extracts the height variable
         Input: xarray DataSet
         Output: xarray DataArray of mean height in km """
    height = dset.height
    height_pdf = dset.height_pdf
    tmp = height_pdf * height
    mean_height = tmp.integrate('z')
    mean_height.attrs['units'] = 'km'
    return mean_height


def get_medianhgt(dset):
    ##IN PROGRESS##
    height = dset.height
    zscore_initial = dset.zscore_initial
    detected = (zscore_initial >= 5.)
    height_cdf = np.concatenate(
        (0.*height_pdf[:, :, 0:1], integrate.cumtrapz(height_pdf, height, axis=-1)), axis=-1)
    median_height = (1. * detected) * dist_percentile(height, height_cdf, 50., axis=- 1)


def dist_percentile(X, CDF_X, percentile, axis=-1):
    # NEED THIS FUNCTION TO INTERPOLATE
    # THE MEDIAN HEIGHT ON EACH CELL
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
#
# POSITIVE DETECTIONS (TRUE / FALSE)
#
# ORBITAL COVERAGE (TRUE / FALSE)
# orbital_coverage = np.isfinite(zscore_initial)
# WHERE IS THE CLOUD (TRUE / FALSE)
# cloud = np.logical_and(detected, orbital_coverage)
# MEDIAN HEIGHT
#


#


# SET UP BASEMAP PROJECTION
# proj = 'laea'  # LAMBERT EQUAL AREA
# wide = 13E6  # meters wide
# high = 13E6  # meters high
# center_lat = 90  # projection center latitude
# center_lon = 153.2540  # projection center longitude (Raikoke)
# m = Basemap(resolution='l', projection=proj, width=wide,height=high, lat_0=center_lat, lon_0=center_lon)

# TRANSFORM GRID ACCORDING TO PROJECTION
# X_grid, Y_grid = m(lon, lat)
# MEAN GRID RESOLUTION
# dX = np.diff(X_grid, axis=-1).mean()
# dY = np.diff(Y_grid, axis=0).mean()

# PLOTTING
# fig = plt.figure(figsize=(12, 8))

# PLOT COASTLINES, LINES OF LONGITUDE,LATITUDE
# coastid = m.drawcoastlines()
# parallels = np.arange(-90., 91, 5)
# meridians = np.arange(-180, 181, 5)
# parid = m.drawparallels(parallels, color='lightgrey', labels=[1, 1, 0, 0])
# merid = m.drawmeridians(meridians, labels=[0, 0, 1, 1], fmt='%3d', color='lightgrey')

# PLOT BACKGROUND
# (IE: WHERE THERE IS ORBITAL COVERAGE IN THE INTERVAL)
# bgid = plt.pcolormesh(X_grid-dX/2., Y_grid-dY/2., np.ma.masked_where(~orbital_coverage,0*X_grid), vmin=-2, vmax=1, cmap='gray')

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
# pltid = m.pcolormesh(X_grid-dX/2., Y_grid-dY/2., np.ma.masked_where(~cloud, median_height), vmin=0, vmax=20, norm=None)
# cbarid = plt.colorbar()
# plt.show()

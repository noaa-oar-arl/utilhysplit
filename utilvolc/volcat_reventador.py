# volcat.py
# A reader for VOLCAT data using xarray
# For use with MONET
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import numpy as np
import numpy.ma as ma
import pandas as pd

"""
This script contains routines that open/read VOLCAT data in xarray format, 
manipulate arrays as necessary, and plots desirable variables.
-------------
Functions:
-------------
open_dataset: opens single VOLCAT file
open_mfdataset: opens multiple VOLCAT files
bbox: finds bounding box around data - used for trimming arrays
_get_time: set time dimension for VOLCAT data
_get_latlon: rename lat/lon, set coordinates of VOLCAT data
get_height: returns array of ash top height from VOLCAT
get_radius: returns array of ash effective radius from VOLCAT
get_mass:  returns array of ash mass loading from VOLCAT
get_atherr: returns array of ash top height error from VOLCAT
plot_height: plots ash top height from VOLCAT
plot_radius: plots ash effective radius from VOLCAT
plot_mass: plots ash mass loading from VOLCAT
------------
"""


def open_dataset(fname):
    """Opens single VOLCAT file"""
    print(fname)
    dset = xr.open_dataset(fname, engine='pynio',mask_and_scale=False, decode_times=False)
    dset = dset.rename({"lines":'y',"elements":'x'})
    #dset = dset.rename({"Dim1": 'y', "Dim0": 'x'})
    dset = _get_latlon(dset)
    dset = _get_time(dset)
    return dset


def open_dataset2(fname):
    """Opens single VOLCAT file in reventador format """
    print(fname)
    dset = xr.open_dataset(fname, mask_and_scale=False, decode_times=False)
    #dset = dset.rename({"Dim1":'y',"Dim0":'x'})
    #dset = _get_latlon(dset)
    #dset = _get_time(dset)
    return dset


def open_mfdataset(fname):
    """Opens multiple VOLCAT files"""
    print(fname)
    dset = xr.open_mfdataset(fname, concat_dim='time', decode_times=False, mask_and_scale=False)
    from glob import glob
    from numpy import sort
    files = sort(glob(fname))
    das = []
    for i in files:
        das.append(open_dataset(i))
    dset = xr.concat(das, dim='time')
    dset = _get_latlon(dset)
    dset = dset.rename({"lines": 'y', "elements": 'x'})
    return dset


def bbox(darray):
    """Returns bounding box around data
    Input: Must be dataarray
    Outupt: Lower left corner, upper right corner of bounding box
    around data"""

    import numpy as np
    arr = darray[0, :, :].values
    a = np.where(arr != -999.)
    box = (np.min(a[0]-3), np.min(a[1])-3, np.max(a[0]+3), np.max(a[1])+3)
    tmp = list(box)
    tmp2 = [0 if i < 0. else i for i in tmp]
    bbox = tuple(([tmp2[0], tmp2[1]], [tmp2[2], tmp2[3]]))
    return bbox


#def _get_latlon(dset):
#    dset = dset.set_coords(['latitude', 'longitude'])
#    return dset

def _get_latlon(dset):
      dset = dset.rename({'pixel_latitude':'latitude'})
      dset = dset.rename({'pixel_longitude':'longitude'})
      dset = dset.set_coords(['latitude','longitude'])
      #test = str(dset.Image_Date)[1:] + str(dset.Image_Time)
      #time = pd.to_datetime(test,format='%y%j%H%M%S')
      #dset['time'] = time
      lat = dset.latitude[:,:] * dset.latitude.scale_factor
      lon = dset.longitude[:,:] * dset.longitude.scale_factor
      dset['latitude'] = lat
      dset['longitude'] = lon
      return dset

def _get_time(dset):
      import pandas as pd
      test = str(dset.Image_Date)[2:-1] + str(dset.Image_Time)[1:5]
      time = pd.to_datetime(test,format='%y%j%H%M')
      dset = dset.assign_coords({'time':time})
      #dset = dset.set_coords(['time'])
      #dset.dims['time'] = time
      dset = dset.expand_dims(dim='time')
      return dset

def get_height(dset):
    """Returns array with retrieved height of the highest layer of ash."""
    """Default units are km above sea-level"""
    height = dset.vcat_ashprop_13_15_16_ash_top_height
    box = bbox(height)
    height = height[:, box[0][0]:box[1][0], box[0][1]:box[1][1]]
    height = height.where(height != height._FillValue)
    return height


def get_radius(dset):
    """Returns 2d array of ash effective radius"""
    """Default units are micrometer"""
    radius = dset.vcat_ashprop_13_15_16_ash_effective_radius
    box = bbox(radius)
    radius = radius[:, box[0][0]:box[1][0], box[0][1]:box[1][1]]
    radius = radius.where(radius != radius._FillValue)
    return radius


def get_mass(dset):
    """Returns 2d array of ash mass loading"""
    """Default units are grams / meter^2"""
    mass = dset.vcat_ashprop_13_15_16_ash_mass_loading
    box = bbox(mass)
    mass = mass[:, box[0][0]:box[1][0], box[0][1]:box[1][1]]
    mass = mass.where(mass != mass._FillValue)
    return mass


def mass_sum(dset):
    mass = get_mass(dset)
    mass2 = mass.where(mass > 0., 0.0).values
    mass_sum = np.sum(mass2)
    return mass_sum


def get_time(dset):
    time = dset.time_coverage_start
    return time


def get_atherr(dset):
    """Returns array with uncertainty in ash top height from VOLCAT."""
    """Default units are km above sea-level"""
    height_err = dset.vcat_ashprop_13_15_16_Ash_Height_Error
    height_err = height_err.where(height_err != height_err._FillValue, drop=True)
    return height_err

# Trim VOLCAT array


def trim_arrray(dset):
    """Trim the VOLCAT array around data
    Make smaller for comparison to HYSPLIT """

# Plotting variables


def plot_height(dset):
    """Plots ash top height from VOLCAT
    Does not save figure - quick image creation"""
    fig = plt.figure('Ash_Top_Height')
    title = 'Ash Top Height (km)'
    ax = fig.add_subplot(1, 1, 1)
    plot_gen(dset, ax, val='height', time=None, plotmap=True,
             title=title)


def plot_radius(dset):
    """Plots ash effective radius from VOLCAT
    Does not save figure - quick image creation"""
    fig = plt.figure('Ash_Effective_Radius')
    title = 'Ash effective radius ($\mu$m)'
    ax = fig.add_subplot(1, 1, 1)
    plot_gen(dset, ax, val='radius', time=None, plotmap=True,
             title=title)


def plot_mass(dset):
    fig = plt.figure('Ash_Mass_Loading')
    ax = fig.add_subplot(1, 1, 1)
    plot_gen(dset, ax, val='mass', time=None, plotmap=True,
             title='Ash_Mass_Loading')


def plot_gen(dset, ax,  val='mass', time=None, plotmap=True,
              title=None):
      """Plot ash mass loading from VOLCAT
      Does not save figure - quick image creation"""
      #lat=dset.latitude
      #lon=dset.longitude
      if val == 'mass':
          mass=get_mass(dset)
      elif val == 'radius':
          mass=get_radius(dset)
      elif val == 'height':
          mass=get_height(dset)
      if time and 'time' in mass.coords:
         mass = mass.sel(time=time)
      elif 'time' in mass.coords:
         mass = mass.isel(time=0)
      lat=mass.latitude
      lon=mass.longitude
      if plotmap:
          m=plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
          m.add_feature(cfeat.LAND)
          m.add_feature(cfeat.COASTLINE)
          m.add_feature(cfeat.BORDERS)
          plt.pcolormesh( lon, lat, mass, transform=ccrs.PlateCarree())
      else:
          plt.pcolormesh( lon, lat, mass )
      plt.colorbar()
      plt.title(title)
      plt.show()

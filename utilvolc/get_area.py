import numpy as np
import xarray as xr
from utilvolc import volcat

def get_area(
    dset, write=False, clip=True
    ):
    """Calculates the area (km^2) of each volcat grid cell
    Converts degress to meters using a radius of 6378.137km.
    Input:
    dset : xarray DataArray with coordinates latitude (y,x) and longitude (y,x)
    write: boolean (default: False) Write area to file
    clip: boolean (default: True) Use clipped array around data, reduces domain
    output:
    area: xarray containing gridded area values
    """
    d2r = np.pi / 180.0  # convert degress to radians
    d2km = 6378.137 * d2r  # convert degree latitude to kilometers
 

    #if clip == True:
    #    mass = volcat.get_mass(dset)
    #    mass = mass.isel(time=0)
    #    # mass = volcat.get_mass(dset)
    #else:
    #    mass = dset.ash_mass_loading.isel(time=0)
    if 'time' in dset.coords and 'time' in dset.dims:
        dset = dset.isel(time=0)
    # latitude (y,x)
    # longitude (y,x)    
    # make sure dimensions are in the right order.
    dset = dset.transpose('y','x')
    lat =  dset.latitude.values
    lon =  dset.longitude.values
    latrad = lat * d2r  # Creating latitude array in radians
    coslat = np.cos(latrad) * d2km * d2km
    shape = np.shape(dset)
    # Make shifted lat and shifted lon arrays to use for calculations
    lat_shift = lat[1:, :]
    lon_shift = lon[:, 1:]
    # Adding row/column of nans to shifted lat/lon arrays
    to_add_lon = np.empty([shape[0]]) * np.nan
    to_add_lat = np.empty([shape[1]]) * np.nan
    # Back to xarray for calculations
    lat2 = xr.DataArray(np.vstack((lat_shift, to_add_lat)), dims=["y", "x"])
    lon2 = xr.DataArray(np.column_stack((lon_shift, to_add_lon)), dims=["y", "x"])

    # area calculation
    area = abs(lat - lat2) * abs(abs(lon) - abs(lon2)) * coslat
    area.name = "area"
    area.attrs["long_name"] = "area of each lat/lon grid box"
    area.attrs["units"] = "km^2"
    # Reformatting array attributes
    if write == True:
        directory = wdir + "area/"
        match=''
        #if correct_parallax == True:
        #    areafname = "area_" + match + "_pc.nc"
        #else:
        areafname = "area_" + match + ".nc"
        #print(directory + areafname)
        area.to_netcdf(directory + areafname)
    dset.close()
    return area


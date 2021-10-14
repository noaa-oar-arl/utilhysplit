# tropomi.py
# A reader/converter for tropomi data

import xarray as xr
import numpy as np
import pandas as pd
from glob import glob


def make_nc(filelist, latfile, lonfile, varname='SO2', netdir=None, write=False):
    """ To convert HyunCheol's regridded text files to a netcdf.
    Inputs:
    filelist: list of full data filenames (including directory)
    latfile: full filename containing latitudes (string)
    lonfile: full filename containing longitudes (string)
    varname: name of data variable (string)
    netdir: default=None - netcdf file directory (string) if writing netcdf file
    write: boolean
    Outputs:
    list of xarray dataarrays 
    If write=True, netcdf file is made
    """
    # Read in lat/lon files
    lat = pd.read_csv(latfile, delimiter='\s+', header=None)
    lon = pd.read_csv(lonfile, delimiter='\s+', header=None)
    # Convert to numpy
    latnp = lat.to_numpy()
    lonnp = lon.to_numpy()
    # Convert to xarray
    latxr = xr.DataArray(data=latnp, dims=('x', 'y'), name='latitude')
    lonxr = xr.DataArray(data=lonnp, dims=('x', 'y'), name='longitude')
    # Creating original array
    latlonxr = xr.merge([latxr, lonxr])
    so2xr_list = []
    for f in filelist:
        data = pd.read_csv(f, delimiter='\s+', header=None)
        datanp = data.to_numpy()
        dataxr = xr.DataArray(data=datanp, dims=('x', 'y'), name=varname)

        so2xr = xr.merge([latlonxr, dataxr])
        so2xr = so2xr.set_coords('latitude')
        so2xr = so2xr.set_coords('longitude')
        start = f.rfind('/')+1
        end = f.find('.txt')
        so2xr.attrs['original file'] = f[start:end]

        if write:
            netfile = netdir+f[start:end]+'.nc'
            so2xr.to_netcdf(netfile)
            print(netfile+' created!')
        so2xr_list.append(so2xr)

    return so2xr_list

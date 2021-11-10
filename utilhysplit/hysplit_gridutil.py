import os
import datetime
import numpy as np
import xarray as xr
from monetio import hysplit

def attr_check(dset):
    for key in ['llcrnr latitude','llcrnr longitude','Latitude Spacing','Longitude Spacing']:
        if key not in dset.attrs: 
           return False
    return True

def align_grids(grid1,grid2):
    # takes care of making sure grids are aligned.
    # can use when grid definition not in the attributes.
    grida, gridb = xr.align(grid1, grid2, join='outer')

    # always update the attributes?
    #if not attr_check(grid1) or not attr_check(grid2):     
    attrs = calc_grids(grid1, grid2, verbose=False)
    grida.attrs.update(attrs)
    gridb.attrs.update(attrs)

    grida = hysplit.reset_latlon_coords(grida)
    gridb = hysplit.reset_latlon_coords(gridb)

    return grida, gridb

def compare_grids(c1,c2,verbose=False, tolerance=1e-5):
    """
    Returns:
    True if grids are the same
    False if grids are not the same
    """
    grid1 = find_grid_specs(c1)
    grid2 = find_grid_specs(c2)
    check = [] 
    for key in ['llcrnr latitude','llcrnr longitude','Latitude Spacing','Longitude Spacing']:
        check.append(np.abs(grid1[key]-grid2[key]))
    if not check_grids(check) and verbose:
       print('Grids do not match')
       print('check values', check)
    return check_grids(check)

def check_grids(check, tolerance=1e-5, verbose=False):
    for val in check:
        if val > tolerance: 
           return False
    return True

def convert_lon(grid):
    # if corner longitude is negative
    # convert to positive and then reset the coordinates
    attrs = find_grid_specs(grid)
    if attrs['llcrnr longitude'] < 0:
       attrs['llcrnr longitude'] += 360
       grid.attrs.update(attrs)
       grid2 = hysplit.reset_latlon_coords(grid)
       return grid2
    else:
       return grid

def find_grid_specs(grid,verbose=False):
    """
    grid : xarray DataSet or DataArray with regular lat-lon grid.
    Returns:
    attrs : dictionary with attributes specifying the grid.
    Note that the extent of the grid is just calculated from the
    maximum value of the latitude and longitude in the file and may
    not reflect the true extent of the original grid.

    Note : this function may not work  if there are nans in the
    latitude longitude field.
    """
    xv = grid.x.values
    iii = len(xv)-2
    lon1 = grid.isel(x=iii).longitude.values[0]
    lon2 = grid.isel(x=iii+1).longitude.values[0]
    dlon = np.abs(lon2 - lon1)
    xval = grid.isel(x=iii).x.values-1
    corner_lon = lon1 - xval*dlon
         
    yv = grid.y.values
    iii = len(yv)-2
    lat1 = grid.isel(y=iii).latitude.values[0]
    lat2 = grid.isel(y=iii+1).latitude.values[0]
    dlat = np.abs(lat2 - lat1)
  
    yval = grid.isel(y=iii).y.values-1
    corner_lat = lat1 - yval*dlat

    maxlon = np.max(grid.longitude.values)
    maxlat = np.max(grid.latitude.values)

    # add 5 to extent for buffer.
    nlat = (maxlat - corner_lat) / dlat + 5
    nlon = (maxlon - corner_lon) / dlon + 5

    # round dlon and dlat
    dlon = np.round(dlon*10000)/10000.0
    dlat = np.round(dlon*10000)/10000.0

    attrs = {'llcrnr latitude':  corner_lat,
             'llcrnr longitude': corner_lon,
             'Latitude Spacing': dlat,
             'Longitude Spacing' : dlon,
             'Number Lat Points' : nlat,
             'Number Lon Points' : nlon}

    return attrs

def calc_grids(c1,c2,verbose=False):
    """
    Returns grid specs that will cover two grids which are matching but
    the extent may be offset in space.
    """
    grid1 = find_grid_specs(c1)
    grid2 = find_grid_specs(c2)
    #check = [] 
    #for key in ['llcrnr latitude','llcrnr longitude','Latitude Spacing','Longitude Spacing']:
    #    check.append(np.abs(grid1[key]-grid2[key]))

    if not compare_grids(c1,c2,verbose=True):
       print('Warning: calc_grids : grids are not the same')
       print('grid1' , grid1)
       print('grid2' , grid1)
       return {}
    attrs = grid1
    nlat = np.max([grid1['Number Lat Points'],grid2['Number Lat Points']])
    nlon = np.max([grid1['Number Lon Points'],grid2['Number Lon Points']])

    attrs.update({'Number Lat Points' : nlat,
                  'Number Lon Points' : nlon})
    return attrs


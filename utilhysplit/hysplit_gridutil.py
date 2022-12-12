import datetime
import os

import numpy as np
import xarray as xr

from monetio.models import hysplit

"""
Change log

2022 5  Dec AMC updated find_grid_specs to be more robust.
2022 12 Dec AMC updated find_grid_specs to make sure dlat, dlon are positive

"""


def findclose(lat, latval):
    try:
        rval = np.where(np.isclose(lat, latval,atol=1e-3))[0][0]
    except Exception as eee:
        print(eee)
        print(lat)
        print(latval)
        print(np.isclose(lat,latval))
        rval = None
    return rval

def change_grid(dset,newattr):
    """
    Latitude and longitude values stay the same.
    x and y coordinate values change.
    Used when a new corner point for the grid is needed. 
    """
    latra = dset.isel(x=0).latitude.values
    lonra = dset.isel(y=0).longitude.values
    ynew, xnew = get_new_indices(latra,lonra,newattr)
    dset = dset.assign_coords(y=ynew)
    dset = dset.assign_coords(x=xnew)
    return dset

def get_new_indices(latra, lonra, attr):
    """
    return integer indices for latitude and longitude arrays.
    """
    #llcrnr_lat = attr['llcrnr latitude']
    #llcrnr_lon = attr['llcrnr longitude']
    #nlat = attr['Number Lat Points']
    #nlon = attr['Number Lon Points']
    #dlat = attr['Latitude Spacing']
    #dlon = attr['Longitude Spacing']

    # These are arrays for the entire grid.
    lat, lon = hysplit.getlatlon(attr)

    # lat is the latitude values for the entire grid
    # latra is the array for the subset of the grid in the array.
    # Add one because convention is that index starts at 1 not 0.
    ilat2 = [findclose(lat, latra[x])+1 for x in np.arange(0, len(latra))]
    ilon2 = [findclose(lon, lonra[x])+1 for x in np.arange(0, len(lonra))]
    return ilat2, ilon2

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
    # if grids are not compatible then calc_grids returns empty dictionary.
    attrs = calc_grids(grid1, grid2, verbose=False)
    grida.attrs.update(attrs)
    gridb.attrs.update(attrs)
    grida = hysplit.reset_latlon_coords(grida)
    gridb = hysplit.reset_latlon_coords(gridb)

    return grida, gridb

def compare_grids(c1,c2,verbose=False, tolerance=1e-3):
    """
    Returns:
    True if grids are the same
    False if grids are not the same
    """
    grid1 = find_grid_specs(c1)
    grid2 = find_grid_specs(c2)
    check = [] 
    keylist = ['llcrnr latitude','llcrnr longitude','Latitude Spacing','Longitude Spacing']
    for key in keylist:
        check.append(np.abs(grid1[key]-grid2[key]))
    if not check_grids(check,tolerance) and verbose:
       print('Grids do not match {:0.2e}'.format(tolerance))
       for val in zip(keylist,check):
           if val[1]>=tolerance:
               print('{} : {:0.2e} : {} : {}'.format(val[0],val[1],grid1[val[0]],grid2[val[0]]))
    return check_grids(check,tolerance)

def check_grids(check, tolerance=1e-5, verbose=False):
    for val in check:
        if np.abs(val) > tolerance: 
           return False
    return True


## DON'T use anymore. 
## grids should always be -180 to 180.
#def convert_lon(gridin):
#    """
#
#    """

    # -180 -179 -178 -177 -176
    #   1    2    3    4    5

    # 360   359  358  357  356

    # if corner longitude is negative
    # convert to positive and then reset the coordinates

    # first change the lat lon points.

#    grid = gridin.copy()

#    latra = gridin.sel(x=1).latitude.values
#    lonra = gridin.sel(x=1).longitude.values
    # all lon values get 180.0 added to them.
#    lonra = lonra + 180.0


#    attrs = find_grid_specs(grid)
#    if attrs['llcrnr longitude'] < 0:
#       attrs['llcrnr longitude'] += 360
#       attrs['Number Lon Points'] += 360
#       grid.attrs.update(attrs)
#       grid2 = hysplit.reset_latlon_coords(grid)
#       return grid2
#    else:
#       return grid

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
    maxlon = np.max(grid.longitude.values)
    maxlat = np.max(grid.latitude.values)

    xv = grid.x.values
    lon1 = grid.sel(x=xv[0]).longitude.values[0]
    lon2 = grid.sel(x=xv[1]).longitude.values[0]
    dlon = np.abs((lon2-lon1)/(xv[1]-xv[0]))
    xval = xv[0]
    corner_lon = lon1 - (xval-1)*dlon
    #print('HERE', lon1, xval-1, dlon, corner_lon)

         
    yv = grid.y.values
    lat1 = grid.sel(y=yv[0]).latitude.values[0]
    lat2 = grid.sel(y=yv[1]).latitude.values[0]
    dlat = np.abs((lat2 - lat1)/(yv[1]-yv[0]))
    yval = yv[0]
    corner_lat = lat1 - (yval-1)*dlat
    if verbose: print(lat1, (yval-1), dlat) 

    maxlon = np.max(grid.longitude.values)
    maxlat = np.max(grid.latitude.values)

    nlat = yv[-1]
    nlon = xv[-1]

    #corner_lon = np.round(corner_lon*1000)/1000.0
    #corner_lat = np.round(corner_lat*1000)/1000.0
    # hysplit always uses -180 to 180 grid.
    if corner_lon < -180.0001:
       corner_lon += 360
    if corner_lon >= 179.999:
       corner_lon += -360 
    # round dlon and dlat
    dlon = np.round(dlon*1000)/1000.0
    dlat = np.round(dlon*1000)/1000.0

    attrs = {'llcrnr latitude':  corner_lat,
             'llcrnr longitude': corner_lon,
             'Latitude Spacing': dlat,
             'Longitude Spacing' : dlon,
             'Number Lat Points' : nlat,
             'Number Lon Points' : nlon}

    return attrs

def calc_grids(c1,c2,tolerance=1e-3,verbose=False):
    """
    Returns grid specs that will cover two grids which are matching but
    the extent may be offset in space.
    The two grids must have the same attributes.
    """
    grid1 = find_grid_specs(c1)
    grid2 = find_grid_specs(c2)
    #check = [] 
    #for key in ['llcrnr latitude','llcrnr longitude','Latitude Spacing','Longitude Spacing']:
    #    check.append(np.abs(grid1[key]-grid2[key]))

    if not compare_grids(c1,c2,tolerance=tolerance,verbose=True):
       print('Warning: calc_grids : grids are not the same')
       print('grid1' , grid1)
       print('grid2' , grid2)
       return {}
    attrs = grid1
    nlat = np.max([grid1['Number Lat Points'],grid2['Number Lat Points']])
    nlon = np.max([grid1['Number Lon Points'],grid2['Number Lon Points']])

    attrs.update({'Number Lat Points' : nlat,
                  'Number Lon Points' : nlon})
    return attrs


def find_longitude_range(grid):
    """
    grid : xarray DataSet or DataArray with longitude and latitude values
    returns a corner longitude based on the longitude range.
    No longitude greater than 180 then return -180.0
    If longitude greater than 180 then return 0.
    Else return the minimum longitude value.
    """
    lonra = grid.longitude.values
    minlon = np.nanmin(lonra)
    maxlon = np.nanmax(lonra)
    if maxlon <=180.0:
       crnrlon = -180
    elif minlon >=0 and maxlon <=360:
       crnrlon = 0.0
    else:
       print('grid has non valid longitude values')
       print('{} to {}'.format(minlon,maxlon))
       crnrlon = minlon
    return crnrlon 


def compare_grids_compat(c1in, c2in,tol=1e-4,verbose=False):
    c1 = c1in.copy()
    c2 = c2in.copy()
    grid1 = find_grid_specs(c1)
    grid2 = find_grid_specs(c2)
    # latitude and longitude spacing must be the same 
    keylist = ['Latitude Spacing','Longitude Spacing']
    for key in keylist:
      check = np.abs(grid1[key]-grid2[key])
      if check > tol: 
         print('Spacing not the same for the grids')
         print('{} : grid1 {}: grid2{}'.format(key,grid1[key],grid2[key]))
         return False

    dlat = grid1['Latitude Spacing']
    dlon = grid1['Longitude Spacing']
    # difference in corner latitude / longitude must be multiple of
    # dlat, dlon spacing.
    latkey = 'llcrnr latitude'
    lonkey = 'llcrnr longitude'
    key =  latkey
    check = np.abs(grid1[key]-grid2[key])
    #check1 = check%dlat
    check1 = np.round(check/dlat) - check/dlat
    check2 = np.abs(check1-dlat)
    if check1 > tol and check2 > tol: 
       print('Cannot alter corners to get grids to match')
       print('{} : grid1 {}: grid2 {}'.format(key, grid1[key],grid2[key]))
       print('check1 : {}: check2 {} : tol {}'.format(check1, check2, tol))
       return False
    key =  lonkey
    check = np.abs(grid1[key]-grid2[key])
    #check1 = check%dlon
    check1 = np.round(check/dlon) - check/dlon
    check2 = np.abs(check2-dlon)
    if check1 > tol and check2 > tol: 
       print('Cannot alter corners to get grids to match')
       print('{} : grid1 {}: grid2 {}'.format(key, grid1[key],grid2[key]))
       print('check: {}, check1 : {}: check2 {} : tol {}'.format(check,check1,check2,tol))
       return False
    return True

def create_global_grid(dlat,dlon,crnrlon):
    # create a global grid.
    nlat = np.ceil(180.0/dlat) + 1
    nlon = np.ceil(360.0/dlon) + 1
    attrs = {'llcrnr latitude'   : -90,
             'llcrnr longitude'  : crnrlon,
             'Latitude Spacing'  : dlat,
             'Longitude Spacing' : dlon,
             'Number Lat Points' : nlat,
             'Number Lon Points' : nlon}
    return attrs

import os
import datetime
import numpy as np
import xarray as xr
from monetio import hysplit

def align_grids(grid1,grid2):
    # takes care of making sure grids are aligned.
    # can use when grid definition not in the attributes.
    checks, gridhash = calc_grids(grid1, grid2, verbose=False)
    grid1, grid2 = xr.align(grid1, grid2, join='outer')
    nlat = (gridhash['maxlat'] - gridhash['corner_lat']) / gridhash['dlat'] + 5
    nlon = (gridhash['maxlon'] - gridhash['corner_lon']) / gridhash['dlon'] + 5

    attrs = {'llcrnr latitude': gridhash['corner_lat'],
             'llcrnr longitude': gridhash['corner_lon'],
             'Latitude Spacing': gridhash['dlat'],
             'Longitude Spacing' : gridhash['dlon'],
             'Number Lat Points' : nlat,
             'Number Lon Points' : nlon}
    grid1.attrs.update(attrs)
    grid2.attrs.update(attrs)
    grid1 = hysplit.reset_latlon_coords(grid1)
    grid2 = hysplit.reset_latlon_coords(grid2)
    return grid1, grid2

def compare_grids(c1,c2,verbose=False, tolerance=1e-5):
    check, junk1 = calc_grids(c1,c2,verbose=verbose)
    for val in check:
        if val > tolerance: 
           print('Grids do not match')
           print(grid1)
           print(grid2)
           return False
    return True

def compare_grids_old(c1,c2,verbose=False):
    check, dlat, dlon, minlat, minlon, maxlat, maxlon  = calc_grids(c1,c2,verbose=verbose)
    if check[0] < 1e-5 and check[1] < 1e-5 and check[2]: return True
    else: return False

def find_grid_specs(grid,verbose=False):
    """
    grid : xarray DataSet or DataArray with regular lat-lon grid.
    """
    xv = grid.x.values
    iii = len(xv)-2
    lon1 = grid.isel(x=iii).longitude.values[0]
    lon2 = grid.isel(x=iii+1).longitude.values[0]
    dlon = np.abs(lon2 - lon1)
    xval = grid.isel(x=iii).x.values-1
    corner_lon = lon1 - xval*dlon
    #print(lon2,lon1,iii,xval,corner_lon)
 
    yv = grid.y.values
    iii = len(yv)-2
    lat1 = grid.isel(y=iii).latitude.values[0]
    lat2 = grid.isel(y=iii+1).latitude.values[0]
    dlat = np.abs(lat2 - lat1)
  
    yval = grid.isel(y=iii).y.values-1
    corner_lat = lat1 - yval*dlat

    maxlon = np.max(grid.longitude.values)
    minlon = np.min(grid.longitude.values)
    maxlat = np.max(grid.latitude.values)
    minlat = np.min(grid.latitude.values)
    rhash = {'dlon':dlon,
             'dlat':dlat,
             'corner_lat':corner_lat,
             'corner_lon':corner_lon,
             'minlon':minlon,
             'maxlon':maxlon,
             'minlat':minlat,
             'maxlat':maxlat}      
    return rhash

def calc_grids(c1,c2,verbose=False):
    grid1 = find_grid_specs(c1)
    grid2 = find_grid_specs(c2)
    check = [] 
    for key in ['dlat','dlon','corner_lat','corner_lon']:
        check.append(np.abs(grid1[key]-grid2[key]))

    dlon = grid1['dlon']
    dlat = grid1['dlat']
    dlon = np.round(dlon*10000)/10000.0
    dlat = np.round(dlon*10000)/10000.0
    minlat = np.nanmin([grid1['minlat'],grid2['minlat']])
    minlon = np.nanmin([grid1['minlon'],grid2['minlon']])
    maxlat = np.nanmax([grid1['maxlat'],grid2['maxlat']])
    maxlon = np.nanmax([grid1['maxlon'],grid2['maxlon']])
    corner_lat = grid1['corner_lat']
    corner_lon = grid1['corner_lon']
    rhash = {'dlon':dlon,
             'dlat':dlat,
             'corner_lat':corner_lat,
             'corner_lon':corner_lon,
             'minlon':minlon,
             'maxlon':maxlon,
             'minlat':minlat,
             'maxlat':maxlat}
    
    return check, rhash

def calc_grids_old(c1,c2,verbose=False):
    # c1 xarray 
    # c2 xarray
    # check to see if latitude and longitude on two different
    # xarrays have the same index number.
    # This is needed to use the align function in xarray.
    # if returns True then xarray align function can be used.
    check3 = True
    xv = c1.x.values
    xv2 = c2.x.values

    t1 = np.array([xv[i] - xv[i-1] for i in np.arange(1,len(xv))])
    t2 = np.array([xv2[i] - xv2[i-1] for i in np.arange(1,len(xv2))])
  
    if np.any(t1!=1):
       print('Warning. x values not sequential for grid 1')
       check3 = False
    if np.any(t2!=1):
       print('Warning. x values not sequential for grid 2')
       check3 = False

    overlap = [x for x in xv if x in xv2]
    if not overlap:
       print('WARNING: no overlap between grids')
       print('Min values', np.min(xv), np.min(xv2)) 
       print('Max values', np.max(xv), np.max(xv2)) 
    # compare minimum index values and get the largest
    miniii = np.max([np.min(xv),np.min(xv2)])
    # compare maximum index values and get the smallest
    maxiii = np.min([np.max(xv),np.max(xv2)])
    # get value where they overlap.
    if overlap: iii = overlap[0]
    #iii = miniii + np.ceil((maxiii - miniii)/2.0)
    x1 = c1.sel(x=iii).longitude.values[0]
    x2 = c2.sel(x=iii).longitude.values[0]
    x3 = c1.sel(x=iii+1).longitude.values[0]
    check1 = x2-x1

    yv = c1.y.values
    yv2 = c2.y.values
    t1 = np.array([yv[i] - yv[i-1] for i in np.arange(1,len(yv))])
    t2 = np.array([yv2[i] - yv2[i-1] for i in np.arange(1,len(yv2))])
  
    if np.any(t1!=1):
       print('Warning. y values not sequential for grid 1')
       check3 = False
    if np.any(t2!=1):
       print('Warning. y values not sequential for grid 2')
       check3 = False
    # compare minimum index values and get the largest
    miniii = np.max([np.min(yv),np.min(yv2)])
    # compare maximum index values and get the smallest
    maxiii = np.min([np.max(yv),np.max(yv2)])
    # get value halfway between.
    jjj = miniii + np.ceil((maxiii - miniii)/2.0)


    y1 = c1.sel(y=jjj).latitude.values[0]

    y2 = c2.sel(y=jjj).latitude.values[0]
    y3 = c2.sel(y=jjj+1).latitude.values[0]

    check2 = y2-y1
    checktot =  check1 < 1e-5 and check2 < 1e-5

    if not checktot: verbose=True

    if verbose:
        print('x range {} to {}'.format(xv[0],xv[-1]))
        print('x range {} to {}'.format(xv2[0],xv2[-1]))
        print('index {} values {} {} diff {}'.format(iii,x1,x2,x2-x1))
        print('lon spacing is {}'.format(x3-x1))

        print('y range {} to {}'.format(yv[0],yv[-1]))
        print('y range {} to {}'.format(yv2[0],yv2[-1]))
        print('index {} values {} {} diff {}'.format(iii,y1,y2,y2-y1))
        print('lat spacing is {}'.format(y3-y1))
        print('min lat {} {}'.format(np.min(c1.latitude.values),np.min(c2.latitude.values)))
        print('min lon {} {}'.format(np.min(c1.longitude.values),np.min(c2.longitude.values)))
        print('max lat {} {}'.format(np.max(c1.latitude.values),np.max(c2.latitude.values)))
        print('max lon {} {}'.format(np.max(c1.longitude.values),np.max(c2.longitude.values)))
    dlat = y3-y1
    dlon = x3-x1 
    minlon = np.min([np.min(c1.longitude.values), np.min(c2.longitude.values)])
    minlat = np.min([np.min(c1.latitude.values), np.min(c2.latitude.values)])
    maxlon = np.max([np.max(c1.longitude.values), np.max(c2.longitude.values)])
    maxlat = np.max([np.max(c1.latitude.values), np.max(c2.latitude.values)])
    return [check1, check2, check3], dlat, dlon, minlat,minlon,maxlat,maxlon

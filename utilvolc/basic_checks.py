import os
import datetime
import numpy as np
import xarray as xr

def compare_grids(c1,c2,verbose=False):
    check, dlat, dlon, minlat, minlon, maxlat, maxlon  = calc_grids(c1,c2,verbose=verbose)
    if check[0] < 1e-5 and check[1] < 1e-5 and check[2]: return True
    else: return False


def calc_grids(c1,c2,verbose=False):
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
 
    # compare minimum index values and get the largest
    miniii = np.max([np.min(xv),np.min(xv2)])
    # compare maximum index values and get the smallest
    maxiii = np.min([np.max(xv),np.max(xv2)])
    # get value halfway between.
    iii = miniii + np.ceil((maxiii - miniii)/2.0)
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

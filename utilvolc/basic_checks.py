import os
import datetime
import numpy as np
import xarray as xr


def compare_grids(c1,c2,verbose=False):
    # c1 xarray 
    # c2 xarray
    # check to see if latitude and longitude on two different
    # xarrays have the same index number.
    # This is needed to use the align function in xarray.
    # if returns True then xarray align function can be used.
    xv = c1.x.values
    xv2 = c2.x.values
    iii = xv[10]
    x1 = c1.sel(x=iii).longitude.values[0]
    x2 = c2.sel(x=iii).longitude.values[0]
    x3 = c1.sel(x=iii+1).longitude.values[0]
    check1 = x2-x1

    yv = c1.y.values
    yv2 = c2.y.values
    jjj = yv[10]
    y1 = c1.sel(y=jjj).latitude.values[0]
    y2 = c2.sel(y=jjj).latitude.values[0]
    y3 = c2.sel(y=jjj+1).latitude.values[0]
    check2 = y2-y1

    if verbose:
        print('x range {} to {}'.format(xv[0],xv[-1]))
        print('x range {} to {}'.format(xv2[0],xv2[-1]))
        print('index {} values {} {} diff {}'.format(iii,x1,x2,x2-x1))
        print('lon spacing is {}'.format(x3-x1))

        print('y range {} to {}'.format(yv[0],yv[-1]))
        print('y range {} to {}'.format(yv2[0],yv2[-1]))
        print('index {} values {} {} diff {}'.format(iii,y1,y2,y2-y1))
        print('lat spacing is {}'.format(y3-y1))
    if check1 < 1e-5 and check2 < 1e-5: return True
    else: return False

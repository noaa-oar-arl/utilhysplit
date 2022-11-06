#!/n-home/alicec/anaconda/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
from math import *
import sys
import numpy as np
import datetime
import pandas as pd
#import monet


def coarsen(small, large):
    """
    small : xarray
    large : xarray
    new : xarray
    """

    boundary = 'trim'
    d1 = np.abs(small.latitude.values[2][0] - small.latitude.values[1][0])
    d2 = np.abs(large.latitude.values[2][0] - large.latitude.values[1][0])
    num = int(np.round(d2 / d1))
    print('lat NUM', d2, d1, num)
    new = small.coarsen(y=num, boundary=boundary).mean()

    d1 = np.abs(small.longitude.values[0][1] - small.longitude.values[0][0])
    d2 = np.abs(large.longitude.values[0][1] - large.longitude.values[0][0])
    num = int(np.round(d2 / d1))
    print('lon NUM', d2, d1, num)
    new = new.coarsen(x=num, boundary=boundary).mean()

    #d1 = np.abs(small.z.values[2] - small.z.values[1])
    #d2 = np.abs(large.z.values[2] - large.z.values[1])
    #num = int(np.round(d2 / d1))
    #print('NUMz', d2, d1,num)
    #new = new.coarsen(z=num, boundary=boundary).mean()
    return new


def xslice(dra, longitude, volcat=False):
    # dra is xrarray
    latitude = dra.latitude.values[0][0]
    if volcat:
        # Since VOLCAT is on an irregular grid, making the mid point of the slice be
        # at the designated longitude rather than the start of the slice
        tmp = np.shape(dra)
        tmp1 = round(tmp[0]/2)
        tmp2 = round(tmp[1]/2)
        latitude = dra.latitude.values[tmp1][tmp2]
    xi, yi = dra.monet.nearest_ij(lat=latitude, lon=longitude)
    return dra.isel(x=xi)


def yslice(dra, latitude, volcat=False):
    # dra is xrarray
    longitude = dra.longitude.values[0][0]
    if volcat:
        # Since VOLCAT is on an irregular grid, making the mid point of the slice be
        # at the designated longitude rather than the start of the slice
        tmp = np.shape(dra)
        tmp1 = round(tmp[0]/2)
        tmp2 = round(tmp[1]/2)
        longitude = dra.longitude.values[tmp1][tmp2]
    xi, yi = dra.monet.nearest_ij(lat=latitude, lon=longitude)
    return dra.isel(y=yi)


def zslice(dra, height):
    # dra is xrarray
    import math
    iii = 0
    for val in dra.z.values:
        if height < val:
            break
        if math.isclose(height, val):
            break
        iii += 1
        print(iii, val, height)
    return dra.isel(z=iii)


def threshold(cra, tval=3, tp='log', fillna=True):
    # cra is xrarray
    """
    Apply a threshold to a DataArray
    if tp=='log' then tval indicates how many orders of magnitude
    shoudl be retained.
    """

    if tp == 'log':
        maxlog = np.max(np.log10(cra))
        minlog = np.min(np.log10(cra))
        thresh = 10**(maxlog - tval)
    else:
        thresh = tval
    cra = cra.where(cra >= thresh)
    if fillna:
        cra = cra.fillna(0)
    return cra

#!/n-home/alicec/anaconda/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
from math import *
import sys 
import numpy as np
import datetime
import pandas as pd

def xslice(dra, longitude ):
    # dra is xrarray 
    latitude = dra.latitude.values[0][0]
    xi, yi = dra.monet.nearest_ij(lat=latitude, lon=longitude) 
    return dra.isel(x=xi)

def yslice(dra, latitude ):
    # dra is xrarray 
    longitude = dra.longitude.values[0][0]
    xi, yi = dra.monet.nearest_ij(lat=latitude, lon=longitude) 
    return dra.isel(y=yi)

def zslice(dra, height ):
    # dra is xrarray 
    import math
    iii=0
    for val in dra.z.values:
        if height < val: break
        if math.isclose(height,val): break
        iii+=1 
        print(iii, val, height)
    return dra.isel(z=iii)

def threshold(cra, tval=3, tp='log',fillna=True):
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
    cra = cra.where(cra>=thresh)
    if fillna: cra = cra.fillna(0)
    return cra



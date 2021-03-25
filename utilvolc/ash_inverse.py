import sys
import os
import datetime
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import numpy as np
import numpy.ma as ma
import pandas as pd
from utilvolc import volcat
from monetio.models import hysplit

def testcase(tdir,vdir):
    fname = 'xrfile.445.nc'
    bezyloc = [160.587,55.978]
    vid = 'v300250'
    inverse = InverseAsh(tdir,fname,vdir,vid)
    return inverse

# get cdump output.
# pick time period.
# cdump output for that time period. 

class InverseAsh:

    def __init__(self, tdir, fname,vdir,vid):
        self.tdir = tdir   # directory for hysplit output
        self.fname = fname # name of hysplit output

        self.vdir = vdir   # directory for volcat data
        self.vid = vid   # volcano id.

        # keep volcat arrays for different averaging times. 
        self.volcat_hash = {}
        self.volcat_avg_hash = {}
        self.cdump_hash = {}

        # hysplit output. xarray. 
        cdump = xr.open_dataset(os.path.join(tdir,fname))
        temp = list(cdump.keys())
        cdump = cdump[temp[0]]
        cdump = cdump.isel(source=0)
        # get the mass loading.
        # TODO later may want to exclude certain levels.
        self.cdump = hysplit.hysp_massload(cdump)

    def get_volcat(self, daterange, verbose=False):
        vdir = self.vdir
        vid = self.vid
        tii = self.time_index(daterange[0])
        if tii not in self.volcat_hash.keys(): 
            das = volcat.get_volcat_list(vdir,daterange=daterange,vid=vid) 
            self.volcat_hash[tii] = das
        else:
            das = self.volcat_hash[tii]
        return das

    def clip(self,dummy):
        # clip ra according to where dummy has 0's.
        aaa = np.where(dummy != 0)
        a1 = np.min(aaa[0])
        a2 = np.max(aaa[0])
        b1 = np.min(aaa[1])
        b2 = np.max(aaa[1])
        return a1,a2,b1,b2

    def time_index(self,time):
        timelist = [pd.to_datetime(x) for x in self.cdump.time.values]
        try:
           iii = timelist.index(time)
        except:
           iii = None
        return iii 

    def prepare_one_time(self, daterange):
        vdir = self.vdir
        vid = self.vid
        tii = self.time_index(daterange[0])
        # assume cdump has times at beginning of averaging time.
        cdump_a = self.cdump.sel(time=daterange[0]) 
        
        # need to clip cdump and volcat.
        dummy = cdump_a.sum(dim='ens')
        a1,a2,b1,b2 = self.clip(dummy)
        dummy = dummy[a1:a2,b1:b2]
        #cdump_a = cdump_a[a1:a2,b1:b2]      
        # get list of volcat data arrays. 
        das = self.get_volcat(daterange)
        # regrid volcat to dummy grid.
        regrid_volcat = volcat.average_volcat(das,dummy)
        # average volcat values.
        avg = regrid_volcat.mean(dim='time')
        cdump_a = cdump_a[:,a1:a2,b1:b2]

        self.cdump_hash[tii] = cdump_a
        self.volcat_avg_hash[tii] = avg

        return  cdump_a, avg
        




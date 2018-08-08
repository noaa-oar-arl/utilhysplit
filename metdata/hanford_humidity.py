import os
import subprocess
import pandas as pd
import numpy as np
import pickle as pickle
from optparse import OptionParser
import datetime
import sys
from monet.obs import ish_mod
import monet.obs.obs_util as obs_util
import matplotlib.pyplot as plt
#from monet import MONET

"""
verify NWP met data with met station measurements.

INPUTS: Dates to run
        Areas to look at

"""

def find_range(ds):
    d1 = ds.index[0]
    de = ds.index[-1]
    dt = datetime.timedelta(hours=24)
    done=False
    rlist = []
    jjj=0
    while not done:
       ix = (ds.index > d1) & (ds.index < d1+dt) 
       dsub = ds[ix]
       rrr = dsub.max() - dsub.min()
       rlist.append(rrr)
       d1 = d1 + dt
       if d1 > de: done=True       
       jjj+=1
       if jjj> 400: done=True
    print(dsub)
    return rlist

def relh(x):
    #temp should be in Celsius
    #dewpoint should be in Celsius 
    dewpoint = x['dpt']
    temp = x['t']
    nnn = 7.5 * dewpoint / (237.3 + dewpoint) 
    mmm = 7.5 * temp / (237.3 + temp) 
    vap = 6.11 * 10**(nnn)
    satvap = 6.11 * 10**(mmm)
    rh = vap / satvap * 100
    return rh


def rplot(df):
    fig = plt.figure(1)
    ax = fig.add_subplot(3,1,1)
    ax.plot(df['time'], df['relh'], '-b.')
    ax2 = fig.add_subplot(3,1,2)
    ax2.plot(df['time'], df['dpt'], '-b.')
    ax3 = fig.add_subplot(3,1,3)
    ax3.plot(df['time'], df['t'], '-b.')
    plt.show()

class Mverify(object):

    def __init__(self, dates, area):
        """
        self.sources: pandas dataframe
            sources are created from a CEMSEmissions class object in the get_sources method.
        """
        ##dates to consider.
        self.d1 = dates[0]
        self.d2 = dates[1]
        self.dates = dates
        #area to consider
        self.area = area

        self.metdir = '/pub/archives/wrf27km/'
        self.hdir = '/n-home/alicec/Ahysplit/trunk/exec/'
        self.tdir = '/pub/Scratch/alicec/SO2/'

        self.fignum = 1

    def find_obs(self, isd=True):
        mdata = ish_mod.ISH()
        self.obs = mdata.add_data(self.dates, country=None, box=self.area, resample=False)
        self.obs['relh'] = self.obs.apply(relh, axis=1) 
        self.obs['latlon'] = str(self.obs['latitude']) + ' ' + str(self.obs['longitude'])
        #print(self.obs['latitude'].unique())
        #print(self.obs['longitude'].unique())
        #print(self.obs.columns.values)
        #rplot(self.obs)
        return self.obs

area = [46.3,-120,46.9, -119]
d1 = datetime.datetime(2012,1,1,0)
d2 = datetime.datetime(2012,12,1,0)
sv = Mverify([d1,d2], area)
obs = sv.find_obs()
obs.set_index('time', inplace=True)
relh = obs['relh']
rlist = find_range(relh)
plt.plot(rlist)
plt.show()




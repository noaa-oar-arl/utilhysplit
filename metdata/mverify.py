import os
import subprocess
import pandas as pd
import numpy as np
import pickle as pickle
from optparse import OptionParser
import datetime
import sys
import monet
#from  monet.obs import *
from monet.obs import ish_mod
#import monet.obs.obs_util as obs_util
import matplotlib.pyplot as plt
import seaborn as sns
#from monet import MONET

"""
verify NWP met data with met station measurements.

INPUTS: Dates to run
        Areas to look at

"""

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

def r2plot(df):
    sns.set()
    fig = plt.figure(1)
    ax2 = fig.add_subplot(2,1,1)
    ax2.plot(df['time'], df['dpt'], '-b.')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Dew Point')
    ax3 = fig.add_subplot(2,1,2)
    ax3.plot(df['time'], df['t'], '-b.')
    ax3.set_xlabel('Time')
    ax3.set_ylabel('Temperature')
    ax3 = fig.add_subplot(2,1,2)
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
        #self.obs = monet.obs.aqs.add_data(self.dates)
        self.obs = mdata.add_data(self.dates, country=None, box=self.area, resample=False)
        self.obs['relh'] = self.obs.apply(relh, axis=1) 
        self.obs['latlon'] = str(self.obs['latitude']) + ' ' + str(self.obs['longitude'])
        print(self.obs['latitude'].unique())
        print(self.obs['longitude'].unique())
        print(self.obs.columns.values)
        r2plot(self.obs)

parser = OptionParser()

parser.add_option('--run', action="store_true", dest="runh", default=False)
parser.add_option('--map', action="store_true", dest="emap", default=False)
parser.add_option('--ploto', action="store_true", dest="oplot", default=False)
parser.add_option('--plume', action="store_true", dest="plume", default=False)
parser.add_option('--plote', action="store_true", dest="eplot", default=False)
parser.add_option('--datem', action="store_true", dest="datem", default=False)
parser.add_option('--rundatem', action="store_true", dest="rundatem", default=False)
parser.add_option('--pickle', action="store_true", dest="pickle", default=False)
parser.add_option('--tcm', action="store_true", dest="tcm", default=False)
parser.add_option('--obs', action="store_true", dest="findobs", default=False)
parser.add_option('-x', action="store_false", dest="opkl", default=True)
parser.add_option('-d', type="string", dest="drange", default="2106:1:1:2016:2:1")
parser.add_option('-a', type="string", dest="area", default="HANFORD")
(options, args) = parser.parse_args()

opkl = options.opkl

temp = options.drange.split(':')
try:
    d1 = datetime.datetime(int(temp[0]), int(temp[1]), int(temp[2]), 0)
except:
    print('daterange is not correct ' + options.drange)
try:
    d2 = datetime.datetime(int(temp[3]), int(temp[4]), int(temp[5]), 0)
except:
    print('daterange is not correct ' + temp)

if options.area.strip() == 'ND':
    area = [-105.0, -97.0, 44.5, 49.5]
elif options.area.strip() == 'HANFORD':
    area = [46.3,-120,46.9, -119]
else:
    area = None

d1 = datetime.datetime(2012,1,1,0)
d2 = datetime.datetime(2012,12,1,0)
sv = Mverify([d1,d2], area)
sv.find_obs()


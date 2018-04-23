#!/n-home/alicec/anaconda/bin/python
from __future__ import print_function
import os
import pandas as pd
import numpy as np
from optparse import OptionParser
import datetime
from monet.obs import cems
import sys
from monet import monet
from monet.obs import aqs
import matplotlib.pyplot as plt
#from monet import MONET
from arlhysplit import emittimes
from arlhysplit import runh
import monet.plots as mp
from monet.verification.obsdata import StationData

"""
INPUTS: Dates to run
        Area to consider emissions from
        top level directories for
            1. code executables
            2. where output should be located.
            3. where any csv files are with emissions data.


STEPS

A. Preliminary.
1. Find Emission sources.
2. Fine measurement stations nearby.
3. Produce map of sources and measurements.
4. Create plots of emissions vs. time
5. Create plots of measurements vs. time

B. Create HYSPLIT input files and run HYSPLIT
1. Use emissions source to create HYSPLIT sources
2. Use measurement stations to create concenctration grids
   or just use a large concentration grid that covers alls stations.
3. Run HYSPLIT.
4. create dfile (datem file) for the measurements in the HYSPLIT output directories.
5. Run c2datem and collect output from cfiles in each directory.
6. Create a modeled measurements vs. time using the output from c2datem
7. Compare the modeled and measured concentrations.

"""
class SO2Verify(object):
    """This class can have
       self.cems is a CEMS object
 
    """
    def __init__(self, dates, area):
        self.d1 = dates[0]
        self.d2 = dates[1]
        self.area = area
        self.metdir = '/pub/archives/wrf27km/'
        self.hdir = '/n-home/alicec/hysplit/trunk/exec/'
        self.tdir = '/pub/Scratch/alicec/SO2/'
        self.tmap = None
        self.fignum = 1

    def find_emissions(self):
        mmm = monet.MONET()
        self.cems= mmm.add_obs(obs='cems')
        area = self.area
        efile = 'emission_02-28-2018_103721604.csv'
        self.cems.load(efile, verbose=True)
        self.cems.latlonfilter((area[2], area[0]), (area[3], area[2]))
        self.ehash = self.cems.create_location_dictionary()

   
    def sources(self):
        sources = self.cems.get_var(('so2','lbs'), loc=None, daterange=[self.d1, self.d2])
         

    def plot_emissions(self):
        for loc in self.ehash:
            data = self.cems.get_var(('so2','lbs'), loc=loc, daterange=[self.d1, self.d2])
            data = data * 0.453592
            data = data.rename("so2_kg")
            plt.plot(data, '--b.')   


    def map_emissions(self):
        if not self.tmap: self.create_map()
        fig = plt.figure(self.fignum)
        for loc in self.ehash:
            latlon = self.cems.get_location(loc)
            x, y = self.tmap(latlon[1], latlon[0])
            self.tmap.plot(x,y,'bo') 

    def create_map(self):
        from mpl_toolkits.basemap import Basemap 
        self.tmap = Basemap(llcrnrlon=area[0], llcrnrlat=area[2], urcrnrlon=area[1], urcrnrlat=area[3], projection='cyl', resolution='h')
        self.tmap.drawcoastlines()
        self.tmap.drawmapboundary()
        self.tmap.drawstates()

    def runHYSPLIT(self):
        sources = self.cems.get_var(('so2','lbs'), loc=None, daterange=[self.d1, self.d2])
        sources = sources * 0.453592  #convert from lbs to kg.
        runh.mult_run(1, self.d1, self.d2, sources, hysplitdir=self.hdir, topdirpath=self.tdir, metdirpath=self.metdir)


    def find_obs(self, verbose=False):
        aq = aqs.AQS()
        area = self.area
        aq.add_data([self.d1, self.d2], param=['SO2'], download=True)
        self.aaa = StationData(df=aq.df)
        self.aaa.latlonfilter((area[2], area[0]), (area[3], area[2]))
        if verbose: self.aaa.summarize()   
        self.ohash = self.aaa.get_lhash('siteid')
 
    def map_obs(self):
        if not self.tmap: self.create_map()
        for key in self.ohash:
            latlon = self.ohash[key]
            x, y = self.tmap(latlon[1], latlon[0])
            plt.text(x, y, str(key), fontsize=7, color='red')
        return 1


parser = OptionParser()

parser.add_option('--run', action="store_true", dest="runh", default=False)
parser.add_option('--map', action="store_true", dest="emap", default=False)

(options, args) = parser.parse_args()

d1 = datetime.datetime(2016,1,1,0)
d2 = datetime.datetime(2016,1,7,0)

#lon lon lat lat
area = [-110, -100, 40, 50]
sv = SO2Verify([d1,d2], area)
if options.runh:
    sv.find_emissions()
    sv.runHYSPLIT()

if options.emap:
    #sv.find_emissions()
    #sv.map_emissions()
    sv.find_obs()
    sv.map_obs()
    plt.show()

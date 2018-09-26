import os
import subprocess
import pandas as pd
import numpy as np
import pickle as pickle
from optparse import OptionParser
import datetime
import sys
import seaborn as sns
from monet.obs import cems_mod
from monet.obs import aqs_mod
from monet.obs import airnow
from monet.obs import ish_mod
import monet.obs.obs_util as obs_util
import matplotlib.pyplot as plt
from arlhysplit import runh
from arlhysplit.runh import date2dir
from arlhysplit.runh import source_generator
from arlhysplit.runh import create_plume
from arlhysplit.tcm import TCM
from arlhysplit.models import emittimes
from arlhysplit.models.datem import writedatem_sh
from arlhysplit.models.datem import frame2datem
from arlhysplit.models.datem import mk_datem_pkl
from monet.obs.epa_util import convert_epa_unit

"""
SEmissions class
"""


class SEmissions(object):
    """This class for running the SO2 HYSPLIT verification.
       self.cems is a CEMS object

       methods
       find_emissions
       get_sources
       plot - plots time series of emissions
       map  - plots locations of power plants on map
    """
    def __init__(self, dates, area, states=['nd']):
        """
        self.sources: pandas dataframe
            sources are created from a CEMSEmissions class object in the get_sources method.
        """
        ##dates to consider.
        self.d1 = dates[0]
        self.d2 = dates[1]
        self.states=states
        #area to consider
        self.area = area
        self.pfile = './' + 'obs' + self.d1.strftime("%Y%m%d.") + self.d2.strftime("%Y%m%d.") + 'pkl'
        self.tdir = '/pub/Scratch/alicec/SO2/'
        self.fignum = 1
        ##self.sources is a DataFrame returned by the CEMS class.
        self.cems= cems_mod.CEMS()
        self.sources = pd.DataFrame()

    def find(self, testcase=False, byunit=False, verbose=False):
        """find emissions using the CEMSEmissions class
           
           prints out list of emissions soures with information about them.

           self.ehash : self.ehash used in 
           self.stackhash : not used anywhere else?
        """
        print("FIND EMISSIONS")
        area = self.area
        if testcase:
            efile = 'emission_02-28-2018_103721604.csv'
            self.cems.load(efile, verbose=verbose)
        else:
            self.cems.add_data([self.d1, self.d2], states=self.states,
                                download=True, verbose=True)
        if area:
            self.cems.df = obs_util.latlonfilter(self.cems.df, (area[2], area[0]), (area[3], area[1]))
        self.ehash = self.cems.create_location_dictionary()
        
        self.stackhash = cems_mod.get_stack_dict(cems_mod.read_stack_height(), orispl=self.ehash.keys())

        #key is the orispl_code and value is (latitude, longitude)
        print('List of emissions sources found\n', self.ehash)
        print('----------------')
        namehash = self.cems.create_name_dictionary()
        self.meanhash={}    
        ##This gets a pivot table with rows time.
        ##columns are (orisp, unit_id) or just (orisp)
        data1 = self.cems.cemspivot(('so2_lbs'), daterange=[self.d1, self.d2],
                verbose=True, unitid=byunit)
        print('done with pivot----------------')
        print(self.ehash)
        
        for oris in self.ehash.keys():
            print('----------------')
            try:
                #data = data1[loc].sum(axis=1)
                data = data1[oris]
            except: 
                data = pd.Series()
            qmean = data.mean(axis=0)
            qmax = data.max(axis=0)
            print(namehash[oris])
            print('ORISPL ' + str(oris))
            print(self.ehash[oris])
            if not np.isnan(qmean):
                self.meanhash[oris] = qmean * 0.453592
            else: 
                self.meanhash[oris] = 0
            print('Mean emission (lbs)', qmean)
            print('Maximum emission (lbs)', qmax)
            print('Stack id, Stack height (meters)')
            for val in self.stackhash[oris]:
                print(str(val[0]) + ',    ' + str(val[1]*0.3048))


    def get_so2(self):
        sources = self.get_sources(stype='so2_lbs') 
        sources = sources * 0.453592  #convert from lbs to kg.
        return sources

    def get_heat(self):
        sources = self.get_sources(stype='heat_input (mmbtu)') 
        mult = 1.055e9 / 3600.0  
        sources = sources * mult  #convert from mmbtu to watts
        return sources
 
    def get_sources(self, stype='so2_lbs'):
        """ 
        Returns a dataframe with rows indexed by date.
        column has info about lat, lon, stackheight, orisp code
        values are
        if stype=='so2_lbs'  so2 emissions
        if stype='

        self.ehash is constructed in find. 
        """
        print("GET SOURCES")
        if self.cems.df.empty: self.find()
        sources = self.cems.cemspivot((stype), daterange=[self.d1, self.d2],
                  verbose=False, unitid=False)
        ehash = self.cems.create_location_dictionary()
        stackhash = cems_mod.get_stack_dict(cems_mod.read_stack_height(), orispl=self.ehash.keys())
        ##This block replaces ORISP code with (lat, lon) tuple as headers
        cnew = []
        columns=list(sources.columns.values)
        for oris in columns:
            sid, ht = zip(*stackhash[oris])
            ##puts the maximum stack height associated with that orispl code.
            newcolumn = (ehash[oris][0], ehash[oris][1], np.max(ht), oris)
            cnew.append(newcolumn)
        sources.columns=cnew
        #print(sources[0:20])
        return sources

    def new_create_emittimes(self, edate, schunks=-99, tdir='./'):
        """
        create emittimes file for CEMS emissions.
        edate is the date to start the file on.
        Currently, 24 hour cycles are hard-wired.
        """
        df = self.get_so2()
        dfheat = self.get_heat()
        locs=df.columns.values
        done = False
        while not done
            d1 = edate
            d2 = edate + datetime.timedelta(hours=schunks)
            dftemp = df.loc[d1:d2]
            hdf = dfheat[d1:d2]
            if dftemp.empty(): 
               break
            self.emit_subroutine(dftemp, hdf, tdir)       
            d1 = d2


    #def emit_subroutine(self, df, dfheat):


    def create_emittimes(self, edate, schunks=-99, tdir='./'):
        """
        create emittimes file for CEMS emissions.
        edate is the date to start the file on.
        Currently, 24 hour cycles are hard-wired.
        """
        df = self.get_so2()
        dfheat = self.get_heat()
        locs=df.columns.values
        for hdr in locs:
            print('HEADER', hdr)
            d1 = edate  #date to start emittimes file.
            dftemp = df[hdr]
            dfh = dfheat[hdr]
            
            oris = hdr[3]
            height = hdr[2]
            lat = hdr[0]
            lon = hdr[1]
            ##hardwire 1 hr duraton of emissions.
            record_duration='0100'
            area=1
            odir =  date2dir(tdir, edate, dhour=schunks, chkdir=True)
            ename = odir + 'EMIT' + str(oris) + '.txt'
            efile = emittimes.EmitTimes(filename=ename)
            ##hardwire 24 hour cycle length
            dt = datetime.timedelta(hours=24)
            efile.add_cycle(d1, "0024")
            for date, rate in dftemp.iteritems():
                if date >= edate:
                    heat=dfh[date]
                    check= efile.add_record(date, record_duration, lat, lon, height,
                                     rate, area, heat)
                    if not check: 
                       d1 = d1 + dt
                       efile.add_cycle(d1, "0024")
                       check2= efile.add_record(date, record_duration, lat, lon, height,
                                     rate, area, heat)
                       if not check2: 
                           print('sverify WARNING: record not added to EmitTimes')
                           print(date.strftime("%Y %m %d %H:%M"))
                           print(str(lat), str(lon), str(rate), str(heat))
                           break
            efile.write_new(ename)
  
    def plot(self):
        """plot time series of emissions"""
        if self.cems.df.empty: self.find()
        sns.set()
        namehash = self.cems.create_name_dictionary()
        data1 = self.cems.cemspivot(('so2_lbs'), daterange=[self.d1, self.d2],
                verbose=False, unitid=False)
        for loc in self.ehash:
            fig = plt.figure(self.fignum)
            ax = fig.add_subplot(1,1,1)
            data = data1[loc] * 0.453592
            ax.plot(data, '--b.')   
            plt.ylabel('SO2 mass kg')
            plt.title(str(loc) + ' ' + namehash[loc])
            self.fignum+=1

    def map(self, ax):
        """plot location of emission sources"""
        if self.cems.df.empty: self.find()
        plt.sca(ax)
        fig = plt.figure(self.fignum)
        for loc in self.ehash:
            lat = self.ehash[loc][0]
            lon = self.ehash[loc][1]
            #print('PLOT', str(lat), str(lon))
            #plt.text(lon, lat, (str(loc) + ' ' + str(self.meanhash[loc])), fontsize=12, color='red')
            pstr = str(loc) + ' \n' + str(int(self.meanhash[loc])) + 'kg'
            if self.meanhash[loc] > 1:
                ax.text(lon, lat, pstr, fontsize=12, color='red')
                ax.plot(lon, lat,  'ko')

    def testsources(self):
        if self.sources.empty: self.get_sources()
        sgenerator = source_generator(self.sources)
        for src in sgenerator:
            print(src)


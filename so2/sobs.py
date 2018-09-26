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
"""

def get_tseries(df, siteid, var='obs', svar='siteid', convert=False):
    qqq = df['siteid'].unique()
    df = df[df[svar] == siteid]
    df.set_index('time', inplace=True)
    mult=1
    if convert: mult=1/2.6178
    series = df[var] * mult
    return series


class SObs(object):
    """This class for running the SO2 HYSPLIT verification.
       self.cems is a CEMS object

       methods
       find_emissions

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

        self.metdir = '/pub/archives/wrf27km/'
        self.hdir = '/n-home/alicec/Ahysplit/trunk/exec/'
        self.tdir = '/pub/Scratch/alicec/SO2/'

        self.tmap = None
        self.fignum = 1
        self.rnum = 1  #run number for HYSPLIT runs.

        ##self.sources is a DataFrame returned by the CEMS class.
        self.cems= cems_mod.CEMS()
        self.sources = pd.DataFrame()

        ##self obs is a Dataframe returned by either the aqs or airnow MONET class.
        self.obs = pd.DataFrame()
        ##siteidlist is list of siteid's of measurement stations that we want to look at.
        self.siteidlist = []


    def plot(self):
        """plot time series of observations"""
        sra = self.obs['siteid'].unique()
        print('PLOT OBS')
        print(sra)
        sns.set()
        dist = []
        ptrue=True
        for sid in sra:
            ts = get_tseries(self.obs, sid, var='obs', svar='siteid', convert=True)
            ms = get_tseries(self.obs, sid, var='mdl', svar='siteid')
            dist.extend(ts.tolist())
            if ptrue:
                fig = plt.figure(self.fignum)
                #nickname = nickmapping(sid)
                ax = fig.add_subplot(1,1,1)
                #plt.title(str(sid) + '  (' + str(nickname) + ')' )
                plt.title(str(sid))
                ax.set_xlim(self.d1, self.d2)
                ts.plot()
                ms.plot()
                self.fignum +=1
        #sns.distplot(dist, kde=False)
        #plt.show()           
        #sns.distplot(np.array(dist)/2.6178, kde=False, hist_kws={'log':True})
        #plt.show()
        #sns.distplot(np.array(dist)/2.6178, kde=False, norm_hist=True, hist_kws={'log':False, 'cumulative':True})
        #plt.show()
 
    def save(self):
        pickle.dump(self.obs, open(self.pfile, "wb"))

    def find(self, verbose=False, pload=True, getairnow=False):
        """
        """
        if pload:
           try:
            print('trying to load file')
            print(self.pfile)
            self.obs = pickle.load(open(self.pfile,"rb"))
           except:
            pload=False
            print('Failed to load')
           if pload: print('loaded file')
           if not pload: print('Not loading pkl file' + self.pfile + "\n")
        if not pload:
            if getairnow: 
               aq = airnow.AirNow()
               aq.add_data([self.d1, self.d2], download=True)
            else:
               aq = aqs_mod.AQS()
               aq.add_data([self.d1, self.d2], param=['SO2'], download=True)
            self.obs = aq.df.copy()
            #print(aq.df['qualifier'].unique())
            self.obs = convert_epa_unit(self.obs, obscolumn='obs', unit='UG/M3')
            ##TO DO - something happens when converting units of mdl column.
            #self.obs = epa_util.convert_epa_unit(self.obs, obscolumn='mdl', unit='UG/M3')
        area = self.area
        if area: self.obs = obs_util.latlonfilter(self.obs, (area[2], area[0]), (area[3], area[1]))
        rt = datetime.timedelta(hours=72)
        self.obs = obs_util.timefilter(self.obs, [self.d1, self.d2+rt])
        siteidlist= np.array(self.siteidlist)
        if siteidlist.size: 
            print('HERE -------------------')
            self.obs = self.obs[self.obs['siteid'].isin(siteidlist)]
            iii = 0
            sym =[['r.','bd'],['r.','c.']]
            #for sid in siteidlist:
            #    print('HERE -------------------' + str(sid))
            #    #self.obs.check_mdl(sid, sym=sym[iii])
            #    iii+=1
            #    plt.show()
        if verbose: obs_util.summarize(self.obs)   
        self.ohash = obs_util.get_lhash(self.obs, 'siteid')
        if not pload: self.save_obs()          
        #self.obs.writecsv()

    def mkpkl(self):
        tdir = self.tdir + 'run' + str(self.rnum) + '/'
        mk_datem_pkl(self.rnum, self.d1, self.d2, tdir)


    def obs2datem(self):
        """
        write datemfile.txt. observations in datem format
        """
        sdate = self.d1
        edate = self.d2
        obs_util.write_datem(self.obs, sitename='siteid', drange=[sdate, edate])
 
    def map(self, ax):
        plt.sca(ax) 
        clr=sns.xkcd_rgb["cerulean"]
        #sns.set()
        #if not self.tmap: self.create_map()
        for key in self.ohash:
            latlon = self.ohash[key]
            #x, y = self.tmap(latlon[1], latlon[0])
            plt.text(latlon[1], latlon[0], str(key), fontsize=7, color='red')
            plt.plot(latlon[1], latlon[0],  color=clr, marker='*')
        return 1



#!/n-home/alicec/anaconda/bin/python
from __future__ import print_function
import os
import subprocess
import pandas as pd
import numpy as np
import cPickle as pickle
from optparse import OptionParser
import datetime
from monet.obs import cems
import sys
from monet import monet
from monet.obs import aqs
from monet.obs import airnow
import matplotlib.pyplot as plt
#from monet import MONET
from arlhysplit import emittimes
from arlhysplit import runh
from arlhysplit.datem import writedatem_sh
from arlhysplit.datem import frame2datem
from arlhysplit.runh import source_generator
from arlhysplit.runh import date2dir
from arlhysplit.datem import mk_datem_pkl
from arlhysplit.invert import TCM
#import monet.plots as mp
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
2. Find measurement stations nearby.
3. Produce map of sources and measurements.
4. Create plots of emissions vs. time
5. Create plots of measurements vs. time

B. Create HYSPLIT input files and run HYSPLIT
1. Use emissions source to create HYSPLIT sources
TODO - need to be able to specify stack height for each source.
TODO - resolve ambiguity in local / daylight savings time in the CEMS.

2. Use measurement stations to create concenctration grids
   or just use a large concentration grid that covers alls stations.
3. Run HYSPLIT.
4. create dfile (datem file) for the measurements in the HYSPLIT output directories.
TODO check that averaging time of measurements and HYSPLIT matches 
5. Run c2datem and collect output from cfiles in each directory.
6. Create a modeled measurements vs. time using the output from c2datem
7. Compare the modeled and measured concentrations.

EPA method code 60  = EQSA-0486-060

"""
def concmerge(hdir, files2merge):
    mfile = 'mfiles.txt'
    ofile = 'merged.bin'
    with open(mfile, 'w') as fid:
         for fl in files2merge:
             fid.write(fl + '\n')
    callstr =  hdir + '/conmerge -i' + mfile + ' -o' + ofile 
    subprocess.call(callstr, shell=True)
    return ofile


class PairedSites(object):

    def __init__(self, dates, sidlist, df):
        self.sidlist = sidlist


class SO2Verify(object):
    """This class for running the SO2 HYSPLIT verification.
       self.cems is a CEMS object
 
    """
    def __init__(self, dates, area):
        self.d1 = dates[0]
        self.d2 = dates[1]
        self.pfile = './' + 'obs' + self.d1.strftime("%Y%m%d.") + self.d2.strftime("%Y%m%d.") + 'pkl'
        self.area = area
        self.metdir = '/pub/archives/wrf27km/'
        self.hdir = '/n-home/alicec/hysplit/trunk/exec/'
        self.tdir = '/pub/Scratch/alicec/SO2/'
        self.tmap = None
        self.fignum = 1
        self.sources = pd.DataFrame()
        mmm = monet.MONET()
        self.cems= mmm.add_obs(obs='cems')
        #self.obs = pd.DataFrame()
        self.rnum = 1

    def find_emissions(self):
        #mmm = monet.MONET()
        #self.cems= mmm.add_obs(obs='cems')
        area = self.area
        efile = 'emission_02-28-2018_103721604.csv'
        self.cems.load(efile, verbose=False)
        if area:
            self.cems.latlonfilter((area[2], area[0]), (area[3], area[2]))
        self.ehash = self.cems.create_location_dictionary()

    def get_sources(self):
        if self.cems.df.empty: self.find_emissions()
        sources = self.cems.get_var(('so2','lbs'), loc=None, daterange=[self.d1, self.d2])
        self.sources = sources * 0.453592  #convert from lbs to kg.
        #print(self.sources)
        #sys.exit()

    def plot_emissions(self):
        for loc in self.ehash:
            data = self.cems.get_var(('so2','lbs'), loc=loc, daterange=[self.d1, self.d2])
            data = data * 0.453592
            data = data.rename("so2_kg")
            plt.plot(data, '--b.')   

    def plot_obs(self):
        sra = self.obs.df['siteid'].unique()
        print(sra)
        for sid in sra:
            fig = plt.figure(self.fignum)
            ax = fig.add_subplot(1,1,1)
            self.obs.plotloc(sid, svar='siteid')    
            ax.set_xlim(self.d1, self.d2)
            plt.show()           

    def map_emissions(self):
        if not self.tmap: self.create_map()
        fig = plt.figure(self.fignum)
        for loc in self.ehash:
            latlon = self.cems.get_location(loc)
            x, y = self.tmap(latlon[1], latlon[0])
            self.tmap.plot(x,y,'bo') 

    def create_map(self):
        from mpl_toolkits.basemap import Basemap 
        bfr = 0.5
        area = self.area
        self.tmap = Basemap(llcrnrlon=area[0]-bfr, llcrnrlat=area[2]-bfr, urcrnrlon=area[1]+bfr, urcrnrlat=area[3]+bfr, projection='cyl', resolution='h')
        self.tmap.drawcoastlines()
        self.tmap.drawmapboundary()
        self.tmap.drawstates()

    def runHYSPLIT(self):
        if self.sources.empty:
           self.get_sources() 
        #selfsources = self.cems.get_var(('so2','lbs'), loc=None, daterange=[self.d1, self.d2])
        #sources = sources * 0.453592  #convert from lbs to kg.
        runh.mult_run(self.rnum, self.d1, self.d2, self.sources, hysplitdir=self.hdir, topdirpath=self.tdir, metdirpath=self.metdir)
              
    def save_obs(self):
        pickle.dump(self.obs, open(self.pfile, "wb"))

    def find_obs(self, verbose=False, pload=True, airnow=True):

           
        if pload:
           try:
            #self.pfile = './' + 'obs' + self.d1.strftime("%Y%m%d.") + self.d2.strftime("%Y%m%d.") + 'pkl'
            print('trying to load file')
            self.obs = pickle.load(open(self.pfile,"rb"))
            #self.obs.latlonfilter((area[2], area[0]), (area[3], area[2]))
            #rt = datetime.timedelta(hours=72)
            #self.obs.timefilter([self.d1, self.d2+rt])
            #siteidlist = [380570124, 380570123, 380570118, 380570102, 380570004, 380650002]
            #self.obs.idfilter('siteid', siteidlist)
           except:
            pload=False
           if pload: print('loading obs.pkl file')
           if not pload: print('Not loading pkl file' + self.d1.strftime("%Y/%m/%d ") + self.d2.strftime("%Y/%m/%d\n"))
        #   sys.exit()
        if not pload:
            if airnow: 
               aq = airnow.AirNow()
            else:
               aq = aqs.AQS()
            aq.add_data([self.d1, self.d2], param=['SO2'], download=True)
            self.obs = StationData(df=aq.df)
            print('HERE A') 
            print(aq.df['qualifier'].unique())
            #self.obs.convert_units(unit='UG/M3')  #convert from ppb to ug/m3
            ##TO DO - something happens when converting units of mdl column.
            #self.obs.convert_units(obscolumn='mdl', unit='UG/M3')  #convert from ppb to ug/m3
        area = self.area
        if area: self.obs.latlonfilter((area[2], area[0]), (area[3], area[1]))
        rt = datetime.timedelta(hours=72)
        #self.obs.timefilter([self.d1, self.d2+rt])
       
        siteidlist=None 
        #siteidlist = [380570124, 380570123, 380570118, 380570102, 380570004, 380650002]
        #siteidlist = [381050103, 381050105, 380930101]
        #siteidlist = [381050103, 381050105] 
        siteidlist = [380930101]
        if siteidlist: 
            print('HERE -------------------')
            self.obs.idfilter('siteid', siteidlist)
            iii = 0
            sym =[['r.','bd'],['r.','c.']]
            for sid in siteidlist:
                print('HERE -------------------' + str(sid))
                self.obs.check_mdl(sid, sym=sym[iii])
                iii+=1
                plt.show()
            #self.obs.convert_units(unit='UG/M3')  #convert from ppb to ug/m3
        if verbose: self.obs.summarize()   
        self.ohash = self.obs.get_lhash('siteid')
        if not pload: self.save_obs()          
        self.obs.writecsv()

    def mkpkl(self):
        tdir = self.tdir + 'run' + str(self.rnum) + '/'
        mk_datem_pkl(self.rnum, self.d1, self.d2, tdir)

    def mktcm(self):
        fig = plt.figure(self.fignum)
        ax = fig.add_subplot(1,1,1)
        tdir = self.tdir + 'run' + str(self.rnum) + '/'
        tcm = TCM(self.d1, self.d2, snum=1, edir=self.hdir, pid=1)
        tcm.make_tcm(topdirpath=tdir)
        #tcm.print_rows()
        #tcm.image_tcm()
        #plt.show()
        tcm.get_conc(emissions=1)
        self.tcm = tcm

    def plot_tcm(self):
        clrs = ['--r*', '--b*', '--g.', '--k.', '--c.']
        iii=0
        for stn, time, conc in self.tcm.yield_conc():
            plt.plot(time, conc, clrs[iii], label=str(stn))
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles, labels)
            #dstr = time[0].strftime("%Y%m%d") + '-' + time[-1].strftime("%Y%m%d")
            iii+=1
            dstr = 't1'
        print('SAVING FIGURE')
        plt.savefig('tcm' + str(stn) + '-' + dstr + '.jpg')
        plt.show()

    def plot_both(self): 
        clrs = ['--r*', '--b*', '--g.', '--b.', '--c.']
        iii=0
        sra = self.obs.df['siteid'].unique()
        for stn, time, conc in self.tcm.yield_conc():
            print(stn)
            fig = plt.figure(self.fignum)
            ax = fig.add_subplot(1,1,1)
            plt.plot(time, conc, clrs[1], label=str(stn))
            handles, labels = ax.get_legend_handles_labels()

            sid = int(stn)
            obs = self.obs.plotloc(sid, svar='siteid', get=True)    
            plt.plot(obs, '--k', label='observations')
            ax.set_xlim(self.d1, self.d2)
            ax.legend(handles, labels)
            #iii+=1
            dstr='t1'
            #dstr = time[0].strftime("%Y%m%d") + '-' + time[-1].strftime("%Y%m%d")
            print('SAVING FIGURE')
            plt.savefig(str(stn) + '-' + dstr + '.jpg')
            plt.show()
            df1 = pd.DataFrame(zip(time, conc), columns=['date','vals']) 
            df1.set_index('date', inplace=True)
            df3 = pd.concat([df1, obs], axis=1, join='inner')
            loc = self.tcm.locations[self.tcm.locations['stationid'] == stn]
            print('LOCATION')
            print(loc)
            print(loc['meas_lat'])
            print(loc['meas_lon'])
            df3['lat'] = loc['meas_lat'].tolist()[0]
            df3['lon'] = loc['meas_lon'].tolist()[0]
            df3['sid'] = sid
            df3['duration'] = '0100'
            df3['altitude'] = 20
            df3.reset_index(inplace=True)
            df3.rename(index=str, columns={'index':'date'}, inplace=True)
            #print(df3[0:10])
            writeover=False
            if iii==0: writeover=True
            frame2datem('outfile.txt', df3, header_str='SO2', writeover=writeover)
            iii+=1 

        #plt.show()
        #sra = self.obs.df['siteid'].unique()
        #print(sra)
        #for sid in sra:
        #    self.obs.plotloc(sid, svar='siteid')    

    def rundatem(self):
        if self.sources.empty: self.get_sources()
        sgenerator = source_generator(self.sources)
        pdir = 'none'
        cfilera = []
        nidra = []
        iii=0
        for src in sgenerator:
            sdate = src.sdate
            newdir = date2dir(self.tdir + 'run'+str(self.rnum) + '/', sdate, chkdir=False) 
            if iii==0: 
               pdir = newdir
               os.chdir(newdir)
            #print('directory ' + str(iii) + ' ' + pdir + ' ' + src.nid)
            if newdir != pdir:  
               #self.obs.write_datem(sitename='siteid')
               print('running c2datem ' + newdir, cfilera)
               ##emissions are in kg (1e3 grams). so mult by 1e9 will put them in ug (1e-6 grams).
               writedatem_sh(cfilera, mult='1e9', mdl=self.hdir, add2ra=nidra)
               callstr = 'chmod u+x datem.sh'
               subprocess.call(callstr, shell=True)
               callstr = './datem.sh'
               subprocess.call(callstr, shell=True)
               nidra=[]
               cfilera=[]
            cfilera.append('cdump.' + src.nid)
            nidra.append(src.nid)
            os.chdir(newdir)
            pdir  = newdir
            iii+=1

    def obs2datem(self):
        #if self.obs.empty: self.find_obs()
        self.obs.write_datem()
        if self.sources.empty: self.get_sources()
        sgenerator = source_generator(self.sources)
        pdir = 'none'
        cfilera = []
        nidra = []
        rt = datetime.timedelta(hours=100)
        for src in sgenerator:
            sdate = src.sdate
            newdir = date2dir(self.tdir + 'run'+str(self.rnum) + '/', sdate, chkdir=False) 
            os.chdir(newdir)
            cfilera.append('cdump.' + src.nid)
            nidra.append(src.nid)
            if newdir != pdir and pdir != 'none':  
               self.obs.write_datem(sitename='siteid', drange=[sdate, sdate + rt])
               #self.obs.write_datem(sitename='siteid')
            pdir  = newdir
 
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
parser.add_option('--ploto', action="store_true", dest="oplot", default=False)
parser.add_option('--plote', action="store_true", dest="eplot", default=False)
parser.add_option('--datem', action="store_true", dest="datem", default=False)
parser.add_option('--rundatem', action="store_true", dest="rundatem", default=False)
parser.add_option('--pickle', action="store_true", dest="pickle", default=False)
parser.add_option('--tcm', action="store_true", dest="tcm", default=False)
parser.add_option('--obs', action="store_true", dest="findobs", default=False)
parser.add_option('-x', action="store_false", dest="opkl", default=True)
parser.add_option('-d', type="string", dest="drange", default="2106:1:1:2016:2:1")
parser.add_option('-a', type="string", dest="area", default="ND")
(options, args) = parser.parse_args()

opkl = options.opkl

temp = options.drange.split(':')
try:
    d1 = datetime.datetime(int(temp[0]), int(temp[1]), int(temp[2]), 0)
except:
    print('daterange is not correct ' + temp)
try:
    d2 = datetime.datetime(int(temp[3]), int(temp[4]), int(temp[5]), 0)
except:
    print('daterange is not correct ' + temp)
     

#d1 = datetime.datetime(2016,7,1,0)
#d2 = datetime.datetime(2016,12,31,23)
#d1 = datetime.datetime(2016,1,1,0)
#d2 = datetime.datetime(2016,1,31,23)
#d1 = datetime.datetime(2016,2,1,0)
#d2 = datetime.datetime(2016,4,30,23)

#d1 = datetime.datetime(2018,1,1,0)
#d2 = datetime.datetime(2018,1,31,23)

#d1 = datetime.datetime(2016,5,1,0)
#d2 = datetime.datetime(2016,7,31,23)
#d1 = datetime.datetime(2016,8,1,0)
#d2 = datetime.datetime(2016,10,31,23)
#d2 = datetime.datetime(2016,12,31,23)

#d1 = datetime.datetime(2016,1,1,0)
#d2 = datetime.datetime(2016,1,2,0)

#d1 = datetime.datetime(2016,2,14,0)
#d2 = datetime.datetime(2016,3,31,23)

#lon lon lat lat
#area = [-103.5, -98.5, 44.5, 49.5]
if options.area.strip() == 'ND':
    area = [-105.0, -97.0, 44.5, 49.5]
else:
    area = None
sv = SO2Verify([d1,d2], area)

##emissions are on order of 1,000-2,000 lbs (about 1,000 kg)
##10,000 particles - each particle would be 0.1 kg or 100g.
##0.05 x 0.05 degree area is about 30.25 km^2. 30.25e6 m^2.
##50 meter in the vertical gives 1.5e9 m^3.
##So 1 particle in a 0.05 x 0.05 degree area is 0.067 ug/m3.
##Need a 100 particles to get to 6.7 ug/m3.
##This seems reasonable.

if options.findobs:
    sv.find_obs(pload=opkl)

if options.runh:
    sv.find_emissions()
    sv.runHYSPLIT()

if options.emap:
    sv.find_emissions()
    sv.map_emissions()
    sv.find_obs(pload=opkl)
    sv.map_obs()
    #sv.obs2datem() 
    plt.show()

if options.eplot:
    sv.find_emissions()
    sv.plot_emissions()

if options.oplot:
   sv.find_obs(pload=opkl)
   sv.plot_obs()

if options.datem:
    sv.find_obs(pload=opkl)
    sv.obs2datem() 
    #sv.rundatem() 

if options.rundatem:
    sv.rundatem() 

if options.pickle:
   sv.mkpkl()

if options.tcm:
   sv.mktcm()
   sv.find_obs(pload=opkl)
   sv.plot_both()

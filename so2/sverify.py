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
#from monet import MONET
from arlhysplit import emittimes
from arlhysplit import runh
from arlhysplit.datem import writedatem_sh
from arlhysplit.datem import frame2datem
from arlhysplit.runh import source_generator
from arlhysplit.runh import date2dir
from arlhysplit.runh import create_plume
from arlhysplit.datem import mk_datem_pkl
from arlhysplit.tcm import TCM
from monet.obs.epa_util import convert_epa_unit

"""
INPUTS: Dates to run
        Area to consider emissions from
        list of states to consider emissions from
        top level directories for
            1. code executables
            2. where output should be located.
            3. (optional) where any csv files are with emissions data.

STEPS
A. Preliminary.
1. Find Emission sources.
2. Find measurement stations in area 
3. Produce map of sources and measurements.
4. Create plots of emissions vs. time
5. Create plots of measurements vs. time

B.1 trajectory runs

B.2 Dispersion runs 
Create HYSPLIT input files and run HYSPLIT
1. Use emissions source to create HYSPLIT sources
TODO - need to be able to specify stack height for each source.
       Stack height for sources not available from current data sources.
       These must be looked up individually.
TODO - resolve ambiguity in local / daylight savings time in the CEMS.
2. Use measurement stations to create concenctration grids
   OR just use a large concentration grid that covers alls stations.
3. Run HYSPLIT.
TODO - should these be unit runs so a TCM can be creatd later or do runs in
some kind of chunks. 

C. Evaluate HYSPLIT output and measurements.
1. create dfile (datem file) for the measurements in the HYSPLIT output directories.
TODO check that averaging time of measurements and HYSPLIT matches 
2. Run c2datem and collect output from cfiles in each directory.
3. Create a modeled measurements vs. time using the output from c2datem
4. Compare the modeled and measured concentrations.
TODO what kinds of statistics?
5. create plume plots (TO DO)

EPA method code 60  = EQSA-0486-060

"""



def create_map(fignum):
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    fig = plt.figure(fignum)
    proj = ccrs.PlateCarree()
    ax = plt.axes(projection=proj)
    gl = ax.gridlines(draw_labels=True, linewidth=2, color='gray')
    gl.ylabels_right=False
    gl.xlabels_top = False
    states = cfeature.NaturalEarthFeature(category='cultural',
             name='admin_1_states_provinces_lines', scale='50m',
             facecolor='none')
    ax.add_feature(states, edgecolor='gray')
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.LAKES)
    ax.add_feature(cfeature.RIVERS)
    ax.add_feature(cfeature.COASTLINE)
    return(ax)

def nickmapping(sid):
    """returns the number or letter
     Nick is using to identify the stations"""
    nhash = {}
    sid = int(sid)
    nhash[381050105] = 1
    nhash[381050103] = 2
    nhash[380070002] = 3
    nhash[380530111] = 4
    nhash[380650002] = 5
    nhash[380570004] = 6
    nhash[300530001] = 7
    nhash[380530104] = 8
    nhash[380150003] = 9
    nhash[380570124] = 'A'
    nhash[380570118] = 'B'
    nhash[380250003] = 'C'
    nhash[380130004] = 'D'
    nhash[380930101] = 'E'
    nhash[380570102] = 'F'
    nhash[380530002] = 'G'
    nhash[300830001] = 'H'
    nhash[380570123] = 'J'
    try:
        nickname = nhash[sid]
    except:
        nickname = ''
    return nickname

def get_stack_heights(sid):
     """sid can either be the orispl code (integer)
     or the station name.
     """
     ##coal creek: two units, same height
     ##chimney height 650 feet = 198 m
     ##https://swce.coop/wp-content/uploads/2016/07/Coal-Creek-Information.pdf
     hdict[6030] = 198
     ndict['coal creek'] = 6030


     rt = True
     try: 
         ht == hdict[sid]
     except:
         rt = False
     if rt: 
         return ht
     else: 
         try:
            orispl = ndict[sid.lower().strip()]
         except:
             return -999
         return hdict[orispl]           

def concmerge_old(hdir, files2merge):
    """
    hdir: string
    files2merge: list of strings
    """
    mfile = 'mfiles.txt'
    ofile = 'merged.bin'
    with open(mfile, 'w') as fid:
         for fl in files2merge:
             fid.write(fl + '\n')
    callstr =  hdir + '/conmerge -i' + mfile + ' -o' + ofile 
    subprocess.call(callstr, shell=True)
    return ofile

def get_tseries(df, siteid, var='obs', svar='siteid', convert=False):
    qqq = df['siteid'].unique()
    df = df[df[svar] == siteid]
    df.set_index('time', inplace=True)
    mult=1
    if convert: mult=1/2.6178
    series = df[var] * mult
    return series



class SO2Verify(object):
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

    #def add_siteidlist(self, slist):
    #    self.siteidlist.extend(slist)

    def find_emissions(self, testcase=False):
        """find emissions using the CEMSEmissions class
        """
        print("FIND EMISSIONS")
        area = self.area
        if testcase:
            efile = 'emission_02-28-2018_103721604.csv'
            self.cems.load(efile, verbose=True)
        else:
            self.cems.add_data([self.d1, self.d2], states=self.states, download=True, verbose=True)
        if area:
            self.cems.df = obs_util.latlonfilter(self.cems.df, (area[2], area[0]), (area[3], area[1]))
        self.ehash = self.cems.create_location_dictionary()
        print('List of emissions sources found\n', self.ehash)
        print('----------------')
        namehash = self.cems.create_name_dictionary()
        self.meanhash={}    
        ##This gets a pivot table with rows time.
        ##columns are (orisp, unit_id) or just (orisp)
        data1 = self.cems.cemspivot(('so2_lbs'), daterange=[self.d1, self.d2], verbose=True, unitid=False)
        for loc in self.ehash:
            print('----------------')
            try:
                #data = data1[loc].sum(axis=1)
                data = data1[loc]
            except: 
                data = pd.Series()
            ##print statements show info about emissions point.
            print(str(loc))
            print(namehash[loc])
            print(self.ehash[loc])
            if not np.isnan(data.mean(axis=0)):
                self.meanhash[loc] = data.mean(axis=0) * 0.453592
            else: 
                self.meanhash[loc] = 0
            print('Mean emission (lbs)', data.mean(axis=0))
            print('Maximum emission (lbs)', data.max(axis=0))
   
    def get_sources(self):
        """ 
        load self.sources attribute from the cems class dataframe so2 column.
        """
        print("GET SOURCES")
        if self.cems.df.empty: self.find_emissions()
        sources = self.cems.cemspivot(('so2_lbs'), daterange=[self.d1, self.d2], verbose=True, unitid=False)
        ehash = self.cems.create_location_dictionary()
        ##This block replaces ORISP code with (lat, lon) tuple as headers
        ##TO DO, add stack height and heat realease to tuple. 
        cnew = []
        columns=list(sources.columns.values)
        for ccc in columns:
            cnew.append(ehash[ccc])
        sources.columns=cnew
        #print(sources[0:20])
        self.sources = sources * 0.453592  #convert from lbs to kg.

    def plot_emissions(self):
        """plot time series of emissions"""
        sns.set()
        namehash = self.cems.create_name_dictionary()
        data1 = self.cems.cemspivot(('so2_lbs'), daterange=[self.d1, self.d2], verbose=True, unitid=False)
        for loc in self.ehash:
            fig = plt.figure(self.fignum)
            ax = fig.add_subplot(1,1,1)
            #data = self.cems.get_var(('so2_lbs'), loc=[loc], daterange=[self.d1, self.d2],  verbose=False)
            #data = self.cems.get_var(('so2','lbs'), loc=[loc])
            data = data1[loc] * 0.453592
            #data = data.rename("so2_kg")
            #print(type(data)) #ERASE
            #print(data)       #ERASE
            ax.plot(data, '--b.')   
            plt.ylabel('SO2 mass kg')
            plt.title(str(loc) + ' ' + namehash[loc])
            self.fignum+=1
       
    def plot_obs(self):
        """plot time series of observations"""
        sra = self.obs['siteid'].unique()
        print('PLOT OBS')
        print(sra)
        sns.set()
        fignum=1
        dist = []
        ptrue=True
        for sid in sra:
            ts = get_tseries(self.obs, sid, var='obs', svar='siteid', convert=True)
            ms = get_tseries(self.obs, sid, var='mdl', svar='siteid')
            dist.extend(ts.tolist())
            if ptrue:
                fig = plt.figure(fignum)
                nickname = nickmapping(sid)
                ax = fig.add_subplot(1,1,1)
                plt.title(str(sid) + '  (' + str(nickname) + ')' )
                ax.set_xlim(self.d1, self.d2)
                ts.plot()
                ms.plot()
                fignum +=1
        #sns.distplot(dist, kde=False)
        plt.show()           
        sns.distplot(np.array(dist)/2.6178, kde=False, hist_kws={'log':True})
        plt.show()
        sns.distplot(np.array(dist)/2.6178, kde=False, norm_hist=True, hist_kws={'log':False, 'cumulative':True})
        plt.show()

    def map_emissions(self, ax):
        """plot location of emission sources"""
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
            #latlon = self.cems.get_location(loc)
            #x, y = self.tmap(latlon[1], latlon[0])
            #self.tmap.plot(x,y,'bo') 


    def runHYSPLIT(self):
        """create a HYSPLIT run for each source in self.sources
           self.sources is a pandas dataframe with index being release data and columns headers showing location of releases.
           values are emission rates.
        """
        if self.sources.empty:
           self.get_sources() 
        #selfsources = self.cems.get_var(('so2','lbs'), loc=None, daterange=[self.d1, self.d2])
        #sources = sources * 0.453592  #convert from lbs to kg.
        runh.mult_run(self.rnum, self.d1, self.d2, self.sources, hysplitdir=self.hdir, topdirpath=self.tdir, metdirpath=self.metdir)

    def plotplume(self):
        runh.create_plume(self.rnum, self.d1, self.d2, self.sources, hysplitdir=self.hdir, topdirpath=self.tdir, verbose=True)
 
    def save_obs(self):
        pickle.dump(self.obs, open(self.pfile, "wb"))

    def find_obs(self, verbose=False, pload=True, getairnow=False):
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
        #print('column values of observation' , self.obs.columns.values)
        if area: self.obs = obs_util.latlonfilter(self.obs, (area[2], area[0]), (area[3], area[1]))
        rt = datetime.timedelta(hours=72)
        self.obs = obs_util.timefilter(self.obs, [self.d1, self.d2+rt])
        #print('Site ids in obs dataset', self.obs['siteid'].unique()) 
        siteidlist= np.array(self.siteidlist)
        #siteidlist = [380570124, 380570123, 380570118, 380570102, 380570004, 380650002]
        #siteidlist = [381050103, 381050105, 380930101]
        #siteidlist = [381050103, 381050105] 
        #siteidlist = [380930101]
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
        sns.set()
        clrs = ['--r*', '--b*', '--g.', '--b.', '--c.']
        clrs[0] = sns.xkcd_rgb["sky blue"] #model data with zeros filled in.
        clrs[1] = sns.xkcd_rgb["blue"]  #model data
        clrs[2] = sns.xkcd_rgb["bright green"]  #observation series
        clrs[3] = sns.xkcd_rgb["forest green"]  #obs series again
        iii=0
        print( self.obs.columns.values)
        sra = self.obs['siteid'].unique()
        for stn, time, conc in self.tcm.yield_conc():
            print('HERE', stn, type(time), type(conc))
            fig = plt.figure(self.fignum)
            ax = fig.add_subplot(1,1,1)

            #------------------------------------------------------------
            ##PLOTTING Model data
            plt.plot(time, conc, clrs[1], label=str(stn))
            handles, labels = ax.get_legend_handles_labels()
            #plt.plot(obs)
          
            #------------------------------------------------------------
            ##PLOTTING OBSERVATIONS 
            sid = str(stn)
            nickname = nickmapping(stn)
            plt.title(sid + ' (' + str(nickname) + ')') 
            if stn in sra:
               print('FOUND')
            else:
               print('NOT FOUND', type(sid), type(stn), stn, sra)
            ts = get_tseries(self.obs, sid, var='obs', svar='siteid')
            #obs = self.obs.plotloc(sid, svar='siteid', get=True)    
            #obs = obs_util.plotloc(self.obs, sid, svar='siteid', get=True)    
            ts.plot(color=clrs[2])
            #------------------------------------------------------------

            #plt.plot(ts, '--k', label='observations')
            ax.set_xlim(self.d1, self.d2)
            #ax.legend(handles, labels)
            #iii+=1
            dstr='t1'
            #dstr = time[0].strftime("%Y%m%d") + '-' + time[-1].strftime("%Y%m%d")
            #print('SAVING FIGURE')
            #plt.savefig(str(stn) + '-' + dstr + '.jpg')
            #plt.show()

            ##Done plotting-----------------------------------------------

            #df1 = pd.DataFrame(zip(time, conc), columns=['date','vals']) 
            tdict = {'date': time, 'vals' : conc}
            df1 = pd.DataFrame(tdict)
            df1.set_index('date') 
            #print('obs', ts[0:10])
            df1.set_index('date', inplace=True)
            #print('DF1', df1[0:10])

            ##TODO How is this join being done?
            ##is this only keeping point where there are both measurements and model data?
            #df3 = pd.concat([df1, ts], axis=1, join='inner')
            df3 = pd.concat([df1, ts], axis=1)
            df3.fillna(0, inplace=True)

            loc = self.tcm.locations[self.tcm.locations['stationid'] == stn]
            #df4 = df3.set_index('date')
            plt.plot(df3['obs'], clrs[3])
            plt.plot(df3['vals'], clrs[0])
            #print(loc['meas_lat'])
            #print(loc['meas_lon'])
            df3['lat'] = loc['meas_lat'].tolist()[0]
            df3['lon'] = loc['meas_lon'].tolist()[0]
            df3['sid'] = sid
            df3['duration'] = '0100'
            df3['altitude'] = 20
            df3.reset_index(inplace=True)
            df3.rename(index=str, columns={'index':'date'}, inplace=True)
            #print('DF3', df3[0:10])
            writeover=False
            if iii==0: writeover=True
            frame2datem('outfile.txt', df3, header_str='SO2', writeover=writeover)
            plt.show()
            iii+=1 

        #plt.show()
        #sra = self.obs.df['siteid'].unique()
        #print(sra)
        #for sid in sra:
        #    self.obs.plotloc(sid, svar='siteid')    


    def testsources(self):
        if self.sources.empty: self.get_sources()
        sgenerator = source_generator(self.sources)
        for src in sgenerator:
            print(src)

    def rundatem(self):
        """
        datemfile.txt which contains observations in datem format should already be written.
        Runs c2datem in each directory.  creates datem.sh which will run c2datem and append
        results in one file.
        """
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
               print('running c2datem ' + newdir, cfilera, self.hdir)
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
        """
        write datemfile.txt. observations in datem format
        """
        #if self.obs.empty: self.find_obs()
        obs_util.write_datem(self.obs)
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
               obs_util.write_datem(self.obs, sitename='siteid', drange=[sdate, sdate + rt])
               #self.obs.write_datem(sitename='siteid')
            pdir  = newdir
 
    def map_obs(self, ax):
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


parser = OptionParser()

parser.add_option('--run', action="store_true", dest="runh", default=False)
parser.add_option('--map', action="store_true", dest="emap", default=False)
parser.add_option('--ploto', action="store_true", dest="oplot", default=False)
parser.add_option('--plume', action="store_true", dest="plume", default=False)
parser.add_option('--plote', action="store_true", dest="eplot", default=False, \
                  help='plot emissions')
parser.add_option('--datem', action="store_true", dest="datem", default=False)
parser.add_option('--rundatem', action="store_true", dest="rundatem", default=False)
parser.add_option('--pickle', action="store_true", dest="pickle", default=False)
parser.add_option('--tcm', action="store_true", dest="tcm", default=False)
parser.add_option('--test', action="store_true", dest="runtest", default=False)
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
    print('daterange is not correct ' + options.drange)
try:
    d2 = datetime.datetime(int(temp[3]), int(temp[4]), int(temp[5]), 0)
except:
    print('daterange is not correct ' + temp)

if options.area.lower().strip() == 'nd':
    area = [-105.0, -97.0, 44.5, 49.5]
    state=['nd']
else:
    area = None
    state=[options.area.strip()]

sv = SO2Verify([d1,d2], area, state)


##emissions are on order of 1,000-2,000 lbs (about 1,000 kg)
##10,000 particles - each particle would be 0.1 kg or 100g.
##0.05 x 0.05 degree area is about 30.25 km^2. 30.25e6 m^2.
##50 meter in the vertical gives 1.5e9 m^3.
##So 1 particle in a 0.05 x 0.05 degree area is 0.067 ug/m3.
##Need a 100 particles to get to 6.7 ug/m3.
##This seems reasonable.

if options.runtest:
   sv.testsources()

if options.plume:
   sv.plotplume()

if options.findobs:
    sv.find_obs(pload=opkl)

if options.runh:
    sv.find_emissions()
    sv.runHYSPLIT()

if options.emap:
    ax = create_map(fignum=1)
    sv.find_emissions()
    sns.set()
    sv.map_emissions(ax)
    sv.find_obs(pload=opkl)
    sv.map_obs(ax)
    sv.obs2datem() 
    plt.show()

if options.eplot:
    sv.find_emissions()
    sv.plot_emissions()
    fignum = sv.fignum + 1
    fig = plt.figure(fignum)
    ax = create_map(fignum)
    sv.map_emissions(ax)
    plt.show()

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

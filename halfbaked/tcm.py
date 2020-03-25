# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import sys 
import numpy as np
#import matplotlib.pyplot as plt
import datetime
#from subprocess import call
#from os import path
import os
import pickle as pickle
import pandas as pd
import seaborn as sns

from arlhysplit.runh import RunParams
from arlhysplit.runh import date2dir

from arlhysplit.runh import stations_info
from arlhysplit.runh import source_generator
#from read_data import getdata
import matplotlib.pyplot as plt
import subprocess
from scipy.optimize import minimize
from scipy.stats import linregress
from pylab import rcParams

"""
NAME: invert.py
UID: p102
PGRMMR: Alice Crawford ORG: ARL
This code written at the NOAA Air Resources Laboratory
ABSTRACT: This code contains functions and classes to create concentrations as a function of time at a location from database
CTYPE: source code

-----------------------------------------------------------------------------------------------------
MAIN function


invert - opens pickled pandas dataframes 
helper functions:

cfunc: cost function which is used by the minimization routine.
sfunc: an alternative cost function.

apply_zlevs    :  averages concentration over given vertical levels.
apply_massfrac :  adjusts the multx column for the mass fraction of each particle size.
-----------------------------------------------------------------------------------------------------
NOTES:
works best when number of rows (measurements) about an order of magnitude more than number of columns (source points)
#This may be tricky since we have only 4 measurement points and looking at 14 source points.

Get rid of rows and columns that are all zeros.
Consider getting rid of columns (sources) which don't contribute to many sources.
Consider getting rid of rows (measurements) which don't 
Hourly time resolution of measurements may not be helpful since contain duplicate information.
"""

##------------------------------------------------------------------------------------##
def sourceid2str(sourceid):
    temp = sourceid.split('_')
    nnn = len(temp[0])-2
    latstr = temp[0][0:nnn] + '.' + temp[0][nnn:] 
    nnn = len(temp[1])-2
    lonstr = temp[1][0:nnn] + '.' + temp[1][nnn:nnn+2] 
    return latstr + ' ' + lonstr


def sourceid2latlon(sourceidlist):
    latlonlist = []
    for si in sourceidlist:
        temp = si.split('_')
        lat = int(temp[0]) / 100.0 
        lon = int(temp[1]) / -100.0 
        latlonlist.append((lat,lon))
    return latlonlist

def get_dri_info(ashtype, rh):
    """
    Returns parameters measured by DRI.
    """
    massflux = {}
    ##key is type of ash and relative humdity
    ##value is  a , b, uncertainties in a and b
    ashtype = ashtype.lower()

    ##Values of [a, b, err in a, err in b]
    ##LOG E = aLOG u* + b   (base 10 log)
    ##E = K u*^a
    ## where K = 10**b
    massflux['vtts',25] = [4.53,2.18,0.32,0.08] 
    massflux['vtts',50] = [4.30,2.21,0.11,0.03] 
    massflux['vtts',75] = [3.94,1.66,0.11,0.03] 
    ##10**2.18 = 151
    ##10**2.21 = 162
    ##10**1.66 = 45.7
    massflux['msh',25] = [4.92,3.19,0.45,0.11] 
    massflux['msh',50] = [5.53,3.13,0.34,0.09] 
    massflux['msh',75] = [5.48,2.85,0.19,0.05] 
    ##10**3.19 = 1548
    ##10**3.13 = 1349
    ##10**2.85 = 708

    ##Log E (mg m-2 s-1) = a Log u* (ms-1) + b
    try:
        values = massflux[ashtype, rh]
    except:
        print('WARNING input into get_dir_file in build_conc module not valid', ashtype, rh)
        sys.exit()
    return  massflux[ashtype, rh]




def apply_zlevs(df, zlevs, clevels, verbose=False):
    df = df[df['level'].isin(zlevs)]
    if len(zlevs)==1:
       return df
    else:  #need to add levels together.
        df['massL'] = df['vals'] * df['thickness']  #calculate mass loading for each level.
        df2 = df.groupby(['date','psize','sourceid','ustar','multx','sourcedate']).sum()  #add mass loadings together.
        #make sure all thicknesses the same - they won't be if there is zero concentration
        #for a certain layer. but want to add in the zero concentrations.
        topht = np.max(zlevs)
        #print topht
        tdict = {}  #dictionary with key being top level height and value being level thickness
        prevl = 0
        thickness = 0
        for cl in clevels:
            tdict[cl] = cl-prevl
            prevl = cl 
        for zl in zlevs:
            thickness += tdict[zl]
        df2['thickness'] = thickness
        df2['level'] = topht
        df2['vals'] = df2['massL'] / df2['thickness']  #convert mass loading to concentration.
        if verbose: print(df2[:50])
        df2.reset_index(inplace=True)
        df2.drop(['massL'], axis=1, inplace=True)      #drop the mass loading column. no longer needed.
        #print df2[:10]
        return df2



def apply_massfrac(df, mult, psizelist, massfrac, verbose=False):
    """ adjust the multx column for the mass fraction of each particle size.
    """
    df['multx'] = mult
    df['psize'] = df['psize'] * 10
    df['psize'] = df.psize.astype(int)
    #print 'apply mass frac'
    if verbose==1: print('A dfday -----------')
    if verbose==1: print(df[:5])
    psize = list(map(int, np.array(psizelist) * 10 )) #will use this as an index. Needs to be an integer.
    if len(psize) != len(massfrac):
        print('WARNING: build_conc2 - massfrac list not same size as psizelist')
        print(psize)
        print(massfrac)
        massfrac[0] = -99
        sys.exit()
    ##change multx to account for mass fraction of each particle size.
    if massfrac[0] >= 0:
        df.set_index('psize', inplace=True)
        for particle in zip(psize, massfrac):
         ##getting warning "value is trying to be set on a copy.
         ##however, tested and this seems to work to change value of mulx column in df dataframe
         ##use try, except because if no concentrations for that day then this will result in an error.
            try:
                df['multx'][particle[0]] *=  particle[1]
            except:
               pass
        df.reset_index(inplace=True)
    return df   

def apply_time_ave(vals, tm):
    """
    Takes output from build_conc and creates time average.
    tm should be an integer number of hours to average over.
    vals is output from build_conc  which is a pandas series with the date as an index.
    """
    ##fill in any missing times with 0 values
    vals = vals.asfreq('1H')
    vals.fillna(value=0, inplace=True) 

    if tm != 1:
        ##Resample to get desired time average.
        freqstr = str(tm) + 'H' 
        vals = vals.resample(freqstr).mean()
    return vals


def roundtime(dt):
    return datetime.datetime(dt.year, dt.month, dt.day, 0, 0)

def sfunc(xval, tcm, apr, meas, verbose=3):
    ##xval is the array of emissions. (should be 1d numpy array)
    ##apr is the first guess (apriori) (should be same length as xval and numpy array)
    ##meas are the measured values (should be as many as rows in tcm)
    apriori_err = 1000
    meas_err = 100
    mmm = tcm.fillna(value=0)
    mmm = mmm.as_matrix()
    nrow = mmm.shape[0]
    #print 'NROWS', nrow, mmm.shape
    #print len(xval)
    diffterm = 0
    for iii in range(0,nrow):
        forecast = (mmm[iii]*xval).sum() 
        if verbose < 2: print('FORECAST', forecast)
        diff = (forecast - meas[iii])**2   #difference between forecast and measurement squared
        if verbose < 2: print('DIFF' ,  diff, meas[iii])
        diffterm += diff / meas_err**2
    aterm = (xval - apr)/apriori_err**2
    aterm = aterm.sum()
    #return 0.5 * aterm + 0.5 * diffterm
    if verbose < 3: print(aterm, diffterm)
    return  diffterm


def cfunc(xval, tcm, apr, meas, merr= 100, verbose=3):
    ##xval is the array of emissions. (should be 1d numpy array)
    ##apr is the first guess (apriori) (should be same length as xval and numpy array)
    ##meas are the measured values (should be as many as rows in tcm)
    apriori_err = 100  #supposdely making this smaller will keep solution closer to apriori term.
    meas_err = merr  
    mmm = tcm.fillna(value=0)
    mmm = mmm.as_matrix()
    nrow = mmm.shape[0]
    #print 'NROWS', nrow, mmm.shape
    #print len(xval)
    diffterm = 0
    lval=100
    for iii in range(0,nrow):
        forecast = (mmm[iii]*xval).sum() 
        tmeas = meas[iii] 
        ##if both forecast and model value below value lval, then no difference
        ## between them. 
        if lval > 0:
           if forecast < lval and tmeas < lval: forecast = tmeas
           
        if verbose < 2: print('FORECAST', forecast)
        diff = (forecast - meas[iii])**2   #difference between forecast and measurement squared
        if verbose < 2: print('DIFF' ,  diff, meas[iii])
        diffterm += diff / meas_err**2
    aterm = (xval - apr)**2/apriori_err**2
    aterm = aterm.sum()
    #return 0.5 * aterm + 0.5 * diffterm
    if verbose < 3: print(aterm, diffterm)
    return  diffterm + aterm



class TCM(object):
    """
    creates the TCM for measurement dates between sdate and edate
    Methods

    set_defaults
    weed_sources - look at how concentration from individual sources correlate with measurements. 
                

    """
    #from ashfall_base_iceland import stations_info

    def __init__(self, sdate, edate, snum=5, edir="/n-home/alicec/amchysplit/trunk/exec/", pid=1):
       ##pid is just an identificaton number.
       self.sdate = sdate
       self.edate = edate
       self.pid = pid
       self.set_defaults()
       self.sources = True
       self.fignum=1  #set initial figure number to 1.
       self.edir = edir
       print('DATES', self.sdate, self.edate)
       rcParams['figure.figsize'] = 20, 5
       self.opc_sizerange = [0,50]  #this includes all sizes- OPC can detect up to 32 um.

    def set_psizes(self, psizelist):
        self.psizes = psizelist

    def reset_fignum(self, fignum=1):
        self.fignum=fignum

    def set_defaults(self):
        self.zlevs=[50]   #use concentration grid up to 50m.
        self.psizes = [5] ##look at 5 micron particles??
        self.tdir = '/pub/Scratch/alicec/DOE_ASHFALL/iceland/'
        self.rnum = 1
        #self.stations = ['heim','hval', 'hvol', 'gren']
        #self.stations = ['heim', 'hvol']
        self.mult = 1 

        self.pname = 'tcm.' + self.sdate.strftime("%Y%m%d") + '_' + self.edate.strftime("%Y%m%d") + '.' + str(self.pid).zfill(3)  + '.pkl'

    def plot_pickle(self):
        self.fignum=1
        xval, xbnds  = self.get_firstguess(ftype=self.ftype, imult=self.imult)
        conc = self.get_conc(self.res.x)
        self.compare(conc)
        self.fignum = self.fignum-4
        conc = self.get_conc(xval)     #plot concentrations with a-priori
        self.compare(conc, sym='-.b.')

        self.evsu(self.res.x, var='wind')   #plots emissions vs. ustar.
        self.evsu(self.res.x, var='ustar')  #plots emissions vs. ustar.
        self.fignum = self.fignum-1
        self.evsu(xval, sym='r.')
        self.evsuhexbin(self.res.x)
        self.write_emissions(self.res.x, fname='tcsolve.' + self.method+ '.txt')
        self.plot_sources(self.res.x, cplt=False)

    def image_tcm(self):
        """plots the tcm"""
        sns.set()
        fig = plt.figure(self.fignum)
        ax = sns.heatmap(self.tcm) 
        plt.title('TCM')
        #self.fignum +=1
        #fig = plt.figure(self.fignum)
        #ax = sns.heatmap(self.aug) 
        #plt.title('Augmented Matrix')
        self.fignum +=1

    def make_tcm(self, 
                run_num=-99, verbose=0,   
                psizes=[-99], zlevs=[-99],
                topdirpath = 'none',
                snum=3):
        """one of the main methods.
           creates the tcm matrix 
        """
        sdate = self.sdate
        edate = self.edate

        ##use defaults if parameters not input
        if run_num == -99: run_num = self.rnum 
        else: self.run_num=run_num
        if zlevs[0] == -99: zlevs=self.zlevs
        if psizes[0] == -99: psizes = self.psizes
        #if topdirpath == 'none': topdirpath = self.tdir

        rhash = RunParams(run_num, met_type='WRF', topdirpath=topdirpath)
        #mult=1e-20 ##this is what is needed to create unit mass.
        mult=1 ##use this so that values are closer to 
        #topdirpath = rhash.topdirpath
        done = False
        dt = datetime.timedelta(hours=24)
     
        ##sdate and edate are the measurement dates of interest.
        ##should look at sources up to 72 hours beforehand
        ##cdate controls which sourcedates are used.
        cdate = roundtime(sdate-datetime.timedelta(hours=0))
        iii=0
        self.locations = pd.DataFrame()
        while not done:
            outdir = date2dir(topdirpath,  cdate, dhour=24)  
            try:
                dfday = pickle.load(open(outdir + 'conc_daily.pkl', "rb"))   #load the pickle file.
            except:
                print('WARNING, file not loading ',  outdir + 'conc_daily.pkl')
                dfday = pickle.load(open(outdir + 'conc_daily.pkl', "rb"))   #load the pickle file.
            #dfday = dfday[dfday['sourceid'].isin(self.sourceidlist)]
            if verbose ==1 : print('using pickle file from', outdir)
            print(('working on', cdate, outdir))
            ##pick only sources which contribute to measurement dates of interest.
            #dfday = dfday[dfday['sourcedate'] >= sdate]
            #dfday = dfday[dfday['sourcedate'] <= edate]

            ##replace meas lat and meas lon with station name, 'station' column
            ##drop rows for stations not in the self.stations list.
            #dfday['station'] = 'none'
            print((dfday['meas_lon'].unique()))
            print((dfday['meas_lat'].unique()))
            print((dfday['stationid'].unique()))
            #for stn in self.stations:
            #for stn in dfday['stationid'].values:
                #mloc = stations_info(stn).location
                #print 'working on station ', stn 
                #at = 0.01
                #vvv = np.logical_and(np.isclose(dfday['meas_lat'],mloc[0], atol=at),np.isclose(dfday['meas_lon'],mloc[1], atol=at))
                #vvv = np.logical_and(np.isclose(dfday['meas_lat'],mloc[0]),np.isclose(dfday['meas_lon'],mloc[1]))
                #if np.any(vvv): print 'TRUE'
                #else: print 'NO VALUES'
                #dfday['station'][vvv] = stn

            #dfday = dfday[dfday['stationid'].isin(self.stations)]
            locations = dfday[['meas_lon', 'meas_lat', 'stationid']]
            locations.drop_duplicates(inplace=True)
            self.locations  = pd.concat([self.locations, locations], axis=0)
            dfday.drop(['meas_lat','meas_lon'], axis=1, inplace=True)

            #print('DF1') 
            #print(dfday[0:10])
            #pick particle size to use
            dfday=dfday[dfday['psize'].isin(['1'])]
            #dfday=dfday[dfday['level'].isin(['50'])]
            #print('level') 
            #print(dfday[0:10])
            ##vertical level processing.
            #if not dfday.empty:
            #    if zlevs[0] != -99:
            #        dfday = apply_zlevs(dfday, zlevs, rhash.clevels)
            #    else:      
            #        dfday = apply_zlevs(dfday, rhash.clevels, rhash.clevels)
            ##instead of picking source release dates.
            ##pick measurement dates
            #dfday2 = dfday[dfday['date'] >= sdate]
            #dfday2 = dfday2[dfday['date'] <= edate]
            dfday2 = dfday.copy()
            ##Now should have left only runs which contributed to measurement at certain location
            ## at a certain date/time.
            #dfday2.drop(['relh','u10m','v10m','uwnd1','vwnd1','uwnd2','vwnd2','pres','prss2','prss1','level','thickness','ustar'], axis=1, inplace=True)
          
            #dfday2['wspd'] = (dfday2['u10m']**2 + dfday2['v10m']**2)**0.5

            dfday2.drop(['level','thickness'], axis=1, inplace=True)
            #print 'here', dfday2 
            #dfday2['vals'] = dfday2['vals'] * mult
            if iii==0:
              dfday3 = dfday2.copy()
            else:
              dfday3 = pd.concat([dfday3, dfday2])
            #print dfday3[:20]
            cdate = cdate + datetime.timedelta(hours=24)
            #print 'TEST', cdate, edate
            if cdate > edate:
               done=True
            #print done
            iii+=1
        
        print('final')
        print('PIVOT ------------------')
        self.locations.drop_duplicates(inplace=True)
        if dfday3.empty: print('WARNING, dfday3 empty')
        self.tcm = pd.pivot_table(dfday3, values='vals', index=['stationid', 'date'], 
                                columns=['sourcedate','sourceid','psize'], aggfunc=np.sum)
         
        ##The TCM will not contain rows for which there are measurements but no sources contributed.

    def print_rows(self):
        print('TCM ROW VALUES', self.tcm.index.values)
        #print 'aug ROW VALUES', self.aug.index.values


    def plot_sources(self, emissions, cplt=True):
        """plots emissions from each source, along with ustar and wind speed."""
        sourcehash = {} #key is the sourceid. Value is list of tuples (date, emission, ustar)
        iii = 0
        sourceids = [] 
        ##initialize list for each sourceid
        for val in list(self.tcm.columns.values):
            sourcehash[val[1]] = []
        ##fill in list for each sourceid
        iii=0
        for val in list(self.tcm.columns.values):
            sourcehash[val[1]].append( (val[0], emissions[iii], val[3], val[4]))
            iii+=1
        for source in list(sourcehash.keys()):
            if cplt:
                fig = plt.figure(self.fignum)
                ax1 = fig.add_subplot(3,1,1)
                ax2 = fig.add_subplot(3,1,2)
                ax3 = fig.add_subplot(3,1,3)
            temp = sourcehash[source]
            temp = list(zip(*temp))  #now this goes into 
            dates = temp[0]
            emissions = temp[1]
            ustar = temp[2]
            wind = temp[3]
            corr1 =   np.corrcoef([emissions, ustar])
            corr2 = np.corrcoef([emissions, ustar])
            if cplt:
                ax1.plot(dates, emissions, '--ko')
                ax2.plot(dates, ustar,  '--ko')
                ax3.plot(dates, wind,  '--bo')
                plt.title(source + ' ' + str(corr1))
                self.fignum +=1
            print(source, 'correlation emissions and ustar' , corr1) 
            #print source, 'correlation emissions and wind' , corr2
 
    def print_columns(self):
        print('COLUMNS TCM', end=' ') 
        for val in list(self.tcm.columns.values):
            print(val)
        print('COLUMNS aug', end=' ') 
        for val in list(self.aug.columns.values):
            print(val)
 
    def info_tcm(self):
        ##This prints out the number of values not NAN in each column
        ##which shows how many measurements that source contributed to.
        #print type(tcm)
        column_cnt=  self.aug.count()  #count returns series with number of non-NA/null obs. over requested axis.
        row_cnt=  self.aug.count(axis=1)  #count returns series with number of non-NA/null obs. over requested axis.
        #print type(cnt)
        print("number of measurements that the source term contributes to")
        print('COLUMNS ------------------')
        for val in column_cnt.items():
            print(val[0], val[1])
        print('ROWS ------------------')
        for val in row_cnt.items():
            print(val[0], val[1])
        #for iii in range(1,10):
        #    print self.aug[iii-1: iii]
       
    def datamaxmin(self, data, trange):
        """given a data series, find max and min values within trange of the time"""

        ##data is a time series.
        data1 = data[self.sdate-trange: self.sdate + trange]
        minconc = np.min(data1)           
        maxconc = np.max(data1)           


    def create_heatmap(self, xx, yy, cc, lhash, sh=True, prnt=True):
        ##lhash is a dictionary which contains information about strings to put on plot.
        ##keys should be
        ##colorbar - label for colorbar
        ##xaxis - label for x axis
        ##yaxis - label for y axis
        ##figname - name to save figure to.
        ##title - title of graph
        data =  pd.DataFrame(data={'p': xx, 't': yy, 'c': cc})
        data = data.pivot(index='p', columns = 't', values= 'c')
        cb = sns.heatmap(data, cbar_kws={'label': lhash['colorbar']})
        plt.title(lhash['title'])
        plt.xlabel(lhash['xaxis'])
        plt.ylabel(lhash['yaxis'])
        plt.savefig(lhash['figname'])
        if sh: plt.show() 

     
    def write_gridsearch_output(self, valra, lhash):
        with open(lhash['fname'], 'w') as fid:
             fid.write(lhash['title']+ '\n')
             pstr = list(map(str, self.psizes))
             
             fid.write('Model Particle size ' +  str.join('um ', pstr) + 'um \n')
             fid.write(lhash['hdr'])
             ##print from lowest to highest value of cc. 
             for val in valra:
                 str2write = ''
                 for ttt in val:
                     if isinstance(ttt, int)  : str2write += str(ttt) + ' '
                     if isinstance(ttt, float)  : str2write += "{:.2e}".format(ttt) + ' '
                     if isinstance(ttt, str)  : str2write += ttt + ' '
                 fid.write(str2write.strip() +  '\n') 
                   

    def try_remove(self, lra, val):
        try:
           return lra.remove(val)
        except:
           return lra

        
    def get_source_conc(self, emissions): 
        """ Loops through each column of the TCM. For
           each column (source) returns concentration time series created by that source.
           returns concentration time series from each individual source.
           returns header on tcm column which descries the source."""
        mmm = self.tcm.fillna(value=0)
        headers = list(self.tcm.columns.values)
        mmm = mmm.as_matrix()
        ##array where each column shows concentration from each source.
        conctemp = emissions * mmm
        conctemp = conctemp.T  ##take transpose so can loop through rows.
        fff = []
        sdict = {}
        sidra = []
        for hdr in headers:
            print('HEADER', hdr)
            sid = hdr[1]
            sidra.append(sid)
        sidra = list(set(sidra))
        for sid in sidra:
            sdict[sid] = np.zeros_like(conctemp[0])
        for iii in range(0, conctemp.shape[0]):
            fofu = conctemp[iii]   #this is now concentration time series due to one source (lat lon, time).
            sid = headers[iii][1]
            sdict[sid]+= fofu  #add up  contributions from differnet times for each  lat lon source
        return sdict

        #    yield headers[iii], fofu

    def get_conc(self, emissions):
        """emissions is a numpy array with same length as columns in TCM"""
        mmm = self.tcm.fillna(value=0)
        mmm = mmm.as_matrix()
        fff = []
        ##multiply each row in the matrix by the emission vector.
        for iii in range(0, mmm.shape[0]):  #loop through rows.
            fofu = mmm[iii] * emissions
            forecast = fofu.sum()   #sum contribution from each source.
            fff.append(forecast)    #add concentration to time series.

        dates = self.tcm.index.tolist()
        temp1 = list(zip(*dates))
        stations = set(temp1[0])
        temp2 = list(zip(temp1[0], temp1[1], fff))  #list of tuples with station, date, value
        self.concdf = pd.DataFrame(temp2, columns = ['station', 'Date', 'value'])
        self.stations = stations
        #return temp2
        #print temp2
        #print '**************************'
        #for stn in stations:
        #    vals = [q for q in temp2 if q[0]==stn]
        #    #print stn, vals
        #    yield vals
        #return fff

    def yield_conc(self):
        for stn in self.stations:
            df = self.concdf[self.concdf['station'] == stn]
            yield stn, df['Date'], df['value']   

    def write_emissions(self, emissions, fname='tcsolve.txt'):
        with open(fname, 'w') as fid:
            fid.write('Date,  Result,\n')
            jjj=0
            for cval in emissions:
                 fid.write('T:'+str(jjj) + ',  ' + "{0:.3e}".format(cval) + '\n')
                 jjj+=1


    def westphal(self):
        """ 
        Uses westphal relationship to convert unit mass to ug.     
        """
        #from product.build_conc import mass_flux_eq
        area = 5e8   ##(111e3 m/degree * 0.30 degree)**2 * cosine(63 degrees)  = 5e8 m^2
        uuu = self.get_ustar(thresh=40)      
        mult = 1e-20 / self.mult 
        #ccc = mult * (uuu/100.0)**4 * 1e-5 * area * 3600 *1e9  #convert from unit mass to mass
        #partition of 0.8
        ccc = mult * (0.8 * uuu/100.0)**4 * 1e-5 * area * 3600 *1e9  #convert from unit mass to mass
        fff = self.get_conc(ccc)
        print('WESTPHAL')
        self.write_emissions(ccc, fname='tcsolve.WE.txt')
        self.compare(fff)
    

    def dri(self):
        """ 
        Uses DRI relationship msh 50% ash.    
        """
        #from product.build_conc import mass_flux_eq
        area = 5e8   ##(111e3 m/degree * 0.30 degree)**2 * cosine(63 degrees)  = 5e8 m^2
        uuu = self.get_ustar(thresh=40)      
        mult = 1e-20 / self.mult 
        #ccc = mult * (uuu/100.0)**5.53 * 1349.2 * area * 3600 *1e9 / 1.0e6  #convert from unit mass to mass
        ccc = mult * (0.6 * uuu/100.0)**5.53 * 1349.2 * area * 3600 *1e9 / 1.0e6  #convert from unit mass to mass
        ##lag covered.
        #ccc = mult * (uuu/100.0)**6.39 * 9.0472 * area * 3600 *1e9 / 1.0e6  #convert from unit mass to mass
        fff = self.get_conc(ccc)
        print('DRI')
        self.write_emissions(ccc, fname='tcsolve.WE.txt')
        self.compare(fff)


    def find_rank(self):
        """
        rank of matrix is dimension of vector space spanned by its columns. (number of linearly independent columns.
        This compares rank of the tcm and the augmented matrix. If the rank of the augmented matrix is less than the
        rank of tcm then the system has infinite number of solutions (under-determined). If the rank of augmented
        matrix is greater than the tcm then the system is inconsistent (over-determined), there is no exact solution. 
        If the ranks are equal then there is one solution. 

        over-determined systems are preferable for inverse modeling applications.
        """
        mmm = self.tcm.fillna(value=0)
        mmm = mmm.as_matrix()

        aug = self.aug.fillna(value=0)
        aug = aug.as_matrix()
        print('Rank TCM', np.linalg.matrix_rank(mmm), 'Rank Augmented', np.linalg.matrix_rank(aug))
        return np.linalg.matrix_rank(mmm)

######following methods create and analyze an obsra. Which is a pandas dataframe with matched observations and forecasts by date.

    def init_obsra(self, obs, fc):
        #initialize an obsra with observations and forecasts matched by date.
        #obs and fc are two time series.
        obsra = self.match_obs(obs,fc)
        return obsra

    def add_obsra(self, obs, fc, obsra):
        ##add more observaton forecast pairs to an obsra.
        ##this is useful when wanting to look at obs /forecasts for multiple measurement stations.
        newobsra = self.match_obs(obs, fc)
        obsra = pd.concat([obsra, newobsra])
        return obsra

    def match_obs(self, obs, fc):
        ##input series with dates and times of observations and forecasts.
        ##put them together into a pandas data frame matched by dates.
        ## remove points for which either and observation or a forecast is not available.
        data =  pd.DataFrame(data={'obs': obs, 'fc': fc})
        data.dropna(axis=0, inplace=True)  #remove nan values
        return data

    def resample_obsra(self, obsra, tm=24):
        freqstr = str(tm) + 'H' 
        obsra = obsra.resample(freqstr).mean()
        return obsra

    def find_corr(self, obsra):
        corr1 = obsra.corr()
        return corr1.at['fc','obs']

    def stats(self, s1, s2, thresh=-50, tvalue = 10, plotdata=False):
        ## To do. when s1 or s2 below a threshold, set the value to tvalue.
        ## This is so we don't take into account correlation with values below the
        ## threshold. May not be the best way to do this.
        if thresh > 0:
           below1 = s1 < thresh
           below2 = s2 < thresh
           s1[below1] = tvalue 
           s2[below2] = tvalue  
        compare = pd.concat([s1,s2], axis=1)  #creates a dataframe with columns of data want to compute correlation for.
        compare.dropna(axis=0, inplace=True)  #remove nan values
        corr = compare.corr()                 #calculate correlation.
        if plotdata:
           s1.plot(kind='line')
           s2.plot(kind='line')
           plt.show()
        return corr

    def find_stats(self, obsra, scale=1):
        ##computes fractional bias and rmse of an obsra.
        tvalue = 20
        thresh=100
        bg = 20
        hit = 0   #keep track of events both forecast and observed
        miss= 0   #keep track of events observed but not forecast
        alarm =0  #keep track of events forecast but not observed
        nonevent =0 
        #if thresh > 0:
        #   obsra.loc[obsra.fc < thresh, 'fc'] = tvalue
        #   obsra.loc[obsra.obs < thresh, 'obs'] = tvalue
        #data = newra
        obs = obsra.obs
        fc = obsra.fc * scale
        mse = 0
        nnn = 0
        obsave =0
        fcave = 0
        for pair in zip(obs, fc):
            mse += (pair[0] - pair[1])**2
            nnn+=1
            obsave += pair[0]
            fcave += pair[1]
            if (pair[0]-bg >= thresh and pair[1] >= thresh):
               hit +=1
            elif (pair[0]-bg >= thresh and pair[1] < thresh):
               miss +=1
            elif (pair[0]-bg < thresh and pair[1] > thresh):
               alarm +=1
            elif (pair[0]-bg < thresh and pair[1] < thresh):
               nonevent +=1
            
        if nnn>0:
            obsave = obsave / float(nnn)
            fcave = fcave / float(nnn)
            fb = 2.0 * (fcave - obsave) / (fcave + obsave)
            rmse =  (mse / float(nnn))**0.5
            if (alarm+hit)> 0:
                far =  alarm/ float(alarm + hit)      #fraction of forecasts which are wrong.
            else:
                far = -1
            if(hit+miss) > 0:
                pod = hit /  float(hit + miss)       #fraction of observations which are detected         
            else:
                pod = -1
            #css = hit/float(hit+alarm) - miss/float(miss+nonevent)
            #tss = (hit*nonevent - alarm*miss) / float(hit+miss) / float(alarm+nonevent)
        else:
            print('WARNING, no pairs to find statistics for')
            rmse = -99
            fb = -99
            far =-1
            pod = -1

        dhash = {}
        dhash['rmse'] = rmse
        dhash['fb']= fb
        dhash['css'] = -99
        dhash['tss'] = -99
        return dhash


    def cdf(self, data):
        sdata = np.sort(data)
        y = np.arange(sdata.size)/  float(sdata.size)
        return sdata, y 

    def scale_results(self, obsra,  plotdata='none', method=3):
        ##finds scaling factor for forecast.
        ##Method 1 simply match max observation to max forecast.
        ##Method 2 return slope of line of observations vs. forecast.
        ##plots observations vs. forecast and returns slope
        ##Method 3 uses CDF matching.
        if method==0:
            return 1, [0], obsra      
        elif method==1:
            #scale simply by matching the maxminum observation to the maximum forecast.
            max_obs = np.max(obsra.obs)
            max_fc = np.max(obsra.fc)
            scale =  max_obs / float(max_fc)
        elif method==2 or method ==3:
           #newra = obsra[obsra.obs >= 50]
           #newra = newra[newra.fc > 0]
           newra = obsra
           poly = np.polyfit(newra.fc, newra.obs, 1) #fit line to data.
           scale = poly[0] #slope of the line.
           if plotdata!='none':
              fig = plt.figure(self.fignum)
              ax = fig.add_subplot(1,1,1)
              plt.plot(obsra.fc, obsra.obs, 'b.')
              mval = np.max([np.max(obsra.fc), np.max(obsra.obs)]) 
              xlin = np.arange(1, mval,10)
              ylin = poly[0] * xlin  + poly[1]
              plt.plot(xlin, ylin, '-r')
              ax.set_ylabel('Observations')
              ax.set_xlabel('Forecast')
              plt.title(plotdata)
              plt.savefig(plotdata + '.jpg')
              plt.show()
           if method==2: 
              obs = obsra['obs']
              fc  = obsra['fc']
              fc = fc * scale
              #poly1 = np.polyfit(robs, diff, 2)
              #y2 = poly1[0]*obs**2 + poly1[1]*obs + poly1[2]
              scaled_obsra = self.init_obsra(obs, fc)
              return scale, poly , scaled_obsra

        if method==3:   
           ##CDF matching of the observations and model forecsat.
           robs = np.sort(newra.obs)
           rfc = np.sort(newra.fc)* scale   #use scaling from linear regression.
           diff = rfc - robs                
           pfit = 1
           poly = np.polyfit(robs, diff, pfit)

           obs = obsra['obs']
           fc  = obsra['fc']
           fc = fc * scale
           if pfit==2:
              fc = fc - (poly[0]*fc**2 + poly[1]*fc + poly[2])
           elif pfit==3:
              fc = fc - (poly[0]*fc**3 + poly[1]*fc**2 + poly[2]*fc + poly[3])
           elif pfit==0:
              fc = fc - (poly[0])
           scaled_obsra = self.init_obsra(obs, fc)
           sp=False
           if plotdata == 'all': sp==True
           sp=True
           ##This block plots differences between cdfs and fit. 
           if sp:
               fig = plt.figure(1)
               ax1 = fig.add_subplot(1,1,1)
               ax1.plot(robs, diff,'-b')
               poly1 = np.polyfit(robs, diff, 2)
               y2 = poly1[0]*robs**2 + poly1[1]*robs + poly1[2]
               ax1.plot(robs, y2, '-r')
               poly1 = np.polyfit(robs, diff, 3)
               y3 = poly1[0]*robs**3 + poly1[1]*robs**2 + poly1[2]*robs + poly1[3]
               ax1.plot(robs, y3, '-g')
               plt.show()
     
               #yn = robs - y2
               #plt.plot(robs, yn, '-r')
               yn = robs - y3
               plt.plot(robs, yn, '-g')
               plt.title('Difference between cdfs and fit')
               plt.show() 
           ##Plots forecasts transformed by cdf matching.
           if sp: 
               plt.plot(newra.obs, '-r')  #measurments are red.
               plt.plot(newra.fc*scale, '-b') #original forecasts in blue
               #newfc = newra.fc *scale
               #newfc2 = newfc - (poly[0]*newfc**2 + poly[1]*newfc + poly[2])
               #plt.plot(newfc2, '--g')  #new forecasts in green.
               plt.plot(scaled_obsra['fc'], '--g')  #new forecasts in green.
               plt.title('Green (forecast scaled), red (obs), blue (unscaled forecasts)')
               plt.show()
           ##This block plots the cdf's of the observed and forecast.
           if sp:
               x1, y1 = self.cdf(robs)
               x2, y2 = self.cdf(rfc)
               x3, y3 = self.cdf(fc)
               plt.step(x1, y1, '-r')
               plt.step(x2, y2, '-b')  #blue shows cdf of forecast
               plt.step(x3, y3, '--g') #green shows cdf of scaled forecast.
               plt.title('red(obs), blue(forecast), green(scaled forecast)')
               plt.show()
           ##This block plots relationship between .
           return scale, poly, scaled_obsra


    #def compare_simple(self, fff, sym='--g.', pstats=True, tave=-1):


    def compare(self, fff, sym='--g.', pstats=True, tave=-1, scalemethod=0):
        """create a plot for each station comparing
           forecast concentrations with observed concentrations.
           fff are the forecast concentrations.
 
           returns list of tuples with
           (station name,  correlation coefficient by station, scaling factor, cor coeff for all stations)
           the scaling factor is done by fitting obsra at all the stations,.
           the correlation coefficient is done by station.
        """
        sns.set()
        sns.set_context('talk')
        corrs = []
        c1 = []
        c2 = []
        iii=0

        obsdates = []
        fcdates = []
        fcra = []
        stra = []  #array of array of observations at each station.
        thash={} 
        for stndata in self.get_dates(fff):
            temp1 = list(zip(*stndata))
            dates = temp1[1]
            fff2 = temp1[2]
            stn_name = list(set(temp1[0]))[0]

            mdates, mdata = self.get_data2plot(stn_name)
            obs = pd.Series(mdata.tolist(), index=mdates)
            fc =  pd.Series(fff2, index=dates)
            fc  = fc.asfreq('1H')
            fc.fillna(value=0, inplace=True) 

            if tave > 0:
                ##first fill in missing measurements by forward filling.
                obs  = obs.asfreq('1H')
                obs.fillna(method='ffill', inplace=True) 
                freqstr = str(tave) + 'H' 
                obs = obs.resample(freqstr).mean()

                ##for the forecast data, nans should be zero.
                #fc  = fc.asfreq('1H')
                #fc.fillna(value=0, inplace=True) 
                freqstr = str(tave) + 'H' 
                fc = fc.resample(freqstr).mean()
            print('Scaling' + stn_name)
            stnobsra = self.init_obsra(obs, fc)
            scale, poly, scaled_stnobsra = self.scale_results(stnobsra, method=0)
            print(scale, poly)
            correlation = self.find_corr(scaled_stnobsra)  
            #correlation =  self.stats(obs, fc)
            c1.append(stn_name)
            #c2.append(correlation.at[0,1])
            c2.append(correlation)
            
            
            #obsdates.append(obs.index.tolist())
            #fcdates.append(fc.index.tolist())
            #fcra.append(fc)
            #stra.append(obs)

            obsdates.append(scaled_stnobsra.index.tolist())
            fcdates.append(scaled_stnobsra.index.tolist())
            fcra.append(scaled_stnobsra.fc)
            stra.append(scaled_stnobsra.obs)
            method=2

            if iii==0: 
               self.obsra = self.init_obsra(obs, fc)
            else:
               self.obsra = self.add_obsra(obs, fc, self.obsra)
            iii+=1
        ##scaling factor based on data at all measurement stations.
        scale, poly, scaled_obsra = self.scale_results(self.obsra, method=scalemethod) 
        print('all', scale, poly)
        allcorr = self.find_corr(self.obsra)      ##correlation for all stations.
        scaled_allcorr = self.find_corr(scaled_obsra)      ##correlation for all stations.
        shash = self.find_stats(self.obsra, scale=scale )
        shash2 = self.find_stats(scaled_obsra, scale=1)
        shash['scale'] = scale
        shash['poly'] = poly
        shash['corr'] = allcorr
        shash['scaledcorr'] = scaled_allcorr
        for ky in list(shash2.keys()):
            shash['scaled' + ky] = shash2[ky]  
        thash['all'] = shash

        station_hash={}
        jjj=0
        if pstats:
            for stn_name in c1:
                fig= plt.figure(self.fignum)
                #fig.set_size_inches(18,10)
                plt.plot(obsdates[jjj], stra[jjj]/1000.0, '-k', linewidth=2)
                plt.plot(fcdates[jjj], np.array(fcra[jjj]) / 1000.0, '-b.')
                #aaa = fcra[jjj] * scale
                #newfc = aaa - (poly[0]*aaa**2 + poly[1]*aaa + poly[2])
                #if scalemethod==3:
                #    newfc = aaa - (poly[0]*aaa**3 + poly[1]*aaa**2 + poly[2]*aaa + poly[1])
                #else:
                #    newfc = aaa 
                newfc = fcra[jjj]
                sobs = pd.Series(stra[jjj], index=obsdates[jjj])
                sfc =  pd.Series(fcra[jjj], index=fcdates[jjj])
                stnobsra = self.init_obsra(sobs, sfc)
                station_hash = self.find_stats(stnobsra, scale=scale)       #find statistics for each station. 
                corr = self.find_corr(stnobsra)      ##correlation for all stations.
                station_hash['corr'] = corr

                sobs = pd.Series(stra[jjj], index=obsdates[jjj])
                sfc =  pd.Series(newfc, index=fcdates[jjj])
                stnobsra = self.init_obsra(sobs, sfc)
                temp_hash = self.find_stats(stnobsra, scale=1)       #find statistics for each station. 

                for key in list(temp_hash.keys()):
                    station_hash['scaled'+key] = temp_hash[key] 
                corr = self.find_corr(stnobsra)      ##correlation for all stations.
                station_hash['scaledcorr'] = corr
 

                thash[stn_name] = station_hash
  
                #plt.plot(fcdates[jjj], newfc, '-r')
                ax = plt.gca()
                ax.set_ylabel('Concentration (mg m$^{-3}$)')
                #ax.set_ylim([1,30000])
                #ax.set_ylim([1,30000])
                #ax.set_yscale("log")
                #plt.title(stn_name)
                self.fignum+=1
                #plt.savefig(stn_name + 'WE.jpg')
                plt.savefig(stn_name + 'msh06.jpg')
                #plt.show()
                jjj+=1
        return  thash  #dictionary key is station name or 'all', value is dictionary with key - statistic description, value statistic.
 
    def evsuhexbin(self, emissions, verbose=False,sym='k*', var='ustar'):
        print('plotting hexbin figure ', self.fignum)
        if var=='ustar':
           ustar = self.get_ustar()
        elif var=='wind':
           ustar = self.get_wind()
        vvv = np.where(emissions <=0)
        #emissions[vvv] = np.min(emissions)/10.0
        emissions[vvv] = 1e-2
        logemissions = np.log(emissions)
        fig = plt.figure(self.fignum)
        cb = sns.jointplot(x=ustar, y=logemissions, kind="hex")
        #plt.colorbar()
        self.fignum+=1
           

    def get_data2plot(self, stn):
        """This returns observed data for one station.
           Will include points not in the TCM for
           which model output forecast 0. Returns list of dates and values."""
        data = self.data.reset_index()
        #print data
        data = data[data['station'].isin([stn])]
        return data['date'], data['mdata']


    def get_dates(self, fff):
        """input array, fff with same length as number of rows in tcm.
           Return list of tuples with station, date, fff value"""
        dates = self.tcm.index.tolist()
        temp1 = list(zip(*dates))
        stations = set(temp1[0])
        temp2 = list(zip(temp1[0], temp1[1], fff))  #list of tuples with station, date, value
        #print temp2
        #print '**************************'
        for stn in stations:
            vals = [q for q in temp2 if q[0]==stn]
            #print stn, vals
            #print '**************************'
            yield vals

 
    def tcsolve(self, fname='tcsolve.txt', verbose=False):
        """reads the output of tcsolve. 
           converts emission to concentration.


        The TCM matrix has the following format
           | S1 | S2 | S3 | S4 |...| SN |
        --------------------------------
        M1 | A11 A12   A13  A14 ... A1N 
        M2 | A21 A22   A23  A24 ... A2N 
        M3 | A31 A32   A33  A34 ... A3N
        M4 | A41 A42   A43  A44 ... A4N
        .. | 
        MJ | AJ1 AJ2   AJ3  AJ4 ....AJN


        There is a row for each measurement location / time (total of J rows)
        There is a column for each source location /time (total of N columns).
        Values in the matrix show unit mass each soure contributed to each measurement.
      
        If we have a vector M representing the  measurements and another vector E representing
        the emission from each source, while A is the TCM matrix shown above then we have
         A * E = M

        M1 = E1*A11  + E2*A12 + E3*A13 + ... EN*A1N
        M2 = E2*A21  + E2*A22 + E2*A23 + ... EN*A2N
        
        MJ = E2*AJ1  + E2*AJ2 + E2*AJ3 + ... EN*AJN

        If we know the vector M and matrix A, then we can solve for E. 
        Because the system is usually over or under-determined and the measurements and model
        contain errors, there is usually no exact solution.

        The output of the solver program gives a list which represents the vector, E. 
        To get the model forecast values we multiply this vector by each row in the TCM.

        """
        emissions=[]
        with open(fname, 'r') as fid:
             fid.readline()
             for line in fid:
                 temp = line.split(',')
                 emissions.append(float(temp[1])) 
        self.emissions=np.array(emissions)
        mmm = self.tcm.fillna(value=0)
        mmm = mmm.as_matrix()
        fff = []
        nrow = mmm.shape[0]
        ncol = mmm.shape[1]
        #print 'emissions', self.emissions.shape, mmm[1].shape
        ##go through each row and multiply by emission vector and sum to get concentration
        ##at that measurement source/time.
        for iii in range(0, nrow):
            forecast = (mmm[iii]*self.emissions).sum()
            fff.append(forecast)
            nnn=0
            if verbose:
                print('solver----------------')
                for mval in mmm[iii]:
                    if mval!=0: print(nnn, mval, self.emissions[nnn], mval* self.emissions[nnn])
                    nnn+=1
            if verbose: print('Forecast, data: ', forecast, self.data['mdata'][iii])
        #print self.data.shape
        self.compare(fff)


    def tcmprint(self, fname='tcm.txt', data=[-99]):
        print('IN TCMPRINT')
        mmm = self.aug.fillna(value=0)
        mmm = mmm.as_matrix()
        thash = {}
        ndx = self.aug.columns
        #print 'COLUMNS'
        #for nd in ndx:
        #    print nd
        #print '--------------------'
        #print 'SHAPE', mmm.shape
        with open(fname, 'w') as fid:
            zzz=0
            for ciii in range(0,mmm.shape[1]-1):
                nd = ', '
                #if ciii == mmm.shape[1]-2: nd = ''
                tstr = 'T:' + str(zzz)
                fid.write(tstr + nd)
                thash[tstr] = (ndx[ciii][0], ndx[ciii][1])
                zzz+=1
             
            fid.write('  M:0 \n') 
            
            for riii in range(0,mmm.shape[0]):
                for ciii in range(0,mmm.shape[1]):
                    nd = ', '
                    if ciii == mmm.shape[1]-1: nd = ''
                    if mmm[riii,ciii]==0:
                        fid.write(str(int(mmm[riii,ciii]))+ nd)
                    else:
                        fid.write("{0:.3e}".format(mmm[riii,ciii])+ nd)
                fid.write('\n') 
        for key in list(thash.keys()):
            #print key, thash[key]
            try:
               temp =  self.uhash[thash[key]]
            except:
               print('no key for uhash')
        #print '----------------------------------------'
        #print self.uhash.keys()
        self.thash = thash
        print('OUT TCMPRINT')



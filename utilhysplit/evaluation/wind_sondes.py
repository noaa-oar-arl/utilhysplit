import datetime
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from utilhysplit.evaluation import reliability

# 01/05/2022
# read files for wind sonde comparison from Binyu.

def test1(wmax):
    ddd = 5
    wmin = 320
    wlist = np.arange(wmin,wmax,ddd)
    wlist = [x if x<=360 else x-360 for x in wlist]
    return wlist


def readfile(fname):
    dstr = "%Y%m%d_%H%M00"
    df = pd.read_csv(fname, sep='\s+',parse_dates=['DATE'])
    df['date'] = df.apply(lambda row: datetime.datetime.strptime(row['DATE'],dstr),axis=1)
    return df 

def ccw_diff(start,stop):
    """
    counter clockwise difference
    start : starting angle
    stop : stopping angle
    Returns :
        how many degrees need to move counter-clockwise to get from small to large.
    """
    if start > stop:
       diff = start-stop
    else:
       diff = start + 360 - stop
    return diff

def cw_diff(start, stop):
    """
    clockwise difference
    start : starting angle
    stop : stopping angle
    Returns :
        how many degrees need to move clocwise to get from small to large.
    """
    dres = 0.001
    if stop >= start:
       diff = stop - start
    else:
       diff = stop + 360-start
    if np.abs(diff-360) < dres:
       diff = 0
    return diff    

def diff_wdir(wdir1, wdir2):
    if wdir1 > wdir2:
       big = wdir1
       small = wdir2
    else:
       big = wdir2
       small = wdir1
    diff1 = big - small
    diff2 = small + (360-big)
    if diff1 > diff2: 
       diff = diff2
       # move counter clockwise from smallest to largest.
       direction = 'ccw' #counter clockwise
    else: 
       diff = diff1
       # move clockwise from smallest to largest.
       direction = 'cw'  #clockwise
    return direction, diff

def spread(wlistin,verbose=False):
    """
    wlist: list of wind directions.
    order from least to greatest (0 to 360)
    Take clock-wise difference between each subsequent one.
    Find largest clock-wise gap.
    Set start degree to direction that occurs after this gap.
    Set end degree to direction that occurs before the gap.
    Return:
        start degree
        end degree
        spread 
    """
    wlist = wlistin.copy()
    if isinstance(wlist, pd.Series):
       wlist = wlist.values
    wlist.sort()
    dlist = []
    dmax = 0
    imax = 0
    if verbose: print(wlist)
    for iii, www in enumerate(wlist[1:]):
        diff = cw_diff(wlist[iii],www)
        if verbose: print(www,wlist[iii],diff)
        dlist.append(diff)
        if diff > dmax:
           dmax = diff
           imax = iii
    # check diff between last and first.
    diff = cw_diff(wlist[-1],wlist[0])
    if verbose: print(wlist[-1],wlist[0], diff)
    if verbose: print('current', dmax, imax)
    dlist.append(diff)
    if diff > dmax:
       dmax = diff
       imax = len(wlist)-1
       imin = 0
       if verbose: print('new', dmax, imax)
    else:
       imin = imax + 1
    spread = cw_diff(wlist[imin], wlist[imax])
    return wlist[imin], wlist[imax],  spread

def process_files(flist, addobs=None,wspd=None):
    """
    flist : list of files
    addobs : str. 'ODIR' or 'O_SPEED'
    wspd : float. If not None then removes rows where observed wind speed is below this value.
    RETURNS
    dfall : pandas dataframe.
    """
    dflist = []

    for fname in flist:
        df = readfile(fname)
        ens = fname.split('.')[0]
        if '/' in ens:
            ens  = ens.split('/')[-1]
        df['ens'] = ens
        #df = df[['date','FDIR','ens']].copy()
        dflist.append(df)

    dfall = pd.concat(dflist)
    if wspd:
       dfall = dfall[dfall['O_SPEED']>wspd]

    # creates dataframe with only three columns, date, speed or direction, ensid.
    if isinstance(addobs,str):
        if not addobs in ['ODIR','O_SPEED']:
           print('addobs not recognized {}'.format(addobs))
        else:   
           dfnew = readfile(flist[0])
           obs = dfnew[['date',addobs]].copy()
           obs = obs.drop_duplicates()
           obs['ens'] = 'obs'

           if addobs=='ODIR': new='FDIR'
           elif addobs=='O_SPEED': new = 'F_SPEED'

           obs.columns = ['date',new,'ens']
        dfall = dfall[['date',new,'ens']]
        dfall = pd.concat([dfall,obs])

    return dfall

def get_wspd_df(flist, addobs=True,addspread=False,wspd=None):
    if addobs: aaa = 'O_SPEED'
    else: aaa=None
    dfall = process_files(flist,addobs=aaa,wspd=wspd)
    dfwdir = pd.pivot(dfall, columns = 'ens', values='F_SPEED', index='date')
    return dfwdir

def get_wdir_df(flist, addobs=True,addspread=True,wspd=None):
    """
    flist : list of files
    addobs : boolean
    addspread : boolean
    wspd : float. If not None then removes rows where observed wind speed is below this value.
    RETURNS
    dfwdir : pandas dataframe.
             index is date
             values are wind direction.
             columns are ensemble member, obs
             if addspread then last three columns give starting angle, ending angle and spread.  
    """
    if addobs: aaa = 'ODIR'
    else: aaa = None
    dfall = process_files(flist,addobs=aaa,wspd=wspd)
    dfwdir = pd.pivot(dfall, columns = 'ens', values='FDIR', index='date')
    if addspread:
        dfwdir = add_spread(dfwdir)
    return dfwdir

def add_spread(dfwdirin):
    dfwdir = dfwdirin.copy()
    dftemp = dfwdir.apply(lambda row: spread(row), axis=1)
    dftemp = dftemp.to_frame().reset_index()
    dftemp.columns = ['date','temp']
    dftemp['spread'] = dftemp.apply(lambda row: row['temp'][2],axis=1)
    dftemp['start']  = dftemp.apply(lambda row: row['temp'][0],axis=1)
    dftemp['stop']  = dftemp.apply(lambda row: row['temp'][1],axis=1)
    dftemp = dftemp[['date','spread','start','stop']]
    dfwdir = pd.concat([dfwdirin, dftemp.set_index('date')], axis=1)
    return dfwdir


def recalc_wdir(dfwdir):
    """
    dfwdir : pandas dataframe with wdir for obs, ensemble members.
             Also needs columns start, stop and spread.
    wdir values are recalculated as degrees clock-wise from the starting value.
    """
    df = dfwdir.copy()
    cnames = dfwdir.columns
    cnames = [x for x in cnames if x not in ['spread','start','stop']]

    def calc(xval,start):
        # calculate clockwise distance from start. 
        newval = cw_diff(start,xval) 
        return newval

    for ens in cnames:
        df[ens] = df.apply(lambda row: calc(row[ens],row['start']),axis=1)
    return df 


class Comparison:

    def __init__(self,tdir,tag='A',plist =['P200','P300','P500']):
        """
        tdir : directory where files are located
        tag  : site identifier
        plist : filenames should be ensid.com.p*.txt where p* indicates the pressure                level. plist gives the pressure levels to use.
        """
        self.tag = tag
        self.plist = plist
        self.filehash = {}
        self.dfhash = {}
        self.dfhash2 = {}
        self.rankhash = {}
        for pressure in self.plist:
            fff = '/*.com.{}.txt'.format(pressure)
            self.filehash[pressure] = glob.glob(tdir + fff)    
            self.dfhash[pressure] = get_wdir_df(self.filehash[pressure],
                                             addobs=True,
                                             addspread=True) 
            self.dfhash2[pressure] = get_wdir_df(self.filehash[pressure],
                                             addobs=False,
                                             addspread=True) 


    def _return_plist(self,pressure=None):
        if not isinstance(pressure,(list,str)):
           pressure = self.plist
        elif isinstance(pressure,str):
           pressure = [pressure]
        return pressure

    def plotspread(self,pressure=None):
        pressure = self._return_plist(pressure)
        for prss in pressure:
            sns.set_style('whitegrid')
            if not prss in self.dfhash.keys():
               print('pressure not found {}'.format(prss)) 
               continue
            plt.plot(self.dfhash2[prss].spread, '--b.', linewidth=0.2, label='with obs')
            plt.plot(self.dfhash[prss].spread, 'r+', markersize=4, label='No obs')
            plt.title(self.tag + ' ' + prss)
            plt.xticks(rotation=45)
            plt.show()

    def calc_rank_diagram(self,pressure=None):
        pressure = self._return_plist(pressure)
        # TODO adjust number of axis depending on length of pressure list.
        fig = plt.figure(1,figsize=[12,2.5])
        ax1 = fig.add_subplot(1,3,1)
        ax2 = fig.add_subplot(1,3,2)
        ax3 = fig.add_subplot(1,3,3)
        axlist = [ax1,ax2,ax3]
        for iii, prss in enumerate(pressure):
            # recalculates wind direction with reference to starting angle.
            new = recalc_wdir(self.dfhash[prss])
            new = new.drop(['spread','start','stop'],axis=1)
            if prss not in self.rankhash.keys():
                rank = reliability.Talagrand(thresh=0,nbins=32)
                rank.add_data(new)
                self.rankhash[prss] = rank
            print(prss)
            self.rankhash[prss].plotrank(ax=axlist[iii])
        plt.sca(ax2)
        plt.title(self.tag)
        return fig, axlist

            


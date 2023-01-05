import datetime
import pandas as pd
import numpy as np

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

def spread(wlist,verbose=False):
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

def process_files(flist, addobs=True):
    dflist = []
    for fname in flist:
        df = readfile(fname)
        ens = fname.split('.')[0]
        if '/' in ens:
            ens  = ens.split('/')[-1]
        df['ens'] = ens
        df = df[['date','FDIR','ens']].copy()
        dflist.append(df)
    if addobs:
        dfnew = readfile(flist[0])
        obs = dfnew[['date','ODIR']].copy()
        obs = obs.drop_duplicates()
        obs['ens'] = 'obs'
        obs.columns = ['date','FDIR','ens']
        dflist.append(obs)
    dfall = pd.concat(dflist)
    dfwdir = pd.pivot(dfall, columns = 'ens', values='FDIR', index='date')
    dfwdir = add_spread(dfwdir)
    return dfwdir

def add_spread(dfwdir):
    dftemp = dfwdir.apply(lambda row: spread(row), axis=1)
    dftemp = dftemp.to_frame().reset_index()
    dftemp.columns = ['date','temp']
    dftemp['spread'] = dftemp.apply(lambda row: row['temp'][0],axis=1)
    dftemp['start']  = dftemp.apply(lambda row: row['temp'][1],axis=1)
    dftemp['stop']  = dftemp.apply(lambda row: row['temp'][2],axis=1)
    dftemp = dftemp[['date','spread','start','stop']]
    dfwdir = pd.concat([dfwdir, dftemp.set_index('date')], axis=1)
    return dfwdir

class Comparison:

    def __init__(self):
        return -1

    def readfile(self,fname):
        df = pd.read_csv(fname, sep='\s+',parse_dates=['DATE'])
        return df 

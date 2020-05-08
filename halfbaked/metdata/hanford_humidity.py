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
import seaborn as sns

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
    alist = []
    daylist = []
    jjj=0
    minlist = []
    maxlist = []
    while not done:
       ix = (ds.index > d1) & (ds.index < d1+dt) 
       dsub = ds[ix]
       rrr = dsub.max() - dsub.min()
       #t1 = dsub[dsub == dsub.max()].index
       #t2 = dsub[dsub == dsub.min()].index
       #rdt1 = t1.max() - t2.min() 
       #rdt2 = t1.min() - t2.max() 
       #print(dsub)
       if len(dsub) > 1 : (rlist.append(rrr))
       #print(dsub)
       #print(dsub.mean())
       alist.append(dsub.mean())
       daylist.append(datetime.datetime(d1.year, d1.month, d1.day, 0))
       rlist.append(rrr)
       d1 = d1 + dt
       if d1 > de: done=True       
       jjj+=1
       if jjj> 400: done=True
    pa = pd.Series(alist, index=daylist)
    return rlist, pa

def calc_relh(x):
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
        obs = mdata.add_data(self.dates, country=None, box=self.area, resample=False)
        obs['relh'] = obs.apply(calc_relh, axis=1) 
        obs['latlon'] = str(obs['latitude']) + ' ' + str(obs['longitude'])
        #print(self.obs['latitude'].unique())
        #print(self.obs['longitude'].unique())
        #print(self.obs.columns.values)
        #rplot(self.obs)
        return obs.copy()

logname = 'hanford_relh_log.txt'

area = [46.3,-120,46.9, -119]
year = 1980
endyear = 2010
done = False
iii=0
verbose=False
rangelist = []
while not done:
    area = [46.3,-120,46.9, -119]
    d1 = datetime.datetime(year,1,1,0)
    d2 = datetime.datetime(year,12,31,23)
    sv = Mverify([d1,d2], area)
    obsb = sv.find_obs()
    obsb.set_index('time', inplace=True)
    obsb = obsb[obsb['station name'].isin(['HANFORD', 'HANFORD AIRPORT'])]
    print(obsb.columns.values)
    stationid = obsb['station_id'].unique()
    stationid2 = obsb['station name'].unique()
    lat = obsb['latitude'].unique()
    lon = obsb['longitude'].unique()
    rstr=''
    with open(logname, 'a') as fid:
         for si in stationid:
             rstr = str(year) + ' ' + str(si) + ' '
         for si in stationid2:
             rstr += str(si) + ' '
         for si in lat:
             rstr += str(si) + ' '
         for si in lon:
             rstr += str(si) + ' '
         rstr += '\n'
         fid.write(rstr)
 
    relh = obsb['relh']
    if verbose: print('relh-----')
    if verbose: print(relh[-10:])
    rlist, pa = find_range(relh)
    rangelist.extend(rlist)
    if verbose: print('PA-----')
    if verbose: print(pa[-10:])
    pw = pa.resample('W').mean()
    pm = pa.resample('M').mean()
    if verbose: print(pw[-10:])
    ##pivot table with column as the year and rows being indvidual months
    pmf = pm.reset_index()
    pmf['month'] = pmf['index'].apply(lambda x: x.month)
    pmf['year'] = pmf['index'].apply(lambda x: x.year)
    pmf = pd.pivot_table(pmf, values=0, index=['month'], columns=['year'])
    if verbose: print(pmf)

    ##pivot table with column as the year and rows being indvidual weeks
    pwf = pw.reset_index()
    pwf['week'] = pwf['index'].apply(lambda x: x.week)
    pwf['year'] = pwf['index'].apply(lambda x: x.year)
    pwf = pwf[pwf['year'] == year]  #remove weeks that end in the next year
    pwf = pd.pivot_table(pwf, values=0, index=['week'], columns=['year'])
    if verbose: print(pwf)
    

    if iii==0:
        relhmonth = pmf
        relhweek = pwf
    else:
        relhmonth = pd.concat([relhmonth,pmf], axis=1)
        relhweek = pd.concat([relhweek,pwf], axis=1)
    year = year + 1
    iii += 1
    if year > endyear: done=True
print('-------------------------HERE')
print(relhmonth)
sns.set()
#plt.plot(relhmonth)
relhmonth.plot(legend=False)
ax = plt.gca()
plt.ylabel('Monthly Average Relative Humidity')
plt.savefig('relh_monthly.jpg')
plt.show()

relhweek.plot(legend=False)
ax = plt.gca()
plt.ylabel('Weekly Average Relative Humidity')
plt.savefig('relh_weekly.jpg')
plt.show()


mname = 'monthly_relh.csv'
wname = 'weekly_relh.csv'

ff= '%5.1f'
relhmonth.to_csv(mname, float_format=ff, header=True)
relhweek.to_csv(wname, float_format=ff, header=True)

#pal = sns.palplot(sns.light_palette("navy", as_cmap=True))
pal = sns.diverging_palette(147,280,s=85,l=25, n=7, as_cmap=True)

sns.heatmap(relhweek.transpose(), center=50, cmap=pal)
plt.savefig('relh_weekly_heatmap.jpg')
plt.show()

sns.heatmap(relhmonth.transpose(), center=50, cmap=pal)
plt.savefig('relh_monthly_heatmap.jpg')
plt.show()

rangelist = np.array(rangelist)
rangelist = rangelist[np.logical_not(np.isnan(rangelist))]
sns.distplot(rangelist)
plt.savefig('relh_daily_differences.jpg')
plt.show()
print(rangelist.mean())





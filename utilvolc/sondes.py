# sonde data for Turrialba from Gary Morris
import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import glob
from scipy.signal import find_peaks
from utilhysplit.evaluation import statmain



def get_peaks(olist,plot=True):
    vlist = olist.dropna()
    rra = vlist.values.flatten()
    print(type(rra),rra.size,rra.shape)
    height = 5
    distance = 50
    prominence = 5
    width = [1,500]
    peaks = find_peaks(rra,height=height, distance=distance,prominence=prominence,width=width)
    #vlist = vlist.reset_index()
    tsp = vlist.iloc[peaks[0]]
    if plot:
       plt.plot(olist.values, olist.index.tolist(),'--k',linewidth=1)
       plt.plot(tsp.values, tsp.index.tolist(), 'ro',MarkerSize=10,alpha=0.5)
    return peaks
    




def get_peaksB(vlist, pvalues=[0.99,1], plot=True):
    """
    Creates a CDF from the values in the input series and
    returns a series with same length in which only the values with
    CDF between pvalue[0] and pvalue[1] are not NaN.
    vlist is a pandas series.
    Returns
    tsp: pandas series
         all values which are below the
 
    """
    olist = vlist.copy()
    valA = statmain.probof(vlist.values,pvalues[0])
    valB = statmain.probof(vlist.values,pvalues[1])
    print(valA, valB)
    tsp = vlist[vlist>=valA]
    tsp = tsp[tsp<=valB]
    if plot:
       plt.plot(vlist.index.tolist(),vlist.values,'--k',linewidth=1)
       plt.plot(tsp.index.tolist(), tsp.values,'b.',MarkerSize=5)
    return tsp 





def get_sondes(tdir):
    fnames =  glob.glob(tdir + '/*dat')
    return fnames


class SondeList:

    def __init__(self,name='Sondes'):
        self.name = name
        self.slist = []  # List of Sonde objects
        self.tlist = []  # List of launch dates
        self.totalSO2 = [] # list of integrated SO2 up to 15 km
        self.totalSO2_label =  ''


    def get_peaks(self, plot=True):
        for sss in self.slist:
            sss.get_peaks(plot=plot)
            if plot: 
               plt.show() 
               
                      


    def plot_totalSO2(self):
        d1 = self.tlist[0]
        d2 = self.tlist[-1]
        fig = plt.figure(figsize=[20,5])
        ax = fig.add_subplot(1,1,1)
        plt.plot(self.tlist,self.totalSO2, '--k.')
        plt.plot([d1,d2],[0,0],'--r')
        ax.set_ylabel(self.totalSO2_label)
 
    def add_data(self,sonde):
        self.slist.append(sonde)
        self.tlist.append(sonde.launchdate)
        ky=None
        for key in sonde.ihash.keys():
            if 'integrated so2' in key.lower(): ky = key
        if ky: 
           try: stotal = float(sonde.ihash[ky])
           except: stotal = -10
           self.totalSO2.append(stotal)
           self.totalSO2_label = ky 

class Sonde:
    def __init__(self,fname):
        self.fname = fname
        self.ihash = {}
        self.units = []    # list of strings with units for each column.
        self.df = pd.DataFrame() # dataframe
        self.columns = ['Time','Press','Alt','Temp','RH','O3a','O3b','O3c','wdir','wspd','Tpump','I O3','SO2','GPS Lon','GPS Lat', 'GPS Alt']
        self.launchdate = None
        self.latitude = None
        self.longitude = None
        iii = self.read_header()
        try:
            self.read_data(iii+1)
        except:
            print('Warning could not read data', self.fname)
        self.read_header()
      

    def get_peaks(self,val='SO2',plot=True):
        stemp = self.df
        stemp = stemp[['Alt',val]]
        stemp = stemp.set_index('Alt')
        self.peaks = get_peaks(stemp, plot=plot)
        if plot:
           ax = plt.gca()
           ax.set_ylabel('Altitude (km a.s.l.)')
           unit = 'ppbv'           
           ax.set_xlabel('{} ({})'.format(val, unit))

           tropo = float(self.ihash['Tropopause height (km)'])
           xi = np.min(stemp.dropna().values)
           xf = np.max(stemp.dropna().values)
           print('Tropopause at', tropo)
           if tropo < 15:
               ax.plot([xi,xf],[tropo,tropo],'--b',label='Tropopause')
           mlh = float(self.ihash['Mixed layer height (km)'])
           print('{} {} MLH {}'.format(xi,xf, mlh))
           ax.plot([xi,xf],[mlh,mlh],'-y',label='MLH', linewidth=5, alpha=0.5)
           ax.plot([0,0],[0,15],'-g',linewidth=5,alpha=0.5)
            

 
    def read_header(self):
        ihash = {}
        with open(self.fname,'r') as fid:
             iii=0
             for line in fid.readlines():
                 temp = line.split(':')
                 if 'Time' in line and 'Press' in line and 'Alt' in line: break  
                 ihash[temp[0].strip()] = temp[-1].strip()
                 if 'latitude' in temp[0].lower(): self.latitude = temp[1]
                 if 'longitude' in temp[0].lower(): self.longitude = temp[1]
                 if 'Launch Date' in temp[0]:
                     dstr = temp[1].strip()
                     self.launchdate = datetime.datetime.strptime(dstr, "%Y%m%d") 
                 iii+=1
        self.ihash = ihash
        return iii 

    def read_data(self,skr):
        fname = self.fname
        df = pd.read_csv(fname, header=0, skiprows=skr, sep = '\s+')
        units = df.columns
        units = list(zip(units,self.columns))
        if len(df.columns)==len(self.columns):
            df.columns = self.columns
        else:
            iii = len(df.columns) 
            df.columns = self.columns[0:iii]
            #print('could not replace columns')
            #print(len(df.columns), df.columns)
        # replace missing or bad values with NaN
        df = df.where(df!=9000,np.NaN)
        self.df = df

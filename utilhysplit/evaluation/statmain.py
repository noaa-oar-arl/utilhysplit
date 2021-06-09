# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import sys 
import numpy as np
#import matplotlib.pyplot as plt
import datetime
#from subprocess import call
#from os import path
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import subprocess
import seaborn as sns

"""
NAME: statmain.py
PGRMMR: Alice Crawford ORG: ARL
This code written at the NOAA Air Resources Laboratory

Classes
MatchedData
     methods: tocsv, load, create_obsra, combine, add_obsra, resample_obsra, find_corr, plotseries,
              apply_thresh, find_stats


Functions
find_stats_function

"""
def autocorr(ts1, nlist=None):
    # plots the autocorrelation of both the series.
    alist = [] 
    if not np.any(nlist): 
       nlist = np.arange(0,48)
    for nnn in nlist:
       alist.append(ts1.autocorr(lag=nnn)) 
    #ax.plot(nlist, alist1, 'k.', label='obs')
    #ax.plot(nlist, alist2, 'b.', label='fc')
    return alist 

def get_pixel_matching_threshold(obsra, modelra, threshold=0):
    obsx = np.sort([x for x in obsra.values.flatten() if x > threshold])
    modelx = np.sort(modelra.values.flatten())
    numobs = len(obsx)
    if numobs < len(modelx):
       modelx = modelx[len(modelx)-numobs:]
    else:
       print('Warning Pixel Matching: More observed than modeled values above threshold')
    return modelx[0] 


def degdiff(cc, mm):
    # given two degree measurements, find smallest difference.
    # if moving clockwise gives you the smallest angular difference
    # then difference is positive. Otherwise negative.
    wdiff = (np.min([np.abs(cc-mm), 360-np.abs(cc-mm)]))
    if np.abs(cc-mm) > 360-np.abs(cc-mm):
       wdiff = -1 * wdiff
    return wdiff

def stepfunction(xstep, ystep, iii):
    """
    f(xstep) = ystep
    where xtep, ystep define a step function 
    (for instance output of cdf)
    input iii  (float)
    returns f(iii) (float)
    """
    # this way may be faster.
    zzz=1
    xxx = np.array(xstep)
    yyy = np.array(ystep)
    vpi = np.where(xxx <= iii)
    try:
        zzz = vpi[-1][-1]
    except:
        #print('stepfunction', vpi)
        #print(iii, xstep[0], xstep[-1])
        return yyy[0]
    try:
        mmm = yyy[zzz+1]
    except:
        return yyy[-1]
    return  mmm

    # loop is slow
   
    # zzz = zip(xstep, ystep)
    # jjj=0
    # for val in zzz:
       # find value of x which is closest to input value iii.
    #    if iii >= val[0]: jjj=val[1]
    #    if iii < val[0]: break
    # test to see if they return the same thing.
    #if jjj != mmm:
    #   print('jjj, mmm', jjj, mmm) 
    #return  jjj 


def get_ds(n1, n2):
    # see page 154 of Wilks
    # Note that for the time series data 
    # the measurements are not independent.
    # N1 and N2 should be number of independent measurments

    # Supposedly ds should be increasing with alpha?
    # so that a large value would be needed to reject null
    # hypothesis at 100% than 50%? 
    # but this function seems to be decreasing.
    alpha = np.arange(1,100,1) / 100.0
    ds = 0.5 * (1./n1 + 1./n2) * np.log(alpha/2)
    ds = (-1 * ds) **0.5
    return alpha, ds

def example_ks_test(numpts=100):
    from utilhysplit.evaluation import statmain
    # create two Gaussian distributions which are not the same.
    val1 = np.random.normal(80,10,numpts)
    val2 = np.random.normal(90,10,numpts)
    # get CDF of distributions
    x1,y1 = nancdf(val1,thresh=0) 
    x2,y2 = nancdf(val2,thresh=0) 
    # plot CDFS
    fig = plt.figure()
    ax = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)
    ax.step(x1,y1)
    ax.step(x2,y2)
   
    # plot difference between CDF's 
    xr = kstest_sub(x1,y1,x2,y2)
    ax2.plot(xr[0],xr[1], '-k.')
    ax2.set_ylabel('Difference in CDFs')
    ax.set_ylabel('Probability')
    print('KSP', kstestnan(val1, val2))
    print('KSP', kstest(val1, val2))

def kstestnan(data1, data2,thresh=None):
    print('start kstestnan')
    cx1, cy1 = nancdf(data1,thresh)  
    cx2, cy2 = nancdf(data2,thresh) 
    result =  kstest_sub(cx1,cy1,cx2,cy2)
    return np.max(np.abs(result[1]))

def kstest(data1, data2):
    """
    Kolmogorov-Smirnov test
    data1 : list of data
    daata2 :  list of data

    creates cdf for data1 and data2 and finds different
    between them assuming each are a step function.
   
    """
    # see page 154 of Wilks
    # two sample K-S test. Compare two batches of data to one another
    # under the null hypotheses that they were drawn from the same (but
    # unspecified) distribution. 

    # Ds = max|FN(x1) - Fm(x2|

    # Null hypothesis at the a *100% level is rejected if

    # DS > (-0.5 * (1/n1 + 1/n2) ln(a/2) ^ 0.5

    # n1 = size of sample1
    # n2 = size of sample2

    cx1, cy1 = cdf(data1)  
    cx2, cy2 = cdf(data2) 
    result =  kstest_sub(cx1,cy1,cx2,cy2)
    return np.max(np.abs(result[1]))

def kstest_sub(cx1,cy1,cx2,cy2):
    """
    cx1 : 1D numpy ndarray
    cy1 : 1D numpy ndarray
    cx2 : 1D numpy ndarray
    cy2 : 1D numpy ndarray

    cx1, cy1 can be output from cdf function.
    cx2, cy2 can be output from cdf function.

    Returns 
    [xtot, difflist]
    xtot - value on x axis of CDFS
    difflist - difference between two cdfs

    xtot and difflist have the same length
    """
    difflist = []
    # this is a little slow for functions with
    # a lot of x values. 
    xtot = np.sort(np.append(cx1,cx2))
    for xxx in xtot:
        val1 = stepfunction(cx1, cy1, xxx)
        val2 = stepfunction(cx2, cy2, xxx)
        difflist.append(val2 - val1) 
    return [xtot,difflist]

def probof(data1, probval):
    """
    input 
    data1: list of data.
    probval: number from 0 to 1.
    creates cdf and returns value which has
    cumulative probability equal to probval.
    """
    cx1, cy1 = cdf(data1)  
    return stepfunction(cy1, cx1, probval)

def kstest_answer(data1, data2):
    d1, d2 = kstest(data1, data2)
    return np.max([np.max(d1), np.max(d2)])

def pixel_matched_cdf(modelra,numobs,threshold=0):
    """
    modelra : 1d numpy array of values.
    nans will be removed.
    numobs : number of observations. The modelra will be
             sorted and lower values removed so that length
             matches number of observations.
    """
    modelra = modelra[~np.isnan(modelra)]
    modelx = np.sort(modelra)
    # remove smaller values so array is same size as obs.
    if numobs < len(modelx):
       modelx = modelx[len(modelx)-numobs:]
    else:
       print('WARNING pixel matching. More obs values than modeled')
       print(len(modelx), numobs)
    modely = 1. * np.arange(modelx.size) / float(modelx.size-1)
    return modelx, modely

def nancdf(data,thresh=None):
    """
    remove nans from data before creating cdf.
    """
    data2 = data[~np.isnan(data)]
    if thresh:
       data2 = data2[data2>thresh]
       #vpi = data2 < thresh
       #data2[vpi] = np.nan
       #data2 = data2[~np.isnan(data2)]
    return cdf(data2)

def cdf(data):
    """
    data should be a list of data.
    sdata : sorted list of the data.
    yval  : normalized index of data (from 0 to 1). 

    """
    sdata = np.sort(data)
    # y axis goes from 0 to 1
    yval = 1. * np.arange(sdata.size) / float(sdata.size-1)
    return sdata, yval

def plot_cdf(sdata, y, ax):
    ax.step(sdata, y, '-r')

class EnsembleTimeSeries:

    def __init__(self, obs=pd.Series(), fc=pd.Series(), stn=None):
        self.df = pd.DaataFrame()
        if not obs.empty and not fc.empty:
           self.df = self.create_df(obs,fc)

    def create_df(self,obs,fc):
        return -1

######following methods create and analyze an obsra. Which is a pandas dataframe with matched observations and forecasts by date.
class MatchedData(object):

    #def create_obsra(obs, fc):
    def __init__(self, obs=pd.Series(), fc=pd.Series(), stn=None):
        ##input series with dates and times of observations and forecasts.
        ##put them together into a pandas data frame matched by dates.
        ## remove points for which either and observation or a forecast is not available.
        #self.obsra = self.create_obsra(obs, fc)
        self.stn = stn
        self.obsra = pd.DataFrame()
       
        if not obs.empty and not fc.empty:
           self.obsra = self.create_obsra(obs, fc) 

    def tocsv(self, fname='./mdata.csv', mult=1, verbose=False):
        """
        writes dataframe to csv file
        """
        dfmt = '%Y-%m-%d-%H'
        ff = '%10.4e'
        if verbose: print( 'PRINTING CSV ', fname)
        obsra = self.obsra.copy()
        obsra['fc'] = obsra['fc'] * mult
        obsra.to_csv(fname, float_format= ff, date_format=dfmt, header=True)
 
    def load(self, fname='./mdata.csv', inplace=False):
        """
        loads dataframe from csv file
        """
        obsra = pd.read_csv(fname, index_col=[0], parse_dates=True)
        
        if inplace: 
           if self.obsra.empty: 
              self.obsra = obsra
           else:
              self.obsra = pdf.concat([self.obsra, obsra]) 
        return obsra 

    def create_obsra(self, obs, fc):
        """
        creates dataframe from two series.
        """
        #data = 
        data = pd.concat([obs,fc], axis=1)
        data.columns = ['obs','fc'] 
        #data =  pd.DataFrame(data={'obs': obs, 'fc': fc})
        data.dropna(axis=0, inplace=True)  #remove nan values
        return data

    def combine(self, obsra, inplace=True):
        """
        concatenates a new obsra onto the end of the current one.
        """
        newobsra = pd.concat([self.obsra, obsra])
        if inplace: self.obsra = newobsra
        return newobsra 

    def add_obsra(self, obs, fc, inplace=True):
        """
        """
        ##add more observaton forecast pairs to an obsra.
        ##this is useful when wanting to look at obs /forecasts for multiple measurement stations.
        newobsra = self.create_obsra(obs, fc)
        obsra = pd.concat([self.obsra, newobsra])
        if inplace: self.obsra = obsra
        return obsra

    def resample_obsra(self, tm=24, inplace=True, stype='mean'):
        """
        resamples the array to a new time period.
        """
        freqstr = str(tm) + 'H' 
        if stype=='sum':
           obsra = self.obsra.resample(freqstr).sum()
        else:
           obsra = self.obsra.resample(freqstr).mean()
        if inplace: self.obsra = obsra
        return obsra

    def find_corr(self):
        """
        returns pearson correlation coefficient
        """
       
        corr1 = self.obsra.corr()
        return corr1.at['fc','obs']

    def compare_cdfs(self):
        x1, y1  = cdf(self.obsra['obs'])
        x2, y2 = cdf(self.fcra['fc'])
         
    def plotseries(self, ax, clrs, obs=True, lbl='forecast', mult=1):
        """
        plots the obsra
        """
        time = self.obsra.index.tolist()
        if obs: ax.plot(time, self.obsra['obs'], clrs[0], label='observations', linewidth=2)
        ax.plot(time, self.obsra['fc']*mult, clrs[1], label=lbl)


    def computediff_degrees(self):
        # finds absolute and smallest difference between two angles.
        temp = self.obsra.copy()

        def diff(cc, mm):
            wdiff = (np.min([np.abs(cc-mm), 360-np.abs(cc-mm)]))
            # if moving clockwise gives you the smallest angular difference
            # then difference is positive. Otherwise negative.
            if np.abs(cc-mm) > 360-np.abs(cc-mm):
               wdiff = -1 * wdiff
            return wdiff
        temp['diff'] =  temp.apply(lambda row: diff(row['obs'],row['fc']),
                                    axis=1)

        return temp

    def computediff(self):
        # finds absolute and smallest difference between two angles.
        temp = self.obsra.copy()
        def diff(cc, mm):
            wdiff  = cc - mm
            return wdiff
        temp['diff'] =  temp.apply(lambda row: diff(row['obs'],row['fc']),
                                   axis=1)
        return temp


    def plotdiff(self, ax, wdir=False, ptype='scatter'):
        if wdir:
            temp = self.computediff_degrees()
        else:
            temp = self.computediff()
        #temp.set_index('time', inplace=True)
        temp = temp['diff']
        if ptype=='scatter': 
            ax.plot(temp) 

    def autocorr(self, ax, nlist=None):
        # plots the autocorrelation of both the series.
        ts1 = self.obsra['obs']
        ts2 = self.obsra['fc']       
        alist1 = [] 
        alist2 = [] 
        if not nlist: 
           nlist = np.arange(0,48)
        for nnn in nlist:
           alist1.append(ts1.autocorr(lag=nnn)) 
           alist2.append(ts2.autocorr(lag=nnn))
        ax.plot(nlist, alist1, '-k.', label='obs')
        ax.plot(nlist, alist2, '-b.', label='fc')

    def plotscatter(self, ax):
        """
        plot obsra on x and fc array on y
        """
        ax.plot(self.obsra['obs'], self.obsra['fc'],'.')      

    def apply_thresh(self, thresh1, thresh2):
        if not self.obsra.empty:
            r1 =  self.obsra[self.obsra['obs']>=thresh1]
            r2 = r1[r1['fc'] >= thresh2]
            return r2
        else:
            print('EMPTY')
            return self.obsra

    def stat_time(self, stype='sum'):
        times=[1,6,12,24,48,72,96]
        thash = {}
        for tm in times:
            mra = self.resample_obsra(tm=tm, inplace=False, stype=stype)
            dhash = find_stats_function(mra)            
            thash[tm] = dhash['rmse'] 
        return thash

    def find_stats(self, scale=1, thresh=0, bg=0):
        return find_stats_function(self.obsra)

    def find_basics(self, scale=1, thresh=0, bg=0):
        mean = self.obsra['obs'].mean()
        print('Mean of obs', mean)        
        mean = self.obsra['fc'].mean()
        print('Mean of model', mean)        



def find_stats_function(mra, scale=1, thresh=0, bg=0):
    """        
    computes fractional bias, rmse and other statistics of obsra from
    Matched Data class.
    returns a dictionary 

    scale : value to scale the modeled values by.
    bg : background value for the observations
    thresh: 

    """ 
    #thresh=100
    #bg = 20

    hit = 0   #keep track of events both forecast and observed
    miss= 0   #keep track of events observed but not forecast
    alarm =0  #keep track of events forecast but not observed
    nonevent =0  #keep track of events neither forecast nor observed.

    #number of forecast events is hit + alarm
    #number of observed events is miss + hit

    #if thresh > 0:
    #   obsra.loc[obsra.fc < thresh, 'fc'] = tvalue
    #   obsra.loc[obsra.obs < thresh, 'obs'] = tvalue
    #data = newra

    obs = mra.obs
    fc = mra.fc * scale
    mse = 0
    nnn = 0
    obsave =0
    fcave = 0

    #standard sums
    summ=0
    sumc=0
    sumk=0
    #squared sums
    sqrm=0
    sqrc=0
    sqrb=0
    sqrd=0

    wdiff=0
  
    #for bias
    sumf=0

    for mm, cc in zip(obs, fc):
        summ += mm
        sumc += cc
        sumk += 1
        sqrm += mm*mm
        sqrc += cc*cc
        sqrb += mm*cc
        sqrd += (cc-mm)**2

        # this is for units specifed in degrees (such as wind direction)
        wdiff += (np.min([np.abs(cc-mm), 360-np.abs(cc-mm)]))**2

        sumf += (cc-mm)
        #mse += (pair[0] - pair[1])**2
        #nnn+=1
        #obsave += pair[0]
        #fcave += pair[1]
        #if (pair[0]-bg >= thresh and pair[1] >= thresh):
        #   hit +=1
        #elif (pair[0]-bg >= thresh and pair[1] < thresh):
        #   miss +=1
        #elif (pair[0]-bg < thresh and pair[1] > thresh):
        #   alarm +=1
        #elif (pair[0]-bg < thresh and pair[1] < thresh):
        #   nonevent +=1
        
    if sumk>0:

        mbar = summ/float(sumk)
        cbar = sumc/float(sumk)
        numb = sumk
        rmse = (sqrd/float(sumk))**0.5
        wrmse = (wdiff/float(sumk))**0.5
        bias = sumf/float(sumk)
        if mbar>0 and cbar>0:
           nmse = sqrd/mbar/cbar/float(sumk)
        else:
           nmse = -99
        if mbar>0 or cbar>0:
           frac = 2.0 * bias / (cbar + mbar)
        else:
           frac = -99

        #obsave = obsave / float(nnn)
        #fcave = fcave / float(nnn)
        #fb = 2.0 * (fcave - obsave) / (fcave + obsave)
        #rmse =  (mse / float(nnn))**0.5
        #if (alarm+hit)> 0:
        #    far =  alarm/ float(alarm + hit)      #fraction of forecasts which are wrong.
        #else:
        #    far = -1
        #if(hit+miss) > 0:
        #    pod = hit /  float(hit + miss)       #fraction of observations which are detected         
        #else:
        #    pod = -1
        #css = hit/float(hit+alarm) - miss/float(miss+nonevent)
        #tss = (hit*nonevent - alarm*miss) / float(hit+miss) / float(alarm+nonevent)
    else:
        print('WARNING, no pairs to find statistics for')
        rmse = -99
        fb = -99
        far =-1
        pod = -1
        nmse = -99  
        frac = -99
         
    corr1 = mra.corr()
    corr = corr1.at['fc','obs']

    dhash = {}
    #dhash['hit'] = hit
    #dhash['alarm'] = alarm
    #dhash['miss'] = miss
    #dhash['far'] = far
    #dhash['pod'] = pod
    dhash['rmse'] = rmse
    dhash['nmse'] = nmse
    dhash['fb']= frac
    dhash['css'] = -99
    dhash['tss'] = -99
    dhash['corr'] = corr
    dhash['mbar'] = mbar
    dhash['cbar'] = cbar
    dhash['wrmse'] = wrmse
    return dhash


#!/n-home/alicec/anaconda/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import sys 
import numpy as np
import datetime
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.stats import describe
from utilhysplit.evaluation.statmain import cdf
from utilhysplit.evaluation.statmain import MatchedData
from utilhysplit.evaluation import ensemble_tools

def not_used_cdf(newra, scale=1):
    if method==3:   
       #newra=obsra
       #scale=1
       ##CDF matching of the observations and model forecsat.
       robs = np.sort(newra.obs)
       #rfc = np.sort(newra.fc)* scale   #use scaling from linear regression.
       rfc = np.sort(newra.fc)   #don't use scaling from linear regression.
       rfc = rfc * scale
       diff = rfc - robs                
       poly = np.polyfit(robs, diff, pfit)
       obs = obsra['obs']
       fc  = obsra['fc']
       if usescale: fc = fc * scale
       if pfit==0:
          fc = fc - (poly[0])
       elif pfit==1:
          fc= fc - (poly[0]*fc + poly[1])
       elif pfit==2:
          fc= fc - (poly[0]*fc**2 + poly[1]*fc + poly[2])
       elif pfit==3:
          fc= fc - (poly[0]*fc**3 + poly[1]*fc**2 + poly[2]*fc + poly[3])

def match_max(obsra):
    #scale simply by matching the maxminum observation to the maximum forecast.
    max_obs = np.max(obsra.obs)
    max_fc = np.max(obsra.fc)
    scale =  max_obs / float(max_fc)
    poly = [scale]
    return poly

 
def scale_results(obsraobject,  plotdata='none', method=3, pfit=1, fignum=1):
    """
    INPUT
    obsraobject: MatchedData object

    finds scaling factor for forecast.
    ##Method 0 no scaling is done.
    ##Method 1 simply match max observation to max forecast.

    ##Method 2 return slope of line of observations vs. forecast.
    ##plots observations vs. forecast and returns slope
    ##in this case pfit is a threshold value to use- only fit for obs and fc greater than this value.
    ##Method 3 uses CDF matching.
    ##if method 3 is picked, then pfit gives the order of polynomial that is used.

    ##method 4 returns slope of line of sorted observations vs. sorted
    ##forecasts. Forecasts and observations are sorted from highest to lowest

    """
    #robs  
    #rfc

    obsra = obsraobject.obsra
    if method==0:
        return 1, [0], obsraobject
    elif method==1:
        poly = match_max(obsra)
    else:
       return -1


def line_fitA(obsraobject,pfit):
           if pfit > 0:
              newra = obsraobject.apply_thresh(thresh1=pfit, thresh2=pfit)
              print('USING PFIT for linear scaling', pfit)
           else:
              newra = obsraobject.obsra.copy()
           robs = newra.obs
           rfc = newra.fc 
           return apply_polyfit(rfc,robs)

def line_fitB(obsraobject,pfit=None):
           ##this fits a line to the sorted data
           ##so smallest values matched to smallest values, no matter location in space/time.
           newra = obsraobject.obsra.copy()
           robs = np.sort(newra.obs)
           rfc = np.sort(newra.fc)   
           if pfit:
              ##look only at values larger than pfit index.
              robs = robs[int(pfit):]
              rfc = rfc[int(pfit):]
           return apply_polyfit(rfc,robs)

def linefit_adjust(forecast, obs):
    poly = np.polyfit(rfc, robs, 1) #fit line to data.
    slope = poly[0]
    intercept = poly[1]
    adjusted_forecast = forecast * slope + intercept
    return poly 


def apply_polyfit(rfc,robs):
       try: 
           poly = np.polyfit(rfc, robs, 1) #fit line to data.
       except:
           poly  = [1] 
       return poly

def plot_fit(poly,obsra):
       #sys.exit()
       #slope, intercept, rval = linregress(rfc, robs) #fit line to data.
       #poly = [slope, intercept]
       #scale = poly[0] #slope of the line.
       if plotdata!='none':
           fig = plt.figure(fignum)
           ax = fig.add_subplot(1,1,1)
           plt.plot(rfc, robs, 'b.')
           mval = np.max([np.max(obsra.fc), np.max(obsra.obs)]) 
           xlin = np.arange(1, mval,mval/20)
           ylin = poly[0] * xlin  + poly[1]
           plt.plot(xlin, ylin, '-r')
           ax.set_ylabel('Observations')
           ax.set_xlabel('Forecast')
           #ax.set_yscale('log')
           #ax.set_xscale('log')
           plt.title('method' + str(method) +' ' +  plotdata)
           plt.savefig(plotdata + '.jpg')
           plt.show()

#    if method==2 or method ==1 or method==4: 
#      obs = obsra['obs']
#      fc  = obsra['fc']
#      fc = fc * scale
#      print('Linear scaling', scale)
#      #poly1 = np.polyfit(robs, diff, 2)
#      #y2 = poly1[0]*obs**2 + poly1[1]*obs + poly1[2]
#      print('AAAA')
#      scaled_obsra = MatchedData(obs, fc)
#      return scale, poly , scaled_obsra


def cdf_match_volcat(forecast, volcat, 
                     thresh=0.01, ens=0,
                     scale=1,pfit=1, 
                     makeplots=False,
                     figname=None,
                     figsize=(10,5)):
 
    # use pixel matching to determine what forecast values to use.
    try:
        fc = forecast.isel(ens=ens)
    except:
        fc = forecast.sel(ens=ens)
        #print('ens {}'.format(ens))
        #print(forecast.ens.values)
        #return 0, [0,0], [0]
    thresh2, match = ensemble_tools.get_pixel_match(fc,volcat,thresh) 
    #print(thresh2)
    if np.isnan(thresh2):
       thresh2=0
    if 'ens' in match.dims:
        rfc = match.isel(ens=ens).values.flatten()
    else:      
        rfc = match.values.flatten()
    #print('pixel match threshold', thresh2)
    robs = volcat.values.flatten()
    #if thresh2 == 0:
    #    rfc2 = [x for x in rfc if x >0]
    #else:
    if thresh2!=0: rfc2 = [x for x in rfc if x >=thresh2]
    else: rfc2 = rfc
    robs2 = [x for x in robs if x >=thresh]
    # make sure arrays are the same shape. Two methods.
 
    # pixel matching should take care of making sure
    # the arrays are same length unless there are more
    # obs than forecasts.
    #print(ens, len(rfc2), len(robs2), thresh2)
    # This block chops the longer arrays.
    # If the forecast has less pixels then chop the observations.
    if len(rfc2) < len(robs2):
       print('chop the obs', ens, thresh, thresh2)
       robs2 = np.sort(robs2)
       robs2 = robs2[len(robs2)-len(rfc2):]

    if len(rfc2) > len(robs2):
       print('chop the forecast', ens, thresh, thresh2)
       rfc2 = np.sort(rfc2)
       rfc2 = rfc2[len(rfc2)-len(robs2):]


    # This block adds zeros to shorter arrays.
    #while len(rfc2) > len(robs2):
    #      robs2.append(0)
    #while len(rfc2) < len(robs2):
    #      rfc2.append(0)

    outputs = cdf_match(rfc2, robs2,scale=scale,pfit=pfit,
                        makeplots=makeplots,figname=figname,
                        figsize=figsize)
    return thresh2, outputs[0], outputs[1], outputs[2]

def cdf_match(rfc, robs, scale=None,pfit=1,
              makeplots=False,
              figname=None,
              figsize=(8,5)):
   """
   rfc : list of forecast values
   robs : list of observations
   Currently length of rfc and robs needs to match.
   pfit : int (0,1,2,3) determines order of polynomial to use
          -1 linear fit that must go through 0.
   # RETURNS
   poly :  polynomial to use for scaling.
   fc :    the scaled values.
   rhash : dictionary
   """
   #obs = robs
   #fc  = rfc
   rhash = {}
   sorted_obs = np.sort(robs)
   if scale:
       sorted_fc = np.sort(rfc)* scale    #use scaling from linear regression.
   else:
       sorted_fc = np.sort(rfc)           #don't use scaling from linear regression.
   #diff = sorted_fc -  sorted_obs      
   diff = (sorted_fc -  sorted_obs)   
   #plt.plot(sorted_obs, diff, 'k.')
   #plt.show()         
   fco = np.array(sorted_fc)
   #print('fco', fco)
   #print('diff', diff)
   poly, covmatrix = np.polyfit(fco, diff, deg=np.abs(pfit),rcond=None,cov=True)
   poly, rhash['residuals'], rank, svs, rcond = np.polyfit(fco, diff, deg=np.abs(pfit),rcond=None,full=True)
   if pfit==0:
      y2 = poly[0]*np.ones_like(fco)
   elif pfit==1 or pfit==-1:
      y2 = poly[0]*fco + poly[1]
      if pfit==-1:
          aaa, b,c,d = np.linalg.lstsq(fco[:,np.newaxis],diff,rcond=None)
          y3 = fco*aaa
          poly = np.append(poly,aaa)
   elif pfit==2:
      y2 = poly[0]*fco**2 + poly[1]*fco + poly[2]
   elif pfit==3:
      y2 = poly[0]*fco**3 + poly[1]*fco**2 + poly[2]*fco + poly[3]
   fc = fco - y2
   try:
      fcp = fco - y3
   except:
      fcp = None
   if makeplots:
       fig = plt.figure(1,figsize=figsize)
       ax1 = fig.add_subplot(1,1,1)
       if pfit==-1:
           plot_cdf_match_fit(fco,y2,y3,diff,ax1)
       else: 
           plot_cdf_match_fit(fco,y2,None,diff,ax1)
       xpos=0.1
       ypos=0.8

       if len(poly) == 1:
           dstr = "y = {:0.1f}".format(poly[0])

       elif len(poly)== 2 or pfit==-1:
           dstr = "y = {:0.2f}x + ({:0.1f})".format(poly[0],poly[1])
           xpos=0.1
           ypos=0.6
       elif len(poly)==3 and pfit == -1:
           dstr = "y = {:0.2f}x ".format(poly[2])
       elif len(poly)==3:
           dstr = "y = {:0.2f}$x^2$ + ({:0.1f})x + ({:0.1f})".format(poly[0],poly[1],poly[2])
       elif len(poly)==4:
           dstr = "y = {:0.2f}$x^3$ + ({:0.1f})$x^2$ + ({:0.1f})x + ({:0.1f})".format(poly[0],poly[1],poly[2],poly[3])
       ax1.text(xpos,ypos,dstr,
                transform=ax1.transAxes,
                bbox=dict(facecolor='white'),size=16)

       if pfit == -1:
           dstr = "y = {:0.2f}x ".format(poly[2])
           xpos=0.1
           ypos=0.8
           ax1.text(xpos,ypos,dstr,
                transform=ax1.transAxes,
                bbox=dict(facecolor='white'),size=16)

       ax1.text(xpos,ypos,dstr,
                transform=ax1.transAxes,
                bbox=dict(facecolor='white'),size=16)
 
       
       if figname: plt.savefig('{}_fit.{}'.format(figname,'png'))


       plt.show()
       fig = plt.figure(2,figsize=figsize)
       ax2 = fig.add_subplot(1,1,1)
       plot_cdf_match(robs,fco,fc,fcp)
       dobs = describe(robs)
       dfco = describe(fco)
       dfc =  describe(fc)
       if isinstance(fcp, np.ndarray):
           dfcp = describe(fcp)
       else:
           dfcp = describe(fc)
       mean = [dobs.mean, dfco.mean,dfc.mean,dfcp.mean]
       var = [dobs.variance, dfco.variance,dfc.variance,dfcp.variance]
       skew = [dobs.skewness, dfco.skewness,dfc.skewness,dfcp.skewness]
       kurt = [dobs.kurtosis, dfco.kurtosis,dfc.kurtosis,dfcp.kurtosis]
       small = [dobs.minmax[0], dfco.minmax[0],dfc.minmax[0],dfcp.minmax[0]]
       big = [dobs.minmax[1], dfco.minmax[1],dfc.minmax[1],dfcp.minmax[1]]
       num = [dobs.nobs, dfco.nobs, dfc.nobs, dfcp.nobs]
       name = ['obs','fc','corrected','corrected 0']
       dfout = pd.DataFrame(zip(name,mean,var,skew,kurt,num,small,big))
       dfout.columns=['name','mean','variance','skewness','kurtosis','N','min','max']
       print(dfout)

       if figname: plt.savefig('{}_corrected.{}'.format(figname,'png'))
       plt.show()
   return poly, fc, rhash

def plot_cdf_match_fit(fco,y2,y3,diff,ax1):
    fs=18    
    ax1.plot(fco, diff,'-k',linewidth=4,label='Difference' )
    ax1.plot(fco, y2,'-c', label='fit')
    if isinstance(y3,np.ndarray): ax1.plot(fco, y3,'-r', label='fit')
    ax1.set_ylabel('Difference (g m$^{-2}$)', fontsize=fs)
    ax1.set_xlabel('Forecast value (g m$^{-2}$)', fontsize=fs)
    ax1.tick_params(labelsize=fs)
    plt.tight_layout()

def plot_cdf_match(robs, fco, fc,fcp):
       fs = 18
       x1, y1 = cdf(robs)
       x2, y2 = cdf(fco)
       x3, y3 = cdf(fc)
       if isinstance(fcp, np.ndarray): x4, y4 = cdf(fcp)
       plt.step(x1, y1, '-k', label='Obs',linewidth=4)
       plt.step(x2, y2, '-g', label='Forecast')  #blue shows cdf of forecast
       plt.step(x3, y3, '--c', label='Corrected',linewidth=3) #green shows cdf of scaled forecast.
       if isinstance(fcp,np.ndarray):  plt.step(x4, y4, '--r', label='Corrected',linewidth=3) #green shows cdf of scaled forecast.
       #plt.title(plotdata + ' red(obs), blue(forecast), green(scaled forecast)')
       ax = plt.gca()
       ax.set_ylabel('CDF',fontsize=fs)
       ax.set_xlabel('Mass Loading (g m$^{-2}$)', fontsize=fs)
       handles,labels = ax.get_legend_handles_labels()
       ax.legend(handles,labels,loc='lower right',fontsize=fs)
       ax.tick_params(labelsize=fs)
       plt.tight_layout()

def old_function(): 
   #scaled_obsra = MatchedData(obs, fc)
   sp= True
   ##This block plots differences between cdfs and fit. 
   if sp:
       xxx = fco
       fig = plt.figure(1)
       ax1 = fig.add_subplot(1,1,1)
       ax1.plot(xxx, diff,'-b')
       ##blue is the differences 
         
       poly1 = np.polyfit(xxx, diff, 0)
       y2 = xxx*0 +  poly1[0]
       ax1.plot(xxx, y2, '-k')
       ##black is a line fit to the difference
       fc0 = fco - (poly1[0])

       poly1 = np.polyfit(xxx, diff, 1)
       y2 = poly1[0]*xxx + poly1[1]
       ax1.plot(xxx,y2, '-c')
       ##cyan is a first order fit to the differences 
       fc1 = fco - (poly1[0]*fco + poly1[1])

       poly1 = np.polyfit(robs, diff, 2)
       y2 = poly1[0]*robs**2 + poly1[1]*robs + poly1[2]
       #ax1.plot(robs, y2, '-r')
       fc2 = fco - (poly1[0]*fco**2 + poly1[1]*fco + poly1[2])
       ##red is a second order fit to the differences 

       poly1 = np.polyfit(robs, diff, 3)
       y3 = poly1[0]*robs**3 + poly1[1]*robs**2 + poly1[2]*robs + poly1[3]
       ##green is a third order fit to the differences 
       #ax1.plot(robs, y3, '-g')
       fc3 = fco - (poly1[0]*fco**3 + poly1[1]*fco**2 + poly1[2]*fco + poly1[3])
       plt.show()

       #yn = robs - y2
       #plt.plot(robs, yn, '-r')
       #yn = robs - y3
       #plt.plot(robs, yn, '-g')
       #plt.title('Difference between cdfs and fit')
       #plt.title(plotdata)
       #plt.show() 
   ##Plots forecasts transformed by cdf matching.
   sp=False
   if sp: 
       plt.plot(obs, '-r')  #measurments are red.
       plt.plot(fc*scale, '-b') #original forecasts in blue
       #plt.plot(newra.fc, '-c') #original forecasts in blue
       plt.plot(fc0, '-k')
       plt.plot(fc1, '-c')
       plt.plot(fc2, '-y')
       plt.plot(fc3, '-g')
       #newfc = newra.fc *scale
       #newfc2 = newfc - (poly[0]*newfc**2 + poly[1]*newfc + poly[2])
       #plt.plot(newfc2, '--g')  #new forecasts in green.
       #plt.plot(scaled_obsra.obsra['fc'], '--m')  #new forecasts in green.
       #plt.title(plotdata + 'Green (forecast scaled), red (obs), blue (unscaled forecasts)')
       plt.show()
   ##This block plots the cdf's of the observed and forecast.
   sp=True
   if sp:
       x1, y1 = cdf(robs)
       x2, y2 = cdf(fco)
       x3, y3 = cdf(fc)
       plt.step(x1, y1, '-r')
       plt.step(x2, y2, '-b')  #blue shows cdf of forecast
       plt.step(x3, y3, '--g') #green shows cdf of scaled forecast.
       #plt.title(plotdata + ' red(obs), blue(forecast), green(scaled forecast)')
       plt.show()
   return scale, poly, obs, fc

#----------------------------------------------------------------------------------------------------------------------------------------


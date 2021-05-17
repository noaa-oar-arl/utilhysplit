#!/n-home/alicec/anaconda/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import sys 
import numpy as np
import datetime
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

def cdf(newra, scale=1):
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
       print 'Scaling method 3 with pfit ' , pfit
       if pfit==0:
          fc = fc - (poly[0])
       elif pfit==1:
          fc= fc - (poly[0]*fc + poly[1])
       elif pfit==2:
          fc= fc - (poly[0]*fc**2 + poly[1]*fc + poly[2])
       elif pfit==3:
          fc= fc - (poly[0]*fc**3 + poly[1]*fc**2 + poly[2]*fc + poly[3])


 
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
        #scale simply by matching the maxminum observation to the maximum forecast.
        max_obs = np.max(obsra.obs)
        max_fc = np.max(obsra.fc)
        scale =  max_obs / float(max_fc)
        poly = [scale]
    ##this block fits a line to plot of obs vs forecast.
    ##if method 3, will then go on to do cdf matching.
    ##if method 2, then will return this fit.
    else:
       ##may want to fit line only to larger events
       #newra = obsra[obsra.obs >= 50]
       #newra = newra[newra.fc > 0]
       newra = obsra.copy()

       if method==2 or method==30:
           #newra = obsra.copy()
           if pfit > 0:
              newra = obsraobject.apply_thresh(thresh1=pfit, thresh2=pfit)
              print 'USING PFIT for linear scaling', pfit 
           robs = newra.obs
           rfc = newra.fc 

       elif method==4 or method==3:
           ##this fits a line to the sorted data
           ##so smallest values matched to smallest values, no matter location in time.
           robs = np.sort(newra.obs)
           rfc = np.sort(newra.fc)   #use scaling from linear regression.
           if method==4:
              ##look only at 20 highest values
              print 'METHOD 4 pfit' , pfit
              robs = robs[int(pfit):]
              rfc = rfc[int(pfit):]
       #background = find_background(robs)

       try: 
           poly = np.polyfit(rfc, robs, 1) #fit line to data.
           print 'POLYFIT SUCCEEDED'
           print rfc[0:10]
           print robs[0:10]
           print('---------------------------')
           print rfc[-10:]
           print robs[-10:]
       except:
           poly  = [1] 
           print 'POLYFIT FAILED'
           print rfc[0:10]
           print robs[0:10]
       #sys.exit()
       #slope, intercept, rval = linregress(rfc, robs) #fit line to data.
       #poly = [slope, intercept]
       scale = poly[0] #slope of the line.
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


    if method==2 or method ==1 or method==4: 
      obs = obsra['obs']
      fc  = obsra['fc']
      fc = fc * scale
      print('Linear scaling', scale)
      #poly1 = np.polyfit(robs, diff, 2)
      #y2 = poly1[0]*obs**2 + poly1[1]*obs + poly1[2]
      print('AAAA')
      scaled_obsra = MatchedData(obs, fc)
      return scale, poly , scaled_obsra
    usescale=False
    if method==3:   
       #newra=obsra
       #scale=1
       ##CDF matching of the observations and model forecsat.
       robs = np.sort(newra.obs)
       #rfc = np.sort(newra.fc)* scale   #use scaling from linear regression.
       rfc = np.sort(newra.fc)   #don't use scaling from linear regression.
       if usescale: rfc = rfc * scale
       diff = rfc - robs                
       poly = np.polyfit(robs, diff, pfit)
       obs = obsra['obs']
       fc  = obsra['fc']
       if usescale: fc = fc * scale
       print 'Scaling method 3 with pfit ' , pfit
       if pfit==0:
          fc = fc - (poly[0])
       elif pfit==1:
          fc= fc - (poly[0]*fc + poly[1])
       elif pfit==2:
          fc= fc - (poly[0]*fc**2 + poly[1]*fc + poly[2])
       elif pfit==3:
          fc= fc - (poly[0]*fc**3 + poly[1]*fc**2 + poly[2]*fc + poly[3])
       print('BBBB')
       scaled_obsra = MatchedData(obs, fc)
       #sp=False
       if plotdata == 'all': sp=True
       sp= True
       ##This block plots differences between cdfs and fit. 
       if sp:
           fig = plt.figure(1)
           ax1 = fig.add_subplot(1,1,1)
           ax1.plot(robs, diff,'-b')
           ##blue is the differences 
             
           poly1 = np.polyfit(robs, diff, 0)
           y2 = robs*0 +  poly1[0]
           ax1.plot(robs, y2, '-k')
           ##black is a line fit to the difference
           fc0 = fc - (poly1[0])
 
           poly1 = np.polyfit(robs, diff, 1)
           y2 = poly1[0]*robs + poly1[1]
           ax1.plot(robs, y2, '-c')
           ##cyan is a first order fit to the differences 
           fc1 = fc - (poly1[0]*fc + poly1[1])

           poly1 = np.polyfit(robs, diff, 2)
           y2 = poly1[0]*robs**2 + poly1[1]*robs + poly1[2]
           ax1.plot(robs, y2, '-r')
           fc2 = fc - (poly1[0]*fc**2 + poly1[1]*fc + poly1[2])
           ##red is a second order fit to the differences 

           poly1 = np.polyfit(robs, diff, 3)
           y3 = poly1[0]*robs**3 + poly1[1]*robs**2 + poly1[2]*robs + poly1[3]
           ##green is a third order fit to the differences 
           ax1.plot(robs, y3, '-g')
           fc3 = fc - (poly1[0]*fc**3 + poly1[1]*fc**2 + poly1[2]*fc + poly1[3])
           plt.show()
 
           #yn = robs - y2
           #plt.plot(robs, yn, '-r')
           #yn = robs - y3
           #plt.plot(robs, yn, '-g')
           #plt.title('Difference between cdfs and fit')
           plt.title(plotdata)
           #plt.show() 
       ##Plots forecasts transformed by cdf matching.
       if sp: 
           plt.plot(newra.obs, '-r')  #measurments are red.
           plt.plot(newra.fc*scale, '-b') #original forecasts in blue
           #plt.plot(newra.fc, '-c') #original forecasts in blue
           plt.plot(fc0, '-k')
           plt.plot(fc1, '-c')
           plt.plot(fc2, '-y')
           plt.plot(fc3, '-g')
           #newfc = newra.fc *scale
           #newfc2 = newfc - (poly[0]*newfc**2 + poly[1]*newfc + poly[2])
           #plt.plot(newfc2, '--g')  #new forecasts in green.
           plt.plot(scaled_obsra.obsra['fc'], '--m')  #new forecasts in green.
           plt.title(plotdata + 'Green (forecast scaled), red (obs), blue (unscaled forecasts)')
           plt.show()
       ##This block plots the cdf's of the observed and forecast.
       if sp:
           print 'zzz scale poly', scale, poly
           x1, y1 = cdf(robs)
           x2, y2 = cdf(rfc)
           x3, y3 = cdf(fc)
           plt.step(x1, y1, '-r')
           plt.step(x2, y2, '-b')  #blue shows cdf of forecast
           plt.step(x3, y3, '--g') #green shows cdf of scaled forecast.
           plt.title(plotdata + ' red(obs), blue(forecast), green(scaled forecast)')
           plt.show()
       ##This block plots relationship between .
       return scale, poly, scaled_obsra


def cdf(data):
    sdata = np.sort(data)
    y = np.arange(sdata.size)/  float(sdata.size)
    return sdata, y 

#----------------------------------------------------------------------------------------------------------------------------------------


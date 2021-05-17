import sys
import os
import datetime
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import numpy as np
import numpy.ma as ma
import pandas as pd
import seaborn as sns
from utilvolc import volcat
from monetio.models import hysplit
from utilhysplit import hcontrol
from utilhysplit.evaluation import plume_stat
from utilvolc.ash_inverse import InverseAsh
import utilvolc.ash_inverse as ainv


class AshEval(InverseAsh):

    def __init__(self, tdir, fname,vdir,vid,configdir=None,configfile=None):
        super().__init__(tdir,fname,vdir,vid,configdir,configfile)
        
    def get_cdump(self,tdir,fname):
        cdump = hysplit.open_dataset(os.path.join(tdir,fname))
        self.cdump = cdump 

    def calc_fss(self,forecast,thresh=None):
        volcat = self.match_volcat(forecast)
        evals = forecast.values
        vvals = volcat.values
        scores = plume_stat.CalcScores(volcat, forecast, threshold=0.2)
        return scores

    def compare_forecast_dist(self,forecast,thresh=None):
        from utilhysplit.evaluation import statmain
        volcat = self.match_volcat(forecast)
        evals = forecast.values
        vvals = volcat.values
        exval,eyval =  statmain.nancdf(evals.flatten(),thresh)
        vxval,vyval =  statmain.nancdf(vvals.flatten(),thresh)
        fig = plt.figure(1)
        ax = fig.add_subplot(1,1,1)
        ax.step(exval, eyval, '-r',label='HYSPLIT')
        ax.step(vxval, vyval, '-b',label='Volcat')
        ks1,ks2 = statmain.kstestnan(evals.flatten(),vvals.flatten(),thresh)
        try:
            print('Kolmogorov-Smirnov Parameter {} {}'.format(np.max(np.abs(ks1)), np.max(np.abs(ks2))))
        except:
            return 1
        return np.max(np.abs(ks1)) 

    def compare(self,forecast,thresh=None):
        volcat = self.match_volcat(forecast)
        evals = self.remove_nans(forecast.values.flatten(),thresh)
        vvals = self.remove_nans(volcat.values.flatten(),thresh)
        # number of pixles above threshold in each.
        forecast_parea = len(evals) 
        volcat_parea = len(vvals)
        return forecast_parea, volcat_parea
       
    def compare_forecast(self, forecast,cmap='Blues',ptype='pcolormesh',vloc=None):
        # forecast should be an xarray in mass loading format with no time dimension.
        sns.set()
        sns.set_style('whitegrid')
        fig = plt.figure(figsize=[15,5])
        ax1 = fig.add_subplot(1,3,1)
        ax2 = fig.add_subplot(1,3,2)
        ax3 = fig.add_subplot(1,3,3)

        time = pd.to_datetime(forecast.time.values)
        tii = self.time_index(time)
        print('tii', tii)
        volcat = self.volcat_avg_hash[tii] 
        evals = forecast.values
        vpi = evals < 0.001
        evals[vpi] = np.nan
        norm = self.get_norm(volcat,forecast)
   
        #evals = np.log10(evals)
        #vvals = np.log10(volcat.values)
        vvals = volcat.values
        vpi = vvals < 0.001
        vvals[vpi] =  np.nan
        clevels = [0.2,2,5,10]
        if ptype == 'pcolormesh':
            cb = ax1.pcolormesh(volcat.longitude, volcat.latitude,vvals,norm=norm, cmap=cmap,shading='nearest')
            cb2 = ax2.pcolormesh(forecast.longitude, forecast.latitude,evals,norm=norm, cmap=cmap,shading='nearest')
            #cb2 = ax3.pcolormesh(forecast.longitude, forecast.latitude,evals,norm=norm, cmap=cmap,shading='nearest')
            ax3.contourf(volcat.longitude, volcat.latitude, volcat.values,levels=clevels,cmap='Reds') 
            ax3.contour(forecast.longitude, forecast.latitude, evals,levels=clevels,cmap='viridis') 
        plt.title(time.strftime("%Y %m/%d %H:%M UTC"))
        plt.colorbar(cb, ax=ax1) 
        plt.colorbar(cb2, ax=ax2) 
        ylim = ax1.get_ylim()
        ax2.set_ylim(ylim)
        ax3.set_ylim(ylim)
        xlim = ax1.get_xlim()
        ax2.set_xlim(xlim)
        ax3.set_xlim(xlim)
        if vloc:
           ax1.plot(vloc[0],vloc[1],'y^')
           ax2.plot(vloc[0],vloc[1],'y^')
           ax3.plot(vloc[0],vloc[1],'y^')
        return fig

    def compare_plots(self, daterange,levels=None):
        fig = plt.figure(1,figsize=(10,5))
        ax1 = fig.add_subplot(1,2,1)
        ax2 = fig.add_subplot(1,2,2)
        tii = self.time_index(daterange[0])
        print('tii',tii)
        cdump = self.cdump_hash[tii]
        volcat = self.volcat_avg_hash[tii] 
        csum = cdump.sum(dim='ens')
        volcat.plot.pcolormesh(x='longitude',y='latitude',levels=levels,ax=ax1)
        #cdump.sum(dim='ens').plot.contour(x='longitude',y='latitude',ax=ax2)
        #plt.pcolormesh(csum.longitude, csum.latitude, np.log10(csum),cmap='Reds')
        #plt.pcolormesh(csum.longitude, csum.latitude, csum,cmap='Reds',shading='nearest')
        cb= csum.plot.pcolormesh(x='longitude',y='latitude',cmap='viridis',ax=ax2)
        #cb = plt.pcolormesh(volcat.longitude, volcat.latitude, np.log10(volcat),cmap='Blues',levels=levels)
        #cb = plt.scatter(volcat.longitude, volcat.latitude, c=np.log10(volcat),s=2,cmap='Blues')
        #cb = plt.scatter(volcat.longitude, volcat.latitude, c=volcat.values,s=2,cmap='viridis',levels=levels)
        #cb = plt.contour(volcat.longitude, volcat.latitude, np.log10(volcat),cmap='Blues')
        # plt.colorbar(cb)
        plt.tight_layout()
        return ax1,ax2


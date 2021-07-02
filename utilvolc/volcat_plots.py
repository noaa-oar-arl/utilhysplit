import numpy as np
import pandas as pd
import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import xarray as xr
from utilvolc import volcat

class VolcatPlots:

    def __init__(self,dsetlist):
        self.dset = dsetlist
        self.set_plot_settings()

    def make_arrays(self):
        print('here')
        das = self.dset
 
        #sns.set()
        #sns.set_style('whitegrid')
        masslist = []
        dlist = []
        # from maximum value of retrieval
        htlist = []
        alist = []
        tmasslist =[]
        arealist = []
        maxmass = []
        minmass = []
        mer = []
        massprev = 0
        radius = []
        minradius = []
        maxradius = []
        dtlist = []
        time_elapsed=0

        for iii in np.arange(0,len(das)+1):
            try:
                vmass  = volcat.get_mass(das[iii],clip=True)
            except:
                break
            vht  = volcat.get_height(das[iii],clip=True)
            vrad  = volcat.get_radius(das[iii],clip=True)

            maxmass.append(np.max(vmass))
            minmass.append(np.min(vmass))

            masstg = das[iii].ash_mass_loading_total_mass.values[0]
            tmasslist.append(masstg)

            # date
            dlist.append(das[iii].time.values)

            # max height in detection
            htlist.append(float(np.max(vht)))

            # area of feature
            arealist.append(das[iii].feature_area.values[0])

            # mean effective radius
            radius.append(vrad.mean())
            minradius.append(np.min(vrad))
            maxradius.append(np.max(vrad))

            
            if (pd.to_datetime(das[iii-1].time.values[0]) > pd.to_datetime(das[iii].time.values[0])):
                print('WARNING: time not increasing ')
                print(das[iii-1].time.values)
                print(das[iii].time.values)

            # mer in kg/s. divide by 1e9 to convert from Tg to kg.
            # divide by 10*60 to get from per retrieval to per second.
            if iii>0:
                dt = pd.to_datetime(das[iii].time.values[0]) - pd.to_datetime(das[iii-1].time.values[0])
                mer.append((masstg - massprev)*1e9/(dt.seconds))
            else:
                dt = pd.to_datetime(das[iii+1].time.values[0]) - pd.to_datetime(das[iii].time.values[0])
                mer.append(0)
            massprev = masstg
            dtlist.append(time_elapsed)
            time_elapsed += dt.seconds 

        self.dtlist = dtlist
        self.dlist = dlist
        self.tmasslist = tmasslist
        self.minmass = minmass
        self.maxmass = maxmass
        self.arealist = arealist
        self.htlist = htlist
        self.mer = mer
        self.averad = radius
        self.minrad = minradius
        self.maxrad = maxradius

    def make_spline(self,s=None):
        import scipy.interpolate 
        if not s:
           s = 2.0 / len(self.dtlist)
        #self.spline = scipy.interpolate.CubicSpline(self.dtlist,self.tmasslist)
        self.spline = scipy.interpolate.UnivariateSpline(self.dtlist,self.tmasslist,s=s)

    def set_plot_settings(self):
        self.main_clr = '--b'
        self.spline_clr = '-r'
        self.sub_clrs = ['--r','--c']

    def plot_multiB(self,fignum=1):
        sns.set_style('whitegrid')
        fig = plt.figure(fignum,figsize=[10,5])

        ax1 = fig.add_subplot(2,1,1)
        self.sub_plot_radius(ax1)

        ax3 = fig.add_subplot(2,1,2)
        self.sub_plot_min_mass(ax3,both=True)

        fig.autofmt_xdate()

        return fig

    def plot_multiA(self,fignum=1):
        self.make_spline()

        sns.set_style('whitegrid')
        fig = plt.figure(fignum,figsize=[10,10])

        ax1 = fig.add_subplot(2,2,1)
        self.sub_plot_mass(ax1)

        ax2 = fig.add_subplot(2,2,2)
        self.sub_plot_area(ax2)

        ax3 = fig.add_subplot(2,2,3)
        self.sub_plot_mer(ax3)

        ax4 = fig.add_subplot(2,2,4)
        self.sub_plot_maxht(ax4)
  
        fig.autofmt_xdate()
        plt.tight_layout()
        return fig

    def sub_plot_radius(self,ax):
        yval = self.averad
        yval2 = self.minrad
        yval3 = self.maxrad
        xval = self.dlist
        ax.plot(xval,yval,self.main_clr,label='Average')
        ax.plot(xval,yval2,self.sub_clrs[0],label='Minimum')
        ax.plot(xval,yval3,self.sub_clrs[1],label='Maximum')
        ax.set_ylabel('Effective radius')
        ax.set_xlabel('Time')
        
    def sub_plot_mass(self,ax):
        yval = self.tmasslist
        xval = self.dlist

        ys = self.spline(self.dtlist)
        ax.plot(xval,ys,self.spline_clr, LineWidth=5,alpha=0.5)

        ax.plot(xval,yval,self.main_clr)


        ax.set_ylabel('Total mass (Tg)')
        ax.set_xlabel('Time')

    def sub_plot_max_mass(self,ax):
        yval = self.maxmass
        xval = self.dlist
        ax.plot(xval,yval,self.main_clr)
        ax.set_ylabel('Maximum Mass Loading (g m$^{-2}$)')
        ax.set_xlabel('Time')
     
    def sub_plot_min_mass(self,ax,both=False):
        yval = self.minmass
        yval2 = self.maxmass
        xval = self.dlist
        ax.plot(xval,yval,self.main_clr, label='Minimum')
        if both: ax.plot(xval,yval2,self.sub_clrs[0], label='Maximum')
        ax.set_ylabel('Mass Loading (g m$^{-2}$)')
        ax.set_xlabel('Time')
        ax.set_yscale('log')
        ax.set_ylim([0.001,50])
        
    def sub_plot_area(self,ax,clr=-1):
        yval = self.arealist
        xval = self.dlist
        if clr < 0:
            ax.plot(xval,yval,self.main_clr)
        else:
            ax.plot(xval,yval,self.sub_clrs[clr])
        ax.set_ylabel('Total Area (km$^2$)')
        ax.set_xlabel('Time')

    def sub_plot_mer(self,ax):
        yval = self.mer
        xval = self.dlist
        #ax.plot(xval,yval,self.main_clr)
        mer = self.spline.derivative()
        ys = mer(self.dtlist)
        ax.plot(xval,ys*1e9,self.spline_clr, LineWidth=5,alpha=0.5)
        ax.set_ylabel('kg s$^{-1}$')
        ax.set_xlabel('Time')

    def sub_plot_maxht(self,ax):
        yval = self.htlist
        xval = self.dlist


        ax.plot(xval,yval,self.main_clr)
        ax.set_ylabel('Maximum height km')
        ax.set_xlabel('Time')

def find_area(vmass):
    r2 = vmass.where(vmass>0)
    r2 = r2.fillna(0)
    #place ones where above threshold
    r2 = r2.where(r2<=0)
    r2 = r2.fillna(1)
    return int(r2.sum())


        
        #sns.set()
        #sns.set_style('whitegrid')

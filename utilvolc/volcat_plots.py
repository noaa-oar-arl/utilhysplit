import numpy as np
from scipy.stats import describe
import pandas as pd
import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import xarray as xr
from utilvolc import volcat
from utilhysplit.evaluation import hysplit_boxplots
from utilhysplit.evaluation import statmain



class VolcatPlots:

    def __init__(self,dsetlist):
        """
        dsetlist can be from output of volcat.get_volcat_list function
        """
        # sort dset list by time
        def ftime(x):
            return x.time.values
        dsetlist.sort(key=ftime)
        self.dset = dsetlist
        self.set_plot_settings()

    def make_arrays(self):
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
        self.vmass =  []

        for iii in np.arange(0,len(das)+1):
            try:
                vmass  = volcat.get_mass(das[iii],clip=True)
            except:
                print('cannot get mass',iii)
                break
            vht  = volcat.get_height(das[iii],clip=True)
            vrad  = volcat.get_radius(das[iii],clip=True)

            self.vmass.append(vmass)  
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

    def boxplotdata(self, datelist, vdata):
        self.dj = hysplit_boxplots.prepare_boxplotdata(datelist,vdata)


    def make_boxplot(self, cols=None):
        hysplit_boxplots.make_boxplot(self.dj, cols=cols)

    def volcat_describe_plot(self,threshold=0,nums=None):
        date = []
        ens = []
        mean = []
        var = []
        skew = []
        kurt = []
        num = []
        small = []
        big = []
        vdata = []
        
        clr = ["-m", "-r", "-b", "-c", "-g"]
        if isinstance(nums,list): vlist = self.vmass[nums[0]:nums[1]]
        else: vlist = self.vmass
        for jjj, volc in enumerate(vlist):
            volc = volc.values.flatten()
            volc = [x for x in volc if x > threshold]
            sts = describe(volc)
            vdata.append(volc)
            date.append(self.dlist[jjj][0])
            ens.append('obs')
            mean.append(sts.mean)
            var.append(sts.variance)
            skew.append(sts.skewness)
            kurt.append(sts.kurtosis)
            num.append(sts.nobs)
            small.append(sts.minmax[0])
            big.append(sts.minmax[1])
        self.boxplotdata(date,vdata)
        data = zip(date,ens,mean,var,skew,kurt,num,small,big)
        colnames = ['date','ens','mean','variance','skewness','kurtosis','N','min','max']
        dfout = pd.DataFrame(data)
        dfout.columns = colnames
        self.dfstats = dfout
        return dfout

    def volcat_cdf_plot(self, threshold=0, nums=None):
        #step = 5
        clr = ["-m", "-r", "-b", "-c", "-g"]
        if isinstance(nums,list): vlist = self.vmass[nums[0]:nums[1]]
        else: vlist = self.vmass
        for jjj, volc in enumerate(vlist):
            # print(self.cdump.time.values[iii])
            #volcat = self.volcat_avg_hash[iii]
            volc = volc.values.flatten()
            volc = [x for x in volc if x > threshold]
            try:
                sdata, y = statmain.cdf(volc)
            except:
                print("cannot calculate cdf for {}".format(jjj))
                continue
            ax = plt.gca()
            if jjj % 5 == 0:
                lw = 3
                # print('here')
            else:
                lw = 1
            ax.step(sdata, y, clr[jjj % len(clr)], LineWidth=lw)
        return ax

    def make_spline(self,s=20):
        import scipy.interpolate 
        s = s / len(self.dtlist)
        #self.spline = scipy.interpolate.CubicSpline(self.dtlist,self.tmasslist)
        self.spline = scipy.interpolate.UnivariateSpline(self.dtlist,self.tmasslist,s=s)

    def set_plot_settings(self):
        self.main_clr = '-b.'
        self.spline_clr = '-r'
        self.sub_clrs = ['r.','c.']

    def plot_multiB(self,fignum=1):
        sns.set_style('whitegrid')
        fig = plt.figure(fignum,figsize=[10,5])

        ax1 = fig.add_subplot(2,1,1)
        self.sub_plot_radius(ax1)

        ax3 = fig.add_subplot(2,1,2)
        self.sub_plot_min_mass(ax3,both=True)

        fig.autofmt_xdate()

        return fig

    def plot_dist_stats(self,fignum=2):
        sns.set_style('whitegrid')
        fig = plt.figure(fignum,figsize=[10,10])

        ax1 = fig.add_subplot(2,2,1)
        self.sub_plot_mean(ax1)


        ax2 = fig.add_subplot(2,2,2)
        self.sub_plot_variance(ax2)

        ax3 = fig.add_subplot(2,2,3)
        self.sub_plot_skew(ax3)

        ax4 = fig.add_subplot(2,2,4)
        self.sub_plot_kurt(ax4)

        ax5 = ax4.twinx()
        self.sub_plot_num(ax5)
        ax5.grid(False)

        fig.autofmt_xdate()
        plt.tight_layout()
        return fig

    def plot_multiA(self,fignum=1,smooth=20,yscale='linear'):
        self.make_spline(s=smooth)

        sns.set_style('whitegrid')
        fig = plt.figure(fignum,figsize=[10,10])

        ax1 = fig.add_subplot(2,2,1)
        self.sub_plot_mass(ax1,yscale=yscale)

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
       
    def fit_exp_decay(self,xval,yval):
        ylogval = np.log(yval)
        k, alog = np.polyfit(xval,ylogval,1)
        return k, alog  

    def sub_plot_mean(self,ax):
        xval = self.dfstats['date']
        yval = self.dfstats['mean']
        yval2 = self.dfstats['min']
        yval3 = self.dfstats['max']
        ax.plot(xval,yval,self.main_clr,label='mean')
        ax.plot(xval,yval2,'--k',label='min')
        ax.plot(xval,yval3,'--k',label='max')
        ax.set_ylabel('mass loading (g m$^{-2}$)')
        ax.set_xlabel('Time')
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles,labels)
        ax.set_yscale('log')

    def sub_plot_variance(self,ax):
        xval = self.dfstats['date']
        yval = self.dfstats['variance']
        ax.plot(xval,yval,self.main_clr)
        ax.set_ylabel('variance')
        ax.set_xlabel('Time')

    def sub_plot_kurt(self,ax):
        xval = self.dfstats['date']
        yval = self.dfstats['kurtosis']
        ax.plot(xval,yval,self.main_clr)
        ax.set_ylabel('Kurtosis')
        ax.set_xlabel('Time')

    def sub_plot_skew(self,ax):
        xval = self.dfstats['date']
        yval = self.dfstats['skewness']
        ax.plot(xval,yval,self.main_clr)
        ax.set_ylabel('Skewness')
        ax.set_xlabel('Time')
       
       
    def sub_plot_num(self,ax):
        xval = self.dfstats['date']
        yval = self.dfstats['N']
        ax.plot(xval,yval,'-k',Linewidth=5,alpha=0.5)
        ax.set_ylabel('Number of Observations')
        ax.set_xlabel('Time')
       

 
    def sub_plot_mass(self,ax,yscale='ln'):
        import matplotlib.ticker as mtick
        def ticks(y,pos):
            pstr = '{:.0f}'.format(np.log(y))
            return r'$e^{}$'.format(pstr)

        yval = self.tmasslist
        xval = self.dlist

        #ys = self.spline(self.dtlist)
        #ax.plot(xval,ys,self.spline_clr, LineWidth=3,alpha=0.6)
        eee = 80
        #ax.plot(xval[20:eee],yval[20:eee],self.main_clr)
        ax.plot(xval,yval,self.main_clr)
        #ax.plot(xval,np.log(yval),self.main_clr)
        if yscale == 'ln': 
           ax.semilogy(basey=np.e) 
           ax.yaxis.set_major_formatter(mtick.FuncFormatter(ticks)) 
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

    def sub_plot_mer(self,ax, yscale='linear', spline=False):
        yval = self.mer
        xval = self.dlist
        if not spline:
            ax.plot(xval,yval,self.main_clr)
        else:
            mer = self.spline.derivative()
            ys = mer(self.dtlist)
            ax.plot(xval,ys*1e9,'r.', LineWidth=2,alpha=0.8)
            #ax.plot(xval,-1*ys*1e9,'k.', LineWidth=2,alpha=0.8)
        ax.set_ylabel('kg s$^{-1}$')
        ax.set_xlabel('Time')
        #ax.semilogy(basey=np.e) 
        if yscale=='log': ax.set_yscale('log')

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

import numpy as np
from scipy.stats import describe
import pandas as pd
import datetime
import logging
import matplotlib.pyplot as plt
import os
import seaborn as sns
import xarray as xr
from utilvolc import volcat
from utilhysplit.evaluation import hysplit_boxplots
from utilhysplit.evaluation import statmain
from utilvolc.qvainterface import DFInterface
from utilhysplit.plotutils import colormaker

logger = logging.getLogger(__name__)
#TODO - is find_area still used? also in ash_inverse.py

# 2023 18 July (amc) some of the values such as masstg are no longer being returned as lists.

class VolcatPlotDF(DFInterface):

    def __init__(self,edf):

        self._required = ['time','area','height','mass','radius',
                          'minmass','maxmass','minradius','maxradius',
                          'platform_ID','ID']

        if self.are_required_columns_present(edf,verbose=False):
            self._edf = edf
        else:
            self._edf = pd.DataFrame()

    @property
    def edf(self):
        return self._edf.copy()

    def required_columns(self):
        return self._required

 
    def save(self, name=None,overwrite=False):
        """
        save dataframe as csv file
        """
        if name: self.savename=name
        rval=True
        if not os.path.isfile(self.savename) or overwrite:
            self.edf.to_csv(self.savename, float_format="%.5e",index=False)
        else:
            self.add_csv(self.savename)
            logger.warning('CSV file exists. adding before saving')
            self.edf.to_csv(self.savename, float_format="%.5e", index=False)
            self.edf.to_csv(self.savename, index=False)
        return rval 

    # extra method
    def read(self, cname=None):
        """
        reads csv file that was previously saved
        """
        if not cname:
            cname = os.path.join(ndir, "Volcat.csv")
        dtp = {"time": "datetime64[ns]"}
        df = pd.read_csv(cname, sep=",", parse_dates=["time"])
        return df

    def add_csv(self, cname=None):
        """
        reads csv file that was previously saved and adds it to dataframe.
        """
        dftemp = self.read(cname)
        self.add_df(dftemp)

    @staticmethod
    def calc_mer(df):
        """
        INPUT
        df : dataframe with 'time' and 'mass' in Tg columns
        OUTPUT
        df : dataframe with mer column in kg/s 
        """
        df.sort_values(by='time',inplace=True)
        df['temp'] = df.time.diff()
        def toseconds(x):
            try: rval = x.seconds 
            except: rval=0
            return rval
        df['dt']=df['temp'].apply(toseconds)
        # mer in kg/s. multiply by 1e9 to convert from Tg to kg.
        df['dm']=1e9*df.mass.diff()              
        df['mer'] = df['dm']/df['dt'] 
        df['time_elapsed'] = df.dt.cumsum()
        return df 

    def add_df(self, df):
        """
        df : pandas DataFrame
        """
        # can utilize output of flist2eventdf
        # change all column names to lower case
        complete=True

        # check that required columns are in dataframe
        columns = df.columns
        complete = self.are_required_columns_present(df)

        if complete:
            if self._edf.empty:
                if isinstance(df, pd.core.frame.DataFrame):
                    self._edf = df
            else:
                if isinstance(df, pd.core.frame.DataFrame):
                    self._edf = pd.concat([self._edf, df])
            #self._edf.drop_duplicates()
        self._edf.sort_values(by='time',inplace=True)
        self._edf.drop_duplicates(inplace=True)
        return complete

    # extra method
        
    def are_required_columns_present(self,df,verbose=True):
        answer = True
        for req in self._required:
            if req not in df.columns:
               if verbose: print('WARNING, data frame does not contain required column {}'.format(req))
               answer=False
        return answer
 

class VolcatPlots:

    def __init__(self,dsetlist,volcano_name='Unknown',tdir='./'):
        """
        dsetlist can be from output of volcat.get_volcat_list function
        """
        self._vdf = VolcatPlotDF(edf=pd.DataFrame())
        self.volcano_name = volcano_name
        self.tdir = tdir
        self.set_plot_settings()      
 
    def make_savename(self):
        sname = '{}_vplots.csv'.format(self.volcano_name)
        return sname

    def save(self):
        self._vdf.savename =  self.make_savename()
        self._vdf.save(overwrite=False)

    def read(self):
        sname =  self.make_savename()
        self._vdf.add_csv(sname)

    def add_dsetlist(self,dsetlist):
        return -1
        # sort dset list by time
        #def ftime(x):
        #    tval = x.time.values
        #    if isinstance(tval,(list,np.ndarray)): return tval[0]
        #    return tval
        #dsetlist.sort(key=ftime)
        #self.dset = dsetlist
        #self.set_plot_settings()

    # returns the dataframe
    @property
    def vdf(self):
        return self._vdf.edf

    def empty(self):
        return False

    def make_arrays(self,das):
        masslist = []
        dlist = []
        # from maximum value of retrieval
        htlist = []
        alist = []
        tmasslist =[]
        arealist = []
        maxmass = []
        minmass = []
        radius = []
        minradius = []
        maxradius = []
        sid = []
        satellite = []
 
        # TO DO  - need to save this in some kind of structure?
        self.vmass =  []


        n_added=0
        for iii in np.arange(0,len(das)):
            checkid = das[iii].attrs['id']

            # if already exists in the dataframe then do not add.
            if not self.vdf.empty:
                if checkid in self.vdf.ID.unique(): 
                   continue
            n_added+=1
            try:
                vmass  = volcat.get_mass(das[iii],clip=True)
            except:
                print('cannot get mass',iii)
                continue
            try:
                vht  = volcat.get_height(das[iii],clip=True)
            except:
                print('cannot get height',iii)
                continue
            vrad  = volcat.get_radius(das[iii],clip=True)
            #dsetlist.append(das[iii])
            self.vmass.append(vmass)  
            maxmass.append(float(np.max(vmass).values))
            minmass.append(float(np.min(vmass).values))

            # this is no longer being returned as a list.
            # not sure if this is related to how the expansion of dims
            # were changed in volcat open_dataset and get_data functions.
            massval = das[iii].ash_mass_loading_total_mass.values
            #if isinstance(massval, list): massval = massval[0]
            masstg = float(massval)
            tmasslist.append(masstg)

            # date
            dlist.append(das[iii].time.values)

            # max height in detection
            htlist.append(float(np.max(vht)))

            # area of feature
            area = float(das[iii].feature_area.values)
            if isinstance(area,list): area = area[0]
            arealist.append(area)

            # mean effective radius
            radius.append(float(vrad.mean().values))
            minradius.append(float(np.min(vrad).values))
            maxradius.append(float(np.max(vrad).values))

            satellite.append(das[iii].attrs['platform_ID'])
            sid.append(das[iii].attrs['id'])
        #self.dset = dsetlist
        logger.info('Number added {} out of {}'.format(n_added,iii))
        if n_added >0:
            zzz = zip(dlist,arealist,htlist,tmasslist,radius,
                      minmass,maxmass,minradius,maxradius,satellite,sid)
            df = pd.DataFrame(list(zzz))
            
            df.columns = ['time','area','height','mass','radius',
                          'minmass','maxmass','minradius','maxradius',
                          'platform_ID','ID']
            self._vdf.add_df(df)


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
       
        dlist = self.vdf['time'] 
        clr = ["-m", "-r", "-b", "-c", "-g"]
        if isinstance(nums,list): vlist = self.vmass[nums[0]:nums[1]]
        else: vlist = self.vdf['mass']
        for jjj, volc in enumerate(vlist):
            volc = volc.values.flatten()
            volc = [x for x in volc if x > threshold]
            try:
                sts = describe(volc)
            except:
                print('no mass in file', jjj)
                continue
            vdata.append(volc)
            date.append(dlist[jjj])
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

    def volcat_cdf_plot(self, threshold=0, nums=None, skip=1):
        #step = 5
        clr = ["-m", "-r", "-b", "-c", "-g"]
        if isinstance(nums,list): vlist = self.vmass[nums[0]:nums[1]]
        else: vlist = self.vmass
        for jjj, volc in enumerate(vlist):
            if not jjj%skip==0: continue
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
            ax.step(sdata, y, clr[jjj % len(clr)], linewidth=lw)
        return ax

    def make_spline(self,s=20,vdf=None):
        import scipy.interpolate 
        if not isinstance(vdf,pd.DataFrame): 
            df = self.vdf
        else:
            df = vdf

        df = self._vdf.calc_mer(df)
        # this deals with if there are duplicate times in the file.
        # TODO - is using max appropriate?
        # did this partially to preserve the time column which is used for plotting.
        # sum and mean do not preserve the time column.
        df2 = df.groupby('time_elapsed').max()
        tmasslist = df2['mass']
        dtlist = df2.index
        dtlist = dtlist
        s = s / float(len(dtlist))
        #self.spline = scipy.interpolate.CubicSpline(self.dtlist,self.tmasslist)
        self.spline = scipy.interpolate.UnivariateSpline(dtlist,tmasslist,s=s)
        return df2

    def set_plot_settings(self):
        self.main_clr = 'b'
        self.spline_clr = 'r'
        self.sub_clrs = ['r','c']

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

    def plot_multiA(self,fignum=1,smooth=20,yscale='linear', bysensor=True):

        sns.set_style('whitegrid')
        fig = plt.figure(fignum,figsize=[10,10])

        ax1 = fig.add_subplot(2,2,1)
        ax2 = fig.add_subplot(2,2,2)
        ax3 = fig.add_subplot(2,2,3)
        ax4 = fig.add_subplot(2,2,4)
     
        if bysensor:  
            sensors = self.vdf.platform_ID.unique()
            cmaker = colormaker.ColorMaker('viridis',len(sensors),ctype='hex',transparency=None)
            clist = cmaker()
            for iii, sensor in enumerate(self.vdf.platform_ID.unique()):
                newdf = self.vdf[self.vdf['platform_ID']==sensor]
                self.main_clr='#' + clist[iii]
                self.sub_plot_mass(ax1,vdf=newdf,yscale=yscale)
                self.sub_plot_area(ax2,vdf=newdf)
                self.sub_plot_maxht(ax4,vdf=newdf)
                self.sub_plot_mer(ax3,vdf=newdf,smooth=smooth)
        else:
            self.sub_plot_mass(ax1,yscale=yscale)
            self.sub_plot_area(ax2)
            self.sub_plot_mer(ax3,smooth=smooth)
            self.sub_plot_maxht(ax4)
  
        fig.autofmt_xdate()
        plt.tight_layout()
        return fig

    def sub_plot_radius(self,ax):
        yval = self.vdf['radius']
        yval2 = self.vdf['minradius']
        yval3 = self.vdf['maxradius']
        xval = self.vdf['time']
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
        ax.plot(xval,yval,'-k',linewidth=5,alpha=0.5)
        ax.set_ylabel('Number of Observations')
        ax.set_xlabel('Time')
       

 
    def sub_plot_mass(self,ax,vdf=None,yscale='ln',label=None):
        import matplotlib.ticker as mtick
        def ticks(y,pos):
            pstr = '{:.0f}'.format(np.log(y))
            return r'$e^{}$'.format(pstr)

        if not isinstance(vdf,pd.DataFrame): 
            yval = self.vdf['mass']
            xval = self.vdf['time']
        else:
            yval = vdf['mass']
            xval = vdf['time']

        print(self.main_clr)
        ax.plot(xval,yval,color=self.main_clr,linestyle='-',marker='.')
        #ax.plot(xval,np.log(yval),self.main_clr)
        if yscale == 'ln': 
           ax.semilogy(base=np.e) 
           ax.yaxis.set_major_formatter(mtick.FuncFormatter(ticks)) 
        ax.set_ylabel('Total mass (Tg)')
        ax.set_xlabel('Time')
        plt.xticks(rotation=45)

    def sub_plot_max_mass(self,ax):
        yval = self.vdf['maxmass']
        xval = self.vdf['time']
        ax.plot(xval,yval,self.main_clr)
        ax.set_ylabel('Maximum Mass Loading (g m$^{-2}$)')
        ax.set_xlabel('Time')
     
    def sub_plot_min_mass(self,ax,both=False):
        yval = self.vdf['minmass']
        yval2 = self.vdf['maxmass']
        xval = self.vdf['time']
        ax.plot(xval,yval,self.main_clr, label='Minimum')
        if both: ax.plot(xval,yval2,self.sub_clrs[0], label='Maximum')
        ax.set_ylabel('Mass Loading (g m$^{-2}$)')
        ax.set_xlabel('Time')
        ax.set_yscale('log')
        ax.set_ylim([0.001,50])
        
    def sub_plot_area(self,ax,vdf=None,clr=-1):
        if not isinstance(vdf,pd.DataFrame): 
            yval = self.vdf['area']
            xval = self.vdf['time']
        else:
            yval = vdf['area']
            xval = vdf['time']
        if clr < 0:
            ax.plot(xval,yval,color=self.main_clr,linestyle='-',marker='.')
        else:
            ax.plot(xval,yval,self.sub_clrs[clr])
        ax.set_ylabel('Total Area (km$^2$)')
        ax.set_xlabel('Time')

    def sub_plot_mer(self,ax, vdf=None,yscale='linear', smooth=0):
        if not isinstance(vdf,pd.DataFrame): 
           df = self.vdf
        else:
           df = vdf
        df = self._vdf.calc_mer(df)
        xval = df['time']
        yval = df['mer']
        ax.plot(xval,yval,color=self.main_clr,linestyle='',marker='.')
        if smooth != 0:
            splinedf = self.make_spline(s=smooth,vdf=df)
            mer = self.spline.derivative()
            ys = mer(splinedf.index)
            ax.plot(splinedf.time,1e9*ys,self.main_clr, linestyle='-',linewidth=4,alpha=0.8)
            #ax.plot(xval,-1*ys*1e9,'k.', LineWidth=2,alpha=0.8)
        ax.set_ylabel('kg s$^{-1}$')
        ax.set_xlabel('Time')
        #ax.semilogy(basey=np.e) 
        if yscale=='log': ax.set_yscale('log')

    def sub_plot_maxht(self,ax,vdf=None):
        if not isinstance(vdf,pd.DataFrame): 
            yval = self.vdf['height']
            xval = self.vdf['time']
        else:
            yval = vdf['height']
            xval = vdf['time']

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

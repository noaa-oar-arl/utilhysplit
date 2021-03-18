import os
import datetime
import matplotlib.pyplot as plt
#import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import seaborn as sns
from utilhysplit.profile import callprofile
from utilhysplit.profile import MeteoProfile
from utilhysplit.plotutils import colormaker


def example_use(hdir, tdir):
    from utilhysplit.metfiles import MetFiles
    lat = 40
    lon = -85
    d1 = datetime.datetime(2020,7,22)
    duration = 6*24

    mf = MetFiles('{}/gdas1/gdas1.%b%y.week'.format(tdir))
    gfslist = mf.get_files(d1,duration)

    mf = MetFiles('{}/wrf27km/inst/%Y/wrfout_d01_%Y%m%d.ARL')
    wrflist = mf.get_files(d1,duration)

    mf = MetFiles('{}/gfs0p25/%Y%m%d_gfs0p25'.format(tdir))
    gfs025list = mf.get_files(d1,duration)

    mf = MetFiles('{}/nam12/%Y%m%d_nam12'.format(tdir))
    nam12 = mf.get_files(d1,duration)

    mf = MetFiles('{}/nams/%Y%m%d_hysplit.t00z.namsa'.format(tdir))
    nams = mf.get_files(d1,duration)

    mf = MetFiles('{}/hrrr/%Y%m%d_%H-??_hrrr'.format(tdir))
    hrrr = mf.get_files(d1,duration)

    metlist = [gfs025list,gfslist,nam12,nams,wrflist,hrrr]
    labels = ['gfs0p25','gdas1','nam12','nams','wrf','hrrr']

    d2 = d1 + datetime.timedelta(hours=duration)
    cp = CompareMetProfile(metlist, d1,d2,lat,lon,hdir,labels=labels)
    cp.call_and_read(nstop=200)


    cp.compare_variables()

    # see what surface and 3d variables are availalbe in each data set.
    print('Surface variables -----------------------------------')
    print(cp.df2)
    print('3D variables -----------------------------------')
    print(cp.df3)
  
    # look at plots of all surface variables. 
    cp.standard_surface_plots(plotall=True)

    # look at correlations between datasets  for surface variables
    cp.check_correlation('PBLH',plot=True,verbose=True)
  
    # look at differences between datasets for surface variables. 
    cp.check_diff('PBLH') 



class CompareMetProfile:
    """
    Use output from the profile program to compare different meteorological files.
   
    """

    def __init__(self,metlist,d1,d2,lat,lon,hdir,wdir='./',labels=None):
        """
        metlist : list of lists. 
                  The inner lists can be output from the MetFiles class get_files method.
                  [(metdir, metfilename),(metdir,metfilename),(metdir,metfilename)]
                  Each list corresponds to a different met data set.

        d1 : start date
        lat : latitude to input into profile
        lon : longitude to input into profile
        hdir : directory path to profile executable
        wdir : where to write output of profile
        labels : short names for the different meteorological files.
                 Should have same length as metlist.     
                  
        """

        self.metlist = metlist   # list of list of met filenames
        self.d1 = d1
        self.d2 = d2
        self.lat = lat
        self.lon = lon
        self.wdir = wdir
        if not labels: self.labels=[]
        else:  self.labels = labels
        self.dt = 1
        self.proflist = []
        self.pnamelist = []
        self.hdir = hdir  # directory for hysplit executables
        self.threedlist = []
        self.twodlist = []
        self.clist = [] 


    def get_colors(self,cmap='viridis'):
        nvals = len(self.proflist)
        cm = colormaker.ColorMaker(cmap,nvals,ctype='rgb')
        self.clist = cm()

    def generate_prof(self):
        # each profile has a label and a color associated with it.
        for iii, prof in enumerate(self.proflist):
            yield  prof, self.labels[iii], self.clist[iii] 


    def compare_variables(self):
        """
        Creates
        df3 : dataframe. rows are 3d variables in met files. Colums are met data sets. 
              values indicate whether met da ta set contains that variable.
        df2 : dataframe. rows are 2d variables in met files. Colums are met data sets. 
              values indicate whether met da ta set contains that variable.
        """
        threedlist = []
        twodlist = []
        dflist = []
        for iii, prof in enumerate(self.proflist):
            threedlist.extend(prof.var3d)
            twodlist.extend(prof.var2d)
        self.threedlist = list(set(threedlist))
        self.twodlist = list(set(twodlist))
        for prof, label, color in self.generate_prof():
            temp = []
            for var in self.threedlist:
                if var in prof.var3d:
                   temp.append(True)
                else: 
                   temp.append(False)
            df = pd.DataFrame(data=temp, index=self.threedlist, columns=[label])
            dflist.append(df)
        df3 = pd.concat(dflist, axis=1, join='outer')

        dflist = []
        for prof, label, color in self.generate_prof():
            temp = []
            for var in self.twodlist:
                if var in prof.var2d:
                   temp.append(True)
                else: 
                   temp.append(False)
            df = pd.DataFrame(data=temp, index=self.twodlist, columns=[label])
            dflist.append(df)
        df2 = pd.concat(dflist, axis=1, join='outer')
        self.df3 = df3
        self.df2 = df2

    def call_and_read(self,nstop=24,call=True, overwrite=False):
        dt = self.dt
        lat = self.lat
        lon = self.lon
        d1 = self.d1
        d2 = self.d2
      
        proflist = []
        pnamelist = []
        
        for iii, metnames in enumerate(self.metlist):
            for jjj, met_tuple in enumerate(metnames):
                mdir = met_tuple[0]
                met = met_tuple[1] 
                pname = '{}_profile.txt'.format(met)
                # don't call if file already exitst unless overwrite is true.
                test = overwrite or not os.path.isfile(os.path.join(self.wdir,pname)) 
                if call and test: 
                   print('calling profile')
                   callprofile(self.hdir, mdir,  met, lat, lon, dt, mfname=pname, nstop=nstop)
                pnamelist.append(pname)
                if jjj==0: prof = MeteoProfile(fname=pname, pdir=self.wdir, datesvalid=[d1,d2])
                else: prof.readnew(os.path.join(self.wdir,pname))
            proflist.append(prof)
            if not self.labels: labels.append(pname)    
        self.proflist = proflist
        self.pnamelist= pnamelist
        self.get_colors()

    def read(self):
        self.call_and_read(call=False)

    def plot_ts(self,yvarname):
        sns.set()
        sns.set_style('whitegrid')
        fig = plt.figure(1)
        ax1 = fig.add_subplot(1,1,1)
        for prof,label,color in self.generate_prof():
            print(label, '-----')
            xvar = prof.get_var(yvarname)
            if xvar!= -1: ax1.plot(prof.date_ra, xvar,marker='.', color=color, label=label)
        handles, labels = ax1.get_legend_handles_labels()
        plt.legend(handles,labels,loc='best',prop={'size':10})
        ax1.set_ylabel(yvarname)
        fig.autofmt_xdate()

    def plot(self,xvarname,yvarname):
        fig = plt.figure(1)
        ax1 = fig.add_subplot(1,1,1)
        for prof,label,color in self.generate_prof():
            xvar = prof.get_var(xvarname)
            yvar = prof.get_var(yvarname)
            ax1.plot(xvar, yvar, label=label)
            
    def standard_surface_plots(self,plotall=True):
        varlist = ['PBLH','USTR','SHTF','T02M','TPP1','TPP6','SHTF','DSWF','U10M','V10M']
        if plotall: varlist = self.twodlist
        for var in varlist:
            self.plot_ts(var)
            plt.show()

    def create_frame(self,var):
        dflist = []
        for prof,label,color in self.generate_prof():
            time = prof.date_ra
            xvar = prof.get_var(var)
            if xvar != -1: 
                dflist.append(pd.DataFrame(xvar,index=time,columns=[label]))
        return pd.concat(dflist,axis=1,join='outer')

    def check_correlation(self,var,thresh=0.90,plot=True,verbose=True):
        df = self.create_frame(var)
        if plot: heatmap(df)
        # calculates correlation between each column.
        corr = df.corr()        
        if verbose: print(corr)
        if np.any(corr[corr<thresh]):
           print('WARNING {} {}'.format(var, thresh))
        else:
           print('Passed {} {}'.format(var, thresh)) 

    def check_diff(self,var,thresh=0.90,plot=True,verbose=True):
        df = self.create_frame(var)
        columns = df.columns.values
        iii=1
        fig = plt.figure(1)
        ax1 = fig.add_subplot(1,1,1)
        for col in columns[1:]:
            temp = df[[columns[0], col]]
            diff = temp.diff(axis=1)
            diff = diff[col].dropna(axis=0)
            ax1.plot(diff,color=self.clist[iii], label=self.labels[iii])
            print('{} mean difference with {} is  {}'.format(col, columns[0], diff.mean()))
            iii+=1
        fig.autofmt_xdate()
        handles, labels = ax1.get_legend_handles_labels()
        plt.legend(handles,labels,loc='best',prop={'size':10})
        ax1.set_ylabel('difference for {}'.format(var))
        #plt.show()
        #return df

def heatmap(df):
    fs = 14
    plt.matshow(df.corr(),cmap='viridis')

    plt.xticks(range(df.shape[1]), df.columns, fontsize=14,rotation=90)
    plt.gca().xaxis.tick_bottom()
    plt.yticks(range(df.shape[1]), df.columns, fontsize=14)
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=14)
    plt.title("Feature Correlation Heatmap", fontsize=fs) 

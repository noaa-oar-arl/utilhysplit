import os
import datetime
import matplotlib.pyplot as plt
#import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import seaborn as sns
from utilhysplit.hysplit_profile import callprofile
from utilhysplit.hysplit_profile import MeteoProfile
from utilhysplit.plotutils import colormaker

def test_wind_direction():
    """
    checks that wind direction function is working.
    """

    u = [10]
    v = [0]
    wdir, wspd = wind_direction(v,u)
    print('270 degrees', wdir[0])
    print('10 m/s', wspd[0])

    u = [-10]
    v = [0]
    wdir, wspd = wind_direction(v,u)
    print('90 degrees', wdir[0])
    print('10 m/s', wspd[0])
   
    u = [0]
    v = [5]
    wdir, wspd = wind_direction(v,u)
    print('180 degrees', wdir[0])
    print('5 m/s', wspd[0])

    u = [0]
    v = [-5]
    wdir, wspd = wind_direction(v,u)
    print('360 degrees', wdir[0])
    print('5 m/s', wspd[0])
 
    u = [-5]
    v = [-5]
    wdir, wspd = wind_direction(v,u)
    print('45 degrees', wdir[0])
    print('{} m/s'.format(50**0.5), wspd[0])

    u = [-5]
    v = [5]
    wdir, wspd = wind_direction(v,u)
    print('135 degrees', wdir[0])
    print('{} m/s'.format(50**0.5), wspd[0])

    u = [5]
    v = [-5]
    wdir, wspd = wind_direction(v,u)
    print('315 degrees', wdir[0])
    print('{} m/s'.format(50**0.5), wspd[0])


def wind_direction(vwind, uwind):
        #vwind is magnitude of wind going from south to north
        #uwind is magnitude of wind going from West to East
        #print(len(vwind), len(uwind))
        
        uwind = np.array(uwind)
        vwind = np.array(vwind)
        zrs = np.where(uwind==0)
        uwind[zrs] = 0.001
        wind_dir = np.arctan(vwind/uwind)*180/np.pi
        upos = np.where(uwind >= 0)
        wind_dir[upos]= 270 - wind_dir[upos] 
        uneg = np.where(uwind < 0)
        wind_dir[uneg] = 90 - wind_dir[uneg] 

        wspd = (uwind**2 + vwind**2)**0.5
        #print(wind_dir.shape , wind_dir[2])
        #print('U' , uwind.shape, uwind[2])
        #print('V' , vwind.shape, vwind[2])
        return wind_dir, wspd

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
        """
        creates proflist
        """

        dt = self.dt
        lat = self.lat
        lon = self.lon
        d1 = self.d1
        d2 = self.d2
      
        proflist = []
        pnamelist = []
        
        for iii, metnames in enumerate(self.metlist):
            print(iii, metnames)
            for jjj, met_tuple in enumerate(metnames):
                mdir = met_tuple[0]
                met = met_tuple[1] 
                pname = '{}_profile.txt'.format(met.replace('/',''))
                # don't call if file already exitst unless overwrite is true.
                test = overwrite or not os.path.isfile(os.path.join(self.wdir,pname)) 
                if call and test: 
                   print('calling profile',self.hdir,mdir,met)
                   callprofile(self.hdir, mdir,  met, lat, lon, dt, mfname=pname, nstop=nstop)
                pnamelist.append(pname)
                if jjj==0: prof = MeteoProfile(fname=pname, pdir=self.wdir, datesvalid=[d1,d2])
                else: prof.readnew(os.path.join(self.wdir,pname))
                print(pname, met, mdir)
            proflist.append(prof)
            if not self.labels: labels.append(pname)    
        self.proflist = proflist
        self.pnamelist= pnamelist
        self.get_colors()

    def read(self):
        self.call_and_read(call=False)

    def set_ref(self,ref):
        self.ref = ref

    def plot_ts2(self,yvarname,fignum=1):
        sns.set()
        sns.set_style('whitegrid')
        fig = plt.figure(fignum)
        ax1 = fig.add_subplot(1,1,1)
        for prof,label,color in self.generate_prof():
            print('THIS IS', label)
            fig = plt.figure(fignum)
            ax1 = fig.add_subplot(1,1,1)
            #print(label, '-----')
            xvar = prof.get_var(yvarname)
            if label == self.ref: lw=5
            else: lw=1
            if xvar!= -1: ax1.plot(prof.date_ra, xvar,marker='.', color=color, label=label,LineWidth=lw)
           
            handles, labels = ax1.get_legend_handles_labels()
            plt.legend(handles,labels,loc='best',prop={'size':10})
            ax1.set_ylabel(yvarname)
            fig.autofmt_xdate()
            plt.show()
        return ax1

    def plot_ts(self,yvarname,fignum=1):
        sns.set()
        sns.set_style('whitegrid')
        fig = plt.figure(fignum)
        ax1 = fig.add_subplot(1,1,1)
        for prof,label,color in self.generate_prof():
            #print(label, '-----')
            xvar = prof.get_var(yvarname)
            if label == self.ref: lw=5
            else: lw=1
            if xvar!= -1: ax1.plot(prof.date_ra, xvar,marker='.', color=color, label=label,LineWidth=lw)
        handles, labels = ax1.get_legend_handles_labels()
        plt.legend(handles,labels,loc='best',prop={'size':10})
        ax1.set_ylabel(yvarname)
        fig.autofmt_xdate()
        return ax1

    def plot(self,xvarname,yvarname):
        fig = plt.figure(1)
        ax1 = fig.add_subplot(1,1,1)
        for prof,label,color in self.generate_prof():
            xvar = prof.get_var(xvarname)
            yvar = prof.get_var(yvarname)
            ax1.plot(xvar, yvar, label=label)
            
    def standard_surface_plots(self,plotall=True,varlist=None,fignum=1,together=True):
        print(type(varlist),varlist)
        if isinstance(varlist,list):
            varlist = varlist
        else:
            varlist = ['PBLH','USTR','SHTF','T02M','TPP1','TPP6','SHTF','DSWF','U10M','V10M']
        print(type(varlist))
        #varlist = ['PBLH','USTR','SHTF','T02M','TPP1','TPP6','SHTF','DSWF','U10M']
        if plotall: varlist = self.twodlist
        varlist.sort()
        for iii, var in enumerate(varlist):
            print(var)
            if together:
                ax = self.plot_ts(var,fignum=fignum+iii)
            else:
                ax = self.plot_ts2(var,fignum=fignum+iii)
            #plt.show()
        return ax

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

    def windspd(self,date,legend=False):
        fig = plt.figure(10)
        ax = fig.add_subplot(1,1,1)
        fig1 = plt.figure(11)
        ax1 = fig1.add_subplot(1,1,1)
        #iii=0
        for prof, label, color in self.generate_prof():
            if 'VWND_rot' not in prof.var3d: continue
            df = prof.get_3dvar_df()
            uwnd = df['UWND_rot'].values
            vwnd = df['VWND_rot'].values
            wdir, wspd = wind_direction(vwnd,uwnd)
            df['wdir'] = wdir
            df['wspd'] = wspd
            df = df[df['time'] == date]
            df.set_index('PRES1',inplace=True)
            if label == self.ref: 
               lw=1
               ms=10
               alpha=1
            else: 
               lw=1
               ms=5
               alpha=0.5
            dfdir = df['wdir']
            dfspd = df['wspd']
            try:
               ax.plot(dfdir,df.index,color=color,marker='.',
                       label=label,linewidth=lw,alpha=alpha,
                       MarkerSize=ms)
            except:
                print('failed {}'.format(label))
            try:
               ax1.plot(dfspd,df.index,color=color,marker='.',
                       label=label,linewidth=lw,alpha=alpha,
                       MarkerSize=ms)
            except:
                print('failed {}'.format(label))
            #if iii> 1: break
            #iii+=1
        handles, labels = ax.get_legend_handles_labels()
        if legend: plt.legend(handles,labels,loc='best',prop={'size':10})
        ax.set_ylabel('Level (mb)')
        ax.invert_yaxis()
        ax1.invert_yaxis()
        return ax, ax1

    def check_3d(self,var,date,legend=False):
        fig = plt.figure(10)
        ax = fig.add_subplot(1,1,1)
        #iii=0
        for prof, label, color in self.generate_prof():
            if var not in prof.var3d: continue
            df = prof.get_3dvar_df()
            df = df[df['time'] == date]
            df.set_index('PRES1',inplace=True)
            if label == self.ref: 
               lw=8
               alpha=1
            else: 
               lw=1
               alpha=0.5
            if var == 'WWND' and 'DIFW' in df.columns.values:
               print('adding difw to wwnd', label)
               df['NEW'] = df['WWND'] + df['DIFW']
               df = df['NEW']
            else:
               df = df[var]
            try:
               ax.plot(df,df.index,color=color,marker='.',label=label,linewidth=lw,alpha=alpha)
            except:
                print('failed {}'.format(label))
            #if iii> 1: break
            #iii+=1
        handles, labels = ax.get_legend_handles_labels()
        if legend: plt.legend(handles,labels,loc='best',prop={'size':10})
        ax.set_ylabel('Level (mb)'.format(var),fontsize=20)
        if(var=='VWND_rot'): xlabel='V (m/s)'
        elif(var=='UWND_rot'): xlabel='U (m/s)'
        elif(var=='WWND'): xlabel='W (m/s)'
        elif(var=='RELH'): xlabel='Relative Humidity (%)'
        elif(var=='TEMP'): xlabel='Temperature ($^o$C)'
        else: xlabel=var
        ax.set_xlabel('{}'.format(xlabel),fontsize=20)
        ax.invert_yaxis()
        return ax

def heatmap(df):
    fs = 14
    plt.matshow(df.corr(),cmap='viridis')

    plt.xticks(range(df.shape[1]), df.columns, fontsize=14,rotation=90)
    plt.gca().xaxis.tick_bottom()
    plt.yticks(range(df.shape[1]), df.columns, fontsize=14)
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=14)
    plt.title("Feature Correlation Heatmap", fontsize=fs) 

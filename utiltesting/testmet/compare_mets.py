import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from utilhysplit.profile import callprofile
from utilhysplit.profile import MeteoProfile
from utilhysplit import colormaker

def windspd(x1, x2):
    x1 = np.array(x1)
    x2 = np.array(x2)

    spd = (x1**2 + x2**2)**0.5
    return spd

class CompareMetProfile:

    def __init__(self,metlist,mdirlist,d1,d2,lat,lon,hdir,wdir='./',labels=None):
        self.metlist = metlist   # list of met filenames
        self.mdirlist = mdirlist # list of met directories.
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

    def call_and_read(self,nstop=24,call=True):
        dt = self.dt
        lat = self.lat
        lon = self.lon
        d1 = self.d1
        d2 = self.d2
      
        proflist = []
        pnamelist = []
        
        for iii, met in enumerate(self.metlist):
            pname = '{}_profile.txt'.format(met)
            if call: 
               callprofile(self.hdir, self.mdirlist[iii], met, lat, lon, dt, mfname=pname, nstop=nstop)
            pnamelist.append(pname)
            prof = MeteoProfile(fname=pname, pdir=self.wdir, datesvalid=[d1,d2])
            proflist.append(prof)
            if not self.labels: labels.append(pname)    
        self.proflist = proflist
        self.pnamelist= pnamelist
        self.get_colors()

    def read(self):
        self.call_and_read(call=False)

    def plot_ts(self,yvarname):
        sns.set()
        fig = plt.figure(1)
        ax1 = fig.add_subplot(1,1,1)
        for prof,label,color in self.generate_prof():
            print(label, '-----')
            xvar = prof.get_var(yvarname)
            if xvar!= -1: ax1.plot(prof.date_ra, xvar,marker='.', color=color, label=label)
        handles, labels = ax1.get_legend_handles_labels()
        plt.legend(handles,labels,loc='best',prop={'size':10})
        ax1.set_ylabel(yvarname)

    def plot(self,xvarname,yvarname):
        fig = plt.figure(1)
        ax1 = fig.add_subplot(1,1,1)
        for prof,label,color in self.generate_prof():
            xvar = prof.get_var(xvarname)
            yvar = prof.get_var(yvarname)
            ax1.plot(xvar, yvar, label=label)
            
    def standard_plots(self):
        varlist = ['PBLH','USTR','SHTF']
        for var in self.twodlist:
            print(var)
            self.plot_ts(var)
            plt.show()


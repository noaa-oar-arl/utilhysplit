# sonde data for Turrialba from Gary Morris
import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import glob
from scipy.signal import find_peaks
from utilhysplit.evaluation import statmain
from utilhysplit import hcontrol
from utilhysplit import metfiles
from utilhysplit import hysplit_profile

def get_peaks(olist,plot=True):
    """
    utilizes the scipy.signal find_peaks method.
    more applicable than the get_peaksB method.
    """
    vlist = olist.dropna()
    rra = vlist.values.flatten()
    print(type(rra),rra.size,rra.shape)
    height = 5
    distance = 50
    prominence = 5
    width = [1,500]
    peaks = find_peaks(rra,height=height, distance=distance,prominence=prominence,width=width)
    #vlist = vlist.reset_index()
    tsp = vlist.iloc[peaks[0]]
    if plot:
       plt.plot(olist.values, olist.index.tolist(),'--k',linewidth=1)
       plt.plot(tsp.values, tsp.index.tolist(), 'ro',MarkerSize=10,alpha=0.5)
    return peaks
    




def get_peaksB(vlist, pvalues=[0.99,1], plot=True):
    """
    Creates a CDF from the values in the input series and
    returns a series with same length in which only the values with
    CDF between pvalue[0] and pvalue[1] are not NaN.
    vlist is a pandas series.
    Returns
    tsp: pandas series
         all values which are below the
 
    """
    olist = vlist.copy()
    valA = statmain.probof(vlist.values,pvalues[0])
    valB = statmain.probof(vlist.values,pvalues[1])
    print(valA, valB)
    tsp = vlist[vlist>=valA]
    tsp = tsp[tsp<=valB]
    if plot:
       plt.plot(vlist.index.tolist(),vlist.values,'--k',linewidth=1)
       plt.plot(tsp.index.tolist(), tsp.values,'b.',MarkerSize=5)
    return tsp 





def get_sondes(tdir):
    fnames =  glob.glob(tdir + '/*dat')
    return fnames


class SondeList:

    def __init__(self,name='Sondes'):
        self.name = name
        self.slist = []  # List of Sonde objects
        self.tlist = []  # List of launch dates
        self.totalSO2 = [] # list of integrated SO2 up to 15 km
        self.totalSO2_label =  ''


    def get_peaks(self, plot=True):
        for sss in self.slist:
            sss.get_peaks(plot=plot)
            if plot: 
               plt.show() 
     
    def match_profile2sonde(self,wdir):
        for iii,sss in enumerate(self.slist):
            pname = 'profile.{}'.format(iii)
            p2 = profile.MeteoProfile(pname,wdir)
            df = p2.get_3dvar_df()
            tdiff = []
            tlist = []
            for ttt in df.time.unique():
                ttt = pd.to_datetime(ttt)
                if ttt > sss.launchdate:
                   diff = ttt-sss.launchdate
                else:  
                   diff = sss.launchdate-ttt
                tdiff.append(diff.seconds/3600 + diff.days*24)
                tlist.append(ttt)
            jjj = tdiff.index(np.min(tdiff))
            tvalue = tlist[jjj]
            dft = df[df['time'] == tvalue]
            sss.add_profile(dft) 
            #print(list(zip(tdiff, tlist))) 
            #print('Matched', tvalue, sss.launchdate)            
            #print('\n\n') 

    def profile_script(self,scriptname,hdir):
        rstr = ''
        for iii, sss in enumerate(self.slist):
            rrr = sss.make_profile_str(hdir,suffix=iii)         
            rstr += rrr
            rstr += '\n\n'
        with open(scriptname, 'w') as fid:
             fid.write(rstr)
        return rstr
                      
    def make_control(self,wdir,duration=-24):
        for iii, sss in enumerate(self.slist):
            sss.make_control(wdir,suffix=iii,duration=duration)
         

    def plot_totalSO2(self):
        d1 = self.tlist[0]
        d2 = self.tlist[-1]
        fig = plt.figure(figsize=[20,8])
        ax = fig.add_subplot(1,1,1)
        plt.plot(self.tlist,self.totalSO2, '--k.')
        plt.plot([d1,d2],[0,0],'-r')
        fs = 18
        ax.set_ylabel(self.totalSO2_label, fontsize=fs)
        ax.set_xlabel('Launch Date', fontsize=fs)
 
    def add_data(self,sonde):
        self.slist.append(sonde)
        self.tlist.append(sonde.launchdate)
        ky=None
        for key in sonde.ihash.keys():
            if 'integrated so2' in key.lower(): ky = key
        if ky: 
           try: stotal = float(sonde.ihash[ky])
           except: stotal = -10
           self.totalSO2.append(stotal)
           self.totalSO2_label = ky 

class Sonde:
    def __init__(self,fname):
        self.fname = fname
        self.ihash = {}
        self.units = []    # list of strings with units for each column.
        self.df = pd.DataFrame() # dataframe
        self.columns = ['Time','Press','Alt','Temp','RH','O3a','O3b','O3c','wdir','wspd','Tpump','I O3','SO2','GPS Lon','GPS Lat', 'GPS Alt']
        self.launchdate = None
        self.latitude = None
        self.longitude = None
        iii = self.read_header()
        try:
            self.read_data(iii+1)
        except:
            print('Warning could not read data', self.fname)
        self.read_header()
        self.profile = pd.DataFrame()

        self.get_peaks()


    def compare_profile(self):
        self.get_peaks(yval='Press', plot=False)
        dstr = self.launchdate.strftime("%d%b%Yy%Hh%Mm")        
        
        if 'uwnd' not in self.df.columns:
           self.add_wind_components()

        fig = plt.figure(1,figsize=[8,20])
        ax = fig.add_subplot(3,1,1)
        ax.plot(self.profile.TEMP, self.profile.PRES1, '-c.')
        ax.plot(self.df.Temp, self.df.Press, '-k.',markersize=1,alpha=0.5)
        ax.invert_yaxis()
        ax.set_ylabel('Pressure (hPa)')
        ax.set_xlabel('Temperature ($^o$C')
        #plt.savefig(dstr + 'Temp.png')
        plt.title(self.launchdate.strftime("%d %b %Y %H:%M"))        

        wspd = self.get_ave('wspd')
        ax2 = fig.add_subplot(3,1,2)
        ax2.plot(self.profile.wspd, self.profile.PRES1, '-co')
        ax2.plot(self.df.wspd, self.df.Press, 'k.',markersize=1,alpha=0.5)
        ax2.plot(wspd.ave, wspd.index.to_list(), '-k',linewidth=3,alpha=0.8)
        ax2.invert_yaxis()
        ax2.set_ylabel('Pressure (hPa)')
        ax2.set_xlabel('Wind Speed (m/s)')
        for yval in self.peaklist.index.to_list():
            xval = [np.min(self.df.wspd), np.max(self.df.wspd)] 
            ax2.plot(xval, [yval,yval], '-r', linewidth=5, alpha = 0.5)


        #plt.savefig(dstr + 'WSPD.png')
        #plt.show() 

        uwnd = self.get_ave('uwnd')
        vwnd = self.get_ave('vwnd')
        temp = pd.concat([uwnd,vwnd],axis=1)
        temp.columns = ['uwnd','uave','vwnd','vave']
        temp['wdir'] = temp.apply(lambda row: profile.wind_direction(row.vave, row.uave),axis=1)      
        #fig = plt.figure(1,figsize=[5,10])
        ax3 = fig.add_subplot(3,1,3)
        ax3.plot(self.profile.wdir, self.profile.PRES1, '-co')
        ax3.plot(self.df.wdir, self.df.Press, 'k.',markersize=1,alpha=0.5)
        ax3.plot(temp.wdir, temp.index.to_list(), 'ko',markersize=2,alpha=0.8)
        ax3.invert_yaxis()
        ax3.set_ylabel('Pressure (hPa)')
        ax3.set_xlabel('Wind Direction (degrees)')
        
        ax3.plot([90,90],[880,0],'-g', linewidth=5, alpha=0.5)
 
        for yval in self.peaklist.index.to_list():
            xval = [np.min(temp.wdir), np.max(temp.wdir)] 
            ax3.plot(xval, [yval,yval], '-r', linewidth=5, alpha = 0.5)
        plt.savefig(dstr + 'WDIR.png')
        plt.show() 
   
 

    def get_ave(self, var):
        temp = self.df[['Press',var]]
        temp = temp.set_index('Press')
        temp['ave'] = temp.rolling(window=100,center=True,min_periods=10).mean()
        return temp

 
    def add_wind_components(self):
        self.df['wnd'] = self.df.apply(lambda row: profile.wind_components(row.wdir, row.wspd),axis=1)
        self.df['uwnd'] = self.df.apply(lambda row: row['wnd'][1],axis=1)
        self.df['vwnd'] = self.df.apply(lambda row: row['wnd'][0],axis=1)
        self.df.drop('wnd', axis=1,inplace=True)     
 
    def add_profile(self, profiledf):
        """
        profiledf : pandas DataFrame 
        """
        self.profile = profiledf
        self.profile['wspd'] = self.profile.apply(lambda row: (row['UWND']*row['UWND'] + row['VWND']*row['VWND'])**0.5, axis=1)
        self.profile['wdir'] = self.profile.apply(lambda row: profile.wind_direction(row['VWND_rot'],row['UWND_rot']),axis=1)

    def get_peaks(self,val='SO2',yval='Alt',plot=True):
        stemp = self.df
        stemp = stemp[[yval,val]]
        stemp = stemp.set_index(yval)
        self.peaks = get_peaks(stemp, plot=plot)
        self.peaklist = stemp.iloc[self.peaks[0]]
        if plot:
           ax = plt.gca()
           ax.set_ylabel('Altitude (km a.s.l.)')
           unit = 'ppbv'           
           ax.set_xlabel('{} ({})'.format(val, unit))

           tropo = float(self.ihash['Tropopause height (km)'])
           xi = np.min(stemp.dropna().values)
           xf = np.max(stemp.dropna().values)
           print('Tropopause at', tropo)
           if tropo < 15:
               ax.plot([xi,xf],[tropo,tropo],'--b',label='Tropopause')
           mlh = float(self.ihash['Mixed layer height (km)'])
           print('{} {} MLH {}'.format(xi,xf, mlh))
           ax.plot([xi,xf],[mlh,mlh],'-y',label='MLH', linewidth=5, alpha=0.5)
           ax.plot([0,0],[0,15],'-g',linewidth=5,alpha=0.5)
           plt.title(self.launchdate.strftime("%d %b %Y %H:%M"))        
           dstr = self.launchdate.strftime("%d%b%Yy%Hh%Mm")        
           plt.savefig(dstr + 'Profile.png')
    
    def read_header(self):
        ihash = {}
        with open(self.fname,'r') as fid:
             iii=0
             for line in fid.readlines():
                 temp = line.split(':')
                 if 'Time' in line and 'Press' in line and 'Alt' in line: break  
                 ihash[temp[0].strip()] = temp[-1].strip()
                 if 'latitude' in temp[0].lower(): self.latitude = float(temp[1])
                 if 'longitude' in temp[0].lower(): self.longitude = float(temp[1])
                 if 'Launch Date' in temp[0]:
                     dstr = temp[1].strip()
                     sdate = datetime.datetime.strptime(dstr, "%Y%m%d") 
                 if 'Launch Time' in temp[0]:
                     hour = int(temp[1])
                     minute = int(temp[1])
                     second = int(temp[2])
                     self.launchdate = datetime.datetime(sdate.year, sdate.month, sdate.day, hour,minute,second)  
                 iii+=1
        self.ihash = ihash
        return iii 

    def read_data(self,skr):
        fname = self.fname
        df = pd.read_csv(fname, header=0, skiprows=skr, sep = '\s+')
        units = df.columns
        units = list(zip(units,self.columns))
        if len(df.columns)==len(self.columns):
            df.columns = self.columns
        else:
            iii = len(df.columns) 
            df.columns = self.columns[0:iii]
            #print('could not replace columns')
            #print(len(df.columns), df.columns)
        # replace missing or bad values with NaN
        df = df.where(df!=9000,np.NaN)
        self.df = df

    def make_control(self,wdir,suffix,duration=-24):
        control = hcontrol.HycsControl(rtype='trajectory')
        td = self.launchdate
        newdate = datetime.datetime(td.year,td.month,td.day,td.hour,td.minute)
        #if td.minute > 30:
        #    newdate += datetime.timedelta(hours=1)
        control.add_sdate(newdate)
        control.add_ztop(18000)
        control.add_duration(duration)
        for peak in self.peaklist.index.tolist():
            alt = peak*1000
            control.add_location((self.latitude, self.longitude), alt=alt)
        control.rename('CONTROL.{}'.format(suffix),working_directory=wdir)
        metstr = '/pub/archives/gdas0p5/%Y%m%d_gdas0p5'
        mf = metfiles.MetFiles(metstr,hours=24)
        mlist = mf.get_files(newdate,duration)
        control.outfile = 'tdump.{}'.format(suffix)
        for mmm in mlist:
            control.add_metfile(mmm[0],mmm[1])
        control.write()
        return control


    def make_profile_str(self,hdir, suffix):
        td = self.launchdate
        newdate = datetime.datetime(td.year,td.month,td.day,td.hour,td.minute)
        metstr = '/pub/archives/gdas0p5/%Y%m%d_gdas0p5'
        mf = metfiles.MetFiles(metstr,hours=24)
        mlist = mf.get_files(newdate,1)
        mlist = mlist[0]
        print(mlist)
        lat = self.latitude
        lon = self.longitude
        rstr = hdir + 'profile -d{} -f{} -x{} -y{} -t1'.format(mlist[0],mlist[1],lon,lat)
        rstr += '\n'
        rstr += 'mv profile.txt profile.{}'.format(suffix)
        return rstr
 

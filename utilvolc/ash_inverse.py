import sys
import os
import datetime
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import numpy as np
import numpy.ma as ma
import pandas as pd
import seaborn as sns
from utilvolc import volcat
from monetio.models import hysplit

def testcase(tdir,vdir):
    fname = 'xrfile.444.nc'
    bezyloc = [160.587,55.978]
    vid = 'v300250'
    inverse = InverseAsh(tdir,fname,vdir,vid)
    return inverse

# get cdump output.
# pick time period.
# cdump output for that time period. 


class InverseOutDat:

    def __init__(self,wdir,fname):
        # out.dat has estimated release rates in same order as tcm columns.
        # out2.dat has observed(2nd column) and modeled(3rd column) mass loadings.
        self.wdir = wdir
        
    def read_out(self,name):
        wdir = self.wdir
        df = pd.read_csv(os.path.join(wdir,name),sep='\s+',header=None)
        return df

    def get_conc(self, name='out2.dat'):
        df = self.read_out(name)
        df.columns = ['index', 'observed', 'model']
        plt.plot(df['observed'],df['model'],'k.',MarkerSize=3)
        ax = plt.gca()
        ax.set_xlabel('observed')
        ax.set_ylabel('model')
        nval = np.max(df['observed'])
        # plot 1:1 line
        plt.plot([0,nval],[0,nval],'--b',LineWidth=1)
        return df


def get_sourcehash(wdir,configfile):
    from utilvolc.ashapp import ashinverse
    #from utilvolc.ashapp.runhelper import JobSetUp
    from utilvolc.ashapp.runhelper import make_inputs_from_file
    setup = make_inputs_from_file(wdir,configfile)
    setup.add_inverse_params()
    sourcehash = ashinverse.inverse_get_suffix_list(setup.inp)
    return sourcehash 


class InverseAsh:

    def __init__(self, tdir, fname,vdir,vid,configfile=None):
        """
        configfile : full path to configuration file.
        """

        self.tdir = tdir   # directory for hysplit output
        self.fname = fname # name of hysplit output

        self.vdir = vdir   # directory for volcat data
        self.vid = vid   # volcano id.

        # keep volcat arrays for different averaging times. 
        self.volcat_hash = {}
        self.volcat_avg_hash = {}
        self.cdump_hash = {}

        # hysplit output. xarray. 
        cdump = xr.open_dataset(os.path.join(tdir,fname))
        print('opened')
        # turn dataset into dataarray
        temp = list(cdump.keys())
        cdump = cdump[temp[0]]

        # get rid of source dimension (for now)
        cdump = cdump.isel(source=0)

        # the ens dimension holds is key to what emission source was used.
        # the sourcehash is a dictionary
        # key is the ensemble number
        # values is another dictionary with
        # sdate: begin emission
        # edate: end emission
        # bottom : lower height of emission
        # top : upper height of emission.
        if configfile:
            self.sourcehash = get_sourcehash(configfile)

        # get the mass loading.
        # TODO later may want to exclude certain levels.
        self.cdump = hysplit.hysp_massload(cdump)

    def get_volcat(self, daterange, verbose=False):
        vdir = self.vdir
        vid = self.vid
        tii = self.time_index(daterange[0])
        if tii not in self.volcat_hash.keys(): 
            das = volcat.get_volcat_list(vdir,daterange=daterange,vid=vid) 
            self.volcat_hash[tii] = das
        else:
            das = self.volcat_hash[tii]
        return das

    def clip(self,dummy):
        # clip ra according to where dummy has 0's.
        aaa = np.where(dummy != 0)
        a1 = np.min(aaa[0])
        a2 = np.max(aaa[0])
        b1 = np.min(aaa[1])
        b2 = np.max(aaa[1])
        return a1,a2,b1,b2

    def time_index(self,time):
        timelist = [pd.to_datetime(x) for x in self.cdump.time.values]
        try:
           iii = timelist.index(time)
        except:
           iii = None
        return iii 

    def make_tcm(self,tiilist, remove_cols=True):
        # header line indicates release times for each column.
        # I think header can just be a dummy.
        # one column for each release time/location.
        # one row for each measurement location.

        #         ens1  ens2 ens3 ens4 ens5 ens6 ens7 .... Obs
        # (x1,y1)
        # (x2,y1)
        # (x3,y1)
        # .
        # .
        # .
        # (x1,y2)

        # column will be from the ensemble dimension.
        # measurement is from the lat/lon dimension.

        # last column is the value of the observation.

        # number of rows corresponds to number of points where there is an observation
        # only include observations above 0.
        # or include where either observation OR model above 0.
        # or include only where model above 0.

        # right now only one time period.
        tii = tiilist[0]
        cdump = self.cdump_hash[tii]
        avg = self.volcat_avg_hash[tii] 


        s1 = avg.shape[0]*avg.shape[1] 

        model = cdump.stack(pos=['y','x'])
        model = model.transpose('pos','ens')      

        # remove columns which have no contribution at all. 
        if remove_cols:
            model = model.where(model>0)
            model = model.dropna(dim = 'ens', how='all')
 
        model_lat = model.latitude.values.reshape(s1,1)
        model_lon = model.longitude.values.reshape(s1,1)
        columns = model.ens.values

        model = model.values
       
        volc = avg.values.reshape(s1,1)
        volc_lat = avg.latitude.values.reshape(s1,1)
        volc_lon = avg.longitude.values.reshape(s1,1)

        tcm = np.concatenate([model,volc],axis=1)
        if not np.all(volc_lon == model_lon):
           print('WARNING, model and observed locations in tcm not matching')
        if not np.all(volc_lat == model_lat):
           print('WARNING, model and observed locations in tcm not matching')
       
        self.tcm = tcm
        self.tcm_lat = model_lat
        self.tcm_lon = model_lon
        # this contains the keys that can be matched in the sourcehash attribute.
        self.tcm_columns =  columns
        return tcm, model_lat, model_lon, columns

    def plot_tcm(self):
        plt.pcolormesh(np.log10(self.tcm))


    def make_fake(self):
        tcm_real = self.tcm.copy()
       
        nra = []
        nra.append([1,0,0,0,0,0,2])
        nra.append([0,1,0,0,0,0,3])
        nra.append([0,0,1,0,0,0,10])
        nra.append([0,0,0,1,0,0,6])
        nra.append([0,0,0,0,1,0,5])
        nra.append([0,0,0,0,0,1,4])
        self.tcm = np.array(nra)
        self.write_tcm('test_tcm.txt')
        self.tcm = tcm_real

    def write_tcm(self, name):
        astr = ''
        sep = ' '
        hstr = '' # header string
        print(self.tcm.shape)\
        # this is number of columns minus 1.
        print('N_ctrl {}'.format(self.tcm.shape[1]-1))
        print('N_ctrl {}'.format(self.tcm.shape[1]-1))
            
        for iii, line in enumerate(self.tcm):
            for jjj, val in enumerate(line):
                if iii==0:
                   hstr += '43637.750' + sep
                if not np.isnan(val): astr += '{:1.5e}'.format(val)
                else: astr += '{:1.4e}'.format(0)
                astr += sep
            astr += '\n '
            #print(astr)
            #print(line[-1])
            #sys.exit()
            if iii==0: hstr += '\n'
        with open(name, 'w') as fid:
            fid.write(hstr + astr)
        return hstr + astr 

    def make_outdat(self,df):
        # matches emissions from the out.dat file with
        # the date and time of emission.
        # uses the tcm_columns array which has the key
        # and the sourehash dictionary which contains the information.
        datelist = []
        htlist = []
        valra = []
        for val in zip(self.tcm_columns,df[1]):
            shash = self.sourcehash[val[0]]
            datelist.append(shash['sdate'])
            htlist.append(shash['bottom'])
            valra.append(val[1])
        return list(zip(datelist,htlist,valra))

    def plot_outdat(self,vals,log=False,fignum=1):
        fig = plt.figure(fignum, figsize=(10,5))
        vals = list(zip(*vals))
        sns.set()
        if log:
            cb = plt.scatter(vals[0],vals[1],c=vals[2],s=50,cmap='Blues') 
        else:
            cb = plt.scatter(vals[0],vals[1],c=np.log10(vals[2]),s=50,cmap='Blues') 
        plt.colorbar(cb)

    def plot_out2dat(self,daterange,df2,cmap='viridis'):
        sns.set()
        tii = self.time_index(daterange[0])
        volcat = self.volcat_avg_hash[tii] 
        shape = volcat.shape
        model = df2['model'].values
        model = model.reshape(shape[0],shape[1])
        cb = plt.pcolormesh(volcat.longitude, volcat.latitude,model,cmap=cmap)
        plt.colorbar(cb)
        cb2 = plt.scatter(volcat.longitude, volcat.latitude, c=volcat.values,s=10,cmap=cmap)
        plt.colorbar(cb2)
        return volcat 

    def compare_plots(self, daterange):
        tii = self.time_index(daterange[0])
        cdump = self.cdump_hash[tii]
        volcat = self.volcat_avg_hash[tii] 
        csum = cdump.sum(dim='ens')
        plt.pcolormesh(csum.longitude, csum.latitude, np.log10(csum),cmap='Reds')
        cb = plt.pcolormesh(volcat.longitude, volcat.latitude, np.log10(volcat),cmap='Blues')
        plt.colorbar(cb)
        ax = plt.gca()
        return ax

    def prepare_one_time(self, daterange):
        vdir = self.vdir
        vid = self.vid
        tii = self.time_index(daterange[0])
        # assume cdump has times at beginning of averaging time.
        cdump_a = self.cdump.sel(time=daterange[0]) 
        
        # need to clip cdump and volcat.
        dummy = cdump_a.sum(dim='ens')
        a1,a2,b1,b2 = self.clip(dummy)
        dummy = dummy[a1:a2,b1:b2]
        #cdump_a = cdump_a[a1:a2,b1:b2]      
    
        # get list of volcat data arrays. 
        das = self.get_volcat(daterange)

        # regrid volcat to dummy grid.
        regrid_volcat = volcat.average_volcat(das,dummy)

        # average volcat values.
        avg = regrid_volcat.mean(dim='time')
        cdump_a = cdump_a[:,a1:a2,b1:b2]

        self.cdump_hash[tii] = cdump_a
        self.volcat_avg_hash[tii] = avg
        return  cdump_a, avg
        




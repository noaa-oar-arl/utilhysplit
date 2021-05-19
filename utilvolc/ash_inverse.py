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
from utilvolc.basic_checks import compare_grids
from utilvolc.basic_checks import calc_grids
from utilvolc.runhelper import Helper

def testcase(tdir,vdir):
    fname = 'xrfile.invbezy.nc'
    bezyloc = [160.587,55.978]
    vid = 'v300250'
    inverse = InverseAsh(tdir,fname,vdir,vid)
    return inverse

def vincent(tdir,vdir):
    # 0.25 degrees
    fname = 'xrfile.vincent2.nc'
    # 0.05 degrees
    #fname = 'xrfile.446.nc'
    inverse = InverseAsh(tdir,fname,vdir,vid=None)
    return inverse

def vincent13(tdir,vdir):
    fname = 'xrfile.445.nc'
    inverse = InverseAsh(tdir,fname,vdir,vid=None)
    return inverse

# get cdump output.
# pick time period.
# cdump output for that time period. 

def getL1L2_filenames(d1, nfiles=100):
    dt = datetime.timedelta(minutes=10)
    flist = []
    dlist = []
    fmt1 = "%Y%j_%H%M"
    d2 = d1
    for iii in np.arange(0,nfiles):
        fname = 'geocatL2.GOES-16.Full_Disk.{}30.some_vars.nc'.format(d2.strftime(fmt1))
        flist.append(fname)
        dlist.append(d2)
        d2 = d2 + dt
    return flist, dlist

def get_L1L2(tdir, d1, nfiles=100):
    flist, dlist = getL1L2_filenames(d1,nfiles)
    das = []
    for fname in flist:
        fullfname = os.path.join(tdir, fname) 
        if os.path.isfile(fullfname):
           das.append(volcat.open_dataset(fullfname,decode_times=True,mask_and_scale=False))
        else:
           print('NOT adding {}'.format(fullfname))
    return das 

def find_area(vmass):
    r2 = vmass.where(vmass>0)
    r2 = r2.fillna(0)
    r2 = r2.where(r2<=0)
    r2 = r2.fillna(1)
    return int(r2.sum())

def plot_L1L2(das,nnn=100):
    masslist = []
    volume = []
    dlist = []
    htlist = []
    alist = []
    for iii in np.arange(0,nnn):
        try:
            vmass  = volcat.get_mass(das[iii],clip=True)
        except:
            continue
        vht  = volcat.get_height(das[iii],clip=True)
        #print(vmass.sum())
        area = 0.03*0.03*111e3*111e3*np.cos(56*np.pi/180.0) #area of one pixel in meters.
        mass = vmass.sum() * area
        #print('{} kg'.format(mass/1e3))
        #print('{}------'.format(iii))
        masslist.append(float(mass.values))
        dlist.append(das[iii].time.values)
        htlist.append(float(np.max(vht)))
        alist.append(find_area(vmass))
    #print(masslist)
    fig = plt.figure()
    plt.plot(dlist,np.array(masslist)/1000.0,'--bo')
    ax = plt.gca()
    fig.autofmt_xdate()
    ax.set_ylabel('Total mass (kg)')
    ax.set_xlabel('Time')
    #plt.savefig('bezymass.png')
    plt.show()
    fig = plt.figure()
    plt.plot(dlist,htlist,'--bo')
    ax = plt.gca()
    fig.autofmt_xdate()
    ax.set_ylabel('Maximum height (km)')
    ax.set_xlabel('Time')
    #plt.savefig('bezyheight.png')
    plt.show()
    fig = plt.figue()
    plt.plot(dlist,alist,'--bo')
    ax = plt.gca()
    fig.autofmt_xdate()
    ax.set_ylabel('Number of pixels')
    #plt.savefig('bezypixels.png')
    ax.set_xlabel('Time')


class ParametersIn:

    def __init__(self,fname):
        with open(fname, 'r') as fid:
             self.lines = fid.readlines()

    def change_and_write(self,nx_ctrl,fname):
        self.change(nx_ctrl)
        self.write(fname)

    def change(self, nx_ctrl):
        nx_ctrl = int(nx_ctrl)
        self.lines[1] = ' N_ctrl={}\n'.format(nx_ctrl)
        self.lines[2] = ' Nx_ctrl={}\n'.format(nx_ctrl)
        self.lines[40] = ' lbfgs_nbd={}*0'.format(nx_ctrl) 
 
    def write(self,fname):
        with open(fname,'w') as fid:
            [fid.write(line) for line in self.lines]



class InverseOutDat:

    def __init__(self,wdir,fname='out.dat',fname2='out2.dat'):
        # out.dat has estimated release rates in same order as tcm columns.
        # out2.dat has observed(2nd column) and modeled(3rd column) mass loadings.
        self.wdir = wdir
        self.df = pd.DataFrame()       
        self.df2 = pd.DataFrame()
        self.fname = fname
        self.fname2 = fname2

    def read_out(self,name):
        wdir = self.wdir
        df = pd.read_csv(os.path.join(wdir,name),sep='\s+',header=None)
        return df

    def get_emis(self,name=None):
        if not name: 
           name = self.fname
        else:
           self.fname = name
        df = self.read_out(name)
        self.df = df
        return df 

    def emis_hist(self, units='g/h'):
        if units=='g/h':
             vals = np.log10(self.df[1].values)
        nmin = int(np.floor(np.min(vals)))
        nmax = int(np.ceil(np.max(vals)))
        nbins = len(np.arange(nmin,nmax,1))
        plt.hist(vals,bins=nbins)
        ax = plt.gca()
        ax.set_xlabel('Log emission (g/h)')

    def get_conc(self, name=None):
        if not name: name = self.fname2
        df = self.read_out(name)
        df.columns = ['index', 'observed', 'model']
        self.df2 = df
        return df

    def plot_conc(self,cmap='viridis'):
        sns.set()
        sns.set_style('whitegrid')
        df = self.df2
        #plt.plot(df['observed'],df['model'],'k.',MarkerSize=3)
        cb = plt.hist2d(df['observed'],df['model'],cmap=cmap,norm=mpl.colors.LogNorm(),bins=[20,20])
        cbar = plt.colorbar(cb[3])
        cbar.ax.set_ylabel('Number of Points')
        ax = plt.gca()
        ax.set_xlabel('observed')
        ax.set_ylabel('modeled')
        nval = np.max(df['observed'])
        # plot 1:1 line
        plt.plot([0,nval],[0,nval],'--b',LineWidth=1)
        return ax

    def get_vals(self):
        return -1

def get_inp_hash(wdir,configfile):
    from utilvolc.ashapp import ashinverse
    #from utilvolc.ashapp.runhelper import JobSetUp
    from utilvolc.ashapp.runhelper import make_inputs_from_file
    setup = make_inputs_from_file(wdir,configfile)
    setup.add_inverse_params()
    return setup.inp

def get_sourcehash(wdir,configfile):
    from utilvolc.ashapp import ashinverse
    #from utilvolc.ashapp.runhelper import JobSetUp
    from utilvolc.ashapp.runhelper import make_inputs_from_file
    setup = make_inputs_from_file(wdir,configfile)
    setup.add_inverse_params()
    sourcehash = ashinverse.inverse_get_suffix_list(setup.inp)
    return sourcehash 

# For the inverse modeling work flow goes
# inverse = ai.InverseAsh(...)
# loop through inverse.perpare_one_time(....)
# use compare_plotsA to check coverage
# inverse.make_tcm_mult() to create TCM
# plot_tcm to check.
# inverse.write_tcm
# Run the TCM 

# call to InverseOutDat class to read output from inversion algorithm.
# input dataframe from InverseOutDat into
# inverse.plot_outdat() 

class InverseAshEns:
    """
    Inverse runs from a meteorological ensemble.
    """
    def __init__(self, tdirlist, fnamelist,
                 vdir,vid,
                 configdir='./',configfile=None,
                 verbose=False):
        self.invlist = []  # list of  InverseAsh objects
        self.fnamelist = fnamelist
        self.tcm_names = [] # list of names of tcm files written.
        
        # assume strings of form NAME_gep04.nc
        try:
            self.taglist = [x[-8:-3] for x in fnamelist]
        except:
            print('InverseAshEns string not in expected form {}'.format(fnamelist[0]))
            self.taglist = list(map(str,np.arange(0,len(fnamelist))))
        for hruns in zip(tdirlist, fnamelist):
            self.invlist.append(InverseAsh(hruns[0],hruns[1],vdir,vid,configdir,configfile,verbose=verbose))

    def set_directory(self,wdir,execdir):
        self.wdir = wdir
        self.execdir = execdir

    def set_concmult(self,mult):
        for hrun in self.invlist:
            hrun.set_concmult(mult)   

    def write_tcm(self,tcm_name):
        for hrun in zip(self.invlist,self.taglist): 
            tname = tcm_name.replace('.txt','')
            tname = '{}_{}.txt'.format(tname,hrun[1])   
            hrun[0].write_tcm(tname)
            if tname not in self.tcm_names:
                self.tcm_names.append(tname)

    def plot_tcm(self,ensi=None):
        if ensi:
            self.invlist[ensi].plot_tcm()
            return True
        for hrun in self.invlist: 
            hrun.plot_tcm()
            plt.show()
        return True

    def make_tcm_names(self):
        out_name1 = 'out.dat'
        out_name2 = 'out2.dat'
        name1 = []
        name2 = []
        for tag in self.taglist:
            name1.append('{}_{}'.format(tag,out_name1))
            name2.append('{}_{}'.format(tag,out_name2))
        return name1, name2 

    def run_tcm(self):
        out_name1 = 'out.dat'
        out_name2 = 'out2.dat'
        inp_name = 'TCM_sum.csv'
        cmd = os.path.join(self.execdir,'new_lbfgsb.x')
        new_name1, new_name2 = self.make_tcm_names()
        for iii, tcm in enumerate(self.tcm_names):
            os.chdir(self.wdir)
            Helper.remove(inp_name)
            Helper.copy(tcm,inp_name)
            print(cmd)
            print(os.getcwd())
            #Helper.execute_with_shell(cmd)
            Helper.execute(cmd)
            Helper.move(out_name1, new_name1[iii]) 
            Helper.move(out_name2, new_name2[iii]) 
                   
    def plot_outdat(self):
        ilist = self.read_outdat()
        for iii, io in enumerate(ilist):
            print(self.taglist[iii])
            df2 = io.get_conc()
            #df = io.get_emis()
            try: 
                ax = io.plot_conc()
            except:
                pass
            plt.show()
            #io.emis_hist()
            #ax = plt.gca()
            #ax.set_xlim(-1,10)
            #plt.show()  
 
    def read_outdat(self):
        name1,name2 = self.make_tcm_names()
        ilist = []
        for outdat in zip(name1,name2):
            io = InverseOutDat(self.wdir,outdat[0],outdat[1])
            ilist.append(io)
        return ilist

    def plot_outdat_ts(self):
        sns.set_style('whitegrid')
        fig = plt.figure(1)
        ax = fig.add_subplot(1,1,1)
        ilist = self.read_outdat()
        clr = ['--k','--r','--b','--g','--c','--y','--m']
        jjj=0
        for iii, io in enumerate(ilist):
            df = io.get_emis()
            self.invlist[0].plot_outdat_ts(df,ax=ax,clr=clr[jjj])
            jjj+=1
            if jjj >= len(clr): jjj=0
        fig.autofmt_xdate()
        plt.show()

    def make_tcm_mult(self,tiilist,remove_cols=True,remove_rows=True):
        for hrun in self.invlist: 
            a,b,c, = hrun.make_tcm_mult(tiilist,remove_cols,remove_rows)

    def prepare_one_time(self,daterange):
        for iii, hrun in enumerate(self.invlist):
            #print(self.taglist[iii])
            hrun.prepare_one_time(daterange)    

    def compare_plotsA(self,daterange=None,tii=None,ens=None):
        for iii, hrun in enumerate(self.invlist):
            print(self.taglist[iii])
            
            hrun.compare_plotsA(daterange, tii)
            plt.show()


class InverseAsh:

    def __init__(self, tdir, fname,vdir,vid,
                 configdir='./',configfile=None,
                 verbose=False):
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
        # multiplication factor if more than 1 unit mass released.
        self.concmult = 1
        #
        self.get_cdump(tdir,fname,verbose)
        self.add_config_info(configdir,configfile)

    def get_cdump(self,tdir,fname,verbose=False):
        # hysplit output. xarray. 
        cdump = xr.open_dataset(os.path.join(tdir,fname))
        if verbose: print('opened',tdir,fname)
        # turn dataset into dataarray
        temp = list(cdump.keys())
        cdump = cdump[temp[0]]
        # get rid of source dimension (for now)
        cdump = cdump.isel(source=0)
        self.cdump = cdump

    def add_config_info(self,configdir, configfile):
        # the ens dimension holds is key to what emission source was used.
        # the sourcehash is a dictionary
        # key is the ensemble number
        # values is another dictionary with
        # sdate: begin emission
        # edate: end emission
        # bottom : lower height of emission
        # top : upper height of emission.
        if configfile:
            self.sourcehash = get_sourcehash(configdir, configfile)
            self.inp = get_inp_hash(configdir, configfile)
        else:
            self.sourcehash = {}
            self.inp = {}

    def add_inp(self,configdir, configfile):
        self.inp = get_inp_hash(configdir, configfile)

    def set_concmult(self,concmult):
        self.concmult = concmult

    def get_volcat(self, daterange, verbose=False):
        vdir = self.vdir
        vid = self.vid
        tii = self.time_index(daterange[0])
        #if tii not in self.volcat_hash.keys(): 
        done=True
        if done: 
            das = volcat.get_volcat_list(vdir,daterange=daterange,vid=vid,decode_times=True,
                                         verbose=verbose) 
            self.volcat_hash[tii] = das
        else:
            das = self.volcat_hash[tii]
        # create one dataset with dimension of time.
        if len(das) > 0:
            vset = xr.concat(das,dim='time')
        else:
            print('No volcat files found ')
            return das
        #vra = vset.ash_mass_loading
        #vra = vra.fillna(0)
        #vmean = vra.mean(dim='time')

        return vset

    def clip(self,dummy, buf=0):
        # clip ra according to where dummy has 0's.
        aaa = np.where(dummy > 0)
        a1 = np.min(aaa[0])
        a2 = np.max(aaa[0])
        b1 = np.min(aaa[1])
        b2 = np.max(aaa[1])

        if a2+buf < dummy.y.values[-1]:
           a2 = a2+buf

        if b2+buf < dummy.x.values[-1]:
           b2 = b2+buf

        a1 = a1 - buf
        b1 = b1 - buf
        if a1<0: a1=0
        if b1<0: b1=0

        return a1,a2,b1,b2

    def print_times(self):
        timelist = [pd.to_datetime(x) for x in self.cdump.time.values]
        for time in timelist:
            print(time.strftime("%Y %m %d %H:%Mz"))

    def get_time(self,tii):
        timelist = [pd.to_datetime(x) for x in self.cdump.time.values]
        return timelist[tii] 

    def time_index(self,time):
        timelist = [pd.to_datetime(x) for x in self.cdump.time.values]
        try:
           iii = timelist.index(time)
        except:
           iii = None
        return iii 

    def make_tcm_mult(self,tiilist,remove_cols=True,remove_rows=True):
        # make the tcm for multiple time periods.
        tcmlist = []
        latlist = []
        lonlist = []
        for tii in tiilist:
            print(self.cdump.time.values[tii])
            tcm, model_lat, model_lon, columns = \
                 self.make_tcm(tii,remove_cols=False,remove_rows=remove_rows)
            tcmlist.append(tcm)
            latlist.append(np.array(model_lat))
            lonlist.append(np.array(model_lon))
        t3 = np.concatenate(tcmlist,axis=0)
        lat = np.concatenate(latlist,axis=0)
        lon = np.concatenate(lonlist,axis=0)
        self.latlist = np.array(latlist)
        self.lonlist = np.array(lonlist)
        if remove_cols:
           nmax = t3.shape[1]
           iremove = []
           # very last column is obs. So start with second to last column.
           for nnn in np.arange(nmax-2,0,-1):
               test = t3[:,nnn]
               # remove colum if it is all 0's.
               if np.all(test==0.0): iremove.append(nnn)
               else: break
           t3 = np.delete(t3,iremove,axis=1) 
        self.tcm = t3
        self.tcm_lat = lat
        self.tcm_lon = lon
        self.latlist = np.array(latlist)
        self.lonlist = np.array(lonlist)
        return t3, lat, lon

    def make_tcm(self,tii, remove_cols=True, remove_rows=False):
        # remove rows means remove rows with observations of 0 or nan.

        # header line indicates release times for each column.
        # I think header can just be a dummy and it is added when writing to file.
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
        cdump = self.cdump_hash[tii]
        avg = self.volcat_avg_hash[tii] 
        cdump = cdump * self.concmult 

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

        volc = volc.flatten()
        if remove_rows:
           vpi = np.where(volc>0)
           model = model[vpi]
           volc = volc[vpi]
           model_lat = model_lat[vpi]
           model_lon = model_lon[vpi]
           volc_lon = volc_lon[vpi]
           volc_lat = volc_lat[vpi]
        volc = volc.reshape(volc.shape[0],1)

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
        plt.pcolormesh(np.log10(self.tcm),cmap='tab20')

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
        print('output file {}'.format(name))
            
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

    def make_outdat(self,dfdat):
        """
        dfdat : pandas dataframe output by InverseOutDat class get_emis method.
        Returns
        vals : tuple (date, height, emission mass)
        """
        # matches emissions from the out.dat file with
        # the date and time of emission.
        # uses the tcm_columns array which has the key
        # and the sourehash dictionary which contains the information.
        datelist = []
        htlist = []
        valra = []
        for val in zip(self.tcm_columns,dfdat[1]):
            shash = self.sourcehash[val[0]]
            datelist.append(shash['sdate'])
            htlist.append(shash['bottom'])
            valra.append(val[1])
        vals =  list(zip(datelist,htlist,valra))
        return vals

    def make_outdat_df(self,dfdat):
        #dfdat : pandas dataframe output by InverseOutDat class get_emis method.
        vals = self.make_outdat(dfdat)
        vals = list(zip(*vals))
        ht = vals[1]
        time = vals[0]
        #emit = np.array(vals[2])/1.0e3/3600.0 
        emit = np.array(vals[2])
        dfout = pd.DataFrame(zip(time,ht,emit))
        dfout = dfout.pivot(columns=0,index=1)
        return dfout

    def plot_outdat(self,vals,
                    log=False,
                    fignum=1,
                    cmap='Blues',
                    unit='kg/s',
                    thresh=0):
        """
        vals is output by make_outdat.
        """
        fig = plt.figure(fignum, figsize=(10,5))
        vals = list(zip(*vals))
        sns.set()
        sns.set_style('whitegrid')
        # output in kg/s?/
        if unit == 'kg/s':
            emit = np.array(vals[2])/1.0e3/3600.0 
        elif unit == 'kg/h':
            emit = np.array(vals[2])/1.0e3 
        elif unit == 'g/h':
            emit = np.array(vals[2])/1.0
        vpi = np.where(emit < thresh)
        emit[vpi]=0
        ht = np.array(vals[1])/1e3
        if not log:
            cb = plt.scatter(vals[0],ht,c=emit,s=100,cmap=cmap,marker='s') 
        else:
            
            cb = plt.scatter(vals[0],ht,c=np.log10(emit),s=100,cmap=cmap,marker='s') 
            #cb = plt.pcolormesh(vals[0],ht,emit,cmap=cmap) 
        cbar = plt.colorbar(cb)
        cbar.ax.set_ylabel(unit)
        fig.autofmt_xdate()

    def plot_outdat_ts(self,dfdat,log=False,fignum=1,unit='kg/s',
                       ax=None,clr='--ko'):
        # plots time series of MER. summed along column.
        #dfdat : pandas dataframe output by InverseOutDat class get_emis method.
        if not ax:
            sns.set()
            sns.set_style('whitegrid')
            fig = plt.figure(fignum, figsize=(10,5))
            ax = fig.add_subplot(1,1,1)
        df = self.make_outdat_df(dfdat)     
        sns.set()
        ts = df.sum()
        if unit == 'kg/s':
           yval = ts.values/3.6e6
        elif unit == 'g/h':
           yval = ts.values
        plt.plot([x[1] for x in ts.index.values], yval, clr)
        #fig.autofmt_xdate()
        ax.set_ylabel('MER {}'.format(unit))
        #return ax, df

    def plot_out2dat_times(self,df2,cmap='viridis'):
        return -1 

    def plot_out2dat_scatter(self,tiilist,df2,vloc,cmap='Blues'):
        if isinstance(tiilist,int): tiilist = [tiilist]
        sns.set()
        ppp=0
        #tii = self.time_index(daterange[0])
        modelall = df2['model'].values
        nnn=0
        for tii in tiilist:
            print(self.cdump.time.values[tii])
            fig = plt.figure(figsize=[10,5])
            ax1 = fig.add_subplot(1,2,1)
            ax2 = fig.add_subplot(1,2,2)
            volcat = self.volcat_avg_hash[tii] 
            shape = volcat.shape
            lon = self.lonlist[tii-1]
            lat = self.latlist[tii-1]
            model = modelall[nnn:nnn+len(lon)]
            volcat = self.volcat_avg_hash[tii] 
            r2 = volcat.where(volcat>0)
            cb = ax1.scatter(lon, lat,c=model,cmap=cmap,s=10,marker='o')
            cb2 = ax2.scatter(volcat.longitude, volcat.latitude, c=r2.values,s=10,cmap=cmap,marker='o')
            nnn=len(lon) 
            if vloc:
               ax1.plot(vloc[0],vloc[1],'y^')
               ax2.plot(vloc[0],vloc[1],'y^')
            plt.colorbar(cb, ax= ax1)
            plt.colorbar(cb2, ax=ax2)
            

    def plot_out2dat(self,tiilist,df2,cmap='viridis',vloc=None,ptype='pcolormesh'):
        # tiilist needs to be same order as for tcm.
        # df2 is from InverseOutDat get_conc method.
        # doesn't work when only later time periods are used!
        if isinstance(tiilist,int): tiilist = [tiilist]
        sns.set()
        ppp=0
        #tii = self.time_index(daterange[0])
        for tii in tiilist:
            print(self.cdump.time.values[tii])
            fig = plt.figure(figsize=[10,5])
            ax1 = fig.add_subplot(1,2,1)
            ax2 = fig.add_subplot(1,2,2)
            volcat = self.volcat_avg_hash[tii] 
            shape = volcat.shape
            model = df2['model'].values
            temp = model[ppp:ppp+shape[0]*shape[1]]
            model = temp.reshape(shape[0],shape[1])
            vpi = np.where(model < 0.01)
            model[vpi] = np.nan
            ppp = ppp+ shape[0]*shape[1]
            r2 = volcat.where(volcat>0)
            m_max = np.nanmax(model)
            v_max = np.nanmax(r2.values)
            m_min = np.nanmin(model)
            v_min = np.nanmin(r2.values)
            p_min = np.nanmin([m_min,v_min])
            p_max = np.nanmax([m_max,v_max])
            norm = mpl.colors.Normalize(vmin=p_min, vmax=p_max) 
            print(np.nanmax(model), np.nanmax(r2.values))
            if ptype == 'pcolormesh':
                cb = ax1.pcolormesh(volcat.longitude, volcat.latitude,model,norm=norm, cmap=cmap,shading='nearest')
                cb2 = ax2.pcolormesh(volcat.longitude, volcat.latitude,r2.values,norm=norm, cmap=cmap,shading='nearest')
            #cb = ax1.scatter(volcat.longitude, volcat.latitude,c=np.log10(model),cmap=cmap,s=50,marker='s')
            else:
                #lon = self.lonlist[tii-1][0:14]
                #lat = self.latlist[tii-1][0:21]
                #xv, yv = np.meshgrid(lon,lat)     
                lon = self.lonlist[tii-1][0:14]
                lat = self.latlist[tii-1][0:21]
                xv, yv = np.meshgrid(lon,lat) 
                print('vshape', volcat.longitude.shape)    
                print(xv.shape, yv.shape)        
                print('model shape', model.shape)
                cb = ax1.scatter(xv, yv,c=model,cmap=cmap,s=10,marker='o')
                cb2 = ax2.scatter(volcat.longitude, volcat.latitude, c=r2.values,s=10,cmap=cmap,marker='o')
            plt.colorbar(cb, ax= ax1)
            #cb2 = ax2.scatter(volcat.longitude, volcat.latitude, c=r2.values,s=10,cmap=cmap,marker='o')
            plt.colorbar(cb2, ax=ax2)
            if vloc:
               ax1.plot(vloc[0],vloc[1],'y^')
               ax2.plot(vloc[0],vloc[1],'y^')
            plt.show()
        return volcat 

    def compare_plotsA(self, daterange,tii=None,levels=None):
        fig = plt.figure(1,figsize=(10,5))
        ax1 = fig.add_subplot(1,1,1)
        if not tii:
            tii = self.time_index(daterange[0])
        print('tii',tii)
        cdump = self.cdump_hash[tii]
        volcat = self.volcat_avg_hash[tii] 
        csum = cdump.sum(dim='ens')
        print(cdump.time)
        #volcat.plot.pcolormesh(x='longitude',y='latitude',levels=levels,ax=ax1)
        #cdump.sum(dim='ens').plot.contour(x='longitude',y='latitude',ax=ax2)
        try:
            #plt.pcolormesh(csum.longitude, csum.latitude, np.log10(csum),cmap='Reds',shading='nearest')
            plt.pcolormesh(csum.x, csum.y, np.log10(csum),cmap='Reds',shading='nearest')
        except:
            print('FAILED max value', np.max(csum))
            print('------------------')
        #plt.pcolormesh(csum.longitude, csum.latitude, csum,cmap='Reds',shading='nearest')
        #cb= csum.plot.pcolormesh(x='longitude',y='latitude',cmap='viridis',ax=ax2)
        cb = plt.pcolormesh(volcat.x, volcat.y, np.log10(volcat),cmap='Blues',shading='nearest')
        #cb = plt.scatter(volcat.longitude, volcat.latitude, c=np.log10(volcat),s=2,cmap='Blues')
        #cb = plt.scatter(volcat.longitude, volcat.latitude, c=volcat.values,s=2,cmap='viridis',levels=levels)
        #cb = plt.contour(volcat.longitude, volcat.latitude, np.log10(volcat),cmap='Blues')
        # plt.colorbar(cb)
        plt.tight_layout()


    def get_norm(self,model,r2):
        m_max = np.nanmax(model)
        v_max = np.nanmax(r2.values)
        m_min = np.nanmin(model)
        v_min = np.nanmin(r2.values)
        p_min = np.nanmin([m_min,v_min])
        p_max = np.nanmax([m_max,v_max])
        norm = mpl.colors.Normalize(vmin=p_min, vmax=p_max) 
        return norm

    
    def match_volcat(self,forecast):
        time = pd.to_datetime(forecast.time.values)
        tii = self.time_index(time)
        volcat = self.volcat_avg_hash[tii] 
        return volcat

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

    def remove_nans(self,data,thresh=None):
        # remove nans
        data2 = data[~np.isnan(data)]
        if thresh:   
           data2 = data2[data2>thresh]
           # apply threshold.
           #vpi = data2 < thresh
           #data2[vpi] = np.nan
           #data2 = data2[~np.isnan(data2)]
        return data2
           
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

    def prepare_times(self,daterangelist,das=None):
        for dates in daterangelist:
            subdas = []
            for sd in das:
                if pd.to_datetime(sd.time.values) < dates[1] and \
                   pd.to_datetime(sd.time.values) >= dates[0]:
                   subdas.append(sd)
            print('preparing', dates)
            self.prepare_one_time(dates,das=subdas)

    def add_volcat_hash(self, volcat_hash):
        # allow to update from another instance of the class.
        self.volcat_avg_hash.update(volcat_hash)

    def add_cdump_hash(self, cdump_hash):
        self.cdump_hash.update(cdump_hash)

    def align_grids(self,grid1,grid2):
        # takes care of making sure grids are aligned.
        checks, dlon,dlat,minlat,minlon, maxlat,maxlon = calc_grids(grid1, grid2, verbose=False)
        latra = np.clip(np.arange(minlat,maxlat+dlat,dlat),None,maxlat)
        lonra = np.clip(np.arange(minlon,maxlon+dlon,dlon),None,maxlon)
        if latra[-1] - latra[-2] < 1e-3: latra = latra[0:-1]
        if lonra[-1] - lonra[-2] < 1e-3 : lonra = lonra[0:-1]
        grid1, grid2 = xr.align(grid1, grid2, join='outer')

        # align function will put nans in lat,lon coordinates. need to redo these coords.
        mgrid = np.meshgrid(lonra, latra)
        grid1 = grid1.assign_coords(longitude=(("y","x"),mgrid[0]))
        grid1 = grid1.assign_coords(latitude=(("y","x"),mgrid[1]))
        grid2 = grid2.assign_coords(longitude=(("y","x"),mgrid[0]))
        grid2 = grid2.assign_coords(latitude=(("y","x"),mgrid[1]))
       
        return grid1, grid2

    def prepare_one_time(self, daterange, das=None):
        vdir = self.vdir
        vid = self.vid
        tii = self.time_index(daterange[0])
        # assume cdump has times at beginning of averaging time.
        cdump_a = hysplit.hysp_massload(self.cdump.sel(time=daterange[0])) 
        if not das:
           vset = self.get_volcat(daterange)

        buf = 5
        # clip the volcat array before aligning.
        vra = vset.ash_mass_loading
        a1,a2,b1,b2 = self.clip(vra.sum(dim='time'),buf=buf)
        vra = vra[:,a1:a2,b1:b2]

        # clip the cdump array before aligning.
        if 'ens' in cdump_a.coords:
            dummy = cdump_a.sum(dim='ens')
        else:
            dummy = cdump_a
        a1,a2,b1,b2 = self.clip(dummy,buf=5)
        if 'ens' in cdump_a.coords:
            cdump_a = cdump_a[:,a1:a2,b1:b2]
        else:
            cdump_a = cdump_a[a1:a2,b1:b2]

        if compare_grids(cdump_a, vra):
            cdump_a, vra = self.align_grids(cdump_a, vra)
        #if compare_grids(cdump_a, vra):
        #    cdump_a, vra = xr.align(cdump_a, vra, join='outer')
        else:
            print('prepare_one_time: grids cannot be aligned')
            return False

        # TO DO: may not need to clip again.
        # need to clip cdump and volcat.
        if 'ens' in cdump_a.coords:
            dummy = cdump_a.sum(dim='ens')
        else:
            dummy = cdump_a
        a1,a2,b1,b2 = self.clip(dummy)
        dummy = dummy[a1:a2,b1:b2]

        aa1,aa2,bb1,bb2 = self.clip(vra.sum(dim='time').fillna(0))
        # use bounds which encompass all obs and model data. 
        a1 = np.min([a1,aa1])
        a2 = np.max([a2,aa2])
        b1 = np.min([b1,bb1])
        b2 = np.max([b2,bb2])
        buf = 5
        if a2+buf < cdump_a.y.values[-1] and a2+buf < vra.y.values[-1]:
           a2 = a2+buf
        if b2+buf < cdump_a.x.values[-1] and b2+buf < vra.x.values[-1]:
           b2 = b2+5

        a1 = a1 -5
        b1 = b1 -5
        if a1<0: a1=0
        if b1<0: b1=0

        avra = vra.fillna(0).mean(dim='time')

        if 'ens' in cdump_a.coords:
            self.cdump_hash[tii] = cdump_a[:,a1:a2,b1:b2]
        else:
            self.cdump_hash[tii] = cdump_a[a1:a2,b1:b2]
        self.volcat_avg_hash[tii] = avra[a1:a2,b1:b2]
        #vra = vra[a1:a2,b1:b2]
        #vra = vra.fillna(0)
        #vmean = vra.mean(dim='time')
        return cdump_a, avra

    def compare_time_ave(self, daterange):
        """
        creates plots illustrating how time averaging affects 
        volcat data.
        """
        sns.set()
        vset = self.get_volcat(daterange)
        vset = vset.ash_mass_loading
        for time in vset.time.values:
            temp = vset.sel(time=time)
            a1,a2,b1,b2 = self.clip(temp.fillna(0))
            temp = temp[a1:a2,b1:b2]
            temp.plot.pcolormesh(x='longitude', y='latitude')
            plt.show()
        avra = vset.fillna(0).mean(dim='time')
        avra2 = vset.mean(dim='time')
        a1,a2,b1,b2 = self.clip(avra)
        avra = avra[a1:a2,b1:b2] 
        avra2 = avra2[a1:a2,b1:b2] 
        print('Average with nans set to 0')
        avra.plot.pcolormesh(x='longitude', y='latitude')
        plt.show()
        print('Average of values above 0')
        avra2.plot.pcolormesh(x='longitude', y='latitude')
        plt.show()
        diff = avra2.fillna(0) - avra
        diff.plot.pcolormesh(x='longitude', y='latitude')
 
        
def make_control(efile,
                 tdir='./',
                 sname = 'CONTROL.default',
                 wdir='./'
                 ):
    control = hcontrol.HycsControl(fname=sname,working_directory=tdir)
    control.read()
    # number of locations need to be written.
    nlocs = efile.cycle_list[0].nrecs 
    nspecies = len(efile.sphash.keys())
    stime = efile.sdatelist[0]
    done=False
    while not done:
        if nlocs == control.nlocs: 
           done=True
        else:
           control.add_dummy_location()
    # set to 0 to use emit times file.
    for species in control.species:
        species.rate = 0
        species.duration = 0
    for grid in control.concgrids:
        grid.outfile = 'cdump.emit'
        grid.outdir = wdir
    # simulation start date same as emit-times 
    control.rename('CONTROL.emit',wdir)
    control.add_sdate(stime)
    control.write()
    return control
    

def make_setup( tdir = './',
                sname='SETUP.default',
                wdir = './',
                efilename = 'emit.txt',
               ):
    setup = hcontrol.NameList(sname,tdir)
    setup.read()
    setup.add('efile',efilename)
    setup.add('poutf','PARDUMP.emit')
    setup.rename('SETUP.emit',wdir)
    setup.write()
    return setup

def make_efile(vals,vlat,vlon,
               area=1,
               emis_threshold=50,
               vres=1000, 
               name='emit.txt',
               date_cutoff=None):
    """
    vals : output from make_outdat method in Inverse class.
    vlat : latitude of vent
    vlon : longitude of vent
    area : area to place ash over.
    emis_threshold : do not use emissions below this value.
    vres : vertical resolution. Needed to create last point in vertical line source.
    """
    from utilhysplit import emitimes
    efile = emitimes.EmiTimes()
    duration = "0100"
    ptime = vals[0][0]
    numcycles = 0
    for iii, value in enumerate(vals):
        time = value[0]
        if date_cutoff:
           if time >= date_cutoff: break
        height = value[1]
        emis = value[2]
        if emis < emis_threshold: emis=0
        newcycle = False
        if iii == 0 : newcycle = True
        elif time != ptime: newcycle =True
        else: newcycle=False
        if newcycle:
           numcycles += 1
           efile.add_cycle(time,duration=1)
        #if emis > emis_threshold:
        efile.add_record(time,
                         duration=duration,
                         lat=vlat,
                         lon=vlon,
                         height=height,
                         rate = emis,
                         area = area,  
                         heat = 0,
                         spnum=1,
                         nanvalue=0)
        ptime = time       
        # if there will be a new cycle
        # need to add another record for top height of last point.
        if iii+1 < len(vals): next_time = vals[iii+1][0] 
        else: next_time = vals[0][0]
        if next_time != time and vres > 0:
            efile.add_record(time,
                         duration=duration,
                         lat=vlat,
                         lon=vlon,
                         height=height+vres,
                         rate = 0,
                         area = area,  
                         heat = 0,
                         spnum=1,
                         nanvalue=0)
    efile.write_new(name) 
    return efile







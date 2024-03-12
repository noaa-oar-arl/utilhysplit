import logging
import json
import pandas as pd
import numpy as np
import os
import datetime
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import time
import monet
import sys
import xarray as xr

from utilvolc import volcano_names
from utilvolc.volcano_names import fix_volc_name
from utilvolc.runhelper import Helper
from utilvolc.runhelper import list_dirs
from utilvolc.runhelper import make_dir

from utilvolc.volcat_event import EventDF
from utilvolc.volcat_event import EventStatus

from utilhysplit.runhandler import ProcessList
from utilhysplit.plotutils import map_util
import utilhysplit.evaluation.web_ensemble_plots as wep
from utilhysplit.evaluation import polygon_plots
from utilvolc import make_data_insertion as mdi
from utilhysplit.evaluation import ensemble_tools
from utilhysplit.evaluation.ensemble_tools import ATL, preprocess, topheight 
from utilhysplit.evaluation import ensemble_polygons


#from utilvolc.qvainterface import DFInterface
#from utilvolc.qvainterface import EventInterface

logger = logging.getLogger(__name__)

class ModelEventDF(EventDF):
    required = ["model_name",
                "file_name",
                "start_date",
                "end_date"]

class ModelEventDisplay:

    def __init__(self,eventid):
        pass


class ModelForecast:

    def __init__(self,
                 model_dset,
                 eventid,
                 start_time,               
                 end_time,
                 timedelta=3,
                 attrs={'model':'HYSPLIT'}):

        self.eventid = eventid
        self._model_dset = model_dset  
        self._start_time= start_time   # start time of first forecast
        self._end_time = end_time      # end time of last forecast
        self.timedelta = timedelta    # time between forecasts h
        self.attrs = attrs    
        self.forecast = None

#   properties --------------------------------------------------------------------
    @property
    def start_time(self):
        return self._start_time

    @property
    def end_time(self):
        return self._end_time

    @property
    def model_dset(self):
        return self._model_dset

    @property
    def timedelta(self):
        return datetime.timedelta(hours=self._timedelta)

    @timedelta.setter
    def timedelta(self,hours):
        if isinstance(hours,(int,float)):
           self._timedelta=hours 

    @property
    def timelist(self):
        dt = self.timedelta
        return [self.start_time+n*dt for n in [0,1,2,3,4,5]]

#   --------------------------------------------------------------------
    def make_savename(self):
        fstr = "%Y%m%d%H"
        dstr = self.start_time.strftime(self.start_time,fstr)
        sname = '{}_{}'.format(self.eventid,dstr)
        return sname
 
    def get_sourcelist(self,time_previous=3):
        """
        This is specifically for data insertion files.
        where the sources have name format that can be handled by EmitName class.
        """
        dt = datetime.timedelta(hours=time_previous)
        snames = self.model_dset.source.values
        # use data insertions up to dt h earlier.
        drange = [self.start_time-dt, self.start_time]
        sourcelist,dlist = pick_source(snames, drange)
        # if doesn't return any, then use all.
        if not sourcelist:
           print('NO SOURCES', snames, drange)
           sourcelist = []
        return sourcelist, dlist

    def get_enslist(self):
        enslist = self.model_dset.ens.values
        return enslist

    def make_gridded_qva(self,time_previous=3):
        sourcelist,dlist = self.get_sourcelist(time_previous=time_previous)
        print(sourcelist)
        enslist = self.get_enslist()
        enslist = None #need to fix this.
        #conc = self.model_dset.sel(source=sourcelist)
        conc = self.model_dset.copy()
        qva = ensemble_tools.ATLra(
            conc,
            enslist=enslist,
            sourcelist=sourcelist,
            threshlist=[0.2, 2, 5, 10],
            weights=None,
            **self.attrs
        )
        self.forecast = qva
        return qva    

    def save_gridded_qva(self):
        """
        dset should have units of mg/m3.
        """
        savename = '{}.nc'.format(self.make_savename())
        encoding = {}
        encoding["zlib"] = True
        encoding["complevel"] = 9
        ehash = {}
        ehash["Concentration"] = encoding
        ehash["FrequencyOfExceedance"] = encoding
        self.forecast.to_netcdf(qva_filename, encoding=ehash)


    def get_attrs(self):
        atthash = {}
        atts = self.model_dset.attrs
        for key in ['latitude','longitude','Level top heights (m)','sample time hours','VolcanoName']:
            if key in atts.keys():
               atthash[key] = atts[key]
        return atthash
   
    def plot_mean_mass2(self,vloc):
        from utilhysplit.evaluation import vaa_montage
        from utilhysplit.evaluation.ensemble_tools import preprocess
        dt = datetime.timedelta(hours=3)
        timelist = [self.start_time+n*dt for n in [0,1,2,3,4,5,6]]
        try:
            dset = self.model_dset.sel(time=timelist)
        except:
            print(timelist)
            return
        dset,dim = preprocess(dset) 
        print('HERE',dim)
        #dset = dset.mean(dim=dim)
        dset.attrs.update(self.get_attrs())
        if isinstance(vloc,(tuple,list,np.ndarray)):
           dset.attrs.update({'longitude':vloc[1],'latitude':vloc[0]})
        vm = vaa_montage.VAAMontage(dset,columns=3)
        return vm 
 
    def plot_mean_mass(self,vloc=None,extrasource=None,time_previous=3):
        from utilhysplit.evaluation import vaa_montage
        from utilhysplit.evaluation.ensemble_tools import preprocess
        sourcelist,dlist = self.get_sourcelist(time_previous=time_previous)
        if isinstance(extrasource,list):
           sourcelist.extend(extrasource)
        dt = datetime.timedelta(hours=3)
        templist = [self.start_time+n*dt for n in [0,1,2,3,4,5,6]]
        timelist = []
        for ttt in templist:
            if ttt in self.model_dset.time.values:
               timelist.append(ttt)
            else:
               print('warning time not in dset', ttt)  
        dset = self.model_dset.sel(time=timelist)
        dset,dim = preprocess(dset,sourcelist=sourcelist) 
        print(time_previous, sourcelist, dset.source.values)
        dset = dset.max(dim=dim)
        dset.attrs.update(self.get_attrs())
        if isinstance(vloc,(tuple,list,np.ndarray)):
           dset.attrs.update({'longitude':vloc[1],'latitude':vloc[0]})
        d0 = self.start_time - dt
        dlist.sort()
        dstr = '{} to {}'.format(dlist[0].strftime("%Y %m/%d %H:%M UTC"), dlist[-1].strftime("%m/%d %H:%M UTC"))
        mstr = '\n Number of runs {}. Maximum value from all runs.'.format(len(dlist))
        dset.attrs.update({'source':'Data Insertion {} {}'.format(dstr, mstr)})
        vm = vaa_montage.VAAMontage(dset,columns=3)
        try:
            fig = vm.plotpage()
        except:
            pass
        return vm

    def plot_concentration(self,vloc=None,pagenum=0,time_previous=3):
        from utilhysplit.evaluation import concentration_montage
        from utilhysplit.evaluation.ensemble_tools import preprocess
        sourcelist,dlist = self.get_sourcelist(time_previous=time_previous)
        dt = datetime.timedelta(hours=3)
        timelist = [self.start_time+n*dt for n in [0,1,2,3,4,5,6]]
        dset = self.model_dset.sel(time=timelist)
        print(dset.source.values)
        print(dset)
        dset,dim = preprocess(dset,sourcelist=sourcelist) 
        print(dset.source.values)
        dset = dset.max(dim=dim)
        dset.attrs.update(self.get_attrs())
        if isinstance(vloc,(tuple,list,np.ndarray)):
           dset.attrs.update({'longitude':vloc[1],'latitude':vloc[0]})
        d0 = self.start_time - dt
        dstr = '{} to {}'.format(d0.strftime("%Y %m/%d %H:%M UTC"), self.start_time.strftime("%m/%d %H:%M UTC"))
        dset.attrs.update({'source':'Data Insertion {}'.format(dstr)})
        vm = concentration_montage.ConcentrationMontage(dset,columns=3)
        fig = vm.plotpage(pagenum)
        return vm

    def plot_concentration2(self,vloc):
        from utilhysplit.evaluation import concentration_montage
        from utilhysplit.evaluation.ensemble_tools import preprocess
        dt = datetime.timedelta(hours=3)
        timelist = [self.start_time+n*dt for n in [0,1,2,3,4,5,6]]
        dset = self.model_dset.sel(time=timelist)
        dset,dim = preprocess(dset) 
        print('HERE',dim)
        #dset = dset.mean(dim=dim)
        dset.attrs.update(self.get_attrs())
        if isinstance(vloc,(tuple,list,np.ndarray)):
           dset.attrs.update({'longitude':vloc[1],'latitude':vloc[0]})
        vm = concentration_montage.ConcentrationMontage(dset,columns=3)
        return vm


    def plot_height(self,vloc,thresh=0.2,unit='km'):
        sourcelist = self.get_sourcelist(time_previous=3)
        dt = datetime.timedelta(hours=3)
        timelist = [self.start_time+n*dt for n in [0,1,2,3,4,5]]
        dset = self.model_dset.sel(source=sourcelist).sel(time=timelist)
        wep.height_plot(dset,vlist=vloc,thresh=thresh,unit=unit)

    def get_iwxxm_forecast(self,problev):
        tp = 1
        if problev==50:
            #conc = self.forecast.Concentration
            sourcelist = self.get_sourcelist(time_previous=tp)
            enslist = self.get_enslist()
            print(sourcelist)
            conc = self.model_dset.sel(source=sourcelist,ens=enslist)
            conc = conc.median(dim='source').median(dim='ens') 
        elif problev==90:
            #conc = self.forecast.Concentration90 
            sourcelist = self.get_sourcelist(time_previous=tp)
            enslist = self.get_enslist()
            conc = self.model_dset.sel(source=sourcelist,ens=enslist)
            conc = conc.max(dim='source').max(dim='ens') 
        elif problev==0:
            sourcelist = self.get_sourcelist(time_previous=1)
            print(sourcelist)
            enslist = self.get_enslist()
            conc = self.model_dset.sel(source=sourcelist[-1],ens=enslist)
            print(conc.source.values)
            conc = conc.isel(ens=0) 
        return conc

    def polygon(self,thresh,problev):
        phash = {}
        for time in self.timelist:
            dset = self.get_iwxxm_forecast(problev)
            zlevs = np.arange(0,len(dset.z.values))
            dset = dset.sel(time=time) 
            top,bottom = topheight(dset,time=time,level=zlevs,thresh=thresh)
            top_poly = ensemble_polygons.HeightPolygons(cmap='viridis')
            top_poly.process(top,alpha=0.1)
            phash[time] = top_poly 
        return phash

    def plot_vaa_forecast(self,vloc,thresh=0.2,problev=50,plev=3):
        conc = self.get_iwxxm_forecast(problev)
        #dt = datetime.timedelta(hours=3)
        conc = conc.sel(time=self.timelist)
        vaa = polygon_plots.PlotVAA()
        vaa.setup()
        vaa.model = conc
 
        if plev==0:
            vaa.plot(vloc=vloc,thresh=thresh,ppp='top',plotheight=True,plotmass=False)

        elif plev==1:
            vaa.plot(vloc=vloc,thresh=thresh,ppp='bottom')
     
        elif plev==2:
            vaa.plot(vloc=vloc,thresh=thresh,ppp='top',plotheight=False,plotmass=True)

        # final plot
        elif plev==3:
            vaa.plot_one(vloc=vloc,thresh=thresh)

        return vaa



def pick_source(sss: str, 
                drange: list) -> list:
    """

    """
    new = []
    dlist = []
    for sname in sss:
        #sn = sname.replace('cdump.','cdump_')
        try:
             v = mdi.EmitName(sname)
        except:
             continue
        ddd = v.vhash['observation_date']
        if ddd >= drange[0] and ddd <= drange[1]:
           new.append(sname)
           dlist.append(ddd)
    return new, dlist


class ModelEvent:

    def __init__(self, eventid,volcano_name):

        self.eventdf = ModelEventDF(edf=pd.DataFrame())
        self.status = EventStatus(eventid,"Initialized")
        self.display = ModelEventDisplay(eventid)
        
        self.ndir = "./"  # directory where model output resides.
        self.edir = "./"  # directory where data insertion model output resides.
        self.idir = "./"  # directory where inversion algorithm model output resides.

        self.events ={} # dictionary where list of datasets with model output.

        self.volcano_name = volcano_name
 
        self.eventid = eventid


    @property
    def volcano_name(self):
        return self._volcano_name

    @volcano_name.setter
    def volcano_name(self,vname):
        self._volcano_name = fix_volc_name(vname)
        
    def check_for_DI(self):
        emit_dir = self.edir
        difiles = glob.glob(emit_dir + '/xrfile.*.nc')
        return difiles

    def get_model_runs(self,dilist=None):
        if not isinstance(dilist,list):
            dilist = self.check_for_DI()
        for di in dilist:
            if di not in self.events.keys():
               self.events[di] = xr.open_dataset(di)

    def find_event(self,start,end):
        elist = []
        for key in self.events.keys():
            dset = self.events[key]
            snames = dset.source.values
            # use data insertions up to dt h earlier.
            drange = [start, end]
            sourcelist,dlist = pick_source(snames, drange)
            if sourcelist:
               elist.append(key)
            #if sourcelist:
            #    print('creating forecast for {}'.format(key))
            #    forecast = ModelForcast(dset,self.eventid,start,end)
            #    return forecast
        return elist

    def get_dir(self, inp, verbose=False, make=True):
        """
        inp : dictionary with key VOLCAT_DIR
        set the directory from the dictionary inp.
        """
        # Events class
        if not isinstance(inp, dict):
            return None
        if "VOLCAT_DIR" not in inp.keys():
            logger.warning("get_dir method input does not contain VOLCAT_DIR")
            return None
        tdir = inp["VOLCAT_DIR"]
        if not self.volcano_name:
            self.set_volcano_name()
        if self.volcano_name != "Unknown":
            ndir = os.path.join(inp["VOLCAT_DIR"], self.volcano_name)
        else:
            ndir = inp["VOLCAT_DIR"]
        edir = os.path.join(ndir, "emitimes")
        idir = os.path.join(ndir, "inverse")

        if verbose:
            logger.info("Downloading to {}".format(ndir))
            logger.info("parallax corrected to {}".format(pdir))
            logger.info("emit times files to {}".format(edir))
        self.set_dir(ndir, edir, idir, make=make)
        return ndir

    def set_dir(self, data_dir, emit_dir, inv_dir=None, make=False):
        """
        sets the directories

        make : boolean : if TRUE then checks to see if directory exists and creates it if it doesn't.
        """
        if isinstance(data_dir, str):
            self.ndir = data_dir
            if make and not os.path.isdir(self.ndir):
                make_dir(self.ndir, None, verbose=True)
        if isinstance(emit_dir, str):
            self.edir = emit_dir
            if make and not os.path.isdir(self.edir):
                make_dir(self.edir, None, verbose=True)
        if isinstance(inv_dir, str):
            self.idir = inv_dir
            if make and not os.path.isdir(self.idir):
                make_dir(self.idir, None, verbose=True)

def forecast2QVA(dset, qva_filename, kwargs):
    """
    dset should have units of mg/m3.
    """
    enslist = dset.ens.values
    qva = ensemble_tools.ATLra(
        dset,
        enslist,
        sourcelist=None,
        threshlist=[0.2, 2, 5, 10],
        weights=None,
        **kwargs
    )
    encoding = {}
    encoding["zlib"] = True
    encoding["complevel"] = 9
    ehash = {}
    ehash["Concentration"] = encoding
    ehash["FrequencyOfExceedance"] = encoding
    qva.to_netcdf(qva_filename, encoding=ehash)
    return qva

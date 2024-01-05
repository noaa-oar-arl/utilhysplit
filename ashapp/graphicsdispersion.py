# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# graphicsdispersion.py - 
#
# TODO - THIS IS UNFINISHED AND NEEDS TO COPY GRAPHICS GENERATION IN ashbase.py
# -----------------------------------------------------------------------------
#
#
# -----------------------------------------------------------------------------


# from abc mport ABC, abstractmethod
import datetime
import glob
import logging
import os
import shutil
import subprocess
import sys
import zipfile

import numpy as np
import requests
import xarray as xr

from utilvolc.runhelper import Helper, JobFileNameComposer
from ashapp.ashruninterface import ModelOutputInterface
import ashapp.plotting_functions as plotf
from utilhysplit.evaluation import web_ensemble_plots as wep
from monetio.models import hysplit
from ashapp.cdump2xml import HysplitKml

logger = logging.getLogger(__name__)


# create_maptext
# get_maptext_info
# generate_kmz
# create_zipped_up_file



class GraphicsDispersion(ModelOutputInterface):
    ilist = [("WORK_DIR",'req'),
             ("jobid",'req'),
             ("jobname",'req'),
             ("MAP_DIR",'req'),
             ("HYSPLIT_DIR",'req'),
             ("CONVERT_EXE",'req'),
             ("graphicsResolution",'req'),
             ("mapBackground",'req'),
             ('graphicsResolution','req'),
             ('zip_compression_level','req')]


    def __init__(self,inp):
        self._inputlist = []
        self._outputlist = []
        self._fhash = {}
        self.inp = inp
        self.filelocator = inp
        #self._ncfile = ashnetcdf.HYSPLITAshNetcdf('TEMP')
        self._ncfile = None

    def ingest_model_output(self, modeloutput):
        """ modeloutput - class with ModelOutputInterface as base class
        """
        self.inputlist = modeloutput.outputlist
        self._ncfile = modeloutput._ncfile
    
    def load_ncfile(self,ncfile):
        self._ncfile = ashnetcdf.HYSPLITAshNetcdf(self._fhash['xrfile'])
        #self._ncfile = ncfile

    @property
    def ncfile(self):
        return self._ncfile

    @property
    def inputlist(self):
        return self._inputlist

    @inputlist.setter
    def inputlist(self, inlist):
        self._inputlist = inlist 

    # no setter
    @property
    def outputlist(self):
        outputlist = []
        for fff in self._fhash:
            if os.path.isfile(fff):
               outputlist.append(fff)
        #output.list.append(self._fhash['parxplot'])
        return outputlist


    def _get_filenames(self):
        fhash = {}
        # input files
        fhash['pardump'] = self.filelocator.get_pardump_filename()
        fhash['cdump'] = self.filelocator.get_cdump_filename()
        fhash['xrfile'] = self.filelocator.get_xrfile()

        # output files
        ptype = 'ps'
        # could use html or svg for svg files
        fhash['parxplot'] = self.filelocator.get_parxplot_filename(ptype=ptype)
        fhash['concplot'] = self.filelocator.get_concplot_filename(stage=0,frame=None,ptype=ptype)
        fhash['massplot'] = self.filelocator.get_massloading_filename(stage=0,frame=None,ptype=ptype)
        fhash['kmz'] = self.filelocator.get_kmz_filename()

        fhash['awips'] = self.filelocator.get_awips_filename(stage=0)

        #TODO
        # create montage
        # create allplots.pdf
        # create awips2 netcdf files
        # create massloading kmz file
        # GELABEL files?
        # zip files
        return fhash       
 

    def check(self):
       rval = True
       for fff in self._fhash:
            if not os.path.isfile(fff):
               logger.warning('Missing file {}'.format(fff))
               rval=False
       return rval

    def postprocess(self):
        stage=0

        #create_maptext(inp)

        # massloading plot
        inputname = self._fhash['cdump']
        outputname = self._fhash['massplot']
        plotf.create_massloading_plot(self.inp,inputname, outputname, conc_multiplier=1)

        # create the particle plot.
        pardump_filename = self._fhash['pardump']
        parxplot_outputname = self._fhash['parxplot']
        rval = plotf.create_parxplot(self.inp,pardump_filename, parxplot_outputname)

        # create the concentration plot
        cdump_filename = self._fhash['cdump']
        concplot = self._fhash['concplot']
        plotf.create_concentration_plot(self.inp, cdump_filename, concplot, conc_multiplier=1)

        # kmz file 
        self.generate_kmz() 

        # still need the montage and awips file.
        # also the maptext and gelabel.
        #create_concentration_montage(stage=stage) 
        self.make_awips_netcdf()

        return True

    def load_netcdf(self):
        if os.path.isfile(self._fhash['xrfile']):
           temp = xr.open_dataset(self._fhash['xrfile'])
           self._ncinput = temp[list(temp.data_vars.keys())[0]]
        else:
           logger.warning('could not find {}'.format(self._fhash['xrfile'])) 

    def get_massloading(self):
        cxra = self._ncfile.cxra
        self.mxra = hysplit.hysp_massload(cxra)
        return self.mxra

    def make_awips_netcdf(self):
        import cdump2netcdf
        jobid = self.inp['jobid']
        ghash = {}
        ghash["source_latitude"] = self.inp["latitude"]
        ghash["source_longitude"] = self.inp["longitude"]
        ghash["source_name"] = self.inp["VolcanoName"]
        ghash["emission_start"] = self.inp["start_date"]
        ghash["emission_duration_hours"] = self.inp["emissionHours"]
        mult=1
        # in mg
        #mult = self.get_conc_multiplier()
        mer = mult / 1e6 / 3600  # kg released in a second.
        ghash["MER"] = mer
        ghash["MER_unit"] = "kg/s"
        logger.debug("MULT value for awips {:3e}".format(mult))
        awipsname = self._fhash['awips']
        c2n = cdump2netcdf.Cdump2Awips(
            self._ncfile.cxra, awipsname, munit="mg", jobid=jobid, globalhash=ghash
        )
        awips_files = c2n.create_all_files()
        self._fhash['awipslist'] = awips_files
        compresslevel = self.inp['zip_compression_level']
        plotf.create_zipped_up_file(self._fhash['awips'], awips_files,compresslevel)


        # returns list of awips files that were created.

    def generate_kmz(self):
        jobid = self.inp['jobid']
        units = 'g/m2'
        mxra = self.get_massloading()
        levels = wep.set_levels(mxra)
        # use the mean massloading.
        # this will work for ensemble and deterministic runs.
        mxra = mxra.mean(dim='ens').mean(dim='source')
        legend_label = 'Mass Loading' 
        h2xml = HysplitKml(levels=levels,
                           sourcehash=self.inp,
                           units=units,
                           jobid=jobid,
                           legend_label=legend_label)

        attrs = self._ncfile.cxra.attrs
        kmlname = self._fhash['kmz'].replace('kmz','kml')
        h2xml.create(mxra,attrs,kmlname)


        compresslevel = self.inp['zip_compression_level']
        hdir = self.inp['HYSPLIT_DIR']
        efiles = plotf.get_kmz_files(hdir)
        efiles.append(kmlname)
        h2xml.create_kmz(compresslevel,efiles)
 
    @property
    def filelocator(self):
        #self._filelist = list(set(self._filelist))
        return self._filelocator

    @filelocator.setter
    def filelocator(self, finp):
        """
        input dictionary  with WORK_DIR, jobid, jobname
        or MetFileFinder object.
        """
        if isinstance(finp, dict):
            self._filelocator = JobFileNameComposer(
                finp["WORK_DIR"], finp["jobid"], finp["jobname"]
            )
            # get all the filenames here.
            self._fhash = self._get_filenames()
        else:
            # TO DO? add test for type.
            self._filelocator = finp


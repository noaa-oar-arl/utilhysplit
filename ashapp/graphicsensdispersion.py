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

from utilvolc.runhelper import Helper
from ashapp.filenamer import  JobFileNameComposer
from ashapp.ashruninterface import ModelOutputInterface
import ashapp.plotting_functions as plotf


logger = logging.getLogger(__name__)


# create_maptext
# get_maptext_info
# generate_kmz
# create_zipped_up_file



class GraphicsEnsembleDispersion(ModelOutputInterface):
    ilist = [("WORK_DIR",'req'),
             ("jobid",'req'),
             ("jobname",'req'),
             ("MAP_DIR",'req'),
             ("HYSPLIT_DIR",'req'),
             ("mapBackground",'req'),
             ("CONVERT_EXE",'req'),
             ('graphicsResolution','req')]


    def __init__(self,inp):
        self._inputlist = []
        self._outputlist = []
        self._fhash = {}
        self.inp = inp
        self.filelocator = inp


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
        output.list.append(self._fhash['parxplot'])
        return self._outputlist


    def _get_filenames(self):
        fhash = {}
        # input files
        fhash['pardump'] = self.filelocator.get_pardump_filename()
        fhash['cdump'] = self.filelocator.get_cdump_filename()

        # output files
        ptype = 'ps'
        # could use html or svg for svg files
        fhash['parxplot'] = self.filelocator.get_parxplot_filename(ptype=ptype)
        fhash['concplot'] = self.filelocator.get_concplot_filename(stage=0,frame=None,ptype=ptype)
        fhash['massplot'] = self.filelocator.get_massloading_filename(stage=0,frame=None,ptype=ptype)

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
       for filename in self._outputlist:
           if not os.path.isfile(filename):
              return False
       return True

    def postprocess(self):
        stage=0

        #create_maptext(inp)
        inputname = self._fhash['cdump']
        outputname = self._fhash['massplot']
        #plotf.create_massloading_plot(self.inp,inputname, outputname, conc_multiplier=1)

        # create the particle plot.
        pardump_filename = self._fhash['pardump']
        parxplot_outputname = self._fhash['parxplot']
        #rval = plotf.create_parxplot(self.inp,pardump_filename, parxplot_outputname)

        # create the concentration plot
        cdump_filename = self._fhash['cdump']
        concplot = self._fhash['concplot']
        print('CONCPLOT', concplot)
        #plotf.create_concentration_plot(self.inp, cdump_filename, concplot, conc_multiplier=1)
     
        kmlnames = [] # kml filenames
        kmzfilename = 'kmzfile.kmz' 
        #plotf.generate_kmz(self.inp['HYSPLIT_DIR'],kmlnames,kmzfilename,self.inp['zip_compression_level'])
        #create_concentration_montage(stage=stage) 
        #make_awips_netcdf()
        return True

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
        # after the filelocator is set, need to redo the filenames
        # self._get_filenames()


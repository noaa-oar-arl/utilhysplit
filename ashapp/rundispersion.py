#!/opt/Tools/anaconda3/envs/hysplit/bin/python
# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# ash_run.py - run HYSPLIT model on web and create plots
#
# 01 JUN 2020 (AMC) - adapted from locusts-run.py
# 09 DEC 2022 (AMC) - changed latspan to 180 and lonspan to 360 for global grid.
# 09 DEC 2022 (AMC) - changed numpar to 10000. this should not be hardwired.
# 28 FEB 2023 (AMC) - change to write_cxra so if netcdf file exists then set it as self.cxra
# 28 FEB 2023 (AMC) - check attributes added as stand alone function but still kept as static method.
#
# -----------------------------------------------------------------------------
# To run in offline mode use python ash_run.py -999
#
#
# -----------------------------------------------------------------------------


# from abc mport ABC, abstractmethod
import datetime
import logging
import os
import subprocess
import sys

import numpy as np

from utilhysplit import hcontrol
import utilhysplit.metfiles as metfile
from utilvolc.runhelper import ConcplotColors, Helper, JobFileNameComposer
from ashapp.ashruninterface import ModelRunInterface

# from runhandler import ProcessList
from utilvolc.volcMER import HT2unit
logger = logging.getLogger(__name__)


def FL2meters(flight_level):
    meters = flight_level * 100 / 3.28084
    return int(meters)
    # return int(np.ceil(flight_level / 10.0) * 10)


def make_chemrate(wdir):
    fname = "CHEMRATE.TXT"
    fstr = "1 2 0.01 1.0"
    fpath = os.path.join(wdir, fname)
    with open(fpath, "w") as fid:
        fid.write(fstr)

class RunDispersion(ModelRunInterface):

    def __init__(self, inp):

        self.JOBID = '999'

        # list of keywords the inp dictionary should contain.
        self._ilist = ['meteorologicalData','forecastDirectory','archivesDirectory',
                 'WORK_DIR','HYSPLIT_DIR','jobname','durationOfSimulation','latitude',
                 'longitude','bottom','top','emissionHours','rate','area','start_date',
                 'samplingIntervalHours','jobid']

        self._inp = {}
        self.inp = inp

        self.default_control = 'CONTROL.default'
        self.default_setup =  'SETUP.default'
        self.set_default_control()
        self.set_default_setup()
        self._control = hcontrol.HycsControl(
            fname=self.default_control,
            working_directory=self.inp["DATA_DIR"],
            rtype='dispersion',
        )
        self._setup = hcontrol.NameList(
            fname=self.default_setup, 
            working_directory=self.inp["DATA_DIR"]
        )

        self.metfilefinder = inp
        self.filelocator = inp
        
        # property with only getter method. 
        self._filelist = []

        self.so2 = False

    # no setter
    @property
    def control(self):
        return self._control

    # no setter
    @property
    def setup(self):
        return self._setup


    @property
    def inp(self):
        return self._inp

    @inp.setter
    def inp(self,inp):
        self._inp.update(inp)
        complete = True
        for iii in self._ilist:
            if iii not in self._inp.keys(): 
               logger.warning('Input does not contain {}'.format(iii))
               complete=False
        if 'jobid' in self._inp.keys():
            self.JOBID = self._inp['jobid'] 
        self.set_default_control()
        self.set_default_setup()
        if complete: logger.info('Input contains all fields')
        #self._inp.update(inp)


    @property
    def filelocator(self):
        self._filelist = list(set(self._filelist))
        return self._filelocator        

    @filelocator.setter
    def filelocator(self,inp):
        self._filelocator = JobFileNameComposer(
            inp["WORK_DIR"], inp['jobid'],  inp["jobname"]
        )

    @property
    def metfilefinder(self):
        return self._metfilefinder         

    @metfilefinder.setter
    def metfilefinder(self,inp):
        self._metfilefinder = metfile.MetFileFinder(inp["meteorologicalData"])
        self._metfilefinder.set_forecast_directory(inp["forecastDirectory"])
        self._metfilefinder.set_archives_directory(inp["archivesDirectory"])

    @property
    def filelist(self):
        return self._filelist


    def set_default_setup(self):
        if "setup" in self.inp.keys():
            self.default_setup = self.inp["setup"]
        #else:
        #    self.default_setup = 'SETUP.default'
        #logger.warning("Using setup {}".format(self.default_setup))

    def set_default_control(self):
        if "control" in self.inp.keys():
            self.default_control = self.inp["control"]
        #else:
        #    self.default_control = 'CONTROL.default'
        #logger.warning("Using control {}".format(self.default_control))
 
    def compose_control(self,stage,rtype):
        self.setup_basic_control(stage)
        self.additional_control_setup(stage)
        cdumpfiles = [x.outfile for x in  self.control.concgrids]
        self.control.write()
        return cdumpfiles

    def compose_setup(self, stage):
        self._setup = self.setup_setup(stage)
        self._setup.write(verbose=False)
        return -1

    def create_run_command(self):
        c = [os.path.join(self.inp["HYSPLIT_DIR"], "exec", "hycs_std"), str(self.JOBID)]
        return c

    def run_model(self,overwrite=False):
        # make control and setup files
        cdumpfiles = self.compose_control(stage=0, rtype="dispersion")
        self.compose_setup(stage=0)
        # start run and wait for it to finish..
        run=False
        for cdump in cdumpfiles:
            if not os.path.isfile(cdump): 
               run=True
               break
            logger.warning('CDUMP file already created {}'.format(cdump))
        if run or overwrite:
            c = self.create_run_command()
            logger.info("Running {} with job id {}".format("hycs_std", c[1]))
            Helper.execute(c)
        return -1

    # Additional methods below.
    def setup_setup(self, stage):
        duration = self.inp["durationOfSimulation"]
        # TO DO - may calculate particle number based on MER.
        setup = hcontrol.NameList(
            fname=self.default_setup, working_directory=self.inp["DATA_DIR"]
        )
        setup.read(case_sensitive=False)
        newname = self.filelocator.get_setup_filename(stage)
        pardumpname = self.filelocator.get_pardump_filename(stage)
        self._filelist.append(newname)
       
        setup.rename(name=newname, working_directory=self.inp["WORK_DIR"])
        self._filelist.append(pardumpname)
        setup.add("poutf", pardumpname)

        # setup.add("numpar", "2000")
        # setup.add("numpar", "20000")
        setup.add("numpar", "10000")
        # setup.add("maxpar", "500000")

        # base frequency of pardump output on length of simulation.
        if duration < 6:
            setup.add("ncycl", "1")
            setup.add("ndump", "1")
        elif duration < 12:
            setup.add("ncycl", "2")
            setup.add("ndump", "2")
        else:
            setup.add("ndump", "3")
            setup.add("ncycl", "3")

        keys = list(setup.nlist.keys())
        self.so2 = False
        if "ichem" in keys:
            if int(setup.nlist["ichem"]) == 2:
                make_chemrate(self.inp["WORK_DIR"])
                logger.warning("creating chemrate file for SO2")
                self.so2 = True
        return setup


    def setup_basic_control(self, stage=0, rtype="dispersion"):
        """
        inp : dictionary with inputs
        used for both dispersion and trajectory runs.
        """
        # logger.debug("Setting up control stage {}".format(stage))
        duration = self.inp["durationOfSimulation"]
        stime = self.inp["start_date"]
        metfiles = self.metfilefinder.find(stime, duration)
        newname = self.filelocator.get_control_filename(stage)

        self._control = hcontrol.HycsControl(
            fname=self.default_control,
            working_directory=self.inp["DATA_DIR"],
            rtype='dispersion',
        )
        self._control.read()
        self._control.rename(newname, self.inp["WORK_DIR"])
        # set start date
        self._control.date = stime

        # add met files
        self._control.remove_metfile(rall=True)
        #metfiles = self.metfilefinder.find(stime, duration)
        for mfile in metfiles:
            logger.debug(os.path.join(mfile[0], mfile[1]))
            self._control.add_metfile(mfile[0], mfile[1])
        if not metfiles:
            self.update_run_status(self.JOBID, "TERMINATED. missing met files")
            sys.exit()
        # add duration of simulation
        self._control.add_duration(duration)
        return self._control

    def additional_control_setup(self,  stage=0):
        # for dispersion control file
        # if setting levels here then need to use the set_levels
        # function and also adjust the labeling for the levels.

        # round to nearest latitude longitude point.
        # this is for better matching with gridded volcat data.
        lat = self.inp["latitude"]
        lon = self.inp["longitude"]
        vent = self.inp["bottom"]
        height = self.inp["top"]
        emission = self.inp["emissionHours"]
        rate = self.inp["rate"]
        area = self.inp["area"]

        # add location of eruption
        # with uniform vertical line source
        self._control.remove_locations()
        self._control.add_location((lat, lon), vent, rate=rate, area=area)
        self._control.add_location((lat, lon), height, rate=0.0001, area=area)
        # rename cdump file
        self._control.concgrids[0].outdir = self.inp["WORK_DIR"]
        self._control.concgrids[0].outfile = self.filelocator.get_cdump_filename(stage)
        self._filelist.append(self._control.concgrids[0].outfile)
        logger.info("{}".format(self._control.concgrids[0].outfile))
        # center concentration grid near the volcano
        self._control.concgrids[0].centerlat = int(lat)
        self._control.concgrids[0].centerlon = int(lon)
        # set a global grid
        self._control.concgrids[0].latspan = 180.0
        self._control.concgrids[0].lonspan = 360.0
        # set the levels
        self.set_levels(self._control.concgrids[0])

        # add emission duration
        rate_test = 0
        if not self.so2:
            for spec in self._control.species:
                spec.duration = emission
                rate_test += spec.rate
        # for so2, the sulfate emission has 0 emissions.
        else:
            spec = self._control.species[0]
            spec.duration = emission
            spec2 = self._control.species[1]
            spec2.duration = 0
            rate_test = spec.rate + spec2.rate

        if np.abs(1 - rate_test) > 0.01:
            logger.warning("Non unit emision rate {}".format(rate_test))

        # TO DO
        # set sample start to occur on the hour.
        stime = self.inp["start_date"]
        # if start past the half hour mark, set sample start at next hour.
        if stime.minute > 30:
            stime = stime + datetime.timedelta(hours=1)
        sample_start = stime.strftime("%y %m %d %H 00")
        # logger.debug("Setting sample start {}".format(sample_start))
        self._control.concgrids[0].sample_start = sample_start

        self._control.concgrids[0].interval = (self.inp["samplingIntervalHours"], 0)
        self._control.concgrids[0].sampletype = -1

    def set_qva_levels(self):
        # every 5000 ft (FL50 chunks)
        # Approx 1.5 km.
        levlist_fl = range(50, 650, 50)
        levlist = [FL2meters(x) for x in levlist_fl]
        rlist = []
        plev = "SFC"
        for lev in levlist_fl:
            nlev = "FL{}".format(lev)
            rlist.append("{} to {}".format(plev, nlev))
            plev = nlev
        logger.debug(rlist)
        return levlist, rlist

    def set_levels(self, cgrid):
        levlist, rlist = self.set_qva_levels()
        cgrid.set_levels(levlist)



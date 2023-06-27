#!/opt/Tools/anaconda3/envs/hysplit/bin/python
# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# rundispersion.py - run HYSPLIT dispersion model
#
# 05 Jun  2023 (AMC) -
#
# -----------------------------------------------------------------------------
# Run specifically for traditional volcanic ash with line source from vent
# to input plume height.
#
# TODO - look at how SO2 functionality is incoporated. This needs to be updated.
# TODO - update how vertical levels are specified.
# -----------------------------------------------------------------------------


# from abc mport ABC, abstractmethod
import datetime
import logging
import os

import numpy as np

from utilhysplit import hcontrol
from utilhysplit.metfiles import MetFileFinder
from utilvolc.runhelper import Helper, JobFileNameComposer
from ashapp.ashruninterface import ModelRunInterface
import ashapp.utildatainsertion as utildi

# from runhandler import ProcessList
# from utilvolc.volcMER import HT2unit

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


class RunEmitTimes(ModelRunInterface):
    def __init__(self, inp):
    """
    Build a run based on an emittimes file and inputs
    """
        # 16 instance attributes  may be too many?

        self.JOBID = "999"
        # list of keywords the inp dictionary should contain.
        self._ilist = [
            "meteorologicalData",
            "forecastDirectory",
            "archivesDirectory",
            "WORK_DIR",
            "HYSPLIT_DIR",
            "jobname",
            "durationOfSimulation",
            "latitude",
            "longitude",
            "emitfile",
            "emissionHours",
            "rate",
            "area",
            "start_date",
            "samplingIntervalHours",
            "jobid",
        ]

        self._inp = {}
        self.inp = inp

        self._status = "Initialized"
        self._history = [self._status]

        self.default_control = "CONTROL.default"
        self.default_setup = "SETUP.default"
        self.set_default_control()
        self.set_default_setup()
        self._control = hcontrol.HycsControl(
            fname=self.default_control,
            working_directory=self.inp["DATA_DIR"],
            rtype="dispersion",
        )
        self._setup = hcontrol.NameList(
            fname=self.default_setup, working_directory=self.inp["DATA_DIR"]
        )
        self._metfilefinder = MetFileFinder(inp["meteorologicalData"])
        self._metfilefinder.set_forecast_directory(inp["forecastDirectory"])
        self._metfilefinder.set_archives_directory(inp["archivesDirectory"])

        # self._filehash and self_filelist are
        # set in the filelocator setter method.
        self._filehash = {}
        self._filelist = []
        self.filelocator = inp
        self.so2 = False

    @staticmethod
    def help():
        helphash = {}

    # no setter
    @property
    def status(self):
        return (self._status, self._history)

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
    def inp(self, inp):
        self._inp.update(inp)
        complete = True

        if "jobid" in self._inp.keys():
            self.JOBID = self._inp["jobid"]

        add_inp={}
        if 'start_date' in inp and 'DI_start_date' not in inp:
           add_inp['DI_start_date'] = inp['start_date'] 
        if 'emissionHours' in inp and 'DIrunHours' not in inp:
           add_inp['DIrunHours'] = inp['emissionHours'] 
        if 'samplingIntervalHours' not in inp: 
           add_inp['samplingIntervalHours'] =  1
        self._inp.update(add_inp)

        for iii in self._ilist:
            if iii not in self._inp.keys():
                logger.warning("Input does not contain {}".format(iii))
                complete = False


        self.set_default_control()
        self.set_default_setup()
        if complete:
            logger.info("Input contains all fields")
        # self._inp.update(inp)

    @property
    def filelocator(self):
        self._filelist = list(set(self._filelist))
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
        else:
            # TO DO? add test for type.
            self._filelocator = finp
        # after the filelocator is set, need to redo the filenames
        self._get_filenames()

    @property
    def metfilefinder(self):
        return self._metfilefinder

    @metfilefinder.setter
    def metfilefinder(self, metinp):
        """
        input dictionary with meteorologicalData, forecastDirectory, archivesDirectory
        or MetFileFinder object.
        """
        if isinstance(metinp, dict):
            self._metfilefinder = MetFileFinder(metinp["meteorologicalData"])
            self._metfilefinder.set_forecast_directory(metinp["forecastDirectory"])
            self._metfilefinder.set_archives_directory(metinp["archivesDirectory"])
        else:
            # TO DO? add test for type.
            self._metfilefinder = metinp

    @property
    def filelist(self):
        cdump = self._filehash["cdump"]
        
        metfile = self._metfilefinder.metid
        if isinstance(self._metfilefinder.suffix,str):
           metfile +=  self._metfilefinder.suffix
        eloc = self.inp['emitfile'].split('/')
        eloc = eloc[-1]
        source = eloc.replace('EMIT','')
        self._filelist = [(cdump,source,metfile)]
        return self._filelist

    @property
    def filehash(self):
        return self._filehash

    def _get_filenames(self):
        fhash = {}
        fhash["control"] = self.filelocator.get_control_filename()
        fhash["setup"] = self.filelocator.get_setup_filename()
        fhash["pardump"] = self.filelocator.get_pardump_filename()
        fhash["cdump"] = self.filelocator.get_cdump_filename()
        self._filehash = fhash

    def set_default_setup(self):
        if "setup" in self.inp.keys():
            self.default_setup = self.inp["setup"]
        # else:
        #    self.default_setup = 'SETUP.default'
        # logger.warning("Using setup {}".format(self.default_setup))

    def set_default_control(self):
        if "control" in self.inp.keys():
            self.default_control = self.inp["control"]
        # else:
        #    self.default_control = 'CONTROL.default'
        # logger.warning("Using control {}".format(self.default_control))

    def compose_control(self, rtype):
        self._setup_basic_control()
        self._additional_control_setup()
        cdumpfiles = [x.outfile for x in self.control.concgrids]
        self.control.write()
        if os.path.isfile(self.filehash["control"]):
            self._history.append("CONTROL exists")
        else:
            self._status = "FAILED to write control file"
            self._history.append(self._status)
        return cdumpfiles

    def compose_setup(self):
        self._setup = self.setup_setup()
        self._setup.write(verbose=False)
        if os.path.isfile(self.filehash["setup"]):
            self._history.append("SETUP exists")
        else:
            self._status = "FAILED to write setup file"
            self._history.append(self._status)
        # return cdumpfiles

    def create_run_command(self):
        command = [
            os.path.join(self.inp["HYSPLIT_DIR"], "exec", "hycs_std"),
            str(self.filelocator.job),
        ]
        return command

    def check_status(self):
        if os.path.isfile(self.filelist["cdump"]):
            complete = "COMPLETE cdump exists"
            self._status = complete
            if self._history[-1] != complete:
                self._history.append(self._status)

    def run(self, overwrite=False):
        command = self.run_model(overwrite=overwrite)
        if isinstance(command,list):
            logger.info("execute {}".format(type(command)))
            Helper.execute(command)
        else:
            logger.info("No run to execture")

    def run_model(self, overwrite=False):
        # make control and setup files
        cdumpfiles = self.compose_control(rtype="dispersion")
        self.compose_setup()
        # start run and wait for it to finish..
        run = False
        command = None
        for cdump in cdumpfiles:
            if not os.path.isfile(cdump):
                run = True
                break
            logger.warning("CDUMP file already created {}".format(cdump))
            if not overwrite:
                self._status = "COMPLETE cdump exists"
                self._history.append(self._status)
        if run or overwrite:
            command = self.create_run_command()
        return command

    # Additional methods below.
    def setup_setup(self):
        duration = self.inp["durationOfSimulation"]
        # TO DO - may calculate particle number based on MER.
        setup = hcontrol.NameList(
            fname=self.default_setup, working_directory=self.inp["DATA_DIR"]
        )
        setup.read(case_sensitive=False)

        setup.rename(
            name=self.filehash["setup"], working_directory=self.inp["WORK_DIR"]
        )
        setup.add("poutf", self.filehash["pardump"])

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

        eloc = self.inp['emitfile'].split('/')
        eloc = eloc[-1]
        setup.add('efile',eloc)

        #keys = list(setup.nlist.keys())
        #self.so2 = False
        #if "ichem" in keys:
        #    if int(setup.nlist["ichem"]) == 2:
        #        make_chemrate(self.inp["WORK_DIR"])
        #        logger.warning("creating chemrate file for SO2")
        #        self.so2 = True

        return setup

    def _setup_basic_control(self, rtype="dispersion"):
        """
        inp : dictionary with inputs
        used for both dispersion and trajectory runs.
        """
        inp = utildi.read_emittimes(self.inp['emitfile'])
        self._inp.update(inp)

        duration = self.inp["durationOfSimulation"]
        stime = self.inp["start_date"]
        metfiles = self.metfilefinder.find(stime, duration)
        self._control = hcontrol.HycsControl(
            fname=self.default_control,
            working_directory=self.inp["DATA_DIR"],
            rtype=rtype,
        )
        self._control.read()
        self._control.rename(self.filehash["control"], self.inp["WORK_DIR"])
        # set start date
        self._control.date = stime

        # add met files
        self._control.remove_metfile(rall=True)
        # metfiles = self.metfilefinder.find(stime, duration)
        for mfile in metfiles:
            logger.debug(os.path.join(mfile[0], mfile[1]))
            self._control.add_metfile(mfile[0], mfile[1])
        if not metfiles:
            logger.warning("TERMINATED. missing met files")
            self._status = "FAILED. missing meteorological files"
            self._history.append(self._status)
            # self.update_run_status(self.JOBID, "TERMINATED. missing met files")
        # add duration of simulation
        self._control.add_duration(duration)
        return self._control

    def _additional_control_setup(self):
        # for dispersion control file
        # if setting levels here then need to use the set_levels
        # function and also adjust the labeling for the levels.

        # add as many dummy locations as in emit-times file
        nlocs = self.inp["nlocs"]
        self._control.remove_locations()
        lat = np.floor(self.inp["latitude"])
        lon = np.floor(self.inp["longitude"])
        vent = self.inp["bottom"]
        for loc in np.arange(0, nlocs):
            self._control.add_location((lat, lon), vent, rate=0, area=0)


        # rename cdump file
        self._control.concgrids[0].outdir = self.inp["WORK_DIR"]
        self._control.concgrids[0].outfile = self.filehash["cdump"]
        logger.info("output file set as {}".format(self._control.concgrids[0].outfile))
        # center concentration grid near the volcano
        self._control.concgrids[0].centerlat = int(lat)
        self._control.concgrids[0].centerlon = int(lon)
        # set a global grid
        self._control.concgrids[0].latspan = 180.0
        self._control.concgrids[0].lonspan = 360.0
        # set the levels
        self.set_levels(self._control.concgrids[0])

        # emission durations and rates in control should be 0
        # all emissions come from the emitimes file.
        for spec in self._control.species:
            spec.duration = 0
            spec.rate = 0.0

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

    @staticmethod
    def set_qva_levels():
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




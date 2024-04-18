#!/opt/Tools/anaconda3/envs/hysplit/bin/python
# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# runtrajectory.py - run HYSPLIT dispersion model
#
# 19 Jul  2023 (AMC) - created from rundispersion
#
# -----------------------------------------------------------------------------
# Run trajectories
# -----------------------------------------------------------------------------

# TODO created from rundispersion. needs to be checked over.

# from abc mport ABC, abstractmethod
import datetime
import logging
import os

import numpy as np

from utilhysplit import hcontrol
from utilhysplit.metfiles import MetFileFinder
from utilvolc.runhelper import Helper, is_input_complete
from ashapp.filenamer import  JobFileNameComposer
from ashapp.ashruninterface import ModelRunInterface

# from runhandler import ProcessList
# from utilvolc.volcMER import HT2unit

logger = logging.getLogger(__name__)


def round_start_time(stime, mround=5):
    if isinstance(stime, np.datetime64):
        # Convert numpy.datetime64 to datetime.datetime
        stime = stime.astype('datetime64[m]').astype(datetime.datetime)

    # round down
    minutes = int(stime.minute - stime.minute%mround)
    newtime = datetime.datetime(stime.year, stime.month, stime.day, stime.hour, minutes)
    return newtime

class RunTrajectory(ModelRunInterface):
    ilist = [
        ("meteorologicalData",'req'),
        ("forecastDirectory",'req'),
        ("archivesDirectory",'req'),
        ("WORK_DIR",'req'),
        ("HYSPLIT_DIR",'req'),
        ("durationOfSimulation",'req'),
        ("latitude",'req'),
        ("longitude",'req'),
        ("start_date",'req'),
        ("jobid",'req'),
        ("jobname",'req')
    ]


    def __init__(self, inp, trajgenerator):
        """
        A trajectory run from inputs.

        durationOfSimulation is number or hours for simulation to run.
                             A negative number will result in backwards trajectories.

        trajgenerator is a generator function which yields latitude, longitude, height in meters.
        """

        # 16 instance attributes  may be too many?

        self.JOBID = "999"
        # list of keywords the inp dictionary should contain.

        self._inp = {}
        self.inp = inp

        self._status = "Initialized"
        self._history = [self._status]

        self.default_control = "CONTROL.traj"
        self.default_setup = "SETUP.traj"
        self.set_default_control()
        self.set_default_setup()
        self._control = hcontrol.HycsControl(
            fname=self.default_control,
            working_directory=self.inp["DATA_DIR"],
            rtype="trajectory",
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
        self.trajgenerator = trajgenerator

    @staticmethod
    def help():
        helphash = {}
        helphash["jobname"] = "str : user input string to identify run"
        helphash["jobid"] = "usually an integer assigned by the system for the job"
        helphash["emissionHours"] = "int : time for emissions to last "

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

        complete = is_input_complete(self.ilist, self._inp)
            
        if "jobid" in self._inp.keys():
            self.JOBID = self._inp["jobid"]
        self.set_default_control()
        self.set_default_setup()
        if complete:
            logger.info("Input contains all fields")
        else:
            logger.warning('Input incomplete')
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
        tdump = self._filehash["tdump"]
        metfile = self._metfilefinder.metid
        if isinstance(self._metfilefinder.suffix, str):
            if self.metfilefinder.suffix not in metfile:
                metfile += self._metfilefinder.suffix
        self._filelist = [(tdump, metfile)]
        return self._filelist

    @property
    def filehash(self):
        return self._filehash

    # no level setter needed for trajectory runs.
    @property
    def levelsetter(self):
        return None

    def _get_filenames(self):
        fhash = {}
        fhash["control"] = self.filelocator.get_control_filename()
        fhash["setup"] = self.filelocator.get_setup_filename()
        fhash["tdump"] = self.filelocator.get_tdump_filename()
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
        self.control.write()
        if os.path.isfile(self.filehash["control"]):
            self._history.append("CONTROL exists")
        else:
            self._status = "FAILED to write control file"
            self._history.append(self._status)
        return os.path.join(self._control.outdir, self._control.outfile)

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
            os.path.join(self.inp["HYSPLIT_DIR"], "exec", "hyts_std"),
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
        if isinstance(command, list):
            logger.info("execute {}".format(type(command)))
            Helper.execute(command)
        else:
            logger.info("No run to execute")

    def run_model(self, overwrite=False):
        # make control and setup files
        self.compose_setup()
        tdump = self.compose_control(rtype="dispersion")
        # start run and wait for it to finish..
        run = False
        command = None
        if not os.path.isfile(tdump):
            run = True
        else:
            logger.warning("TDUMP file already created {}".format(tdump))
        if not overwrite:
           self._status = "COMPLETE tdump exists"
           self._history.append(self._status)
        if run or overwrite:
            command = self.create_run_command()
        return command

    # Additional methods below.
    def setup_setup(self):
        setup = hcontrol.NameList(
            fname=self.default_setup, working_directory=self.inp["DATA_DIR"]
        )
        setup.read(case_sensitive=False)
       
        setup.rename(
            name=self.filehash["setup"], working_directory=self.inp["WORK_DIR"]
        )
        keys = list(setup.nlist.keys())
        return setup


    def _setup_basic_control(self, rtype="trajectory"):
        """
        inp : dictionary with inputs
        used for both dispersion and trajectory runs.
        """
        duration = self.inp["durationOfSimulation"]

        # if time step set then need to make sure start time is compatible
        if 'delt' in self._setup.keys:
           rtime = self._setup.nlist['delt']
           stime = round_start_time(self.inp["start_date"])
        else:
           stime = self.inp["start_date"]
        print('STIME', stime)
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
        self._control.add_duration(duration)
        return self._control

    def _additional_control_setup(self):
        self._control.remove_locations()
        [self._control.add_location((x['latitude'], x['longitude']), x['height']) for x in self.trajgenerator]
        self._control.outdir = self.inp['WORK_DIR']
        self._control.outfile = self.filelocator.get_tdump_filename(0)


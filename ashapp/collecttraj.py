import datetime
import logging
import os

import ashapp.utildatainsertion as udi
from ashapp.ashruninterface import ModelCollectionInterface
from ashapp.runtrajectory import RunTrajectory
from utilhysplit.runhandler import ProcessList
from utilvolc.runhelper import is_input_complete
from ashapp.trajectory_generators import timegenerate_traj_from_obsdf
from ashapp.trajectory_generators import timegenerate_height_traj_from_obsdf
from math import floor
logger = logging.getLogger(__name__)


# 2023 Dec 4 (amc) setup method changed to utilize generate_height_traj_from_obsdf generator
#                  old functionality retained in setup2 method.


class CollectTrajectory(ModelCollectionInterface):
    """
    Runs multiple trajectory runs.
    """

    # list of required (req) and optional (opt) inputs
    ilist = [
        ("WORK_DIR", "req"),
        ("jobname", "req"),
        ("start_date", "req"),
        ("jobid", "req"),
    ]
    ilist.extend(RunTrajectory.ilist)

    def __init__(self, inp, jobid):
        self.JOBID = jobid

        self._inp = {}
        self.inp = inp

        self._filehash = {}
        self._filelist = []

        self._status = {"MAIN": "INITIALIZED"}

    @property
    def filehash(self):
        return self._filehash

    @property
    def filelist(self):
        return self._filelist

    @property
    def inp(self):
        return self._inp

    @inp.setter
    def inp(self, inp):
        self._inp.update(inp)

        if "jobid" in self._inp.keys():
            self.JOBID = self._inp["jobid"]

        complete = is_input_complete(self.ilist, self._inp)
        if not complete:
            self.status = "FAILED inputs incomplete"

    @property
    def status(self):
        return self._status

    @status.setter
    def status(self, status):
        self._status = status

    def setup(self, overwrite):
        from ashapp.runtrajectory import RunTrajectory
        command_list = []
        inp = self._inp.copy()
        csvname = '{}/{}.csv'.format(inp['WORK_DIR'],self.JOBID)
        if not os.path.isfile(csvname):
            edate = inp['start_date'] + datetime.timedelta(hours=1)
            daterange = [inp['start_date'], edate] 
            csvname = udi.find_btraj_file(inp['WORK_DIR'],daterange)
            csvname = csvname[-1]
            print('CSVNAME', csvname, daterange)
        #csvname = '{}/{}.csv'.format(inp['WORK_DIR'],self.JOBID)
        iii=0
        # first level to start trajectories at in km above msl. 
        #  use at or slightly below vent height. 
        minht = floor(inp['bottom'])/1000.0
        maxht = floor(inp['top'])/1000.0
        for time, trajgen in timegenerate_height_traj_from_obsdf(csvname,minht=minht,maxht=maxht):
            inp['start_date'] = time
            inp['jobid'] = '{}.{}'.format(self.JOBID, iii)
            run = RunTrajectory(inp,trajgen)
            command = run.run_model(overwrite=False)
            if isinstance(command,(str,list)): command_list.append(command) 
            self._filelist.extend(run.filelist)
            iii += 1
        return command_list

    def setup2(self, overwrite):
        from ashapp.runtrajectory import RunTrajectory
        command_list = []
        inp = self._inp.copy()
        csvname = '{}/{}.csv'.format(inp['WORK_DIR'],self.JOBID)
        iii=0
        for time, trajgen in timegenerate_traj_from_obsdf(csvname):
            inp['start_date'] = time
            inp['jobid'] = '{}.{}'.format(self.JOBID, iii)
            run = RunTrajectory(inp,trajgen)
            command = run.run_model(overwrite=False)
            if isinstance(command,(str,list)): command_list.append(command) 
            self._filelist.extend(run.filelist)
            iii += 1
        return command_list

    def run(self, overwrite=False):
        import time
        command_list = self.setup(overwrite)
        processhandler = ProcessList()
        processhandler.pipe_stdout()
        processhandler.pipe_stderr()
        # suffix = gefs_suffix_list()
        iii=-1
        for iii, command in enumerate(command_list):
            logger.info("Running {} with job id {}".format("hycs_std", command[1]))
            processhandler.startnew(command, self.inp["WORK_DIR"], descrip=str(iii))
            # wait 5 seconds to avoid runs trying to access ASCDATA.CFG files at the
            # same time.
            time.sleep(5)
            num_proces = processhandler.checkprocs()
            total_time = 0
            logger.info('Number of processes running {}'.format(num_proces))
            while num_proces >= 10:
                  num_proces = processhandler.checkprocs()
                  max_time = 10 * 60
                  seconds_to_wait = 10
                  total_time += seconds_to_wait
                  time.sleep(seconds_to_wait)
                  if total_time > max_time: 
                     logger.info('max time reached')
                     break
        if iii==-1: logger.warning('NO COMMANDS EXECUTED')
        # wait for runs to finish
        done = False
        seconds_to_wait = 30
        total_time = 0
        # 60 minutes.
        max_time = 60 * 60
        # max_time = 0.5*60
        while not done:
            num_proces = processhandler.checkprocs()
            if num_proces == 0:
                done = True
            time.sleep(seconds_to_wait)
            total_time += seconds_to_wait
            if total_time > max_time:
                processhandler.checkprocs()
                processhandler.killall()
                logger.warning("HYSPLIT run Timed out")
                done = True



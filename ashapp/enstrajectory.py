# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# ashtrajectory.py - run HYSPLIT model on web and create plots
#
# 01 JUN 2020 (AMC) - adapted from locusts-run.py
# 28 FEB 2023 (AMC) - default is for 10 evenly spaced trajectories from vent to top height.
# -----------------------------------------------------------------------------
# To run in offline mode use python ash_run.py -777
# -----------------------------------------------------------------------------
import datetime
import logging
import os
import time

import numpy as np

from monetio.models import hytraj

from ashapp.ashbase import AshRun
from utilvolc.runhelper import Helper
from ashapp.ashtrajectory import TrajectoryAshRun
from utilhysplit import metfiles
from utilhysplit.runhandler import ProcessList

logger = logging.getLogger(__name__)


def print_usage():
    print(
        """\
USAGE: to test use python ash_run.py run -777
"""
    )


# Base class is AshRun
class EnsTrajectoryRun(TrajectoryAshRun):

    def __init__(self, JOBID):
        super().__init__(JOBID)
        self.ens_suffix_list = metfiles.gefs_suffix_list()
        self.number_of_members = len(self.ens_suffix_list)

    def get_maptext_info(self):
        maptexthash = {}
        rstr = 'HYSPLIT ensemble trajectory run'
        maptexthash['run_description'] = rstr
        maptexthash['infoc'] = 'GEFS'
        return maptexthash


    def run_model(self):
        logger.info('RUNNING MODEL')
        processhandler = ProcessList()
        processhandler.pipe_stdout()
        processhandler.pipe_stderr()
        for stage, suffix in enumerate(self.ens_suffix_list):
            logger.info('WORKING ON {}'.format(suffix))
            cname = self.filelocator.get_tdump_filename(stage+1)
            if os.path.isfile(cname):
               logger.info('File exists {}'.format(cname)) 
               continue
            self.metfilefinder.set_ens_member("." + suffix)
            # make control and setup files
            self.compose_control(stage=stage+1, rtype="trajectory")
            self.compose_setup(stage=stage+1)
            run_suffix = self.filelocator.get_control_suffix(stage+1)
            # start run and wait for it to finish..
            cproc = [os.path.join(self.inp["HYSPLIT_DIR"], "exec", "hyts_std"), run_suffix]
          
            logger.info("Running {} with job id {}".format("hyts_std", cproc[1]))
            processhandler.startnew(cproc,self.inp["WORK_DIR"],descrip=suffix)
            time.sleep(5)
        done = False
        seconds_to_wait = 20
        total_time = 0
        max_time = 60*60
        while not done:
              num_procs = processhandler.checkprocs()
              if num_procs==0:
                 done=True
              time.sleep(seconds_to_wait)
              if total_time > max_time:
                 processhandler.checkprocs()
                 processhandler.killall()
                 logger.warning('HYSPLIT runs timed out')
                 done=True   

    def write_cxra(self):
        """
        write csv file with all trajectories
        """
        logger.info('Writing csv file with tdump information')
        outname = 'tdump.{}.csv'.format(self.JOBID)
        flist = []
        taglist = []
        for stage, suffix in enumerate(self.ens_suffix_list):
            fname = self.filelocator.get_tdump_filename(stage+1)
            flist.append(fname)
            taglist.append(suffix)
        tdumpdf = hytraj.combine_dataset(flist,taglist,renumber=False)
        tdumpdf.to_csv(outname) 

    def after_run_check(self, update=False):
        # Check for the tdump/cdump file
        for stage, suffix in enumerate(self.ens_suffix_list):
            fname = self.filelocator.get_tdump_filename(stage+1)
            if not os.path.isfile(fname):
               return False
        return True   

    def create_trajectory_plot(self,stage):
        # creates one plot for each member
        for stage, suffix in enumerate(self.ens_suffix_list):
            super().create_trajectory_plot(stage=stage+1)

 

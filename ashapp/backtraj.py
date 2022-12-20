import datetime
import logging
import os
import time

import pandas as pd
from monetio.models import hytraj
from utilhysplit import hcontrol
from utilhysplit.runhandler import ProcessList

from ashapp.ashtrajectory import TrajectoryAshRun

logger = logging.getLogger(__name__)

"""
BackTraj class utilizes a csv file with starting points for the back trajectories.
csv file can be generated from the VOLCAT data.
"""


def combine_traj(fnames, csvfile=None):
    """
    combined trajectories in different files into one dataframe.
    """
    trajlist = []
    if csvfile:
        weightcsv = pd.read_csv(csvfile)
    for iii, fnn in enumerate(fnames):
        df1 = hytraj.open_dataset(fnn)
        temp = fnn.split(".")
        trajnum = int(temp[-1])
        df1["trajnum"] = trajnum
        if csvfile:
            temp = weightcsv.loc[trajnum]
            weight = temp.massI
        else:
            weight = 1
        df1["weight"] = weight
        trajlist.append(df1.copy())
    trajdf = pd.concat(trajlist)
    return trajdf


class BackTraj(TrajectoryAshRun):
    def __init__(self, JOBID):
        super().__init__(JOBID)
        self.obsdf = pd.DataFrame()  # dataframe with obs to initiate trajectories from.

    # def __init__(self, wdir='./', controldefault='CONTROL.0'):
    #    self.wdir = wdir
    #    self.default = controldefault

    #    self.control = hcontrol.HycsControl(fname=controldefault, rtype='trajectory', working_directory=wdir)
    #    self.control.read()

    def additional_control_setup(self, control, suffix):
        control.add_sdate(self.inp["start_date"])
        alt = self.inp["top"]
        lat = self.inp["latitude"]
        lon = self.inp["longitude"]

        control.remove_locations()

        control.add_location(latlon=[lat, lon], alt=alt)

        control.rename("CONTROL.{}".format(suffix))
        control.outdir = self.inp["WORK_DIR"]
        control.outfile = self.filelocator.get_tdump_filename(suffix)
        control.outfile = "tdump.{}".format(suffix)

    def read_obsdf(self):
        dname = self.inp["WORK_DIR"]
        fname = "{}.csv".format(self.JOBID)
        obsdf = pd.read_csv(os.path.join(dname, fname), parse_dates=["time"])
        self.obsdf = obsdf

    def run_model(self):
        processhandler = ProcessList()
        processhandler.pipe_stdout()
        processhandler.pipe_stderr()
        self.read_obsdf()
        max_procs = 10
        for iii, row in enumerate(self.obsdf.iterrows()):
            tdumpfile = "tdump.{}".format(iii)
            if not os.path.isfile(tdumpfile):
                self.inp["top"] = row[1].height * 1000
                self.inp["latitude"] = row[1].lat
                self.inp["longitude"] = row[1].lon
                stime = row[1].time
                # import sys
                # time = datetime.datetime.strptime("%Y-%m-%d %H:%M:%S", row[1].time)
                # print(time)
                self.inp["start_date"] = stime
                self.compose_control(iii, rtype="trajectory")
                self.compose_setup(iii)
                run_suffix = str(iii)
                tdumpfile = "tdump.{}".format(iii)
                # run_suffix = self.filelocator.get_control_suffix(stage)
                cproc = [
                    os.path.join(self.inp["HYSPLIT_DIR"], "exec", "hyts_std"),
                    run_suffix,
                ]
                # logger.info('Running {} with job id {}'.format("hyts_std", cproc[1]))
                processhandler.startnew(cproc, self.inp["WORK_DIR"], descrip=run_suffix)
                # wait 5 seconds between run starts.
                time.sleep(5)
            else:
                print("tdump file already exists {}".format(tdumpfile))
            num_procs = processhandler.checkprocs()
            while num_procs >= max_procs:
                print("in loop {}".format(num_procs))
                time.sleep(5)

        # wait for runs to finish
        done = False
        seconds_to_wait = 30
        total_time = 0
        # 60 minutes.
        max_time = 60 * 60
        while not done:
            num_procs = processhandler.checkprocs()
            if num_procs == 0:
                done = True
            time.sleep(seconds_to_wait)
            total_time += seconds_to_wait
            if total_time > max_time:
                processhandler.checkprocs()
                processhandler.killall()
                logger.warning("HYSPLIT run Timed out")
                done = True

    def create_plots(self, redraw=False, stage=0):
        pass

    def after_run_check(self, update=False):
        # Check for the tdump/cdump file
        rval = True
        # fn = self.filelocator.get_tdump_filename(stage=0)
        fn = "tdump.0"
        logger.debug("Looking for tdump file " + fn)
        if not os.path.exists(fn):
            rval = False
            if update:
                logger.error(
                    "******************************************************************************"
                )
                logger.error(
                    "The model has crashed. Check the HYSPLIT Message file for further information."
                )
                logger.error(
                    "******************************************************************************"
                )
        return rval

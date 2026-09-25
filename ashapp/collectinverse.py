import os
import datetime


import logging
from ashapp.ashruninterface import ModelCollectionInterface
from utilhysplit.metfiles import gefs_suffix_list
from utilhysplit.runhandler import ProcessList
import ashapp.utildatainsertion as udi
from utilvolc import runhelper

from utilvolc.invhelper import inverse_get_suffix_list

logger = logging.getLogger(__name__)


class CollectInverse(ModelCollectionInterface):
    """
    Runs unit source runs for volcanic ash inversion
    """
    ilist = [
        "meteorologicalData",
        "forecastDirectory",
        "archivesDirectory",
        "WORK_DIR",
        "HYSPLIT_DIR",
        "jobname",
        "durationOfSimulation",
        "latitude",
        "longitude",
        "bottom",
        "top",
        "emissionHours",
        "timeres",
        "rate",
        "area",
        "start_date",
        "samplingIntervalHours",
        "jobid",
    ]

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
       
        # inversion has a time resolution 
        if 'timeres' not in self._inp.keys():
            logger.info('setting default time resolution {}h'.format(1))
            self._inp['timeres']=1  

        # inversion also has a vertical resolution
        # this is used as input into the 
        if 'inv_vertical_resolution' not in self._inp.keys():
            logger.info('setting default time resolution {}m'.format(1000))
            self._inp['inv_vertical_resolution']=1000  

        logger.info('setting rate to {}'.format(1))
        self._inp['rate'] = 1

        #inverse_hash = inverse_get_suffix_list(inp)
        complete = runhelper.is_input_complete(self.ilist,self._inp)

        if "jobid" in self._inp.keys():
            self.JOBID = self._inp["jobid"]
        if complete:
            logger.info("Input contains all fields")

    # TO DO make member generator it's own property?:q
    # @property(self):
    # def membergenerator(self):

    @property
    def status(self):
        return self._status

    @status.setter
    def status(self, status):
        self._status = status

    def setup(self, overwrite):
        from ashapp.rundispersion import RunDispersion

        inp = self.inp.copy()

        edir_alt = os.path.join(inp["WORK_DIR"], inp["VolcanoName"], "inverse/")
        if os.path.isdir(edir_alt):
            inp["WORK_DIR"] = edir_alt

        #edate = inp["start_date"] + datetime.timedelta(hours=inp["emissionHours"])
        #drange = [inp["start_date"], edate]
        command_list = []

        suffix_hash = inverse_get_suffix_list(self.inp,suffix_type='date')
        for key in suffix_hash.keys(): print(suffix_hash[key])
        #import sys
        #sys.exit()

        for iii, suffix in enumerate(suffix_hash.keys()):
            inp["jobid"] = "{}_{}".format(self.JOBID, suffix)
            inp.update(suffix_hash[suffix])
            inp['source description'] = suffix

            inp['emissionHours'] = self.inp['timeres']
            inp['start_date']  = suffix_hash[suffix]['sdate']

            run = RunDispersion(inp)
            command = run.run_model(overwrite=False)
            self._filehash[suffix] = run.filehash
            logger.info("ADDING {}".format(run.filelist))
            self._filelist.extend(run.filelist)

            self._status[suffix] = run.status
            if "FAILED" in run.status[0] or "COMPLETE" in run.status[0]:
                logger.warning(run.status)
                continue
            if command:
                command_list.append(command)
            del run
        return command_list

    def run(self, overwrite=False):
        import time

        command_list = self.setup(overwrite=overwrite)
        processhandler = ProcessList()
        processhandler.pipe_stdout()
        processhandler.pipe_stderr()
        for iii, command in enumerate(command_list):
            logger.info("Runnning {} with job id{}".format("hycs_std", command[1]))
            processhandler.startnew(command, self.inp["WORK_DIR"], descrip=str(iii))
            # wait 5 seconds to avoid runs trying to access ASCDATA.CFG files at the
            # same time.
            time.sleep(5)
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




import os
import datetime


import logging
from ashapp.ashruninterface import ModelCollectionInterface
from utilhysplit.metfiles import gefs_suffix_list
from utilhysplit.runhandler import ProcessList
import ashapp.utildatainsertion as udi

logger = logging.getLogger(__name__)


class CollectEmitTimes(ModelCollectionInterface):
    """
    Runs regular dispersion runs with GEFS.
    TODO - will this also handle other types of runs with GEFS?
           This could be done just by replacing the RunDispersion class
           with a different  ModelRunInterface class.
    """

    def __init__(self, inp, jobid):
        self.JOBID = jobid

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
            "bottom",
            "top",
            "emissionHours",
            "rate",
            "area",
            "start_date",
            "samplingIntervalHours",
            "jobid",
        ]

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
        complete = True
        for iii in self._ilist:
            if iii not in self._inp.keys():
                logger.warning("Input does not contain {}".format(iii))
                complete = False
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
        from ashapp.runemittimes import RunEmitTimes

        inp = self.inp.copy()

        edir_alt = os.path.join(inp["WORK_DIR"], inp["VolcanoName"], "emitimes/")
        if os.path.isdir(edir_alt):
            inp["WORK_DIR"] = edir_alt

        edate = inp["start_date"] + datetime.timedelta(hours=inp["emissionHours"])
        drange = [inp["start_date"], edate]
        command_list = []

        emitlist = udi.find_emit_file(inp["WORK_DIR"], drange)
        # emitlist is np.ndarray and using not to test a full np.ndarray
        # gives an error that truth value of array with more than one element is ambiguous.
        # however can test a regular list like this.
        if not list(emitlist):
            logger.warning("No EMITTIMES files found in {}".format(inp["WORK_DIR"]))
            self.status = "FAILED no emittimes files found to create runs from"

        for iii, emitfile in enumerate(emitlist):
            suffix = emitfile.split("/")[-1]
            suffix = suffix.replace("EMIT_", "")
            suffix = suffix.replace("EMIT", "")
            inp["jobid"] = "{}_{}".format(self.JOBID, suffix)
            inp["emitfile"] = emitfile
            run = RunEmitTimes(inp.copy())
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
        command_list = self.setup(overwrite)
        processhandler = ProcessList()
        processhandler.pipe_stdout()
        processhandler.pipe_stderr()
        # suffix = gefs_suffix_list()
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


class GEFSEmitTimes(CollectEmitTimes):

    def setup(self, overwrite):
        from ashapp.runemittimes import RunEmitTimes

        inp = self.inp.copy()

        edir_alt = os.path.join(inp["WORK_DIR"], inp["VolcanoName"], "emitimes/")
        if os.path.isdir(edir_alt):
            inp["WORK_DIR"] = edir_alt

        edate = inp["start_date"] + datetime.timedelta(hours=inp["emissionHours"])
        drange = [inp["start_date"], edate]
        command_list = []

        emitlist = udi.find_emit_file(inp["WORK_DIR"], drange)
        # emitlist is np.ndarray and using not to test a full np.ndarray
        # gives an error that truth value of array with more than one element is ambiguous.
        # however can test a regular list like this.
        if  not list(emitlist):
            logger.warning("No EMITTIMES files found in {}".format(inp["WORK_DIR"]))
            self.status = "FAILED no emittimes files found to create runs from"

        for metsuffix in gefs_suffix_list():

            for iii, emitfile in enumerate(emitlist):
                suffix = emitfile.split("/")[-1]
                suffix = suffix.replace("EMIT_", "")
                suffix = suffix.replace("EMIT", "")
                inp["jobid"] = "{}_{}_{}".format(self.JOBID, suffix, metsuffix)
                inp["emitfile"] = emitfile
                run = RunEmitTimes(inp)
                run.metfilefinder.set_ens_member(metsuffix)
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

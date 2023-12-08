import datetime
import logging
import os

import ashapp.utildatainsertion as udi
from ashapp.ashruninterface import ModelCollectionInterface
from ashapp.runemittimes import RunEmitTimes
from utilhysplit.metfiles import gefs_suffix_list
from utilhysplit.runhandler import ProcessList
from utilvolc.runhelper import is_input_complete

logger = logging.getLogger(__name__)

# 2023 Dec 07 (amc) added emit_file_finder attribute to class.


class CollectEmitTimes(ModelCollectionInterface):
    """
    Creates runs from Emittimes files.
    """

    # list of required (req) and optional (opt) inputs
    ilist = [
        ("WORK_DIR", "req"),
        ("jobname", "req"),
        ("HoursToEnd", "req"),
        ("start_date", "req"),
        ("jobid", "req"),
        ("VolcanoName", "req"),
    ]
    ilist.extend(RunEmitTimes.ilist)
    # these inputs will be created.
    ilist.remove(("emitfile", "req"))

    def __init__(self, inp, jobid):
        self.JOBID = jobid

        self._inp = {}
        self.inp = inp

        self._filehash = {}
        self._filelist = []

        self._emit_file_finder = None

        self._status = {"MAIN": "INITIALIZED"}

    @property
    def emit_file_finder(self):
        return self._emit_file_finder

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
        inp = self._inp

        vals = [x[0] for x in RunEmitTimes.ilist]
        newinp = dict((k, inp[k]) for k in vals if k in inp.keys())

        edir_alt = os.path.join(inp["WORK_DIR"], inp["VolcanoName"], "emitimes/")
        if os.path.isdir(edir_alt):
            inp["WORK_DIR"] = edir_alt

        edate = inp["start_date"] + datetime.timedelta(hours=inp["HoursToEnd"])
        drange = [inp["start_date"], edate]
        command_list = []

        self._emit_file_finder.wdir = inp["WORK_DIR"]
        emitlist = self._emit_file_finder.find(drange)
        # emitlist is np.ndarray and using not to test a full np.ndarray
        # gives an error that truth value of array with more than one element is ambiguous.
        # however can test a regular list like this.
        if not list(emitlist):
            logger.warning("No EMITTIMES files found in {}".format(inp["WORK_DIR"]))
            self.status = "FAILED no emittimes files found to create runs from"

        for emitfile in emitlist:
            suffix = emitfile.split("/")[-1]
            suffix = suffix.replace("EMIT_", "")
            suffix = suffix.replace("EMIT", "")
            if str(self.JOBID) not in suffix:
                newinp["jobid"] = "{}_{}".format(self.JOBID, suffix)
            newinp["emitfile"] = emitfile
            run = RunEmitTimes(newinp)
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
            logger.info("Running {} with job id {}".format("hycs_std", command[1]))
            processhandler.startnew(command, self.inp["WORK_DIR"], descrip=str(iii))
            # wait 5 seconds to avoid runs trying to access ASCDATA.CFG files at the
            # same time.
            time.sleep(5)
            num_proces = processhandler.checkprocs()
            total_time = 0
            logger.info("Number of processes running {}".format(num_proces))
            while num_proces >= 10:
                num_proces = processhandler.checkprocs()
                max_time = 10 * 60
                seconds_to_wait = 10
                total_time += seconds_to_wait
                time.sleep(seconds_to_wait)
                if total_time > max_time:
                    logger.info("max time reached")
                    break
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
        inp = self.inp.copy()

        edir_alt = os.path.join(inp["WORK_DIR"], inp["VolcanoName"], "emitimes/")
        if os.path.isdir(edir_alt):
            inp["WORK_DIR"] = edir_alt

        edate = inp["start_date"] + datetime.timedelta(hours=inp["HoursToEnd"])
        drange = [inp["start_date"], edate]
        command_list = []

        self._emit_file_finder.wdir = inp["WORK_DIR"]
        emitlist = self._emit_file_finder.find(drange)
        #emitlist = udi.find_emit_file(inp["WORK_DIR"], drange)
        # emitlist is np.ndarray and using not to test a full np.ndarray
        # gives an error that truth value of array with more than one element is ambiguous.
        # however can test a regular list like this.
        if not list(emitlist):
            logger.warning("No EMITTIMES files found in {}".format(inp["WORK_DIR"]))
            self.status = "FAILED no emittimes files found to create runs from"

        for metsuffix in gefs_suffix_list():
            for emitfile in emitlist:
                suffix = emitfile.split("/")[-1]
                suffix = suffix.replace("EMIT_", "")
                suffix = suffix.replace("EMIT", "")
                if str(self.JOBID) not in suffix:
                    #newinp["jobid"] = "{}_{}".format(self.JOBID, suffix)
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

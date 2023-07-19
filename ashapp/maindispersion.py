# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# maind_dispersion.py - run HYSPLIT model
#
# 14 JUL 2023 (AMC) updated ilist class attributes.
#
# -----------------------------------------------------------------------------

import datetime
import logging
import os
import shutil

import requests

from ashapp import utils
from ashapp.ashruninterface import MainRunInterface
from ashapp.collectemittimes import CollectEmitTimes, GEFSEmitTimes
from ashapp.collectinverse import CollectInverse
from ashapp.ensembledispersion import EnsembleDispersion
from ashapp.graphicsdispersion import GraphicsDispersion
from ashapp.outputdispersion import OutputDispersion
from ashapp.rundispersion import RunDispersion
from ashapp.runtrajectory import RunTrajectory
from ashapp.outputtrajectory import OutputTrajectory
from ashapp.graphicstrajectory import GraphicsTrajectory
from utilvolc.runhelper import complicated2str, is_input_complete


# from ashapp import  utils

# from utilvolc.volcMER import HT2unit

logger = logging.getLogger(__name__)




class MainDispersion(MainRunInterface):

    ilist = []
    ilist.extend(RunDispersion.ilist)
    ilist.extend(OutputDispersion.ilist)
    # these are set in the main routines.
    ilist.remove(("Use_Mastin_eq", "req"))
    ilist.remove(("fraction_of_fine_ash", "req"))

    def __init__(self, inp, JOBID):
        # 14 instance attributes
        self.JOBID = JOBID  # string

        inp["jobid"] = JOBID
        self._inp = {}
        self.inp = inp  # dictionary from JobSetUP
        self.apistr = None
        self.urlstr = None
        self.headerstr = None

        self.filelocator = None
        #self.maptexthash = {}
        self.awips = True

        self._modelrun = RunDispersion(inp)
        inp["Use_Mastin_eq"] = True
        inp["fraction_of_fine_ash"] = 0.01
        self._modeloutput = OutputDispersion(inp, [])
        self._modelgraphics = GraphicsDispersion(inp)

        utils.setup_logger()

 

    @property
    def JOBID(self):
        return self._JOBID

    @JOBID.setter
    def JOBID(self, JOBID):
        self._JOBID = str(JOBID)

    @property
    def inp(self):
        return self._inp

    @inp.setter
    def inp(self, inp):
        self._inp.update(inp)
        complete = is_input_complete(self.ilist, self._inp)
        if not complete:
            logger.warning("Inputs not complete")
        else:
            logger.debug("Inputs complete")

    @property
    def modelrun(self):
        return self._modelrun

    @modelrun.setter
    def modelrun(self, mrun):
        self._modelrun = mrun

    @property
    def modeloutput(self):
        return self._modeloutput

    @modeloutput.setter
    def modeloutput(self, mout):
        self._modeloutput = mout

    @property
    def modelgraphics(self):
        return self._modelgraphics

    def add_api_info(self, apistr, urlstr, headerstr):
        self.apistr = apistr
        self.urlstr = urlstr
        self.headerstr = headerstr

    def update_run_status(self, jobId, status):
        if self.apistr:
            API_KEY = os.environ[self.apistr]
            RUN_URL = os.environ[self.urlstr]
            statusUrl = "{}/status/{}/{}".format(RUN_URL, jobId, status)
            req = requests.put(statusUrl, headers={self.headerstr: API_KEY})
            logger.debug('Requests put {}'.format(req))
        else:
            logger.info("Running in offline test mode")
        logger.info("Posted status {} for job {}".format(status, jobId))

    def handle_crash(self, stage=0):
        self.update_run_status(self.JOBID, "CRASHED")
        logger.info("The model has crashed for job {}.".format(self.JOBID))
        logger.info(datetime.datetime.now())

    def debug_message(self):
        # debug messages
        logger.debug("HYSPLIT_DIR     = {}".format(self.inp["HYSPLIT_DIR"]))
        logger.debug("MAP_DIR         = {}".format(self.inp["MAP_DIR"]))
        logger.debug("WORK_DIR        = {}".format(self.inp["WORK_DIR"]))
        logger.debug("CONVERT_EXE     = {}".format(self.inp["CONVERT_EXE"]))
        logger.debug("GHOSTSCRIPT_EXE = {}".format(self.inp["GHOSTSCRIPT_EXE"]))
        logger.debug("PYTHON_EXE      = {}".format(self.inp["PYTHON_EXE"]))

    def doit(self):
        """
        Main work flow.
        """
        self.debug_message()
        os.chdir(self.inp["WORK_DIR"])
        if not os.path.exists("ASCDATA.CFG"):
            shutil.copyfile(
                self.inp["HYSPLIT_DIR"] + "/bdyfiles/ASCDATA.CFG", "ASCDATA.CFG"
            )
        logger.info("Please wait for further information....")
        logger.info("Model submitted on {}".format(datetime.datetime.now()))

        # run_model will check if the run has already been done.
        # self.modelrun.run_model(overwrite=False)
        self.modelrun.run(overwrite=False)

        logger.info("STATUS {}".format(complicated2str(self.modelrun.status)))

        if not self.modelrun.filelist:
            logger.warning("No model files produced. exiting")
            self.update_run_status(self.JOBID, "FAILED")
            return

        # make the model output.
        self.modeloutput.inputlist = self.modelrun.filelist
        self.modeloutput.postprocess()

        # make the graphics
        if self.modeloutput.check():
            self.modelgraphics.inputlist = self.modeloutput.outputlist
            self.modelgraphics.postprocess()

        # update the run status
        self.update_run_status(self.JOBID, "COMPLETED")

        # cleanup files
        self.cleanup()

    def cleanup(self):
        return True

    def after_run_check(self, stage=0, update=False):
        return True
        # Check for the tdump/cdump file
        # should this be done here?
        # rval = True
        # fnlist = [x for x in self.modelrun.filelist if 'cdump' in x]
        # for fn in fnlist:
        #    if not os.path.exists(fn):
        #        rval = False
        #        logger.info("NOT found cdump file " + fn)
        #    else:
        #        logger.info("found cdump file " + fn[0])

        # if update and not rval:
        #    logger.error(
        #        "******************************************************************************"
        #    )
        #    logger.error(
        #        "The model has crashed. Check the HYSPLIT Message file for further information."
        #    )
        #    logger.error(
        #        "******************************************************************************"
        #    )
        #    self.handle_crash(stage=0)
        # return rval


class MainEmitTimes(MainDispersion):

    ilist = []
    ilist.extend(CollectEmitTimes.ilist)
    ilist.extend(OutputDispersion.ilist)

    def __init__(self, inp, JOBID):

        """
        modelrun attribute is the EnsembleDispersion class.
        """

        self.JOBID = JOBID  # string

        inp["jobid"] = JOBID
        self._inp = {}
        self.inp = inp  # dictionary from JobSetUP
        self.apistr = None
        self.urlstr = None
        self.headerstr = None

        self.filelocator = None
        #self.maptexthash = {}
        self.awips = True

        if inp["meteorologicalData"].lower() == "gefs":
            self._modelrun = GEFSEmitTimes(inp, self.JOBID)
        else:
            self._modelrun = CollectEmitTimes(inp, self.JOBID)

        inp["Use_Mastin_eq"] = False
        inp["fraction_of_fine_ash"] = 1
        self._modeloutput = OutputDispersion(inp, [])
        self._modelgraphics = GraphicsDispersion(inp)

        utils.setup_logger()


class MainInverse(MainDispersion):

    ilist = []
    ilist.extend(CollectInverse.ilist)
    ilist.extend(OutputDispersion.ilist)


    def __init__(self, inp, JOBID):
        """
        modelrun attribute is the CollectInverse class.
        Run unit source runs for inversion
        """

        self.JOBID = JOBID  # string

        inp["jobid"] = JOBID
        self._inp = {}
        self.inp = inp  # dictionary from JobSetUP
        self.apistr = None
        self.urlstr = None
        self.headerstr = None

        self.filelocator = None
        #self.maptexthash = {}
        self.awips = True

        self._modelrun = CollectInverse(inp, self.JOBID)
        inp["Use_Mastin_eq"] = False
        inp["fraction_of_fine_ash"] = 1
        self._modeloutput = OutputDispersion(inp, [])
        self._modelgraphics = GraphicsDispersion(inp)

        utils.setup_logger()


class MainGEFSInverse(MainInverse):

    # same as MainInverse but over-rides the doit method.

    # needs a setter since reset the model output for each GEFS run.
    # @modeloutput.setter
    # def modeloutput(self, mout):
    #    self._modeloutput = mout

    def doit(self):
        from utilhysplit.metfiles import gefs_suffix_list

        """
        The main workflow is implemented for each GEFS members
        """
        self.debug_message()
        os.chdir(self.inp["WORK_DIR"])
        if not os.path.exists("ASCDATA.CFG"):
            shutil.copyfile(
                self.inp["HYSPLIT_DIR"] + "/bdyfiles/ASCDATA.CFG", "ASCDATA.CFG"
            )
        logger.info("Please wait for further information....")
        logger.info("Model submitted on {}".format(datetime.datetime.now()))

        # run_model will check if the run has already been done.
        # self.modelrun.run_model(overwrite=False)

        # create a separate netcdf file for each member?
        inp = self.inp.copy()

        for metsuffix in gefs_suffix_list():
            inp["jobid"] = "{}_{}".format(self.JOBID, metsuffix)
            inp["meteorologicalData"] = "gefs{}".format(metsuffix)
            self.modelrun = CollectInverse(inp, inp["jobid"])
            self.modelrun.run(overwrite=False)

            logger.info(self.modelrun.status)

            if not self.modelrun.filelist:
                logger.warning("No model files produced. exiting")
                self.update_run_status(self.JOBID, "FAILED")
                return

            inp["Use_Mastin_eq"] = False
            inp["fraction_of_fine_ash"] = 1
            self.modeloutput = OutputDispersion(inp, [])
            # make the model output.
            self.modeloutput.inputlist = self.modelrun.filelist
            self.modeloutput.postprocess()
            inp = inp.remove("Use_Mastin_eq")
            inp = inp.remove("fraction_of_fine_ash")

        # make the graphics
        # if self.modeloutput.check():
        #    self.modelgraphics.inputlist = self.modeloutput.outputlist
        #    self.modelgraphics.postprocess()

        # update the run status
        self.update_run_status(self.JOBID, "COMPLETED")

        # cleanup files
        self.cleanup()


class MainEnsemble(MainDispersion):
 
    ilist = []
    ilist.extend(EnsembleDispersion.ilist)
    ilist.extend(OutputDispersion.ilist)


    def __init__(self, inp, JOBID):
        """
        modelrun attribute is the EnsembleDispersion class.
        """

        self.JOBID = JOBID  # string


        inp["jobid"] = JOBID
        self._inp = {}
        self.inp = inp  # dictionary from JobSetUP
        self.apistr = None
        self.urlstr = None
        self.headerstr = None

        self.filelocator = None
        #self.maptexthash = {}
        #self.awips = True

        inp["Use_Mastin_eq"] = True
        inp["fraction_of_fine_ash"] = 0.05
        self._modelrun = EnsembleDispersion(inp, self.JOBID)
        self._modeloutput = OutputDispersion(inp, [])
        self._modelgraphics = GraphicsDispersion(inp)

        utils.setup_logger()


class MainTrajectory(MainDispersion):

    ilist = []
    ilist.extend(RunTrajectory.ilist)
    ilist.extend(OutputTrajectory.ilist)
    # these are set in the main routines.

    def __init__(self, inp, JOBID):
        # 14 instance attributes
        self.JOBID = JOBID  # string

        inp["jobid"] = JOBID
        self._inp = {}
        self.inp = inp  # dictionary from JobSetUP
        self.apistr = None
        self.urlstr = None
        self.headerstr = None

        self.filelocator = None
        #self.maptexthash = {}

        self._modelrun = RunTrajectory(inp)
        self._modeloutput = OutputTrajectory(inp, [])
        self._modelgraphics = GraphicsTrajectory(inp)


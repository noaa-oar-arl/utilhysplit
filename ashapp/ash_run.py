#!/opt/Tools/anaconda3/envs/hysplit/bin/python
# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# ash_run.py - run HYSPLIT model on web (or offline) and create plots
#
# 01 JUN 2020 (AMC) - adapted from locusts-run.py
# 15 NOV 2022 (AMC) - updated how tests can be run
# -----------------------------------------------------------------------------
# To run in offline mode standard dispersion run  python ash_run.py -999
# To run in offline mode standard trajectory run  python ash_run.py -777
# To run in offline mode ensemble dispersion run  python ash_run.py -888
#
# TODO To run in offline mode ensemble trajectory run  python ash_run.py -555
#
#
# -----------------------------------------------------------------------------

#import json
import logging
import os
import sys
import traceback
#from abc import ABC, abstractmethod

import requests

from utilhysplit.metfiles import gefs_suffix_list
from utilvolc.runhelper import JobSetUp, make_inputs_from_file
from utils import setup_logger

logger = logging.getLogger(__name__)


def print_usage():
    print(
        """\
USAGE: ash_run.py RUNTYPE JOBID

where RUNTYPE is either 'run' or 'redraw' or 'test' (without quotations) and JOBID
is the job identification number.

If not using runtype of 'test' then
The following environment variables must be set prior to calling this script:
    VOLCANICASH_API_KEY         - secret key to access the hysplitash web APIs.
    VOLCANICASH_URL             - URL to the hysplitash web application."""
    )


# Base class is AshRun
# EnsembleAshRun inherits from AshRun
# TrajectoryAshRun inherits from AshRun


def create_run_instance(JOBID, inp):
    """
    create the correct object for the type of run and return it.
    """

    if inp["runflag"] == "dispersion":
        if inp["meteorologicalData"].lower() == "gefs":
            from ashensemble import EnsembleAshRun
            arun = EnsembleAshRun(JOBID)
        else:
            from ashbase import AshRun
            arun = AshRun(JOBID)

    elif inp["runflag"] == "inverse":
        from ashinverse import InverseAshRun
        arun = InverseAshRun(JOBID)

    elif inp["runflag"] == "DataInsertion":
        from ashdatainsertion import DataInsertionRun
        arun = DataInsertionRun(JOBID)

    else: # Trajectory run
        from ashtrajectory import TrajectoryAshRun
        arun = TrajectoryAshRun(JOBID)

    arun.add_inputs(inp)
    return arun

def choosetest(JOBID,setup):
    if JOBID == "-446":
        # test inverse dispersion run.
        # currently
        setup.add_inverse_params()
        inp = setup.inp
        inp["runflag"] = "inverse"
        # inp['meteorologicalData'] = 'GFS'
        arun = create_run_instance("446", inp)
    elif JOBID == "-555":
        # test ensemble trajectory run.
        # currently
        inp["runflag"] = "trajectory"
        inp["meteorologicalData"] = "GEFS"
        arun = create_run_instance("1004", inp)
    elif JOBID == "-888":
        # test ensemble dispersion run.
        inp["meteorologicalData"] = "GEFS"
        arun = create_run_instance("1001", inp)
    elif JOBID == "-777":
        # test trajectory  run.
        inp["runflag"] = "trajectory"
        arun = create_run_instance("1003", inp)
    else:
        # test regular dispersion run.
        inp["runflag"] = "dispersion"
        arun = create_run_instance("1005", inp)
    return arun

if __name__ == "__main__":
    # Configure the logger so that log messages appears in the "Model Status" text box.
    setup_logger(level=logging.DEBUG)

    if len(sys.argv) != 3:
        print_usage()
        sys.exit(1)

    RUNTYPE = sys.argv[1]  # 'run' or 'redraw'
    JOBID = sys.argv[2]  # This is a REDRAWID for redraw.
    setup = JobSetUp()

    # create a test to run without inputs from web page.
    if RUNTYPE == 'test':
        logging.getLogger().setLevel(20)
        logging.basicConfig(stream=sys.stdout)
        configname = "config.{}.txt".format(JOBID)
        logger.info("CONFIG FILE {}".format(configname))
        logger.info("TESTING")
        setup = make_inputs_from_file("./", configname)
        inp = setup.inp
        if inp['runflag'] == 'inverse':
           setup.add_inverse_params()
           logger.info("Inverse run")
        inp = setup.inp
        if inp["meteorologicalData"].lower() == "gefs":
            for suffix in gefs_suffix_list():
                tdir = "/hysplit-users/alicec/utilhysplit/utilvolc/ashapp/"
                setup = make_inputs_from_file(tdir, configname)
                setup.add_inverse_params()
                inp = setup.inp
                logger.info("GEFS WORKING on {}".format(suffix))
                JOBIDens = "{}_{}".format(JOBID, suffix)
                inp["gefsmember"] = suffix
                arun = create_run_instance(JOBIDens, inp)
                arun.doit()
        else:
            arun = create_run_instance(JOBID, inp)
            arun.doit()
        sys.exit(1)

    # create test based on JOBID
    elif JOBID in ["-999", "-888", "-777", "-555", "-446"]:
        # config file should be called config.jobid.txt
        logging.getLogger().setLevel(20)
        logging.basicConfig(stream=sys.stdout)
        # inp = setup.make_test_inputs()
        configname = "config.{}.txt".format(JOBID)
        logger.info("CONFIG FILE {}".format(configname))
        setup = make_inputs_from_file("./", configname)
        inp = setup.inp
        # inp = setup.make_test_inputs()
        arun=choosetest(JOBID,setup,inp)
        arun.doit()
        sys.exit(1)
    
    # actual run from web api inputs
    else:
        # arun = AshRun(JOBID)
        apistr = "VOLCANICASH_API_KEY"
        urlstr = "VOLCANICASH_URL"
        headerstr = "VolcanicAsh-API-Key"
        try:
            API_KEY = os.environ[apistr]
            RUN_URL = os.environ[urlstr]
            if RUNTYPE == "redraw":
                # Update JOBID after obtaining the redraw input.
                REDRAWID = sys.argv[2]
                inputUrl = "{}/redrawinput/{}".format(RUN_URL, REDRAWID)
                r = requests.get(inputUrl, headers={headerstr: API_KEY})
                a = r.json()
                JOBID = a["id"]
            else:
                inputUrl = "{}/runinput/{}".format(RUN_URL, JOBID)
                r = requests.get(inputUrl, headers={headerstr: API_KEY})
                a = r.json()
            # print(json.dumps(a, indent=4))
            inp = setup.parse_inputs(a)
            arun = create_run_instance(JOBID, inp)
            arun.add_api_info(apistr, urlstr, headerstr)
            arun.doit()
        except Exception:
            # TODO:
            # Send the error message to the Application layer so that the user could see it.
            logger.error(
                "Oops! Something went wrong. This is an internal error and we will need to fix it."
            )
            logger.error(traceback.format_exc())
            arun.handle_crash()
    sys.exit(0)

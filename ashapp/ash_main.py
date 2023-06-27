#!/opt/Tools/anaconda3/envs/hysplit/bin/python
# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# ash_run.py - run HYSPLIT model on web (or offline) and create plots
#
# 01 JUN 2020 (AMC) - adapted from locusts-run.py
# 15 NOV 2022 (AMC) - updated how tests can be run
# 23 May 2023 (AMC) - made data insertion runs more flexible
# -----------------------------------------------------------------------------
# To run in offline mode standard dispersion run  python ash_run.py -999
# To run in offline mode standard trajectory run  python ash_run.py -777
# To run in offline mode ensemble dispersion run  python ash_run.py -888
#
# TODO To run in offline mode ensemble trajectory run  python ash_run.py -555
#
#
# -----------------------------------------------------------------------------

# import json
import logging
import os
import sys
import traceback

# from abc import ABC, abstractmethod

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
    VOLCANICASH_URL             - URL to the hysplitash web application.

If using 'test' then a configuration file is read. 
The call should be of the form
python ash_run.py test JOBID

The program will look for a file with the name config.JOBID.txt and read it for inputs.
Examples are provided.
In the configuration file, the runtype is specified. 
Further documentation for some run types is below.

------------------------------------------------------------------------------------
DATAINSERTION
Assumes that emit-times file generated from VOLCAT data (or other data) has been
generated.

Looks for a directory comprised of the following from the configuration file
/wdir/volcanoname/emitimes/

In this directory, look for emit-times files. Naming convention for emitimes files should
be EMITIMES_suffix, or EMIT_suffix.

If naming convention is according to volcat then will also use the dates in the 
configuration file to only create runs for those filenames with dates between start_date
and start_date +  emissionHours.

If different naming convention then will simply create runs for all EMIT files in the directory.


"""
    )


# Base class is AshRun
# EnsembleAshRun inherits from AshRun
# TrajectoryAshRun inherits from AshRun


def create_run_instance(JOBID, inp):
    """
    create the correct object for the type of run and return it.
    The system currently supports the following types of runs.

    dispersion deterministic : creates forward dispersion runs
    dispersion gefs : creates forward dispersion runs for each gefs member
    inverse : creates a set of unit source forward runs for creating TCM
    inverse  gefs: creates a set of unit source forward runs for creating TCM for each gefs member.
    DataInsertion : uses previously created  emit-times files to create runs
    BackTrajectory : runs back trajectories using a csv file which has information about starting points.
    trajectory : deterministic : runs forward trajectories at predetermined set of heights
    trajectory : gefs : runs forward trajectories at predetermined set of heights for each GEFS ensemble member.

    """
    logger.info("Creating run Instance")

    if inp["runflag"] == "dispersion":
        if inp["meteorologicalData"].lower() == "gefs":
            from maindispersion import MainEnsemble

            arun = MainEnsemble(inp, JOBID)
            logger.info("Dispersion GEFS")
        else:
            from maindispersion import MainDispersion

            arun = MainDispersion(inp, JOBID)
            logger.info("Dispersion")

    elif inp["runflag"].lower() == "datainsertion":
        # This handles GEFS as well as deterministic runs.
        from maindispersion import MainEmitTimes

        arun = MainEmitTimes(inp, JOBID)
        logger.info("Use EmitTimes files")

    elif inp["runflag"] == "inverse":
        if inp["meteorologicalData"].lower() == "gefs":
            # this one generates a separate netcdf file
            # for each gefs member
            from maindispersion import MainGEFSInverse
            arun = MainGEFSInverse(inp, JOBID)
        else:
            from maindispersion import MainInverse

            arun = MainInverse(inp, JOBID)
        logger.info("Inversion")

    # elif inp["runflag"] == "DataInsertion":
    #    from ashdatainsertion import DataInsertionRun
    #    arun = DataInsertionRun(JOBID)
    #    logger.info('Data Insertion')

    # elif inp["runflag"] == "BackTrajectory":
    #    from backtraj import BackTraj
    #    arun = BackTraj(JOBID)
    #    logger.info('Back Trajectory')

    # elif inp["runflag"] == "trajectory":
    #    if inp["meteorologicalData"].lower() == "gefs":
    #        from enstrajectory import EnsTrajectoryRun
    #        arun = EnsTrajectoryRun(JOBID)
    #        logger.info('Ensemble Trajectory')
    #    else:
    #        from ashtrajectory import TrajectoryAshRun
    #        arun = TrajectoryAshRun(JOBID)
    #        logger.info('Trajectory')

    # arun.inp = inp
    return arun


if __name__ == "__main__":
    # Configure the logger so that log messages appears in the "Model Status" text box.
    # setup_logger(level=logging.DEBUG)
    setup_logger()

    if len(sys.argv) != 3:
        print_usage()
        sys.exit(1)

    RUNTYPE = sys.argv[1]  # 'run' or 'redraw'
    JOBID = sys.argv[2]  # This is a REDRAWID for redraw.
    setup = JobSetUp()

    # create a test to run without inputs from web page.
    if RUNTYPE == "test":
        logging.getLogger().setLevel(20)
        logging.basicConfig(stream=sys.stdout)
        configname = "config.{}.txt".format(JOBID)
        logger.info("CONFIG FILE {}".format(configname))
        logger.info("TESTING")
        setup = make_inputs_from_file("./", configname)
        inp = setup.inp
        if inp["runflag"] == "inverse":
            setup.add_inverse_params()
            logger.info("Inverse run")
        inp = setup.inp

        arun = create_run_instance(JOBID, inp)
        arun.doit()
        sys.exit(1)

    # actual run from web api inputs
    else:
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
            elif RUNTYEP == "dispersion":
                inputUrl = "{}/runinput/{}".format(RUN_URL, JOBID)
                r = requests.get(inputUrl, headers={headerstr: API_KEY})
                a = r.json()
                a = r.json()
            elif RUNTYEP == "datainsertion":
                inputUrl = "{}/datainsertion/{}".format(RUN_URL, JOBID)
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

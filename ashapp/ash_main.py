#!/opt/Tools/anaconda3/envs/hysplit/bin/python
# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# ash_main.py - run HYSPLIT model on web (or offline) and create plots
#
# 01 JUN 2020 (AMC) - adapted from locusts-run.py
# 15 NOV 2022 (AMC) - updated how tests can be run
# 23 May 2023 (AMC) - made data insertion runs more flexible
# 29 JUN 2023 (AMC) - modified ash_run into ash_main.
#                     refactored code for more object comoposition and less inheritance.
# -----------------------------------------------------------------------------
#
# EMAIL 6/2/2023 from Sonny
# -s good to talk with you yesterday about the changes for the volcanic ash web app.
# After the meeting, it occurred to me that we would need an additional web endpoint
# for a python script to call to create/update a volcano event for the data insertion
# feature. The endpoint will need to store volcano events to a database table which
# will provide information needed for data insertion later. Let's circle back in a
# few weeks
# ----------------------------------------------------------------------------

# import json
import logging
import os
import sys
import traceback

import requests

from utils import setup_logger
from utilvolc.runhelper import JobSetUp, make_inputs_from_file

# from abc import ABC, abstractmethod

#pylint: disable-msg=C0103

logger = logging.getLogger(__name__)


def print_usage():
    print(
        """\
USAGE: ash_main.py RUNTYPE JOBID

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
Runs from EmitTimes files

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

def create_run_instance(jid, runinp):
    """
    create the correct object for the type of run and return it.
    The system currently supports the following types of runs.

    dispersion: creates forward dispersion runs
    inverse : creates a set of unit source forward runs for creating TCM
    DataInsertion : uses previously created  emit-times files to create runs

    BackTrajectory : runs back trajectories using a csv file which has information
                     about starting points.
    trajectory : deterministic : runs forward trajectories at predetermined
                 set of heights
    trajectory : gefs : runs forward trajectories at predetermined set of heights
                 for each GEFS ensemble member.
    """
    logger.info("Creating run Instance")

    if runinp["runflag"] == "dispersion":
        if runinp["meteorologicalData"].lower() == "gefs":
            from maindispersion import MainEnsemble
            crun = MainEnsemble(runinp, jid)
            logger.info("Dispersion GEFS")
        else:
            from maindispersion import MainDispersion

            crun = MainDispersion(runinp, jid)
            logger.info("Dispersion")

    elif runinp["runflag"].lower() == "datainsertion":
        # This handles GEFS as well as deterministic runs.
        from maindispersion import MainEmitTimes

        crun = MainEmitTimes(runinp, jid)
        logger.info("Use EmitTimes files")

    elif runinp["runflag"] == "inverse":
        if runinp["meteorologicalData"].lower() == "gefs":
            # this one generates a separate netcdf file
            # for each gefs member
            from maindispersion import MainGEFSInverse

            crun = MainGEFSInverse(runinp, jid)
        else:
            from maindispersion import MainInverse

            crun = MainInverse(runinp, jid)
        logger.info("Inversion")

    # elif runinp["runflag"] == "DataInsertion":
    #    from ashdatainsertion import DataInsertionRun
    #    crun = DataInsertionRun(jid)
    #    logger.info('Data Insertion')

    # elif runinp["runflag"] == "BackTrajectory":
    #    from backtraj import BackTraj
    #    crun = BackTraj(jid)
    #    logger.info('Back Trajectory')

    # elif runinp["runflag"] == "trajectory":
    #    if runinp["meteorologicalData"].lower() == "gefs":
    #        from enstrajectory import EnsTrajectoryRun
    #        crun = EnsTrajectoryRun(jid)
    #        logger.info('Ensemble Trajectory')
    #    else:
    #        from ashtrajectory import TrajectoryAshRun
    #        crun = TrajectoryAshRun(jid)
    #        logger.info('Trajectory')

    # crun.inp = runinp
    return crun


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
        finputs = make_inputs_from_file("./", configname)
        inp = finputs.inp
        if inp["runflag"] == "inverse":
            finputs.add_inverse_params()
            logger.info("Inverse run")
        inp = finputs.inp

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
                rrr = requests.get(inputUrl, headers={headerstr: API_KEY})
                aaa = rrr.json()
                JOBID = aaa["id"]
            elif RUNTYPE == "dispersion":
                inputUrl = "{}/runinput/{}".format(RUN_URL, JOBID)
                rrr = requests.get(inputUrl, headers={headerstr: API_KEY})
                aaa = rrr.json()
            elif RUNTYPE == "datainsertion":
                inputUrl = "{}/datainsertion/{}".format(RUN_URL, JOBID)
                rrr = requests.get(inputUrl, headers={headerstr: API_KEY})
                aaa = rrr.json()
            # print(json.dumps(a, indent=4))
            inp = setup.parse_inputs(aaa)
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

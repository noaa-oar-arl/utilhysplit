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
#import traceback
from ashapp.utils import setup_logger
#from abc import ABC, abstractmethod
from utilvolc import volcat
import monet

logger = logging.getLogger(__name__)


def print_usage():
    print(
        """\
USAGE: process_volcat.py inputdir inpfilename outdir gridspace

inputdir is the directory where volcat file resides
inpfilename is the name of the volcat file
outdir is the directory where the parallax and regridded file should be written
gridspace is the horizontal resolution of regridded data in degrees.
    """
    )


if __name__=="__main__":
    setup_logger(level=logging.DEBUG)
    logger.info("number of input arguments {}".format(len(sys.argv)))
    if len(sys.argv)!= 5:
       print_usage()
       sys.exit(1)

    data_dir = sys.argv[1] # directory for volcat file
    filename = sys.argv[2]  # filename of volcat data to process
    pc_dir   = sys.argv[3]  #
    gridspace = sys.argv[4] # grid spacing to regrid to

    print(data_dir)
    print(filename)
    print(pc_dir)
    print(gridspace)
    logger.info("Working on it")
    volcat.write_parallax_corrected_files(data_dir,pc_dir,flist=[filename],gridspace=float(gridspace),verbose=True)
    
    sys.exit(0)

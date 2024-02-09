# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# fitgmm.py - run HYSPLIT model on web (or offline) and create plots
#
#                     refactored code for more object comoposition and less inheritance.
# -----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------

import logging
import os
import sys
#import traceback
from ashapp.utils import setup_logger

import numpy as np
from monetio.models import pardump
from utilml import par2conc

# from abc import ABC, abstractmethod

# pylint: disable-msg=C0103

logger = logging.getLogger(__name__)


def print_usage():
    print(
        """\
USAGE: fitgmm.py 


"""
    )


if __name__ == "__main__":
    # Configure the logger so that log messages appears in the "Model Status" text box.
    # setup_logger(level=logging.DEBUG)
    setup_logger()

    #if len(sys.argv) != 3:
    #    print_usage()
    #    sys.exit(1)

    print(sys.argv)

    pname = sys.argv[1]  # pardump name
    nnn = int(sys.argv[2])  # number of gaussians to fit
    oname = sys.argv[3]  # outputfilename

    if os.path.isfile(pname):
        pardf = pardump.open_dataset(pname)
    else: 
        logger.warning('file not found {}'.format(pname))
        sys.exit()

    dlist = pardf.date.unique()
    print('The following dates were found in the file')
    print(dlist)

    mass = pardf.pmass.unique()
    if len(mass) > 1:
       logger.warning('particles have different masses')
       print(mass)

    hts = pardf.ht.unique()
    maxht = np.max(hts)

    mlist = par2conc.fit_timeloop(pardf,nnn=nnn,maxht=maxht,method='gmm')

    for mfit in mlist:
        mfit.save()
 
    sys.exit(0)

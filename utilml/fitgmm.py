# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# fitgmm.py - Fits Gaussian Mixture Models to HYSPLIT particle data
# -----------------------------------------------------------------------------
# This script processes HYSPLIT particle dump files and applies Gaussian Mixture
# Model (GMM) fitting to the particle distributions. The fitted models are then
# saved to output files.
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
USAGE: fitgmm.py PARDUMP_FILE NUM_GAUSSIANS OUTPUT_NAME

Arguments:
  PARDUMP_FILE   : Path to the HYSPLIT particle dump file
  NUM_GAUSSIANS  : Number of Gaussians to fit in the mixture model
  OUTPUT_NAME    : Base name for the output file(s)

Example:
  python fitgmm.py pardump.txt 5 output_model
"""
    )


if __name__ == "__main__":
    # Configure the logger so that log messages appear in the "Model Status" text box.
    # setup_logger(level=logging.DEBUG)
    setup_logger()

    # Check for correct number of arguments
    if len(sys.argv) != 4:
        print_usage()
        sys.exit(1)

    print(sys.argv)

    pname = sys.argv[1]  # pardump name
    nnn = int(sys.argv[2])  # number of gaussians to fit
    oname = sys.argv[3]  # outputfilename

    if os.path.isfile(pname):
        pardf = pardump.open_dataset(pname)
    else: 
        logger.warning('File not found: {}'.format(pname))
        sys.exit(1)

    dlist = pardf.date.unique()
    print('The following dates were found in the file:')
    print(dlist)

    mass = pardf.pmass.unique()
    if len(mass) > 1:
       logger.warning('Particles have different masses:')
       print(mass)

    hts = pardf.ht.unique()
    maxht = np.max(hts)

    # Fit Gaussian Mixture Models (GMM) to particle data for each time step
    logger.info(f"Fitting {nnn} Gaussians to particle distribution with max height {maxht}")
    mlist = par2conc.fit_timeloop(pardf, nnn=nnn, maxht=maxht, method='gmm')

    # Save all fitted models to files using the provided output name as base
    logger.info(f"Saving fitted models with base name: {oname}")
    for i, mfit in enumerate(mlist):
        # Pass the output name to the save method or modify the filename
        mfit.save(filename=f"{oname}_{i}")
 
    logger.info("Processing complete")
    sys.exit(0)

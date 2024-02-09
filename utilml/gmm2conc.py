# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# gmm2conc.py - get netcdf with concentrations from fits previously made to pardump outputs.
# -----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------

import logging
import os
import sys
#import traceback
from ashapp.utils import setup_logger
import xarray as xr
import glob
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

    pname = sys.argv[1]  # partial name
    dd = float(sys.argv[2])  # horizontal resolution
    dz = float(sys.argv[3])  # vertical resolution
    oname = sys.argv[4]

    possible = glob.glob('*_covariances.npy')
    print(possible)
    print(pname)
    matches = [x for x in possible if pname in x]
    print('Matches found ', matches)

    dralist = []
    for match in matches:
        base =  match.replace('_covariances.npy','')
        print('Working on {}'.format(base))
        mfit = par2conc.load_massfit(base) 
        dra = mfit.get_conc(dd,dz)
        dra = par2conc.shift_underground(dra)
        #dra = dset[list(dset.data_vars.keys())[0]] 
        dra = dra.to_dataset(name='CONC')
        dralist.append(dra)
 
    dall = xr.merge(dralist,join='outer')
    dall.to_netcdf(oname ) 
    sys.exit(0)

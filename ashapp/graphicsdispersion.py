#!/opt/Tools/anaconda3/envs/hysplit/bin/python
# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# graphicsdispersion.py - 
#
# 28 FEB 2023 (AMC) - check attributes added as stand alone function but still kept as static method.
#
# -----------------------------------------------------------------------------
# To run in offline mode use python ash_run.py -999
#
#
# -----------------------------------------------------------------------------


# from abc mport ABC, abstractmethod
import datetime
import glob
import logging
import os
import shutil
import subprocess
import sys
import zipfile

import numpy as np
import requests
import xarray as xr

from utilvolc.runhelper import ConcplotColors, Helper, JobFileNameComposer
from ashapp.ashruninterface import ModelOutputInterface


logger = logging.getLogger(__name__)


class GraphicsDispersion(ModelOutputInterface):

    def __init__(self,inp):
        self._inputlist = []
        self._outputlist = []
        self.inp = inp

    @property
    def inputlist(self):
        return self._inputlist

    @inputlist.setter
    def inputlist(self, inlist):
        self._inputlist = inlist 

    # no setter
    @property
    def outputlist(self):
        return self._outputlist

    def check(self):
       rval = True
       for filename in self._outputlist:
           if not os.path.isfile(filename):
              return False
       return True

    def postprocess(self):
        return True



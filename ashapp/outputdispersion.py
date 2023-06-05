#!/opt/Tools/anaconda3/envs/hysplit/bin/python
# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# outputdispersion.py - run HYSPLIT model on web and create plots
#
# 03 JUN 2023 (AMC) - 
#
# -----------------------------------------------------------------------------
# To run in offline mode use python ash_run.py -999
#
#
# -----------------------------------------------------------------------------


# from abc mport ABC, abstractmethod
import datetime
import logging
import os
import sys

import numpy as np
import xarray as xr

from monetio.models import hysplit
from utilvolc.runhelper import  Helper
from ashapp.ashruninterface import ModelOutputInterface
from ashapp import  utils
from ashapp.ashnetcdf import HYSPLITAshNetcdf


# from runhandler import ProcessList
from utilvolc.volcMER import HT2unit
logger = logging.getLogger(__name__)
utils.setup_logger()

class OutputDispersion(ModelOutputInterface):
    from ashapp.ashnetcdf import HYSPLITAshNetcdf

    def __init__(self,inp,filelist): 

        self.JOBID = inp['jobid']
        self.inp = inp
        self._inputlist = filelist
        self._outputlist = self.set_outputlist()

        self._ncfile = HYSPLITAshNetcdf(self.cdumpxraname())

    @property
    def inputlist(self):
        return self._inputlist

    @inputlist.setter
    def inputlist(self, inlist):
        self._inputlist = inlist 

    # no setter. read only
    @property
    def outputlist(self):
        return self._outputlist

    def cdumpxraname(self): 
       return "xrfile.{}.nc".format(self.JOBID)

    def set_outputlist(self):
        outputlist = []
        outputlist.append(self.cdumpxraname())
        return outputlist

    def postprocess(self):
       self.get_cdump_xra()      
       self._ncfile.write_with_compression(overwrite=True)
       # combine all cdump files into an xarray and write a netcdf file.
       return -1 

    def check(self):
       rval = True
       for filename in self._outputlist:
           if not os.path.isfile(filename):
              return False
       return True

## additional methods ----------------------------------------------------
    def get_cdump_xra(self):
        """
        reads from file if it exists.
        Otherwise create from list of cdump names.
        updates the multiplication factor if necessary.
        """
        if  self._ncfile.empty():
            blist = []
            cdumpname = [x for x in self.inputlist if 'cdump' in x]
            cdumpname = cdumpname[0]
            logger.info('Creating xra for cdump files {}'.format(cdumpname))
            source_tag = "Line to {:1.0f} km".format(self.inp["top"] / 1000.0)
            met_tag = self.inp["meteorologicalData"]
            blist = [(cdumpname, source_tag, met_tag)]
            print(blist)
            century = 100 * (int(self.inp["start_date"].year / 100))
            species=None
            ainp = self.inp
            ainp['mult'] = 1
            self._ncfile.make_cdump_xra(blist,century,species=species, inp=ainp)
        mult = self.get_conc_multiplier()
        change = self._ncfile.changemult(mult)
        return self._ncfile.cxra 

    def write_cdump_xra(self):
        # convert to mg /m3
        mult = self.get_conc_multiplier()
        change = self._ncfile.changemult(mult)
        if change:
            ainp = self.inp
            ainp['mult'] = mult
            Helper.remove(os.path.join(self.inp["WORK_DIR"], fname))
            self._ncfile.assign_attrs(ainp)
            self._ncfile.write_with_compression()

        if not self._ncfile.empty():
            logger.info(
                "make_awips_netcdf: cxra empty. cannot create awips\
                         files"
            )
            return []


    def make_awips_netcdf(self):
        import cdump2netcdf

        ghash = {}
        ghash["source_latitude"] = self.inp["latitude"]
        ghash["source_longitude"] = self.inp["longitude"]
        ghash["source_name"] = self.inp["VolcanoName"]
        ghash["emission_start"] = self.inp["start_date"]
        ghash["emission_duration_hours"] = self.inp["emissionHours"]
        # in mg
        mult = self.get_conc_multiplier()
        mer = mult / 1e6 / 3600  # kg released in a second.
        ghash["MER"] = mer
        ghash["MER_unit"] = "kg/s"
        logger.debug("MULT value for awips {:3e}".format(mult))
        awipsname = self.filelocator.get_awips_filename(stage=0)
        c2n = cdump2netcdf.Cdump2Awips(
            self._cxra, awipsname, munit="mg", jobid=self.JOBID, globalhash=ghash
        )
        awips_files = c2n.create_all_files()
        # returns list of awips files that were created.
        return awips_files
 
    def get_conc_multiplier(self):
        """
        factor to convert concentration in unit mass /m3 to mg/m3.
        Model must have been run with 1 unit mass/h released.
        """
        height = (self.inp["top"] - self.inp["bottom"]) / 1000.0
        # HT2unit outputs in grams. Then multiply by 1e3 to get mg
        conc_multiplier = 1e3 * HT2unit(height) * self.get_ash_reduction()
        return conc_multiplier

    def get_ash_reduction(self):
        eflag = float(self.inp["eflag"])
        M63 = 0.01  # fraction of  fine ash
        conc_multiplier = 10 ** (-1 * eflag)
        return M63 * conc_multiplier


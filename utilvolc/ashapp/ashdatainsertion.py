# -----------------------------------------------------------------------------

#from abc import ABC, abstractmethod
import datetime
import logging
import os
import time

# from hysplitdata.traj import model
# from hysplitplot import timezone

# import hysplit
import metfiles as metfile
from monetio.models import hysplit
from runhelper import Helper
from runhandler import ProcessList
from ashbase import AshRun
import ensemble_tools
from cdump2xml import HysplitKml
from emitimes import EmitTimes

logger = logging.getLogger(__name__)


def print_usage():
    print(
        """\
USAGE: ashensemble.py JOBID

The following environment variables must be set prior to calling this script:
    RUN_API_KEY         - secret key to access Locusts web APIs.
    RUN_URL             - URL to the Locusts web application."""
    )

class DataInsertionAshRun(AshRun):

    def __init__(self, JOBID):
        super().__init__(JOBID)

    def get_maptext_info(self):
        maptexthash = {}
        rstr = "HYSPLIT Data Insertion."
        maptexthash["run_description"] = rstr
        maptexthash["infoc"] = ""
        return maptexthash

    def read_emittimes(self, emitfile):
        """
        get start date and number of locations from emit file.
        """
        self.file_not_found_error(emitfile, message=True)
        etf = emitimes.EmitTimes(filename=emitfile)
        etf.read_file(num_species=1)
        # look at first emission cycle 
        ecycle = etf.cycle_list[0]
        # number of locations that need to be in CONTROL file.
        nlocs = ecycle.nlocs 
        # starting date of this cycle      
        sdate = ecycle.sdate
        self.inp['start_date'] = sdate

        # duration of this cycle
        #cduration = ecycle.duration

        # get first line locations
        erecord = ecycle.recordra[0]
        self.inp['latitude'] = erecord.lat
        self.inp['longitude'] = erecord.lon

        # set to 0 since these will be from emit-times
        self.inp['rate'] = 0
        self.inp['area'] = 0
        self.inp['bottom'] = 0
        self.inp['top'] = 0
        return nlocs

    def setup_setup(self,stage):
        super().setup_setup(stage)
        setup.add("efile",self.inp['emitfile'])

    def additional_control_setup(self, control, stage=0):
        emitfile = self.inp['emitfile']
        # get number of locations from emit-times file
        # and set some values in inp needed for the control setup.
        nlocs = self.read_emittimes(emitfile)
        super().additional_control_setup()
        # add as many dummy locations as in emit-times file
        control.remove_locations()
        for loc in nlocs:
            control.add_location((lat,lon),vent,rate=rate,area=area) 

    def file_not_found_error(self, fln, message=False):
        if not os.path.exists(fln):
            if message:
                logger.error(
                    "******************************************************************************"
                )
                logger.error(
                    "This file was not found {}\
                              for job {}.".format(
                        fln, self.JOBID
                    )
                )
                logger.error(
                    "******************************************************************************"
                )
            rval = False
        else:
            logger.debug("file found {}".format(fln))
            rval= True
        return rval



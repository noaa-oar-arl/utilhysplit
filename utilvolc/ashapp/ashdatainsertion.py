# -----------------------------------------------------------------------------

#from abc import ABC, abstractmethod
import datetime
import numpy as np
import logging
import os
import time
from glob import glob
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
from emitimes import EmiTimes


logger = logging.getLogger(__name__)


def print_usage():
    print(
        """\
USAGE: ashensemble.py JOBID

The following environment variables must be set prior to calling this script:
    RUN_API_KEY         - secret key to access Locusts web APIs.
    RUN_URL             - URL to the Locusts web application."""
    )


def find_emit_file(wdir):
    elist = glob(os.path.join(wdir,'EMIT_*'))
    return elist[0]

class DataInsertionRun(AshRun):

    def __init__(self, JOBID):
        super().__init__(JOBID)


    def add_inputs(self, inp):
        inp['emissionHours']=0
        inp['samplingIntervalHours']=1
        inp['WORK_DIR'] = os.path.join(inp['WORK_DIR'],
                                            inp['VolcanoName'],
                                            'emitimes/')
        logger.info('Working directory set {}'.format(inp["WORK_DIR"]))
        super().add_inputs(inp)

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
        etf = EmiTimes(filename=emitfile)
        self.inp['emitfile'] = emitfile
        etf.read_file(num_species=1)
        # look at first emission cycle 
        ecycle = etf.cycle_list[0]
        #print('ecycle', ecycle)
        # number of locations that need to be in CONTROL file.
        nlocs = ecycle.nrecs
        logger.info('{} number of locations'.format(nlocs))
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
        self.inp['nlocs'] = nlocs
        return nlocs

    def setup_setup(self,stage):
        setup = super().setup_setup(stage)
        eloc = self.inp['emitfile'].split('/')
        eloc = eloc[-1]
        setup.add("efile",eloc)
        return setup

    def setup_basic_control(self,stage=0,rtype='dispersion'):
        emitfile = find_emit_file(self.inp['WORK_DIR'])
        nlocs = self.read_emittimes(emitfile)
        control = super().setup_basic_control(stage=stage,rtype=rtype)
        return control

    def additional_control_setup(self, control, stage=0):
        #emitfile = find_emit_file(self.inp['WORK_DIR'])
        # get number of locations from emit-times file
        # and set some values in inp needed for the control setup.
        #nlocs = self.read_emittimes(emitfile)
        nlocs = self.inp['nlocs']
        super().additional_control_setup(control,stage=stage)
        # add as many dummy locations as in emit-times file
        control.remove_locations()
        lat = self.inp['latitude']
        lon = self.inp['longitude']
        vent = self.inp['bottom']
        area = self.inp['area']
        rate = self.inp['rate']
        for loc in np.arange(0,nlocs):
            control.add_location((lat,lon),vent,rate=rate,area=area) 
        #control.concgrids[0].outfile = self.filelocator.get_cdump_filename(stage)
        #control.concgrids[0].outdir = self.inp['WORK_DIR']

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



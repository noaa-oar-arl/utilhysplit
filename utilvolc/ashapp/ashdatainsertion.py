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
from utilvolc import write_emitimes as vwe
from utilvolc.ashapp.runhelper import AshDINameComposer

logger = logging.getLogger(__name__)


def print_usage():
    print(
        """\
USAGE: ashensemble.py JOBID

The following environment variables must be set prior to calling this script:
    RUN_API_KEY         - secret key to access Locusts web APIs.
    RUN_URL             - URL to the Locusts web application."""
    )


def find_emit_file(wdir,daterange):
    print('wdir',wdir)
    edf = vwe.get_emit_name_df(wdir)
    edf = edf[edf['algorithm name'] == 'EMIT']
    edf = edf[edf['idate'] >= daterange[0]]
    print(daterange[0], '----------------')
    print(edf['idate'])
   
    edf = edf[edf['idate'] <= daterange[1]]
    print(daterange[1], '----------------')
    print(edf['idate'])
    #elist = glob(os.path.join(wdir,'EMIT_*'))
    elist = edf['filename'].values
    print('elist', elist)
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
        self.filelocator = AshDINameComposer(self.inp['WORK_DIR'],
                                 self.JOBID, 
                                 self.inp['jobname'])
    

    def get_maptext_info(self):
        maptexthash = {}
        rstr = "HYSPLIT Data Insertion."
        maptexthash["run_description"] = rstr
        maptexthash["infoc"] = ""
        return maptexthash

    def read_emittimes(self, emitfile):
        """
        get information from emit-times file including
        start date, number of locations, latitude, longitude
        
        """
        self.file_not_found_error(emitfile, message=True)
        etf = EmiTimes(filename=emitfile)
        self.inp['emitfile'] = emitfile
        etf.read_file(num_species=1)
        # look at first emission cycle 
        ecycle = etf.cycle_list[0]
        #print('ecycle', ecycle)
        # number of locations that need to be in CONTROL file.
        self.inp['nlocs'] = ecycle.nrecs
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

    def setup_setup(self,stage):
        setup = super().setup_setup(stage=self.emitfile)
        # add the emit times file
        eloc = self.inp['emitfile'].split('/')
        eloc = eloc[-1]
        setup.add("efile",eloc)
        return setup

    def setup_basic_control(self,stage=0,rtype='dispersion'):
        edate = self.inp['start_date'] + datetime.timedelta(hours=self.inp['durationOfSimulation'])
        emitfile = find_emit_file(self.inp['WORK_DIR'], [self.inp['start_date'],edate])
        self.emitfile = emitfile
        self.read_emittimes(emitfile)
        control = super().setup_basic_control(stage=emitfile,rtype=rtype)
        return control

    def additional_control_setup(self, control, stage=0):
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

    def run_model(self):
        # make control and setup files
        self.compose_control(stage=0, rtype="dispersion")
        self.compose_setup(stage=0)
        run_suffix = self.filelocator.get_control_suffix(stage=self.emitfile)
        # start run and wait for it to finish..
        c = [os.path.join(self.inp["HYSPLIT_DIR"], "exec", "hycs_std"), str(run_suffix)]
        logger.info("Running {} with job id {}".format("hycs_std", c[1]))
        Helper.execute(c)


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



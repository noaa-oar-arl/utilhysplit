# -----------------------------------------------------------------------------

# from abc import ABC, abstractmethod
import datetime
import logging
import os
import time

import numpy as np

import xarray as xr
# import hysplit
from monetio.models import hysplit
from utilvolc import write_emitimes as vwe
from utilvolc.ashapp.ashbase import AshRun
from utilvolc.ashapp.emitimes import EmiTimes
from utilvolc.ashapp.runhandler import ProcessList
from utilvolc.ashapp.runhelper import AshDINameComposer

# from hysplitdata.traj import model
# from hysplitplot import timezone


logger = logging.getLogger(__name__)


def print_usage():
    print(
        """\
USAGE: ashensemble.py JOBID

The following environment variables must be set prior to calling this script:
    RUN_API_KEY         - secret key to access Locusts web APIs.
    RUN_URL             - URL to the Locusts web application."""
    )


def find_emit_file(wdir, daterange, retype="fname"):
    return find_di_file(wdir, daterange, "EMIT")


def find_cdump_df(wdir, daterange):
    return find_di_file(wdir, daterange, "cdump", rtype="dataframe")


def find_di_file(wdir, daterange, ftype, rtype="fname"):
    edf = vwe.get_emit_name_df(wdir)
    edf = edf[edf["algorithm name"] == ftype]
    edf = edf[edf["idate"] >= daterange[0]]
    edf = edf[edf["idate"] <= daterange[1]]
    # elist = glob(os.path.join(wdir,'EMIT_*'))
    elist = edf["filename"].values
    if rtype == "fname":
        rval =  elist
    elif rtype == "dataframe":
        rval = edf
    return rval

class DataInsertionRun(AshRun):
    #def __init__(self, JOBID):
    #    super().__init__(JOBID)

    def add_inputs(self, inp):
        inp["emissionHours"] = 0
        inp["samplingIntervalHours"] = 1
        inp["WORK_DIR"] = os.path.join(inp["WORK_DIR"], inp["VolcanoName"], "emitimes/")
        logger.info("Working directory set {}".format(inp["WORK_DIR"]))
        super().add_inputs(inp)
        self.filelocator = AshDINameComposer(
            self.inp["WORK_DIR"], self.JOBID, self.inp["jobname"]
        )

    def get_maptext_info(self):
        maptexthash = {}
        rstr = "HYSPLIT Data Insertion."
        maptexthash["run_description"] = rstr
        maptexthash["infoc"] = ""
        return maptexthash

    def get_cdump_xra(self):
        edate = self.inp["start_date"] + datetime.timedelta(
            hours=self.inp["durationOfSimulation"]
        )
        cdf = find_cdump_df(self.inp["WORK_DIR"], [self.inp["start_date"], edate])
        blist = list(cdf["filename"].values)
        alist = []
        if not blist:
            logger.warning("No cdump files found")
            return xr.Dataset()
        for fname in blist:
            alist.append((fname, fname, self.inp["meteorologicalData"]))
        century = 100 * (int(self.inp["start_date"].year / 100))
        cdumpxra = hysplit.combine_dataset(
            alist, century=century, sample_time_stamp="start"
        )
        cdumpxra.attrs["Volcano ID"] = cdf["volcano id"].unique()[0]
        cdumpxra.attrs["Volcano Name"] = self.inp["VolcanoName"]
        cdumpxra.attrs["layer"] = cdf["layer"].unique()[0]
        cdumpxra.attrs["mult"] = 1
        return cdumpxra

    def read_emittimes(self, emitfile):
        """
        get information from emit-times file including
        start date, number of locations, latitude, longitude

        set rate, area, top and bottom to 0 since these values
        will be from the emit-times file.
        """
        self.file_not_found_error(emitfile, message=True)
        etf = EmiTimes(filename=emitfile)
        self.inp["emitfile"] = emitfile
        etf.read_file(num_species=1)
        # look at first emission cycle
        ecycle = etf.cycle_list[0]
        # print('ecycle', ecycle)
        # number of locations that need to be in CONTROL file.
        self.inp["nlocs"] = ecycle.nrecs

        # starting date of this cycle
        sdate = ecycle.sdate
        self.inp["start_date"] = sdate
        # duration of this cycle
        # cduration = ecycle.duration

        # get first line locations
        # erecord = ecycle.recordra[0]
        # self.inp['latitude'] = erecord.lat
        # self.inp['longitude'] = erecord.lon

        # set to 0 since these will be from emit-times
        self.inp["rate"] = 0
        self.inp["area"] = 0
        self.inp["bottom"] = 0
        self.inp["top"] = 0

    def setup_setup(self, stage):
        setup = super().setup_setup(stage=stage)
        # add the emit times file
        eloc = self.inp["emitfile"].split("/")
        eloc = eloc[-1]
        setup.add("efile", eloc)
        return setup

    def setup_basic_control(self, stage="emitfile", rtype="dispersion"):
        self.read_emittimes(stage)
        control = super().setup_basic_control(stage=stage, rtype=rtype)
        return control

    def additional_control_setup(self, control, stage=0):
        nlocs = self.inp["nlocs"]
        super().additional_control_setup(control, stage=stage)
        # add as many dummy locations as in emit-times file
        control.remove_locations()
        lat = np.floor(self.inp["latitude"])
        lon = np.floor(self.inp["longitude"])
        vent = self.inp["bottom"]
        area = self.inp["area"]
        rate = self.inp["rate"]
        for loc in np.arange(0, nlocs):
            control.add_location((lat, lon), vent, rate=rate, area=area)

    def get_conc_multiplier(self):
        return 1

    def after_run_check(self, update):
        edate = self.inp["start_date"] + datetime.timedelta(
            hours=self.inp["durationOfSimulation"]
        )
        rlist = []
        rval = True
        for emitfile in find_emit_file(
            self.inp["WORK_DIR"], [self.inp["start_date"], edate]
        ):
            rval = super().after_run_check(stage=emitfile, update=update)
            rlist.append(rval)
        return np.all(rlist)

    def run_model(self):
        edate = self.inp["start_date"] + datetime.timedelta(
            hours=self.inp["durationOfSimulation"]
        )
        processhandler = ProcessList()
        processhandler.pipe_stdout()
        processhandler.pipe_stderr()
        for emitfile in find_emit_file(
            self.inp["WORK_DIR"], [self.inp["start_date"], edate]
        ):
            fn = self.filelocator.get_cdump_filename(stage=emitfile)
            if os.path.exists(fn):
                logger.info("cdump exists {} continuing to next run".format(fn))
                continue
            # make control and setup files
            self.compose_control(stage=emitfile, rtype="dispersion")
            self.compose_setup(stage=emitfile)
            run_suffix = self.filelocator.get_control_suffix(emitfile)
            # start run and wait for it to finish..
            cproc = [
                os.path.join(self.inp["HYSPLIT_DIR"], "exec", "hycs_std"),
                str(run_suffix),
            ]
            logger.info("Running {} with job id {}".format("hycs_std", cproc[1]))
            processhandler.startnew(cproc, self.inp["WORK_DIR"], descrip=run_suffix)
            # wait 5 seconds between run starts to avoid
            # runs trying to access ASCDATA.CFG at the same time.
            time.sleep(5)
        # wait for runs to finish
        done = False
        seconds_to_wait = 30
        total_time = 0
        # 60 minutes.
        max_time = 60 * 60
        # max_time = 0.5*60
        while not done:
            num_proces = processhandler.checkprocs()
            if num_proces == 0:
                done = True
            time.sleep(seconds_to_wait)
            total_time += seconds_to_wait
            if total_time > max_time:
                processhandler.checkprocs()
                processhandler.killall()
                logger.warning("HYSPLIT run Timed out")
            # Helper.execute(c)

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
            rval = True
        return rval

# -----------------------------------------------------------------------------

# from abc import ABC, abstractmethod
import datetime
import logging
import os
import time
import numpy as np
import pandas as pd
# import xarray as xr

# import hysplit
from monetio.models import hysplit
from utilhysplit.emitimes import EmiTimes
from utilhysplit.runhandler import ProcessList

from ashapp.ashbase import AshRun
from utilvolc import make_data_insertion as mdi
from utilvolc.runhelper import AshDINameComposer

# from hysplitdata.traj import model
# from hysplitplot import timezone


logger = logging.getLogger(__name__)


def print_usage():
    print(
        """\
Run through the ash_run.py executable.

Assumes that emit-times file generated from VOLCAT data (or other data) has been
generated.

Looks for a directory comprised of the following from the configuration file
/wdir/volcanoname/emitimes/

In this directory, look for emit-times files. Naming convention for emitimes files should
be EMITIMES_suffix, or EMIT_suffix.

If naming convention is according to volcat then will also use the dates in the 
configuration file to only create runs for those filenames with dates between start_date
and start_date +  emissionHours.

If different naming convention then will simply create runs for all EMIT files in the directory.

"""
    )

def find_emit_file_alt(wdir):
    """
    Find all files which start with EMIT in the directory
    """
    import glob
    efile = glob.glob(wdir + '/EMIT*') 
    efile = [x.split('/')[-1] for x in efile]
    return efile

def find_emit_file(wdir, daterange, rtype="fname"):
    # first look for files with volcat naming convention.
    elist =  find_di_file(wdir, daterange, "EMIT", rtype=rtype)
    if not(list(elist)):
       elist = find_emit_file_alt(wdir)      
    return elist

def find_cdump_df_alt(wdir, jobid):
    efile = find_emit_file_alt(wdir)
    filelocator = AshDINameComposer(
                  wdir,jobid,jobid)
    cnames = []

    for eff in efile:
        stage = '{}_{}'.format(eff,jobid)
        cfile = filelocator.get_cdump_filename(stage)
        print('find cdump', eff, cfile)
        cnames.append(cfile)
    df = pd.DataFrame(cnames,columns=['filename'])
    df['file descriptor'] = 'cdump'
    df["volcano id"] = 'unknown'
    df['layer'] = 'unknown'
    return df 

def find_cdump_df(wdir, jobid, daterange):
    ftype = 'cdump'
    # this finds all the cdump files.
    dset = find_di_file(wdir, daterange, ftype, rtype="dataframe")
    if isinstance(dset,list):
       dset = find_cdump_df_alt(wdir,jobid)
    # return only cdump files with the jobid in the name.
    return dset[dset.apply(lambda row: jobid in row['filename'],axis=1)]

def find_di_file(wdir, daterange, ftype, rtype="fname"):
    edf = mdi.get_emit_name_df(wdir,ftype)
    if 'file descriptor' not in edf.columns:
        return []
    edf = edf[edf["file descriptor"] == ftype]
    edf = edf[edf["idate"] >= daterange[0]]
    edf = edf[edf["idate"] <= daterange[1]]
    # elist = glob(os.path.join(wdir,'EMIT_*'))
    elist = edf["filename"].values
    if rtype == "fname":
        rval = elist
    elif rtype == "dataframe":
        rval = edf
    #print('FIND DI FILE', wdir, daterange, ftype, edf, rval)
    return rval


class DataInsertionRun(AshRun):
    """
    INPUTS
    DIrunHours is set to inp['emissionHours']
    The inp['emissionHours'] set to 0. 

    WORK_DIR is set to the input work directory + volcano name + emitimes
             This may need to be changed.

    """


    # def __init__(self, JOBID):
    #    super().__init__(JOBID)

    def add_inputs(self, inp):
        # this start date is the time that the data insertion runs begin.
        # later start_date will be used for the start time of each individual run.
        inp["DI_start_date"] = inp['start_date'] 
        inp["DIrunHours"] = inp["emissionHours"]
        inp["emissionHours"] = 0
        inp["samplingIntervalHours"] = 1
        inp["WORK_DIR"] = os.path.join(inp["WORK_DIR"], inp["VolcanoName"], "emitimes/")
        logger.info("Working directory set {}".format(inp["WORK_DIR"]))
        super().add_inputs(inp)
        if inp["meteorologicalData"].lower() == "gefs":
            logger.info("ens member {}".format(inp["gefsmember"]))
            self.metfilefinder.set_ens_member("." + inp["gefsmember"])
        self.awips = False
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
        import sys
        edate = self.inp["DI_start_date"] + datetime.timedelta(
            hours=self.inp["DIrunHours"]
        )
        logger.warning('DATES {} {} {}'.format(self.inp['DI_start_date'], edate, self.inp['DIrunHours']))
        cdf = find_cdump_df(self.inp["WORK_DIR"], self.JOBID, [self.inp["start_date"], edate])
        blist = list(cdf["filename"].values)
        alist = []
        # if not blist:
        #    logger.warning("No cdump files found")
        #    return xr.Dataset()
        for fname in blist:
            logger.info("Adding to netcdf file {}".format(fname))
            alist.append((fname, fname, self.inp["meteorologicalData"]))
        century = 100 * (int(self.inp["DI_start_date"].year / 100))
        cdumpxra = hysplit.combine_dataset(
            alist, century=century, sample_time_stamp="start"
        )
        cdumpxra.attrs["Volcano ID"] = cdf["volcano id"].unique()[0]
        cdumpxra.attrs["Volcano Name"] = self.inp["VolcanoName"]
        cdumpxra.attrs["layer"] = cdf["layer"].unique()[0]
        cdumpxra.attrs["mult"] = 1
        return cdumpxra

    def create_plots(self, redraw=False, stage=0):
        # TODO currently no plots automatically created.
        return True 

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
        self.read_emittimes(stage.replace('_' + self.JOBID,''))
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
        logger.warning("after run check")
        edate = self.inp["start_date"] + datetime.timedelta(
            hours=self.inp["DIrunHours"]
        )
        logger.warning(
            "after run check {} {} {}".format(
                self.inp["start_date"], edate, self.inp["DIrunHours"]
            )
        )
        rlist = []
        rval = True
        for emitfile in find_emit_file(
            self.inp["WORK_DIR"], [self.inp["start_date"], edate]
        ):
            logger.warning("Looking for {}".format(emitfile))
            rval = super().after_run_check(stage='{}_{}'.format(emitfile,self.JOBID), update=update)
            logger.warning("FOUND {}".format(rval))
            rlist.append(rval)
        return np.all(rlist)

    def run_model(self):
        edate = self.inp["start_date"] + datetime.timedelta(
            hours=self.inp["DIrunHours"]
        )
        processhandler = ProcessList()
        processhandler.pipe_stdout()
        processhandler.pipe_stderr()
        for emitfile in find_emit_file(
            self.inp["WORK_DIR"], [self.inp["start_date"], edate]
        ):
            stage = '{}_{}'.format(emitfile,self.JOBID)
            #stage = stage.replace('.','')
            fnn = self.filelocator.get_cdump_filename(stage=stage)
            logger.info("cdump filename {}".format(stage))
            logger.info("stage {}".format(stage))
            if os.path.exists(fnn):
                logger.info("cdump exists {} continuing to next run".format(fnn))
                continue
            # make control and setup files
            self.compose_control(stage=stage, rtype="dispersion")
            self.compose_setup(stage=stage)
            run_suffix = self.filelocator.get_control_suffix(stage)
            # start run and wait for it to finish..
            cproc = [
                os.path.join(self.inp["HYSPLIT_DIR"], "exec", "hycs_std"),
                str(run_suffix),
            ]
            #import sys
            logger.info("Running {} with job id {}".format("hycs_std", cproc[1]))
            #sys.exit()
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

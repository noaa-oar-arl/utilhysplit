# -----------------------------------------------------------------------------

# from abc import ABC, abstractmethod
import logging
import os
import pandas as pd

# import xarray as xr

# import hysplit
from utilhysplit.emitimes import EmiTimes

# from ashapp.ashbase import AshRun
from utilvolc import make_data_insertion as mdi
from utilvolc.runhelper import AshDINameComposer

# from hysplitdata.traj import model
# from hysplitplot import timezone


logger = logging.getLogger(__name__)


def find_emit_file_alt(wdir):
    """
    Find all files which start with EMIT in the directory
    """
    import glob

    efile = glob.glob(wdir + "/EMIT*")
    efile = [x.split("/")[-1] for x in efile]
    return efile


def find_emit_file(wdir, daterange, rtype="fname"):
    # first look for files with volcat naming convention.
    elist = find_di_file(wdir, daterange, "EMIT", rtype=rtype)
    if not (list(elist)):
        elist = find_emit_file_alt(wdir)
    return elist


def find_cdump_df_alt(wdir, jobid):
    efile = find_emit_file_alt(wdir)
    filelocator = AshDINameComposer(wdir, jobid, jobid)
    cnames = []

    for eff in efile:
        stage = "{}_{}".format(eff, jobid)
        cfile = filelocator.get_cdump_filename(stage)
        print("find cdump", eff, cfile)
        cnames.append(cfile)
    df = pd.DataFrame(cnames, columns=["filename"])
    df["file descriptor"] = "cdump"
    df["volcano id"] = "unknown"
    df["layer"] = "unknown"
    return df


def find_cdump_df(wdir, jobid, daterange):
    ftype = "cdump"
    # this finds all the cdump files.
    dset = find_di_file(wdir, daterange, ftype, rtype="dataframe")
    if isinstance(dset, list):
        dset = find_cdump_df_alt(wdir, jobid)
    # return only cdump files with the jobid in the name.
    d2 = dset[dset.apply(lambda row: jobid in row["filename"], axis=1)]
    return d2


def find_di_file(wdir, daterange, ftype, rtype="fname"):
    edf = mdi.get_emit_name_df(wdir, ftype)
    if "file descriptor" not in edf.columns:
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
    # print('FIND DI FILE', wdir, daterange, ftype, edf, rval)
    return rval


def add_inputs(inp):
    # this start date is the time that the data insertion runs begin.
    # later start_date will be used for the start time of each individual run.
    inp["DI_start_date"] = inp["start_date"]
    inp["DIrunHours"] = inp["emissionHours"]
    inp["emissionHours"] = 0
    inp["samplingIntervalHours"] = 1
    inp["WORK_DIR"] = os.path.join(inp["WORK_DIR"], inp["VolcanoName"], "emitimes/")
    logger.info("Working directory set {}".format(inp["WORK_DIR"]))
    return inp


def get_maptext_info():
    maptexthash = {}
    rstr = "HYSPLIT Data Insertion."
    maptexthash["run_description"] = rstr
    maptexthash["infoc"] = ""
    return maptexthash


def read_emittimes(emitfile):
    """
    get information from emit-times file including
    start date, number of locations, latitude, longitude
    set rate, area, top and bottom to 0 since these values
    will be from the emit-times file.
    """
    inp = {}
    # self.file_not_found_error(emitfile, message=True)
    etf = EmiTimes(filename=emitfile)
    inp["emitfile"] = emitfile
    etf.read_file(num_species=1)
    # look at first emission cycle
    ecycle = etf.cycle_list[0]
    # print('ecycle', ecycle)
    # number of locations that need to be in CONTROL file.
    inp["nlocs"] = ecycle.nrecs

    # starting date of this cycle
    sdate = ecycle.sdate
    inp["start_date"] = sdate
    # duration of this cycle
    # cduration = ecycle.duration

    # get first line locations
    # erecord = ecycle.recordra[0]
    # self.inp['latitude'] = erecord.lat
    # self.inp['longitude'] = erecord.lon

    # set to 0 since these will be from emit-times
    inp["rate"] = 0
    inp["area"] = 0
    inp["bottom"] = 0
    inp["top"] = 0
    return inp

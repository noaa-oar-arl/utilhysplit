import logging
import json
import pandas as pd
import numpy as np
import os
import datetime
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import time
import monet
import sys
import xarray as xr

from utilvolc import volcano_names
from utilvolc.volcano_names import fix_volc_name
from utilvolc.runhelper import Helper
from utilvolc.runhelper import list_dirs
from utilvolc.runhelper import make_dir

from utilvolc import volcat
from utilvolc.volcat import flist2eventdf

from utilhysplit.runhandler import ProcessList
from utilhysplit.plotutils import map_util
import utilhysplit.evaluation.web_ensemble_plots as wep
from utilvolc import make_data_insertion as mdi
from utilhysplit.evaluation import ensemble_tools

from utilvolc.volcat_files import EventFile
from utilvolc.volcat_files import get_summary_file_df
from utilvolc.volcat_files import get_log_files
from utilvolc.volcat_files import check_file

from utilvolc import volcat_plots as vp

from utilvolc.qvainterface import DFInterface
from utilvolc.qvainterface import EventInterface


logger = logging.getLogger(__name__)

"""
 Volcat Event file handling

 Changelog
 2023 Jul 26 AMC  created from material in qva_logic.
 2023 Aug 10 AMC  created OneEventTime to combine volcat files w
                  ith same time but different feature id.
 2023 Sep 14 AMC  started adding average method to Events class
                  started fleshing out EventStatus class.
"""

#TODO 
# integrate the OneEventTime class into the workflow
# complete the average method or other way of creating averaged volcat data.
# 


class EventDisplay:
    """
    Helper class for Events
    """

    def __init__(self, eventid="TEST", volcano_name="Unknown", tdir="./"):

        self.eventid = eventid
        self.volcano_name = fix_volc_name(volcano_name)
        self._vplot = vp.VolcatPlots(
            pd.DataFrame(), volcano_name=self.volcano_name, tdir=tdir
        )
        self.tdir = tdir

    @property
    def tdir(self, tdir):
        return self._tdir

    @tdir.setter
    def tdir(self, tdir):
        # update the directory for VolcatPlots
        # read any csv file that is there.
        self._vplot.tdir = tdir
        print("HERE", tdir)
        self._vplot.read()
        self._tdir = tdir

    def add(self, events):
        print("HERE add")
        self._vplot.make_arrays(events)
        print("HERE add")
        self._vplot.save()


class EventStatus:
    """
    Helper class for Events
    """

    def __init__(self, eventid="TEST", status=None):

        self.eventid = eventid
        self.status = status
      
        self.start_time=None   # start of events
        self.current_time=None # last time updated
        self.end_time=None     # time event was ended.

    @property
    def status(self):
        return self._status

    @status.setter
    def status(self,status):
        if not isinstance(status,str):
           try: status = str(status)
           except Exception as eee:
                logger.warning(eee)
        slist = ['Ended','Failed','Ongoing','Paused']
        if status in slist:
           self._status = status
        else:
           self._status = status          

    @property
    def start_time(self):
        return self._start_time

    @start_time.setter
    def start_time(self,time):
        if isinstance(time,datetime.datetime):
            self._start_time=time


class EventDF(DFInterface):
    # columns which are required.
    required = [
        "sensor_name",
        "feature_id",
        "event_file",
        "volcano_lat",
        "volcano_lon",
        "volcano_name",
        "observation_date",
    ]

    def __init__(self, edf, fid="TEST"):

        if self.are_required_columns_present(edf):
            self._edf = edf
        else:
            self._edf = pd.DataFrame()

        self.savename = "{}.csv".format(fid)

    @property
    def edf(self):
        return self._edf.copy()

    @property
    def required_columns(self):
        return self.required

    def save(self, overwrite=False):
        """
        save dataframe as csv file
        """
        rval = True
        if not os.path.isfile(self.savename) or overwrite:
            self.edf.to_csv(self.savename, float_format="%.1f", index=False)
        else:
            logger.warning("CSV file exists. cannot save")
            rval = False
        return rval

    # extra method
    def read(self, ndir, cname=None):
        """
        reads csv file that was previously saved
        """
        if not cname:
            cname = os.path.join(ndir, "Events.csv")
        else:
            cname = os.path.join(ndir, cname)
        dtp = {"observation_date": "datetime64[ns]"}
        df = pd.read_csv(cname, sep=",", parse_dates=["observation_date"])
        return df

    def add_csv(self, ndir, cname=None):
        """
        reads csv file that was previously saved and adds it to dataframe.
        """
        dftemp = self.read_csv(ndir, cname)
        self.add_df(dftemp)

    def are_required_columns_present(self, df):
        answer = True
        for req in self.required:
            if req not in df.columns:
                logger.info(
                    "EventDF WARNING, data frame does not contain required column {}".format(
                        req
                    )
                )
                answer = False
        return answer

    def add_df(self, df):
        """
        df : pandas DataFrame
        """
        # can utilize output of flist2eventdf
        # change all column names to lower case
        complete = True

        # check that required columns are in dataframe
        columns = df.columns
        newc = [x.lower() for x in columns]
        df.columns = newc
        complete = self.are_required_columns_present(df)

        if complete:
            if self._edf.empty:
                if isinstance(df, pd.core.frame.DataFrame):
                    self._edf = df
            else:
                if isinstance(df, pd.core.frame.DataFrame):
                    self._edf = pd.concat([self._edf, df])
            self._edf.drop_duplicates()
        return complete

    # extra method
    def add_events(self, eventlist):
        """
        eventlist: list of EventFile objects.
        """
        # elist - list of dataframes associated with the EventFile.
        elist = []
        for efile in eventlist:
            df = efile.df.copy()
            # for key in efile.attrs:
            #    df[key.lower()] = efile.attrs[key]
            if self.are_required_columns_present(df):
                elist.append(df)
        dfnew = pd.concat(elist)
        self.add_df(dfnew)


class OneEventTime:
    def __init__(self, das):
        """
        list of volcat files which are valid for the same time but have different feature id.
        """
        self.das = das  # list of volcat files which are valid for the same time
        t1 = das[0].time.values
        for ddd in das:
            if ddd.time.values != t1:
                print("WARNING. list has different times")

        self.mass = xr.DataArray()
        self.height = xr.DataArray()
        self.totalmass = None
        self.totalarea = None
        self.gridspace = None

    # no setter
    @property
    def dset(self):
        dset = xr.Dataset(
            {
                "ash_mass_loading": self.mass,
                "ash_cloud_height": self.height,
                "ash_mass_loading_total_mass": self.totalmass,
                "feature_area": self.totalarea,
            }
        )
        return dset

    def create(self, gridspace):
        self.gridspace = gridspace
        newlist = self.regrid(gridspace)
        self.combine(newlist)
        return self.dset

    def regrid(self, gridspace):
        newlist = []
        for dset in self.das:
            newlist.append(volcat.correct_pc(dset, gridspace=gridspace))
        return newlist

    def combine(self, newlist):
        masslist = []
        hlist = []
        totalmass = 0
        totalarea = 0
        for new in newlist:
            masslist.append(volcat.get_mass(new))
            hlist.append(volcat.get_height(new))
            totalmass += new.ash_mass_loading_total_mass.values
            totalarea += new.feature_area.values
        mass = volcat.combine_regridded(masslist)
        height = volcat.combine_regridded(hlist)

        time = list(set(mass.time.values))
        if len(time) > 1:
            write("warning, multiple times in OneEventTime")
        else:
            time = time[0]

        mass = mass.sum(dim="time")
        mass["time"] = time
        mass = mass.expand_dims(dim="time")

        height = height.max(dim="time")
        height["time"] = time
        height = height.expand_dims(dim="time")

        self.mass = mass
        self.height = height
        self.totalmass = totalmass
        self.totalarea = totalarea


class Events:
    """
    combines information from multiple Event Files.

    Properties
        VolcatPlots


    Attributes
        self.edf :   pandas DataFrame : information on event files
        self.ndir : str : file path.

        self.events     : list of xarray DataSets with volcat data
        self.pcevents   : list of xarray DataSets with volcat data parallax corrected.

        self.volcano_name : str

    WorkFlows
        1. add the self.edf dataframe.
           a. This may be done by using the output of flist2eventdf which creates a dataframe from list of volcat file names.
              this method is suitable for files which have already been downloaded.
           b. This may also be done from the output of get_log_files. This method is suitable for files which need to be downloaded.
           c. read a previously written csv file for this event.
        2. download event files if necessary.

        3. read event files into self.events

        4a. create the emit-times files from the self.events
        4b. create the parallax corrected files from the self.events
        4c. get the information for unit source runs for inverse modeling.

        NO. ONLY handle the VOLCAT data NOT model runs.
        5a. call to ash_run.py to start data insertion runs.
        5b. call to ash_run.py to start inverse modeling runs.

        6. read the parallax corrected files in self.pcevents



    """

    def __init__(self, eventid="TEST", volcano_name=None, eventlist=None, inp=None):
        """

        inp : dictionary with information for setting the directories.
        """
        # Events class

        # need to add using
        #         add_df method
        #         add_csv
        #         add_events

        self.eventdf = EventDF(edf=pd.DataFrame())
        if isinstance(eventlist, (list, np.ndarray)):
            self.eventdf.add_events(eventlist)
        elif isinstance(eventlist, (pd.core.frame.DataFrame)):
            self.eventdf.add_df(eventlist)

        self.set_volcano_name(vname=volcano_name)  # str volcano name

        # need to set using get_dir or set_dir methods.
        self.ndir = "./"  # directory where VOLCAT files are
        self.pdir = "./"  # directory where parallax corrected volcat files are
        self.edir = "./"  # directory where emit-times files are
        self.idir = "./"  # directory where inversion runs are
        # try the get_dir method.
        self.get_dir(inp)

        self.status = EventStatus(eventid, "Initialized")

        # send display output to the parallax corrected directory?
        self.display = EventDisplay(eventid, self.volcano_name, tdir=self.pdir)

        self.events = []  # list of xarray DataSets with volcat data
        self.pcevents = []  # list of xarray DataSets with volcat data
        self.maxi = 0

        self.eventid = eventid

    def plotsummary(self, smooth=20, yscale="linear"):
        if not self.events:
            logger.warning("May need to add events before plotting summary")
        self.display.add(self.events)
        self.display._vplot.plot_multiA(smooth=smooth, yscale=yscale)

    @property
    def edf(self):
        # always returns a copy of the dateframe.
        return self.eventdf.edf

    def check_val(self, val):
        """
        plots val vs. observation date.
        """
        # Events class
        dtemp = self.edf
        sns.set()
        for fid in dtemp[val].unique():
            dtemp2 = dtemp[dtemp[val] == fid]
            plt.plot(dtemp2["observation_date"], dtemp2[val], ".")
        fig = plt.gcf()
        fig.autofmt_xdate()
        ax = plt.gca()
        ax.set_ylabel(val)

    def check_sensor(self):
        self.check_val("sensor_name")

    def check_feature_id(self):
        self.check_val("feature_id")

    def set_dir(self, data_dir, parallax_dir, emit_dir, inv_dir=None, make=False):
        """
        sets the directories

        make : boolean : if TRUE then checks to see if directory exists and creates it if it doesn't.
        """
        if isinstance(data_dir, str):
            self.ndir = data_dir
            if make and not os.path.isdir(self.ndir):
                make_dir(self.ndir, None, verbose=True)
        if isinstance(parallax_dir, str):
            self.pdir = parallax_dir
            if make and not os.path.isdir(self.pdir):
                make_dir(self.pdir, None, verbose=True)
        if isinstance(emit_dir, str):
            self.edir = emit_dir
            if make and not os.path.isdir(self.edir):
                make_dir(self.edir, None, verbose=True)
        if isinstance(inv_dir, str):
            self.idir = inv_dir
            if make and not os.path.isdir(self.idir):
                make_dir(self.idir, None, verbose=True)

    def set_volcano_name(self, vname=None):
        """
        First tries to use volcano name from dataframe.
        If that does not exist will use vname if input as string.
        If that is None then set name to 'Unknown'
        """
        vlist = []
        vstr = [x for x in self.edf.columns if "VOLCANO_NAME" in x.upper()]
        # print(vstr)
        if vstr:
            vstr = vstr[0]
            vlist = self.edf[vstr].unique()
        if len(vlist) > 1:
            print("WARNING: not setup for multiple volcano in same Event class")
        elif len(vlist) == 1:
            vstr = vlist[0]
            self.volcano_name = fix_volc_name(vstr)
        elif isinstance(vname, str):
            self.volcano_name = fix_volc_name(vname)
        else:
            self.volcano_name = "Unknown"
        return self.volcano_name

    def get_dir(self, inp, verbose=False, make=True):
        """
        inp : dictionary with key VOLCAT_DIR
        set the directory from the dictionary inp.
        """
        # Events class
        if not isinstance(inp, dict):
            return None
        if "VOLCAT_DIR" not in inp.keys():
            logger.warning("get_dir method input does not contain VOLCAT_DIR")
            return None
        tdir = inp["VOLCAT_DIR"]
        if not self.volcano_name:
            self.set_volcano_name()
        if self.volcano_name != "Unknown":
            ndir = os.path.join(inp["VOLCAT_DIR"], self.volcano_name)
        else:
            ndir = inp["VOLCAT_DIR"]
        pdir = os.path.join(ndir, "pc_corrected")
        edir = os.path.join(ndir, "emitimes")
        idir = os.path.join(ndir, "inverse")

        if verbose:
            logger.info("Downloading to {}".format(ndir))
            logger.info("parallax corrected to {}".format(pdir))
            logger.info("emit times files to {}".format(edir))
        self.set_dir(ndir, pdir, edir, idir, make=make)
        return ndir

    def download(self, inp, daterange=None, verbose=False, log=None):
        """
        Obtains list of VOLCAT event files to download from self.edf.
        Obtains location to download files to from inp or from self.ndir.
        Checks to see if files exist in location and if they do not already exist,
        downloads the files.

        INPUTS
        inp : dictionary with information on directory
        verbose : boolean

        UTILIZES
        self.ndir : uses this or inp for place to download files to.
        self.edf   : event_url column provides urls for file downloads
        """
        # Events class
        df = self.edf.copy()
        if daterange:
            df = df[
                (df.observation_date > daterange[0])
                & (df.observation_date <= daterange[1])
            ]

        if not self.ndir:
            ndir = self.get_dir(inp)
        ndir = self.ndir
        if verbose:
            print("Downloading to {}".format(ndir))
        if "event_url" not in df.columns:
            logger.warning("event_url not available")
            return -1
        failed_list = []
        for eurl in df["event_url"]:
            file_download = check_file(eurl, ndir, verbose=True)
            if file_download:
                os.system("wget -P" + ndir + " " + eurl)
                file_download = check_file(eurl, ndir, verbose=True)
                if file_download:
                    print("FAILED to create {}".format(eurl))
                    jjj = df[df["event_url"] == eurl]
                    jjj = jjj.observation_date.values[0]
                    print("ZZZZ", jjj)
                    atemp = df[df["observation_date"] == jjj]
                    if len(atemp) > 1:
                        print(atemp[["observation_date", "feature_id", "event_file"]])
                        print("------------")
                    failed_list.append(eurl)
                else:
                    if verbose:
                        print("Downloaded", eurl)
            else:
                if verbose:
                    print("Already Downloaded", eurl)

        if isinstance(log, str) and failed_list:
            print("WRITING TO {}".format(log))
            if os.path.isfile(log):
                with open(log, "r") as lid:
                    current = lid.readlines()
            else:
                current = []
            newfailed = 0
            with open(log, "a") as lid:
                for eurl in failed_list:
                    if eurl + "\n" not in current:
                        lid.write(eurl)
                        lid.write("\n")
                        newfailed += 1
            print("NEW FAILURES in {} {}".format(log, newfailed))

    def get_closest_time(self, target_time, nmatches=1):
        """
        returns row of edf with cloest time.
        """
        df2 = self.edf
        df2 = df2[["observation_date", "event_file"]]
        # create column with time delta object showing time differences
        df2["diff"] = target_time - df2["observation_date"]
        # create column with absolute difference in hours.
        df2["hdiff"] = df2.apply(
            lambda row: np.abs(row["diff"].days * 24 + row["diff"].seconds / 3600.0),
            axis=1,
        )
        # get row with the smallest absolute difference
        best = df2.nsmallest(nmatches, "hdiff")
        return best

    @staticmethod
    def get_flist(df, ndir):
        """
        df : pandas dataframe : needs to have "event_file" as a column
        ndir : str : directoy where event files can be found
        RETURNS
        yeslist : list : event files which exist in ndir
        """
        flist = df["event_file"].unique()
        yeslist = [x for x in flist if os.path.isfile(os.path.join(ndir, x))]
        return yeslist

    def get_missing_flist(self):
        # Events class
        flist = self.edf["event_file"].unique()
        nolist = [x for x in flist if not os.path.isfile(os.path.join(self.ndir, x))]
        return nolist

    def check_for_volcat(self):
        vdir = self.ndir
        difiles = glob.glob(vdir + "/*VOLCAT.*.nc")
        df = flist2eventdf(difiles, {"VOLCANO_NAME": self.volcano_name})
        ## TO DO - read config file to give summary of run.
        return df

    def check_for_DI(self):
        emit_dir = self.edir
        difiles = glob.glob(emit_dir + "/xrfile.*.nc")
        ## TO DO - read config file to give summary of run.
        return difiles

    def check_for_inv(self):
        inv_dir = self.idir
        difiles = glob.glob(inv_dir + "/xrfile.*.nc")
        ## TO DO - read config file to give summary of run.
        return difiles

    def write_emit(self, overwrite=False, verbose=False):
        """
        write emit-times files.
        """
        emit_dir = self.edir
        volc_dir = self.ndir

        pollnum = 1
        pollpercents = [1]
        layer = 0.0

        date_time = None
        clip = False

        if not volc_dir or not emit_dir:
            logger.warning("Directories do not exist. use set_dir to add them")
        flist = self.get_flist(self.edf, volc_dir)
        for iii, flt in enumerate(flist):
            volcemit = mdi.InsertVolcat(
                emit_dir,
                volc_dir,
                date_time,
                fname=flt,
                layer=layer,
                pollnum=pollnum,
                pollpercents=pollpercents,
            )
            if not volcemit.check_for_file() or overwrite:
                oname = volcemit.write_emit(area_file=False, clip=clip, verbose=verbose)
                try:
                    oname = volcemit.write_emit(
                        area_file=False, clip=clip, verbose=verbose
                    )
                except Exception as eee:
                    logger.warning(eee)
                    logger.warning("emit file could not be written")
                    continue
                logger.info(
                    "Emit-times file written {}".format(volcemit.make_emit_filename())
                )
                print("CREATED", oname)
            else:
                logger.info(
                    "Emit-times file exists {}".format(volcemit.make_emit_filename())
                )

    def write_parallax_corrected(
        # Events class
        self,
        gridspace=None,
        daterange=None,
        verbose=False,
    ):

        pdir = self.pdir
        ndir = self.ndir
        if not pdir or not ndir:
            logger.warning("Directories do not exist. use set_dir to add them")

        df = self.edf
        if daterange:
            df = df[
                (df.observation_date > daterange[0])
                & (df.observation_date <= daterange[1])
            ]

        flist = self.get_flist(df, ndir)
        # flist = [flist[0]]
        print("Events class, Number of parallax files to write {}".format(len(flist)))
        # sys.exit()
        for iii, flt in enumerate(flist):
            if iii % 20 == 0:
                logger.info("Working on file {}".format(iii))
                print("Working on file {}".format(iii))
            # daterange only works with this if flist isn't used.
            volcat.write_parallax_corrected_files(
                ndir,
                pdir,
                flist=[flt],
                verbose=verbose,
                daterange=None,
                gridspace=gridspace,
            )

    def get_flistdf(self):
        """
        returns list of files sorted by date and feature id.
        """
        df2 = self.edf
        slist = ["observation_date", "feature_id"]
        alist = slist.copy()
        # alist.append("event_file")
        # alist.append("sensor_name")
        alist.append("event_file")
        alist.append("sensor_name")
        df2 = df2[alist].sort_values(by=slist, ascending=True)
        df2 = df2.reset_index()
        df2 = df2.drop("index", axis=1)
        # flist = df2['event_file'].values
        return df2

    def generate_combined(self, gridspace):
        for event in self.generate_fid:
            newlist = []
            ilist = event.index.values
            newlist = self.events[ilist[0], ilist[-1] + 1]
            oet = OneEventTime(newlist)
            oet.create(gridspace=gridspace)
            yield oet.dset

    def generate_fid(self):
        """
        find files which have same observation date in order to combine them.
        """
        jtemp = self.get_flistdf()
        jdt = jtemp["observation_date"].unique()
        for jjj in jdt:
            atemp = jtemp[jtemp["observation_date"] == jjj]
            if len(atemp) > 1:
                yield atemp

    def check_fid(self):
        """
        Indicates whether different features are contained in different files.
        """
        for fid in self.generate_fid():
            print(fid[["observation_date", "feature_id", "event_file"]])
            print("------------")
        print("done checking fid")

    def get_volcat_events(self, bysensor=None, daterange=None, verbose=False):
        # Events class
        # close any open files first.
        if not self.ndir:
            print("WARNING need directory information")
            return False
        for event in self.events:
            event.close()
        df2 = self.get_flistdf()
        if daterange:
            df2 = df2[
                (df2.observation_date > daterange[0])
                & (df2.observation_date <= daterange[1])
            ]
        if bysensor:
            df2 = df2[df2["sensor_name"] == bysensor]
        flist = self.get_flist(df2, self.ndir)  # df2['event_file'].values
        print("Number of files {}".format(len(flist)))
        das = volcat.get_volcat_list(
            self.ndir,
            flist=flist,
            decode_times=True,
            correct_parallax=False,
            verbose=verbose,
        )
        if verbose:
            print("get_volcat_events {} {}".format(len(das), len(flist)))

        def ftime(x):
            tval = x.time.values
            if isinstance(tval, list):
                return tval[0]
            return tval

        das.sort(key=ftime)
        # self.pcevents = das

        self.events = das
        self.maxi = len(das)
        return flist


    def subset_time(self,d1,d2):
        timelist = [pd.to_datetime(x.time.values[0]) for x in self.pcevents] 
        zzz = zip(timelist,self.pcevents)
        elist = [x for x in zzz if x[0]>=d1 and x[0]<d2]
        elist = list(zip(*elist))[1]
        aset = volcat.combine_regridded(elist) 
        return aset

    def average_mass(self,d1,d2):
        aset = self.subset_time(d1,d2)
        aset = aset.ash_mass_loading
        return aset.mean(dim='time')

    def max_height(self,d1,d2):
        aset = self.subset_time(d1,d2)
        aset = aset.ash_cloud_height
        return aset.max(dim='time')


    def get_volcat_events_pc(self, bysensor=None, verbose=False, daterange=None):
        for event in self.pcevents:
            event.close()

        df2 = self.get_flistdf()
        if daterange:
            df2 = df2[
                (df2.observation_date > daterange[0])
                & (df2.observation_date <= daterange[1])
            ]

        if bysensor:
            df2 = df2[df2["sensor_name"] == bysensor]
        flist = self.get_flist(df2, self.ndir)
        print("Number of files {}".format(len(flist)))
        wdir = os.path.join(self.ndir, "pc_corrected")
        flist = [x.replace(".nc", "_pc.nc") for x in flist]
        # flist = [flist[0]]
        das = volcat.get_volcat_list(
            wdir, flist=flist, correct_parallax=False, verbose=verbose
        )

        def ftime(x):
            return x.time.values[0]
        das.sort(key=ftime)
        self.pcevents = das
        self.maxi = len(das)
        return flist

    @property
    def vloc(self):
        if "volcano_lat" in self.edf.columns and "volcano_lon" in self.edf.columns:
            vloc = [self.edf.volcano_lat.unique()[0], self.edf.volcano_lon.unique()[0]]
        else:
            vloc = [-999, -999]
        self._vloc = vloc
        return self._vloc

    def get_vloc(self):
        if "volcano_lat" in self.edf.columns and "volcano_lon" in self.edf.columns:
            vloc = [self.edf.volcano_lat.unique()[0], self.edf.volcano_lon.unique()[0]]
        else:
            vloc = [-999, -999]
        return vloc

    # def boxplot(self, vplot, bstep=10):
    #    vplot.make_boxplot(np.arange(0, self.maxi - 1, bstep))

    # def vplots(self, clr="-k"):
    # event class
    #    from utilvolc import volcat_plots as vp

    #    vplot = vp.VolcatPlots(None,volcano_name=self.volcano_name)
    # read any data that has been saved in a csv file.
    # vplot.read()
    #     vplot.main_clr = clr
    #     vplot.make_arrays(self.events)
    # vplot.volcat_describe_plot()
    #     fig1 = vplot.plot_multiA(fignum=1, smooth=0.08, yscale="linear")
    # plt.show()
    # fig2 = vplot.plot_multiB(fignum=2)
    #     return vplot

    #def volcat2pol(self, iii):
    #    pcdas = self.pcevents[iii]
    #    vmass = volcat.get_mass(das[iii], clip=True)

    def compare_pc(
        self, pstep, daterange=None, fid=None, central_longitude=0, vlist=None
    ):
        # TODO this is diagnostic only. probably move elsewhere.

        from utilvolc import volcat_2dplots as v2d

        vloc = self.get_vloc()
        fidcheck = fid
        das = self.events
        pcdas = self.pcevents
        if not isinstance(vlist, list):
            vlist = list(np.arange(0, self.maxi, pstep))
        jtemp = self.get_flistdf()
        for iii in vlist:

            fid = das[iii].feature_id.values
            if fidcheck:
                if fid != fidcheck:
                    continue
            vmass = volcat.get_mass(das[iii], clip=True)
            pcmass = volcat.get_mass(pcdas[iii], clip=True)
            if daterange:
                if pd.to_datetime(vmass.time.values) < daterange[0]:
                    continue
                if pd.to_datetime(vmass.time.values) > daterange[1]:
                    continue

            transform = wep.get_transform(central_longitude=-180)
            fig, axarr = plt.subplots(
                nrows=1,
                ncols=2,
                figsize=(10, 5),
                constrained_layout=False,
                subplot_kw={"projection": transform},
            )
            axlist = axarr.flatten()
            ax = axlist[0]
            ax2 = axlist[1]
            print("feature id", das[iii].feature_id.values)
            print("feature id", pcdas[iii].feature_id.values)
            # sns.set()
            # print("sensor", jtemp["sensor_name"].values[iii])
            # print("feature id", jtemp["feature_id"].values[iii])
            print(vmass.time.values, pcmass.time.values)

            checkmass = volcat.check_total_mass(das[iii])
            checkmasspc = volcat.check_total_mass(pcdas[iii])
            print("mass {:0.3e} {:0.3e}".format(checkmass, checkmasspc))
            print("volcat mass {:0.3e}".format(volcat.get_total_mass(das[iii])))

            temp = vmass.isel(time=0)
            # plt.sca(ax)
            cb = ax.pcolormesh(
                temp.longitude.values,
                temp.latitude.values,
                np.log10(temp.values),
            )
            # parallax corrected on the right.
            plt.colorbar(cb)
            if vloc[0] != -999:
                ax.plot(vloc[1], vloc[0], "m^", markersize=5)
            temp = pcmass.isel(time=0)
            plt.sca(ax2)
            cb = ax2.pcolormesh(
                temp.longitude.values,
                temp.latitude.values,
                np.log10(temp.values),
            )
            plt.colorbar(cb)
            if vloc[0] != -999:
                ax2.plot(vloc[1], vloc[0], "m^", markersize=5)
            # ax = plt.gca()
            transform = wep.get_transform()
            wep.format_plot(ax, transform)
            wep.format_plot(ax2, transform)

            plt.tight_layout()
            plt.show()
            # plt.close()

    def plots_with_vaas(self, vaas, pstep=1, pc=True):
        vloc = self.get_vloc()

        def ftime(x):
            return x.time.values[0]

        das = self.pcevents
        das2 = self.events

        das.sort(key=ftime)
        das2.sort(key=ftime)
        vlist = list(np.arange(0, self.maxi, pstep))
        for iii in vlist:
            intime = pd.to_datetime(das[iii].time.values[0])
            ptime = pd.to_datetime(das2[iii].time.values[0])
            print(intime, ptime, intime == ptime)
            matches = vaas.find_time_match(vname=None, intime=intime, forecast=0, dt=1)
            if not matches:
                print("NO VAA for {}".format(intime))
                continue
            fig = plt.figure(1, figsize=(10, 5))
            ax = fig.add_subplot(1, 1, 1)
            vmass = volcat.get_mass(das[iii], clip=True)
            vmass2 = volcat.get_mass(das2[iii], clip=True)
            temp = vmass.isel(time=0)
            temp2 = vmass2.isel(time=0)
            cb = plt.pcolormesh(
                temp.longitude.values,
                temp.latitude.values,
                np.log10(temp.values),
                cmap="Reds",
            )
            cb2 = plt.pcolormesh(
                temp2.longitude.values,
                temp2.latitude.values,
                np.log10(temp2.values),
                cmap="Blues",
            )
            plt.colorbar(cb)
            plt.colorbar(cb2)
            plt.plot(vloc[1], vloc[0], "m^", markersize=5)
            for mmm in matches:
                vaa = vaas.ilist[mmm]
                vaa.plot_vaa(ax=ax)
                plt.title(intime)
            plt.show()

    def plots(
        self,
        pstep,
        pc=False,
        levels=[0.02, 0.2, 0.3, 2, 5, 10, 50],
        vlist=None,
        central_longitude=0,
    ):
        from matplotlib.colors import BoundaryNorm
        import cartopy
        from utilhysplit.plotutils import vtools

        transform = cartopy.crs.PlateCarree(central_longitude=central_longitude)
        volcat_transform = cartopy.crs.PlateCarree(central_longitude=0)

        vloc = self.get_vloc()
        if pc:
            das = self.pcevents
        else:
            das = self.events
        if isinstance(vlist, int):
            vlist = [vlist]
        elif isinstance(vlist, (list, np.ndarray)):
            vlist = vlist
        else:
            vlist = list(np.arange(0, self.maxi, pstep))
        for iii in vlist:
            fig, axrr = plt.subplots(
                nrows=1,
                ncols=2,
                figsize=(10, 5),
                constrained_layout=True,
                subplot_kw={"projection": transform},
            )
            axlist = axrr.flatten()
            ax = axlist[0]
            ax2 = axlist[1]
            # ax = fig.add_subplot(1, 2, 1)
            # ax2 = fig.add_subplot(1, 2, 2)
            vht = volcat.get_height(das[iii], clip=True)
            sns.set()
            print(iii)
            print("total mass", das[iii].ash_mass_loading_total_mass.values)
            print("area", das[iii].feature_area.values)
            try:
                print("instrument---", das[iii].attrs["instrument_ID"])
            except Exception as eee:
                print("EXCEPTION in plots", eee)
                print(das[iii].attrs)
                pass
            try: 
                vht.isel(time=0).plot.pcolormesh(
                    ax=ax,
                    x="longitude",
                    y="latitude",
                    cmap="Reds",
                    transform=volcat_transform,
                )
            except:
                print('Failed at ', iii)
                return das[iii] 

            # ax.set_xlim(xmin,xmax)
            # ax.set_ylim(ymin,ymax)
            # plt.tight_layout()
            # plt.show()
            # plt.savefig('./animations/bezy_volcat_mass{:03d}.png'.format(iii))
            # plt.close()

            # plt.savefig('bezy_volcat_2040_ht.png')
            # print(np.max(vht))
            # plt.show()

            mass = volcat.get_mass(das[iii], clip=True)
            # sns.set()
            # vht.isel(time=0).plot.pcolormesh(ax=ax2, x='longitude',y='latitude',cmap='Reds')
            temp = mass.isel(time=0)
            # cb = plt.pcolormesh(
            #    temp.longitude.values, temp.latitude.values, np.log10(temp.values)
            # )
            cmap = plt.get_cmap("viridis")
            norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
            # mass.isel(time=0).plot.pcolormesh(
            #    ax=ax2, x="longitude", y="latitude", cmap="Reds",transform=volcat_transform
            #    )

            cb = ax2.pcolormesh(
                temp.longitude.values,
                temp.latitude.values,
                temp.values,
                norm=norm,
                transform=volcat_transform,
            )
            plt.colorbar(cb)
            ax2.plot(vloc[1], vloc[0], "m^", markersize=5, transform=volcat_transform)
            # plt.savefig('./animations/bezy_volcat_mass{:03d}.png'.format(iii))
            # print(np.max(vht))
            # plt.tight_layout()
            # plt.close()
            wep.format_plot(ax, transform)
            wep.format_plot(ax2, transform)
            plt.show()

        return ax, ax2, temp


def create_event_from_fnames(inp, vname, volclistfile=None):
    """
    Create Event object  from the filenames of the event netcdf files.
    """

    if not isinstance(volclist, str):
        volclist = "/hysplit-users/alicec/ashapp/data/volclist.txt"
    vprops = volcano_names.VolcList(volclistfile)
    vname = fix_volc_name(vname)
    files = glob.glob("{}/{}/{}".format(inp["VOLCAT_DIR"], vname, "*VOLCAT*nc"))
    record = vprops.get_record(vname)
    record = record.to_dict(orient="records")[0]
    logdf = volcat.flist2eventdf(files, record)
    eve = Events(volcano_name=vname, eventlist=logdf, inp=inp)
    return eve


def create_event(indf, fdir, vlist=None, verbose=False):
    """
    Create Event object from dataframe with logfiles.

    indf : DataFrame with the following columns
          : volcano_name - str
          : log - str - name of EventFile logfile (.json file)
          :
    fdir  : location of log files
          : can be string
          :
    vlist : list of volcanoes to create  Events objects for
    verbose : boolean

    output

    eventhash : dictionary
                key is volcano name. value is Events object.
    """
    required = ["volcano_name", "log"]
    for req in required:
        if req not in indf.columns:
            logger.warning(
                "create_event function required column {} not available in input".format(
                    req
                )
            )

    if isinstance(fdir, dict):
        inp = fdir
        fdir = fdir["VOLCAT_LOGFILES"]
    else:
        inp = {}
    eventhash = {}  # key is volcano name. value is Event class instance

    # determining vlist - list of volcano names
    if not isinstance(vlist, (list, np.ndarray)):
        vlist = indf["volcano_name"].unique()
    else:
        # so vlist doesn't have to be case dependent.
        templist = indf["volcano_name"].unique()
        vlist2 = []
        for vvv in vlist:
            temp = [x for x in templist if vvv.lower() in x.lower()]
            vlist2.extend(temp)
        vlist = vlist2

    # loop through each volcano name.
    for volc in vlist:
        print(volc)
        eventlist = []
        df2 = indf[indf["volcano_name"] == volc]
        efile = df2["log"].unique()

        for efl in efile:
            if not isinstance(efl, str):
                print("NOT string", efl)
                continue
            evo = EventFile(efl, fdir)
            success = evo.open()
            if not success:
                print("cannot open {} {}".format(fdir, efl))
                print("\n")
                continue
            eventlist.append(evo)
        print("EVENTLIST", len(eventlist))
        if eventlist:
            eve = Events(volcano_name=volc, eventlist=eventlist, inp=inp)

            # set the directories for the event.
            if isinstance(fdir, dict):
                eve.get_dir(inp, verbose=verbose, make=True)
            eventhash[volc] = eve
        else:
            print("NO events found")

    return eventhash

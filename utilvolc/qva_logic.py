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

from utilvolc.ashutil import fix_volc_name
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


logger = logging.getLogger(__name__)

"""
 Workflow for creating QVA (quantitative volcanic ash) forecasts

 Changelog
 2022 Dec 8  AMC  added date input to function get_summary_file_df
 2023 Feb 28 AMC  added error checking for empty files to get_summary_file_df
 2023 May 05 AMC  fixed bug in get_summary_file_df
 2023 May 05 AMC  added forecast2QVA
 2023 May 05 AMC  changed write_emitimes to write_di_emitimes
 2023 May 05 AMC  added mapping with cartopy to plots method
 2023 May 05 AMC  added write_emit method to Events class
 2023 May 05 AMC  added check_for_DI, check_for_inv, check_for_volcat  method to Events class
 2023 May 05 AMC  started adding a Region class, DisplayEvent class, match_time function
 2023 May 26 AMC  moved some functions / classes to volcat_files.py
 



 Issue = vl.get_files pauses and asks if you want to overwrite the logfiles.
 TODO  -  figure out what permissions need to be set for files.
 TODO  -  generation of parallax corrected files is slow.
 TODO  - generate forecasts from TCM rather than a new run of HYSPLIT.
"""

# Prata 2018 uses dosage 


def workflow():
    """
     Set the directories.
        # inp = {}
        # inp['JPSS_DIR'] = '/pub/jpsss_upload'
        # inp['VOLCAT_LOGFILES'] = '/pub/ECMWF/JPSS/VOLCAT/LogFiles/'
        # inp['VOLCAT_DIR']= '/pub/ECMWF/JPSS/VOLCAT/Files/'

    # WorkFlow object has following attributes
        # dirhash
        # greenlist
        # sumdf
        # logdf
        # ehash

    # hours = 7*24
    ## Create instance of WorkFlow class
    # work = volcat_logic.WorkFlow(inp,verbose=False)
    ## create the sumdf attribute
    ## download log files
    # work.get_volcat(hours=hours)

    # Check if more than one file for an observation date.
    # This can occur if two files have different feature ID's
      # work.check_fid()

    #------------------------------------    
    # TYPES OF FILES
    #------------------------------------    
    # Event summary files. 
            SummaryFile class. json file. pushed to ARL.
    # Event Files. 
            EventFile class. json file. pulled.
    # Data Files. 
            see classes and methods in volcat.py. netcdf files. pulled.  
    # 
    """
    #------------------------------------    
    # COMPLETED:
    # done - check for event summary files and read
    # use a text log file to keep track.
    # (later?) - decide which file list json files to pull (may not be needed).
    # (done) - pulling all event log files (json format).
    # (done) check if modified on their site. if not modified don't pull them again.
    #  no log file.
    # (done) specify vaac region to pull netcdf files from event log files
    # (done) -  pull all event files (netcdf). Organized by volcano name (folder)
    # checks to see if file already exists. only pulls non-existing files.
    # (later?) check if file meets other requirements (time resolution???) as needed.
    # (done) - parallax corrected files automatically generated
    # (done) - generate plots of total mass, total area, max top height for event (defined by events in event log file). Uses volcplot.py functions.
    # (done) - function generated to list unique eruption event times in each volcano folder so images can be created for individual eruption events
    # (done) Make list of workflow for functions
    # (done) Function to write red_list.txt based on green and yellow list files. Green and
    # yellow list files are written by hand.
    # (done) make emit-times files - might want to modify for more flexibility
    # (done) area calculation is necessary for emitimes files generation
    # no separate area file generated.
    # (done) fixed naming of emit-times file issues - now includes event datetime, imagedatetime, and volcano id in name of file

    # IN PROGRESS:
    # functions can be added
    # TO DO: combine g001 g002 g003 etc. files.
    #        for now only use g001 but will need to add them together later.
    # SANGAY eruption may have some examples of this.
    # a) Do we Need to merge files that have the same timestamp?
    #    Alternative is to just do separate HYSPLIT runs for them.
    #    However need to be careful of combining them for ensemble relative frequency then.
    #    How good a classifier is the event time? (Probably not great?)
    #     (i) CASE 1 is that they have the same image date time.
    #                                  the same event date time.
    #                                  different image identifier g001, g002, g003 etc.
    #         we think that in this case the ash is probably close together and could be
    #         merged for 1 emit-times file.

    #     (ii) CASE 2 is that they have the same image date time.
    #                                  the different event date time.
    #                                  same or different g001, g002, g003 etc.
    #         we think that in this case the ash clouds are more likely to be far apart.
    #         may want seperate emit-times files - what would we do with them?
    #             may be useful for evaluation.
    #             will have to be combined for the forecast.

    #       The problem with the seperate is that you have to be careful when you combine
    #       them into the probabilistic forecast.

    # b) Keep track of the event dates is useful
    # starting dispersion runs from volcat and emitimes files

    # in progress: false alarm logic
    # List volcanoes that are green light - yes it is probably real
    # Volcanoes that are active
    # List volcanoes that are red light - we dont think this is active
    # Delete files in red category after a week or so
    # List volcanoes that are yellow light - this should be checked by a person
    # process but dont post to webpage
    # (skip for now?) some decision logic about when and what kind of HYSPLIT runs we make.

    # NEXT STEPS:
    # automatic runs are they triggered by a request from web.
    # Limited number of runs for active volcano and then link to READY website
    # where they can choose to modify inputs.
    # populate with latitude longitude what satellite retrievals are available.etc

    # run HYSPLIT data insertion
    # write CONTROL file
    # write SETUP file

    # regrid volcat to HYSPLIT grid
    # needed for inverse modeling
    # needed for evaluation
    # write regridded file for every volcat retrieval
    # Perform the regridding at the same time as the parallax correction.

    # run HYSPLIT inverse modeling
    # make unit source runs (maybe triggered by web request before data available?)
    #      (Note number 5 could be an input)
    #      run duration 5 hours (this means you can only use observations out 5 hours)
    #      runs that start at hour 0 - 5 hour duration
    #      runs that start at hour 1 - 4 hour duration
    #      runs that start at hour 2 - 3 hour duration
    #      runs that start at hour 3 - 2 hour duration
    #      runs that start at hour 4 - 1 hour duration
    #      etc. so they all end at the same time.
    # pick observation(s) to use in the run.
    #      at first can only use observations out to hour 5.
    # For eruptions lasting longer than 5 hours.
    #      pick up the pardump files from the last runs and continue them for longer.
    #      runs that start at hour 5 - 5 hour duration
    #      runs that start at hour 6 - 4 hour duration
    #
    # solve the TCM
    # make the forecast using the source term.

    # Postprocessing
    # generate ensemble netcdfs.
    # ensemble weighting.
    # graphics

    # Evaluation against observations.
    return 0



class Region:
    """
    matches different events with the same time
    """

    def __init__(self,target_time,tres,region_name='Region'):
        self.region_name = region_name
        self.eventlist = []

    def add_event(self, event):
        dset = self.match_time(event)
        self.eventlist.extend(dset)
    
def match_time(events,target_time,tres):
    matches = []
    tres = datetime.timedelta(hours=tres/60.0)
    for eee in events:
        stime = pd.to_datetime(eee.time.values[0])
        if stime > target_time:
           if stime-target_time < tres: matches.append(eee)
        else:
           if target_time -stime < tres: matches.append(eee)
    return matches


#def plot_times(events):
#    fig = plt.figure(1)
#    ax = fig.add_subplot(1,1,1)
#    for eve in events:
     
class DisplayEvent:
    """
    Helper class for Events
    """
  
    def __init__self(eventid='TEST'):
        self.eventid = eventid




class EventStatus:
    """
    Helper class for Events
    """
  
    def __init__self(eventid='TEST'):
        self.eventid = eventid


class Events:
    """
    combines information from multiple Event Files.

    Attributes
        self.df :   pandas DataFrame : information on event files
        self.ndir : str : file path.

        self.events     : list of xarray DataSets with volcat data
        self.pcevents   : list of xarray DataSets with volcat data parallax corrected.

        self.volcano_name : str

    WorkFlows
        1. add the self.df dataframe.
           a. This may be done by using the output of flist2eventdf which creates a dataframe from list of volcat file names.
              this method is suitable for files which have already been downloaded.
           b. This may also be done from the output of get_log_files. This method is suitable for files which need to be downloaded.
           c. read a previously written csv file for this event.
        2. download event files if necessary.

        3. read event files into self.events

        4a. create the emit-times files from the self.events
        4b. create the parallax corrected files from the self.events
        4c. get the information for unit source runs for inverse modeling.


        5a. call to ash_run.py to start data insertion runs.
        5b. call to ash_run.py to start inverse modeling runs.
       
        6. read the parallax corrected files in self.pcevents

        

    """

    def __init__(self,eventid='TEST'):
        # Events class

        # need to add using 
        #         add_eventdf method 
        #         add_csv
        #         add_events
        self.df = pd.DataFrame()

        # need to set using get_dir or set_dir methods.
        self.ndir = None # directory where VOLCAT files are
        self.pdir = None # directory where parallax corrected volcat files are
        self.edir = None # directory where emit-times files are
        self.idir = None # directory where inversion runs are

        self.events = []    # list of xarray DataSets with volcat data
        self.pcevents = []  # list of xarray DataSets with volcat data
        self.maxi = 0

        self.eventid = eventid
        self.volcano_name = None # str volcano name

    def save_csv(self, inp=None):
        """
        save dataframe as csv file
        """
        # Events class
        if not self.ndir:
            self.get_dir(inp)
        fname = os.path.join(self.ndir, "Events.csv")
        self.df.to_csv(fname, float_format="%.1f", index=False)

    def read_csv(self, cname=None):
        """
        reads csv file that was previously saved
        """
        if not self.ndir:
           logger.warning('No data directory set.')
        if not cname:
            cname = os.path.join(self.ndir, "Events.csv")
        else:
            cname = os.path.join(self.ndir, cname)
        dtp = {"observation_date": "datetime64[ns]"}
        df = pd.read_csv(cname, sep=",", parse_dates=["observation_date"])
        return df

    def add_csv(self, cname=None):
        """
        reads csv file that was previously saved and adds it to dataframe.
        """
        # Events class
        dftemp = self.read_csv(cname)
        if not self.df.empty:
            dftemp = pd.concat(self.df, dftemp)
        dftemp.drop_duplicates()
        self.df = dftemp


    def add_eventdf(self, df):
        """
        df : pandas DataFrame
        """
        # can utilize output of flist2eventdf
        # change all column names to lower case
        columns = df.columns
        newc = [x.lower() for x in columns]
        df.columns = newc

        if self.df.empty:
            if isinstance(df, pd.core.frame.DataFrame):
                self.df = df
        else:
            if isinstance(df, pd.core.frame.DataFrame):
                self.df = pd.concat([self.df, df])
        if 'volcano_lat' not in self.df.columns:
           print('add_events warning. no volcano_lat column', self.df.columns)

    def add_events(self, eventlist):
        """
        eventlist: list of EventFile objects.
        """
        # Events class
        elist = []
        for eve in eventlist:
            df = eve.df.copy()
            #print('HERE----', df.columns)
            
            for key in eve.attrs:
                df[key] = eve.attrs[key]
            elist.append(df)

        if self.df.empty:
            if isinstance(elist, list):
                self.df = pd.concat(elist)
            elif isinstance(elist, pd.core.frame.DataFrame):
                self.df = elist
        else:
            if isinstance(elist, list):
                df = pd.concat(elist)
                self.df = pd.concat([self.df, df])
            elif isinstance(elist, pd.core.frame.DataFrame):
                self.df = pd.concat([self.df, elist])

        columns = self.df.columns
        newc = [x.lower() for x in columns]
        self.df.columns = newc
        #print('NEWC ---', self.df.columns)

        if 'volcano_lat' not in self.df.columns:
           print('add_events warning. not VOLCANO_LAT column',  self.df.columns)


    def check_val(self, val):
        """
        plots val vs. observation date.
        """
        # Events class
        dtemp = self.df.copy()
        sns.set()
        for fid in dtemp[val].unique():
            dtemp2 = dtemp[dtemp[val] == fid]
            plt.plot(dtemp2["observation_date"], dtemp2[val], ".")
        fig = plt.gcf()
        fig.autofmt_xdate()
        ax = plt.gca()
        ax.set_ylabel(val)

    def check_sensor(self):
        # Events class
        self.check_val("sensor_name")

    #def filter_sensor(self, sensorname):
    #    # Events class
    #    dtemp = self.df.copy()
    #    dtemp = dtemp[dtemp['sensor_name']==sensorname]
    #    self.df = dtemp

    def check_feature_id(self):
        # Events class
        self.check_val("feature_id")
        # dtemp = self.df
        # sns.set()
        # for fid in dtemp["feature_id"].unique():
        #    dtemp2 = dtemp[dtemp["feature_id"] == fid]
        #    plt.plot(dtemp2["observation_date"], dtemp2["feature_id"], ".")
        # fig = plt.gcf()
        # fig.autofmt_xdate()
        # ax = plt.gca()
        # ax.set_ylabel("feature id")

    def set_dir(self, data_dir, parallax_dir, emit_dir,inv_dir=None,make=False):
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


    def set_volcano_name(self,vname=None):
        """
        First tries to use volcano name from dataframe.
        If that does not exist will use vname if input as string.
        If that is None then set name to 'Unknown'
        """

        vlist = [] 
        vstr = [x for x in self.df.columns if 'VOLCANO_NAME' in x.upper()]
        #print(vstr)
        if(vstr): 
           vstr = vstr[0] 
           vlist = self.df[vstr].unique()
        if len(vlist) > 1:
            print("WARNING: not setup for multiple volcano in same Event class")
        elif len(vlist) == 1:
            vstr = vlist[0]
            self.volcano_name = fix_volc_name(vstr)
        elif isinstance(vname,str):
            self.volcano_name = fix_volc_name(vname)        
        else:
            self.volcano_name = 'Unknown'
        return self.volcano_name

    def get_dir(self, inp, verbose=False,make=True):
        """
        inp : dictionary with key VOLCAT_DIR
        set the directory from the dictionary inp.
      
        """
        # Events class
        tdir = inp["VOLCAT_DIR"]
        if not self.volcano_name: self.set_volcano_name() 
        ndir = os.path.join(inp['VOLCAT_DIR'], self.volcano_name)
        pdir = os.path.join(ndir, 'pc_corrected')
        edir = os.path.join(ndir, 'emitimes')
        idir = os.path.join(ndir, 'inverse')

        if verbose:
            logger.info("Downloading to {}".format(ndir))
            logger.info("parallax corrected to {}".format(pdir))
            logger.info("emit times files to {}".format(edir))
        self.set_dir(ndir,pdir,edir,idir,make=make)
        return ndir


    def download(self, inp, daterange=None, verbose=False, log=None):
        """
        Obtains list of VOLCAT event files to download from self.df. 
        Obtains location to download files to from inp or from self.ndir.
        Checks to see if files exist in location and if they do not already exist,
        downloads the files.        

        INPUTS
        inp : dictionary with information on directory
        verbose : boolean
        
        UTILIZES
        self.ndir : uses this or inp for place to download files to.
        self.df   : event_url column provides urls for file downloads 
        """
        # Events class
        df = self.df.copy()
        if daterange:
           df = df[(df.observation_date>daterange[0]) & (df.observation_date<=daterange[1])]

        if not self.ndir:
            ndir = self.get_dir(inp)
        ndir = self.ndir   
        if verbose:
            print("Downloading to {}".format(ndir))
        if 'event_url' not in df.columns:
            logger.warning('event_url not available')
            return -1
        failed_list = []
        for eurl in df["event_url"]:
            file_download = check_file(eurl, ndir, verbose=True)
            if file_download:
                os.system("wget -P" + ndir + " " + eurl)
                file_download = check_file(eurl, ndir, verbose=True)
                if file_download:
                    print("FAILED to create {}".format(eurl))
                    jjj = df[df['event_url']==eurl]
                    jjj = jjj.observation_date.values[0]
                    print('ZZZZ', jjj)                  
                    atemp = df[df["observation_date"] == jjj]
                    if len(atemp) > 1:
                        print(atemp[["observation_date", "feature_id",'event_file']]) 
                        print("------------")
                    failed_list.append(eurl)
                else:
                    if verbose:
                        print("Downloaded", eurl)
            else:
                if verbose:
                    print("Already Downloaded", eurl)

        if isinstance(log,str) and failed_list:
            print('WRITING TO {}'.format(log))
            if os.path.isfile(log):
                with open(log, 'r') as lid:
                     current = lid.readlines()
            else:
                current = []
            newfailed = 0
            with open(log, 'a') as lid:
                for eurl in failed_list:
                   if eurl + '\n' not in current:
                       lid.write(eurl) 
                       lid.write('\n')
                       newfailed += 1
            print('NEW FAILURES in {} {}'.format(log, newfailed)) 

    def get_closest_time(self, target_time, nmatches=1):
        # Events class
        df2 = self.df.copy()
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
    def get_flist(df,ndir):
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
        flist = self.df["event_file"].unique()
        nolist = [x for x in flist if not os.path.isfile(os.path.join(self.ndir, x))]
        return nolist

    def check_for_volcat(self):
        vdir = self.ndir
        difiles = glob.glob(vdir + '/*VOLCAT.*.nc')
        df = flist2eventdf(difiles, {'VOLCANO_NAME':self.volcano_name})
        ## TO DO - read config file to give summary of run.
        return df

    def check_for_DI(self):
        emit_dir = self.edir
        difiles = glob.glob(emit_dir + '/xrfile.*.nc')
        ## TO DO - read config file to give summary of run.
        return difiles

    def check_for_inv(self):
        inv_dir = self.idir
        difiles = glob.glob(inv_dir + '/xrfile.*.nc')
        ## TO DO - read config file to give summary of run.
        return difiles

    def write_emit(self,overwrite=False,verbose=False):
        """
        write emit-times files.
        """
        emit_dir = self.edir
        volc_dir = self.ndir
 
        pollnum=1
        pollpercents = [1]
        layer = 0.0

        date_time=None
        clip = False

        if not volc_dir or not emit_dir:
           logger.warning('Directories do not exist. use set_dir to add them')
        flist = self.get_flist(self.df,volc_dir)
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
                try:
                    oname = volcemit.write_emit(area_file=False, clip=clip, verbose=verbose)
                except Exception as eee:
                    logger.warning(eee)
                    logger.warning('emit file could not be written')
                    continue
                logger.info(
                    "Emit-times file written {}".format(volcemit.make_emit_filename())
                )
                print('CREATED', oname)
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
           logger.warning('Directories do not exist. use set_dir to add them')

        df = self.df.copy()
        if daterange:
           df = df[(df.observation_date>daterange[0]) & (df.observation_date<=daterange[1])]

        flist = self.get_flist(df,ndir)
        #flist = [flist[0]]
        print("Events class, Number of parallax files to write {}".format(len(flist)))
        #sys.exit()
        for iii, flt in enumerate(flist):
            if iii%20 == 0:
               logger.info('Working on file {}'.format(iii))
               print('Working on file {}'.format(iii))
            # daterange only works with this if flist isn't used.      
            volcat.write_parallax_corrected_files(
                ndir,
                pdir,
                flist=[flt],
                verbose=verbose,
                daterange = None,
                gridspace=gridspace,
            )

    def get_flistdf(self):
        # Events class
        """
        returns list of files sorted by date and feature id.
        """
        df2 = self.df.copy()
        slist = ["observation_date", "feature_id"]
        alist = slist.copy()
        #alist.append("event_file")
        #alist.append("sensor_name")
        alist.append("event_file")
        alist.append("sensor_name")
        df2 = df2[alist].sort_values(by=slist, ascending=True)
        # flist = df2['event_file'].values
        return df2

    def check_fid(self):
        # Events class 
        jtemp = self.get_flistdf()
        jdt = jtemp["observation_date"].unique()
        for jjj in jdt:
            atemp = jtemp[jtemp["observation_date"] == jjj]
            if len(atemp) > 1:
                print(atemp[["observation_date", "feature_id",'event_file']]) 
                print("------------")
        print("done checking fid")

    def get_volcat_events(self, bysensor=None, daterange=None, verbose=False):
        # Events class 
        # close any open files first.
        if not self.ndir:
           print('WARNING need directory information')
           return False 
        for event in self.events:
            event.close()
        df2 = self.get_flistdf()
        if daterange:
           df2 = df2[(df2.observation_date>daterange[0]) & (df2.observation_date<=daterange[1])]
        if bysensor:
            df2 = df2[df2["sensor_name"] == bysensor]
        print(df2)
        flist = self.get_flist(df2,self.ndir)  # df2['event_file'].values
        print('Number of files {}'.format(len(flist)))
        das = volcat.get_volcat_list(
            self.ndir, flist=flist, decode_times=False,correct_parallax=False, verbose=verbose
        )
        if verbose:
            print("get_volcat_events {} {}".format(len(das),len(flist)))

        def ftime(x):
            return x.time.values[0]
        das.sort(key=ftime)
        #self.pcevents = das

        self.events = das
        self.maxi = len(das)
        return flist

    def get_volcat_events_pc(self, bysensor=None, verbose=False,daterange=None):
        for event in self.pcevents:
            event.close()


        df2 = self.get_flistdf()
        if daterange:
           df2 = df2[(df2.observation_date>daterange[0]) & (df2.observation_date<=daterange[1])]

        if bysensor:
            df2 = df2[df2["sensor_name"] == bysensor]
        flist = self.get_flist(df2,self.ndir)
        print('Number of files {}'.format(len(flist)))
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

    def get_vloc(self):
        if 'volcano_lat' in self.df.columns and 'volcano_lon' in self.df.columns:
            vloc = [self.df.volcano_lat.unique()[0], self.df.volcano_lon.unique()[0]]
        else:
            vloc = [-999,-999]
        return vloc

    def boxplot(self, vplot, bstep=10):
        vplot.make_boxplot(np.arange(0, self.maxi - 1, bstep))

    def vplots(self, clr="-k"):
        # event class
        from utilvolc import volcat_plots as vp

        vplot = vp.VolcatPlots(self.events)
        vplot.main_clr = clr
        vplot.make_arrays()
        vplot.volcat_describe_plot()
        fig1 = vplot.plot_multiA(fignum=1, smooth=0.08, yscale="linear")
        #plt.show()
        #fig2 = vplot.plot_multiB(fignum=2)
        return vplot

    def volcat2poly(self,iii):
        pcdas = self.pcevents[iii]
        vmass = volcat.get_mass(das[iii], clip=True)

    def compare_pc(self, pstep, daterange=None, fid=None, central_longitude=0,vlist=None):
        vloc = self.get_vloc()
        fidcheck = fid 
        das = self.events
        pcdas = self.pcevents
        if not isinstance(vlist,list):
            vlist = list(np.arange(0, self.maxi, pstep))
        jtemp = self.get_flistdf()
        for iii in vlist:
            #fig = plt.figure(1, figsize=(10, 5))
            #ax = fig.add_subplot(1, 2, 1)
            #ax2 = fig.add_subplot(1, 2, 2)
            #ax = map_util.draw_map(1,ax)
            #ax2 = map_util.draw_map(1,ax2)

            fid = das[iii].feature_id.values[0]
            if fidcheck:
               if fid != fidcheck: continue
            vmass = volcat.get_mass(das[iii], clip=True)
            pcmass = volcat.get_mass(pcdas[iii], clip=True)
            if daterange:
               if pd.to_datetime(vmass.time.values[0]) < daterange[0]: continue
               if pd.to_datetime(vmass.time.values[0]) > daterange[1]: continue

            transform = wep.get_transform(central_longitude=-180)
            fig,axarr = plt.subplots(nrows=1,ncols=2,figsize=(10,5),
                                     constrained_layout=False,
                                     subplot_kw={"projection":transform})
            axlist = axarr.flatten()
            ax = axlist[0]
            ax2 = axlist[1]
            print("feature id", das[iii].feature_id.values) 
            print("feature id", pcdas[iii].feature_id.values) 
            #sns.set()
            #print("sensor", jtemp["sensor_name"].values[iii])
            #print("feature id", jtemp["feature_id"].values[iii])
            print(vmass.time.values, pcmass.time.values)

            checkmass = volcat.check_total_mass(das[iii])
            checkmasspc = volcat.check_total_mass(pcdas[iii])
            print("mass {:0.3e} {:0.3e}".format(checkmass, checkmasspc))
            print("volcat mass {:0.3e}".format(volcat.get_total_mass(das[iii])))

            temp = vmass.isel(time=0)
            #plt.sca(ax)
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
            #ax = plt.gca()
            transform = wep.get_transform()
            wep.format_plot(ax,transform)
            wep.format_plot(ax2,transform)

            plt.tight_layout()
            plt.show()
            #plt.close()

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

    def plots(self, pstep, pc=False, levels=[0.02,0.2,0.3,2,5,10,50],vlist=None,central_longitude=0):
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
        if isinstance(vlist,int):
            vlist = [vlist]
        elif isinstance(vlist,(list,np.ndarray)):
            vlist = vlist
        else:
            vlist = list(np.arange(0, self.maxi, pstep))
        for iii in vlist:
            fig,axrr = plt.subplots(nrows=1,ncols=2,figsize=(10, 5),
                       constrained_layout=True,subplot_kw={'projection':transform})
            axlist = axrr.flatten()
            ax = axlist[0]
            ax2 = axlist[1]
            #ax = fig.add_subplot(1, 2, 1)
            #ax2 = fig.add_subplot(1, 2, 2)
            vht = volcat.get_height(das[iii], clip=True)
            sns.set()
            print(iii)
            print("total mass",  das[iii].ash_mass_loading_total_mass.values[0])
            print("area", das[iii].feature_area.values[0])
            try:
                print("instrument---", das[iii].attrs['instrument_ID'])
            except Exception as eee:
                print('EXCEPTION ', eee)
                print(das[iii].attrs)
                pass 
            vht.isel(time=0).plot.pcolormesh(
                ax=ax, x="longitude", y="latitude", cmap="Reds",transform=volcat_transform
            )

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
            #cb = plt.pcolormesh(
            #    temp.longitude.values, temp.latitude.values, np.log10(temp.values)
            #)
            cmap = plt.get_cmap('viridis')
            norm = BoundaryNorm(levels,ncolors=cmap.N,clip=False)
            #mass.isel(time=0).plot.pcolormesh(
            #    ax=ax2, x="longitude", y="latitude", cmap="Reds",transform=volcat_transform
            #    )

            cb = ax2.pcolormesh(
                temp.longitude.values, temp.latitude.values, temp.values,
                norm=norm,transform=volcat_transform
                )
            plt.colorbar(cb)
            ax2.plot(vloc[1], vloc[0], "m^", markersize=5,transform=volcat_transform)
            # plt.savefig('./animations/bezy_volcat_mass{:03d}.png'.format(iii))
            # print(np.max(vht))
            #plt.tight_layout()
            #plt.close()
            wep.format_plot(ax,transform)
            wep.format_plot(ax2,transform)
            plt.show()

        return ax, ax2,temp
        
    # def plots(self):
    #    flist = self.df['event_file'].values
    #    das = volcat.get_volcat_list(self.ndir, flist=flist, correct_parallax=False, verbose=False)
    #    self.das = das
    #    maxi = len(das)




class WorkFlow:
    """
    Attributes:
        dirhash : dictionary : directory locations
        greenlist : list :
        sumdf : DataFrame
        logdf : DataFrame
        ehash : dictionary
    """

    def __init__(self, inp, verbose=False):
        """
        inp : dictionary
              key = JPSS_DIR
              key = VOLCAT_LOGFILES
              key = VOLCAT_DIR
        """
        # set directories
        self.dirhash = self.set_directories(inp)
        # get list of volcanoes to process data for
        self.greenlist = self.get_greenlist()
        
        # initialize attributes
        # DataFrame with info from summary files
        self.sumdf = pd.DataFrame()
        # DataFrame with info from log files
        # the informaton from the sumdf file is merged into
        # the logdf file, so logdf contains all information in sumdf.
        self.logdf = pd.DataFrame()

        # dictionary. key is volcano name. value is Event object
        self.ehash = {}

    def get_all_volcanos(self, hours=24):
        """
        hours : integer
        """
        # get pandas dataframe with summary file.
        self.sumdf = get_summary_file_df(self.dirhash["JPSS_DIR"], hours=hours)
        # make the greenlist all the volcanoes in the summary files
        greenlist = self.sumdf["volcano_name"].unique()
        return greenlist

    def add_to_greenlist(self, vname):
        self.greenlist.append(vname)
        self.greenlist = list(set(self.greenlist))

    def reset_greenlist(self, hours=24):
        """
        hours : integer
        """
        self.greenlist = self.get_all_volcanos(hours=hours)

    def get_inactive(self, hours=24):
        self.sumdf = get_summary_file_df(self.dirhash["JPSS_DIR"], hours=hours)
        allvolc = self.sumdf["volcano_name"].unique()
        inactive = [x for x in self.greenlist if x not in allvolc]
        return inactive

    def get_redlist(self, hours=24):
        self.sumdf = get_summary_file_df(self.dirhash["JPSS_DIR"], hours=hours)
        allvolc = self.sumdf["VOLCANO_NAME"].unique()
        # red = [x for x in self.greenlist if x not in allvolc]
        red = [x for x in allvolc if x not in self.greenlist]
        return red

    def get_sumdf(self, hours=24):
        return get_summary_file_df(self.dirhash["JPSS_DIR"], hours=hours)

    def get_volcat(self, hours=24,  verbose=False):
        """
        Retrieves log files from ftp.
        inp : dictionary with file locations
        hours : integer indicating how many hours back to find files.
        """
        # get pandas dataframe with summary file.
        # the summary files are pushed to us and have information about the log files.
        self.sumdf = get_summary_file_df(self.dirhash["JPSS_DIR"], hours=hours)

        # download log files for volcanoes in greenlist
        khash = {}

        # we pull the log files which we want.
        khash["volcano_name"] = self.greenlist
        self.logdf = get_log_files(
            self.sumdf, inp=self.dirhash,verbose=verbose,**khash
        )
        logger.info('retrieved log files') 

        # download the event netcdf files
        self.ehash = self.create_ehash()
        #indf = self.logdf[self.logdf['status']==True]
        #self.ehash = create_event(indf, self.dirhash, vlist=self.greenlist)
        for volc in self.ehash.keys():
            print("download {} netcdf files".format(volc))
            logfile = '/hysplit-users/alicec/{}_FAILED.txt'.format(fix_volc_name(volc))
            self.ehash[volc].download(self.dirhash,log=logfile,verbose=verbose)


    def create_ehash(self):
        # download netcdf files
        indf = self.logdf[self.logdf['status']==True]
        ehash = create_event(indf, self.dirhash, vlist=self.greenlist)
        return ehash

    def check_log_files(self):
        dfbad = self.logdf[self.logdf['status']==False]
        logger.info('Number of files that could not be retrieved {}'.format(len(dfbad.log_url.values)))
        return dfbad

    def check_fid(self):
        for volc in self.ehash.keys():
            print(volc)
            self.ehash[volc].check_fid()

    def set_directories(self, inp):
        """
        keys should be
        JPSS_DIR - location of xml files that are pushed to our system
        VOLCAT_LOGFILES - location to download the log files
        VOLCAT_DIR - location for netcdf files to be downloaded.
        """
        return inp

    def get_greenlist(self):
        # WorkFlow class
        gfile = os.path.join(self.dirhash["VOLCAT_DIR"], "green_list.txt")
        with open(gfile, "r") as fid:
            greenlist = fid.readlines()
        greenlist = [x.strip() for x in greenlist]
        return greenlist

    def write_parallax_corrected(self, gridspace=None, elist=None):
        # WorkFlow class

        failed_list = []
        # if not isinstance(elist,list):
        #    elist = self.ehash.keys()

        print("WRITING FOR THESE", elist)
        for volc in elist:
            print("writing parallax corrected for {}".format(volc))
            try:
                self.ehash[volc].write_parallax_corrected(gridspace=gridspace)
            except Exception as eee:
                print("FAILED to create")
                print(eee)
                print("----------------")
                failed_list.append(volc)
        print("parallax correction failed for {}".format(", ".join(failed_list)))
        return failed_list




def file_progression():
    import qva_logic as vl

    # Step 0:
    df = vl.get_summary_file_df(fdir, hours=x)
    # creates a dataframe which contains information from the summary files that have
    # been pushed to fdir over the last x hours.

    # Step 1:

    # get event log files.

    # Step 1:
    vl.get_files()
    # In step 1, you can specify the vaac region (vaac='Washington') for files to be
    # downloaded. If the region is not specified, then all files are download
    # Can also specify verbose=True
    # In this step, all event log files are downloaded, event log files are parsed for
    # the vaac region. If the event log file is for the specified region, then the netcdf
    # files in the event log file are downloaded.
    # The netcdf files are filed by volcano in the data_dir (/pub/ECMWF/JPSS/VOLCAT/Files/)
    # Step 1.5:
    vl.make_red_list(data_dir)
    # In this step, the red_list.txt file can be updated. The green_list.txt and yellow_list.txt
    # files are generated by hand at the moment. Volcanoes not on either list are written
    # to red_list.txt. This way we do not process non-active volcano files.
    # It is not necessary to run every time, and should be monitored to ensure we are not
    # overlooking any eruptions that may be occuring but havent been active for quite
    # some time.
    # Step 2:
    data_dir = "/pub/ECMWF/JPSS/VOLCAT/Files/"
    vl.make_pc_files(data_dir, volcano="volcano name", vlist_file="green_list.txt")
    # In step 2, you MUST specify the data directory, which is the parent directory
    # for all the volcanoes. You can also specify verbose=True if desired
    # You can also specify the volcano if desired - creates pc files for that volcano only
    # You can specify a file which contains a list of volcanoes to process
    # If you do not want to specify either volcano or vlist_file, files for all volcanoes will
    # be processed.
    # In this step, a parallax_corrected folder is created in each volcano directory.
    # Then parallax corrected files are written for all netcdf files in within the volcano
    # directory.
    # Step 3:
    events = vl.list_times(data_dir, volcano="volcano name", pc=True)
    vl.multi_plots(
        data_dir,
        volcano="volcano name",
        eventdf=events.Event_Dates[i],
        pc=True,
        saveas=True,
    )
    # In this step, figures can be generated for each eruption showing the timeseries
    # of area, maximum height, total mass, etc.
    # This step requires first creating a list of the individual eruption events, then this list
    # is used to determine the individual eruption events for each volcano.
    # This is done using vl.list_times(). The figures are generated for the individual events
    # Step 4:
    events = vl.list_times(data_dir, volcano="volcano name", pc=True)
    vl.make_data_insertion(
        data_dir,
        volcano="volcano name",
        event_date=events.Event_Dates[i],
        pc=True,
        verbose=False,
    )
    # In step 4, Emitimes files are generated for the available data. If the volcano is
    # specified, then files are generated for only that volcano. If the eruption event time
    # and volcano are specified, then the files are generated only for the specified
    # eruption event.
    # When you have the events list, you can loop through the list (events.Event_Dates[i])
    # for the value in event_date.
    # I have been running into some errors with files that have values of 0:
    # ValueError: zero-size array to reduction operation minimum which has no identity
    # This error comes from the make_1D function in write_emitimes.py. I don't have
    # a fix for it yet, but if you specify the event_date, you can move past the files that
    # are causing a problem.
    # Step 5:
    vl.setup_runs()  # AMC adapt the ashapp functions to do this.
    #  read the emit-times files
    #  might trigger on the emit-times file runs
    # In this step, control and setup files are generated for data insertion runs.
    # IN PROGRESS


def generate_report(vmin=None, vmax=None, greenlist=None, **kwargs):
    import matplotlib.pyplot as plt

    # get_files()
    if "VOLCAT_DIR" in kwargs.keys():
        data_dir = kwargs["VOLCAT_DIR"]
    else:
        data_dir = "/pub/ECMWF/JPSS/VOLCAT/Files/"

    vnames = os.listdir(data_dir)
    redlist = []
    #print(vnames)
    #vnames = vnames[vmin:vmax]
    print(len(vnames))
    for iii, volc in enumerate(vnames):
        if isinstance(greenlist,list):
           if not volc in greenlist: 
              redlist.append(volc)
              continue
        fig = plt.figure(figsize=[10, 2])
        try:
            events = list_times(os.path.join(data_dir, volc), pc=False)
        except Exception as eee:
            print(eee)
            print("warning in generate report for directory ", volc)
            continue
        if 'Event_Dates' in events.keys():
            sns.set()
            plt.plot(events["Event_Dates"], events["Num_Files"], "k.")
            plt.title(volc)
            ax = plt.gca()
            ax.set_ylabel('Number of files for event')
            fig.autofmt_xdate()
            plt.tight_layout()
            plt.show()
        else:
           print(volc, events.keys())
    return redlist 

def create_event(sumdf, fdir, vlist=None, verbose=False):
    """
    sumdf : DataFrame with information from summary files
    fdir  : location of log files
    vlist : list of volcanoes to create  Events objects for
    verbose : boolean

    output

    eventhash : dictionary
                key is volcano name. value is Events object.
    """
    if isinstance(fdir, dict):
        fdir = fdir["VOLCAT_LOGFILES"]

    eventhash = {}  # key is volcano name. value is Event class instance
    if not isinstance(vlist, (list, np.ndarray)):
        vlist = sumdf["volcano_name"].unique()
    else:
        # so vlist doesn't have to be case dependent.
        templist = sumdf['volcano_name'].unique()
        vlist2 = []
        for vvv in vlist:
            temp = [x for x in templist if vvv.lower() in x.lower()]     
            vlist2.extend(temp)
        vlist = vlist2
    for volc in vlist:
        print(volc)
        eventlist = []
        df2 = sumdf[sumdf["volcano_name"]== volc]
        efile = df2["log"].unique()
        for efl in efile:
            evo = EventFile(efl, fdir)
            #evo.open()
            try:
                evo.open()
            except Exception as eee:
                print('cannot open {} {}'.format(fdir, efl))
                print(eee)
                print('\n')
                continue
            #print('Opened {}'.format(efl))
            eventlist.append(evo)
        if eventlist:
            eve = Events()
            eve.add_events(eventlist)
            eventhash[volc] = eve
    return eventhash

def correct_pc(data_dir, newdir="pc_corrected", daterange=None, verbose=False,gridspace=0.1,dfile_list=None):
    """Create pc_corrected folder if not already there.
    Create pc_corrected netcdf file in pc_corrected folder if not already there
    """
    import sys
    # May want to streamline this more so all files are not checked each time!
    # Create pc_corrected netcdf files if not already created, put in pc_corrected folder
    # Make sure data_dir ends with '/'
    # data_dir = os.path.join(data_dir, "")
    # Create pc_corrected folder if not already there

    #edir = '/hysplit-users/alicec/utilhysplit/utilvolc/'
    processhandler = ProcessList()
    helper = Helper
    #processhandler.pipe_stdout()
    #processhandler.pipe_stderr()
    print(os.getcwd())
    #sys.exit()
    make_dir(data_dir, verbose=verbose)
    pc_dir = os.path.join(data_dir, newdir, "")
    # Create list of files original directory
    if not isinstance(dfile_list,list):
        dfile_list = volcat.find_volcat(
            data_dir, vid=None, daterange=daterange, include_last=True, return_val=3
        )
    # dfile_list = sorted(glob(data_dir+'*.nc'))
    # Create hypothetical list of pc corrected files

    file_list = []
    pcfile_list = []
    jjj=0
    #print("volcat files", dfile_list)
    logger.info('Number of files to write {}'.format(len(dfile_list)))
    for iii, element in enumerate(dfile_list):
        if iii%50 == 0:
           logger.info('Working on file {}'.format(iii))
        #print("element ", element)
        #print("--------------")
        s = element.rfind("/")
        fname = element[s + 1 :]
        pcfname = os.path.splitext(fname)[0] + "_pc.nc"
        make_pcfile = check_file(pcfname, pc_dir, verbose=verbose)
        if make_pcfile:
            # Create pc_corrected file if not in pc directory
            flist = [fname]
            volcat.write_parallax_corrected_files(data_dir,pc_dir,flist=flist,gridspace=gridspace,verbose=False)
            if verbose:
                print(os.path.join(data_dir, fname))
            try:
                volcat.write_parallax_corrected_files(data_dir,pc_dir,flist=flist,gridspace=gridspace,verbose=False)
            except Exception as eee:
                print('warning could not write pc corrected file {}'.format(flist[0]))
                print(eee)
                print('----') 
            #cproc = ['python', 'process_volcat.py', data_dir, fname, pc_dir, str(gridspace)]
            #print(cproc)
            #processhandler.startnew(cproc,edir)
            #Helper.execute_with_shell(cproc)
            #make_pcfile = check_file(pcfname, pc_dir, verbose=verbose)
            #print('file written', not(make_pcfile), pcfname, pc_dir)
            #done=False
            #seconds_to_wait=2
            #total_time = 0
            #max_time = 60*6000
            #while not done:
            #      num_procs = processhandler.checkprocs()
            #      logger.info('{}'.format(processhandler.err))
            #      if num_procs==0: done=True
            #      time.sleep(seconds_to_wait)
            #      total_time += seconds_to_wait
            #jjj+=1
            #if jjj>2: sys.exit()
    return None


# def list_dirs(data_dir):
#    """ Lists subdirectories within give directory
#    Inputs:
#    data_dir: directory path of parent directory (string)
#    Outputs:
#    dirlist: list of subdirectories within data_dir
#    """
#    # scan directory works with python 3.5 and later.
#    dirlist = os.scandir(data_dir)
#    newlist = [volc.path for volc in dirlist if volc.is_dir()]
#    return sorted(newlist)


def make_red_list(data_dir):
    """Makes list of volcanoes not expected to be active - left over from green and
    yellow list. If on the red list, parallax files will not be written, emitimes files will not
    be generated, further processing will not occur.
    Could change this to not even download netcdf files, but this will take some further
    logic generation.
    Inputs:
    data_dir: directory path of parent directory (string)
    Outputs:
    Text file listing volcanoes that are likely false positives for VOLCAT
    """
    dirlist = list_dirs(data_dir)
    redold = []
    rednew = []
    red = []
    yellow = []
    green = []
    # Read green list
    gname = os.path.join(data_dir, "green_list.txt")
    if os.path.isfile(gname):
        with open(gname, "r") as fid:
            green = fid.readlines()
            green = [line.rstrip() for line in green]
    else:
        logger.warning("green list file not found {}".format(gname))
    # Read yellow list
    yname = os.path.join(data_dir, "yellow_list.txt")
    if os.path.isfile(yname):
        with open(yname, "r") as fid:
            yellow = fid.readlines()
            yellow = [line.rstrip() for line in yellow]
    else:
        logger.warning("yellow list file not found {}".format(yname))
    # Compare dirlist to green and yellow lists
    # Read yellow list
    rname = os.path.join(data_dir, "red_list.txt")
    if os.path.isfile(rname):
        with open(rname, "r") as fid:
            redold = fid.readlines()
            redold = [line.rstrip() for line in redold]
    for volc in dirlist:
        if volc not in green:
            if volc not in yellow:
                red.append(volc)
                if volc not in redold:
                    rednew.append(volc)

    logger.warning("New volcanoes added to red list: {}".format(", ".join(rednew)))
    # Write red_list.txt file - will overwrite previous file
    with open(os.path.join(data_dir, "red_list.txt"), "w") as fid:
        for volc in red:
            fid.write(volc + "\n")
    os.chmod(data_dir + "red_list.txt", 0o666)
    return None



def make_pc_files(
    data_dir, volcano=None, vlist_file=None, daterange=None, verbose=False
):
    """Makes corrected pc files.
    Might want to streamline the check process at some point. Not necessary now
    Inputs:
    data_dir: parent directory for volcanoes (string)
    volcano: name of specific volcano (string) None by default
    vlist_file: name of file with list of volcanoes (string) None by default
            Written fo use with green_list.txt, yellow_list.txt
    If volcano and vlist_file are both None, pc_corrected files are written
    for all availabe files.
    verbose: boolean
    Outputs:
    Parallax files are generated for specified volcano, or vlist_file, or all volcanoes
    Depending on inputs.
    """
    # Make list of available directories
    dirlist = list_dirs(data_dir)
    if volcano != None:
        if volcano in dirlist:
            dirlist = [volcano]
            # file_dir = os.path.join(data_dir, volcano, '')
            # correct_pc(file_dir, verbose=verbose)
        # if verbose:
        #    print('Parallax corrected files available in '+volcano+' directory')
    elif vlist_file != None:
        with open(vlist_file) as vfile:
            volclist = vfile.readlines()
            volclist = [line.rstrip() for line in volclist]
            file.close()
        for volcano in volclist:
            newlist = []
            if volcano in dirlist:
                newlist.append(volcano)
            dirlist = newlist
            # file_dir = os.path.join(data_dir, volcano, '')
            # correct_pc(file_dir, verbose=verbose)
            # if verbose:
            #    print('Parallax corrected files available in '+volcano+' directory!')
    # else:
    for direct in dirlist:
        file_dir = os.path.join(data_dir, direct, "")
    if verbose:
        print(
            "Parallax corrected files available in these directories: " + str(dirlist)
        )
    return None


def volcplots(das_list, img_dir, saveas=True):
    """Makes time series plots of total mass, total area, MER, max height.
    Inputs:
    dfile_list: list of volcat xarray (list)
    img_dir: filepath of image directory (string)
    saveas: (boolean) default=True
    Outputs:
    Shows 4-panel figure
    Saves figure image filepath if saveas=True
    """
    from utilvolc import volcat_plots as vp
    import matplotlib.pyplot as plt

    # Initalize VolcatPlots class
    vplot = vp.VolcatPlots(das_list)
    vplot.make_arrays()
    vplot.set_plot_settings()
    fig1 = vplot.plot_multiA(fignum=1, smooth=0.08, yscale="linear")
    fig1.autofmt_xdate()
    volcname = das_list[0].attrs["volcano_name"]
    volcano = fix_volc_name(volcname)
    dset_name = das_list[0].attrs["dataset_name"]
    s = dset_name.find("b")
    e = dset_name.rfind("_")
    begin_time = dset_name[s + 1 : e]
    # TO DO: Use volcat.get_volcat_name_df to get begin_time value
    if saveas:
        if pc:
            figname = (
                volcano + "_" + begin_time + "_mass_area_kgs_maxhgt_pc_corrected.png"
            )
        else:
            figname = volcano + "_" + begin_time + "mass_area_kgs_maxhgt.png"
        fig1.savefig(img_dir + figname)
        plt.close()
        return print("Figure saved: " + img_dir + figname)
    else:
        return fig1.show()


def check_volcano(data_dir, volcano=None, verbose=False):
    """Checks to see if specified volcano is in data directory
    Inputs:
    data_dir: data_directory (string)
    volcano: volcano name (string)
    verbose: (boolean)
    Outputs:
    result: (boolean) True if volcano is in data directory, False if not
    """
    # List directories in data_dir
    dirlist = list_dirs(data_dir)
    if volcano is not None:
        if volcano in dirlist:
            return True
        else:
            if verbose:
                return print(volcano + " not in " + str(dirlist))
            else:
                return False
    else:
        if verbose:
            return print("Need to specify volcano name.")
        else:
            return False


def list_times(data_dir, volcano=None, pc=True):
    """Lists all available volcanic beginning event times in given data directory
    Provides number of files attributed to the given beginning event time.
    Used to determine which time to create images.
    Inputs:
    data_dir: data directory (string)
    volcano:  volcano name (string)
    pc: using pc files - default True (boolean)
        assumes pc corrected file in subdirectory called pc_corrected.
    Outputs:
    events: pandas dataframe of available times, number files for each time

    Example input:
       list_times(inp['VOLCAT_DIR'] + 'Popocatepetl/')


    Example output:

        Event_Dates Num_Files
    0   2022-03-14 15:30:30 1
    1   2022-03-15 09:40:30 1
    2   2022-03-15 19:50:30 33
    3   2022-03-16 01:30:29 1
    4   2022-03-16 10:50:30 1
    ... ... ...
    1168    2023-03-28 15:00:30 1
    1169    2022-10-05 11:30:29 2
    1170    2022-10-12 18:20:30 8
    1171    2023-03-25 00:01:30 2
    1172    2022-05-21 22:50:30 5


    """
    import pandas as pd

    volc_check = check_volcano(data_dir, volcano=volcano)
    if volc_check:
        volc_dir = os.path.join(data_dir, volcano, "")
    else:
        volc_dir = data_dir
    if pc:
        volc_dir = os.path.join(volc_dir, "pc_corrected", "")
    # Creating dataframe of filename information
    dataf = volcat.get_volcat_name_df(volc_dir, include_last=True)
    if dataf.empty:
        return dataf
    event_dates = dataf["edate"].unique()
    eventd = pd.DataFrame(event_dates, columns=["Event_Dates"])
    lens = []
    g = 0
    while g < len(event_dates):
        files = dataf.loc[dataf["edate"] == event_dates[g], "filename"]
        lens.append(len(files))
        g += 1
    lensd = pd.DataFrame(lens, columns=["Num_Files"])
    events = pd.concat([eventd, lensd], axis=1)
    return events


def make_volcat_plots(
    data_dir, volcano=None, event_date=None, pc=True, saveas=True, verbose=False
):
    """Calls functions to create plots of volcat files within designated data directory.
    To add: Make flag for calling different plotting funtions with this function?
    Inputs:
    data_dir: path for data directory (string)
    volcano: name of specific volcano (string)
         If None, function goes through all available volcano subdirectories
    event_date: date/time of volcanic eruption (datetime object,datetime64, or timestamp)
    pc: (boolean) default=True - use parallax corrected files
    saveas: (boolean) default=True
    verbose: (boolean) default=False
    Outputs:
    Figures generated in image directory
    """
    from datetime import datetime

    # List directories in data_dir
    dirlist = list_dirs(data_dir)
    datadirs = []
    # Check to see if given volcano is within list of directories (if not None)
    # Generate list of volcano directories
    if volcano:
        if volcano in dirlist:
            datadirs.append(os.path.join(data_dir, volcano, ""))
        else:
            return print(volcano + " not in " + str(dirlist))
    else:
        for volcano in dirlist:
            datadirs.append(os.path.join(data_dir, volcano, ""))
    # Create image directory within volcano directory if it doesnt exist
    img_dirs = []
    for directory in datadirs:
        newdir = "Images"
        make_dir(directory, newdir=newdir, verbose=verbose)
        image_dir = os.path.join(directory, newdir, "")
        img_dirs.append(image_dir)
        # Generate list of files
        if pc:
            # Check if pc_corrected directory exists
            pcdir = "pc_corrected"
            if not os.path.exists(directory + pcdir):
                return print(
                    "pc_corrected directory does not exist! Make " + directory + pcdir
                )
            else:
                volc_dir = directory + pcdir
        else:
            # Using non-parallax corrected files
            volc_dir = directory
        if event_date:
            # Files only with specific event date
            das_list = volcat.get_volcat_list(volc_dir, fdate=event_date)
        else:
            # All files in directory
            das_list = volcat.get_volcat_list(volc_dir)
        # Generate plots
        volcplots(das_list, image_dir, pc=pc, saveas=saveas)
    return print("Figures generated in " + str(img_dirs))


def multi_plots(
    data_dir, volcano=None, eventdf=None, pc=True, saveas=True, verbose=False
):
    """Creates multiple plots in a volcano directory based on event
    dataframe from vl.list_times().
    Inputs:
    data_dir: data parent directory (string)
    volcano: volcano name (string)
    eventdf: event dataframe (pandas dataframe)
    pc: use parallax corrected files? (boolean)
    saveas: save figures? (boolean)
    verbose: (boolean)
    """
    if eventdf is not None:
        a = 0
        while a < len(eventdf):
            event = eventdf["Event_Dates"][a]
            numfiles = eventdf["Num_Files"][a]
            if numfiles > 4:
                make_volcat_plots(
                    data_dir,
                    volcano=volcano,
                    event_date=event,
                    pc=pc,
                    saveas=saveas,
                    verbose=verbose,
                )
            a += 1
        if verbose:
            return "Image files created for events with more than 4 observation files"
        else:
            return None
    else:
        return "Missing event dataframe: use volcat_logic.list_times()"


def update_vaac_files():
    """Writes three files based on data from Washington vaac webpage.
    This is tabled for now.
    Lists can be generated by hand from the Washington VAAC webpage."""
    from lxml import html
    from bs4 import BeautifulSoup
    import requests

    page = requests.get("https://www.ssd.noaa.gov/VAAC/ARCH21/archive.html")
    tree = html.fromstring(page.content)
    volcanoes = tree.xpath('//span[@id="quicklinks"]/text()')


def write_di_emitimes(
    data_dir,
    volcano=None,
    event_date=None,
    pc=True,
    verbose=True,
    clip=True,
    overwrite=False,
):
    """
    Writes emitimes file for data insertion from volcat netcdf file(s).
    Inputs:
    data_dir: path for data directory (string)
    volcano: name of specific volcano (string)
    event_date: event date (datetime object,datetime64, or timestamp)
    pc: (boolean) default=True - use parallax corrected files
    verbose: (boolean)
    Return:
    None
    """
    from utilvolc import make_data_insertion as mdi

    pollnum = 1
    pollpercents = [1]
    layer = 0.0

    # List directories in data_dir
    dirlist = list_dirs(data_dir)
    datadirs = []
    # Check to see if volcano is in data_directory
    if volcano != None:
        volc_check = check_volcano(data_dir, volcano=volcano, verbose=verbose)
        if volc_check:
            datadirs.append(os.path.join(data_dir, volcano, ""))
    else:
        if verbose:
            return "Must specify volcano"
        else:
            return None
            # Create emitimes directory within volcano directory if it doesnt exist
    emit_dirs = []
    for directory in datadirs:
        newdir = "emitimes"
        make_dir(directory, newdir=newdir, verbose=verbose)
        emit_dir = os.path.join(directory, newdir, "")
        emit_dirs.append(emit_dir)
        # Generate list of files
        if pc:
            # Check if pc_corrected directory exists
            pcdir = "pc_corrected"
            if not os.path.exists(directory + pcdir):
                return print("pc_corrected directory does not exist!")
            else:
                volc_dir = directory + pcdir
        else:
            # Using non-parallax corrected files
            volc_dir = directory
        # Dataframe of files in directory
        dframe = volcat.get_volcat_name_df(volc_dir)
        # sort by image date
        # dframe = dframe.sort_values(by="idate")
        # print(dframe)
        # pick by event date if indicated
        if event_date:
            dframe = dframe[dframe["edate"] == event_date]
        else:
            print('available event dates', dframe["edate"].unique())
            #return

        imgdates = dframe["idate"].unique()
        if verbose: print('available image dates', imgdates)
        # loop through image dates
        # write one emit-times file per image date.
        # sometimes this will utilize more than one volcat file.
        for iii, idate in enumerate(imgdates):
            # get all files with that image date
            filenames = dframe.loc[dframe["idate"] == idate, "filename"].tolist()
            #print('working on ', filenames)
            date_time = pd.to_datetime(idate)
            # Initialize write_emitimes function
            volcemit = mdi.InsertVolcat(
                emit_dir,
                volc_dir,
                date_time,
                fname=filenames,
                layer=layer,
                pollnum=pollnum,
                pollpercents=pollpercents,
            )
            #print('FILENAME', filenames)
            #overwrite=True
            if not volcemit.check_for_file() or overwrite:
                try:
                    oname = volcemit.write_emit(area_file=False, clip=clip, verbose=verbose)
                except:
                    logger.warning('emit file could not be written')
                    continue
                logger.info(
                    "Emit-times file written {}".format(volcemit.make_emit_filename())
                )
                print('CREATED', oname)
            else:
                logger.info(
                    "Emit-times file exists {}".format(volcemit.make_emit_filename())
                )
        return None



def forecast2QVA(dset,qva_filename,kwargs):
    """
    dset should have units of mg/m3.
    """
    enslist = dset.ens.values
    qva = ensemble_tools.ATLra(dset,enslist,sourcelist=None, threshlist=[0.2,2,5,10],weights=None,**kwargs)
    encoding = {}
    encoding['zlib']=True
    encoding['complevel']=9
    ehash={}
    ehash['Concentration'] = encoding
    ehash['FrequencyOfExceedance'] = encoding
    qva.to_netcdf(qva_filename,encoding=ehash) 
    return qva


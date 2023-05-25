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
from utilvolc.runhelper import Helper
from utilvolc.runhelper import list_dirs
from utilvolc.runhelper import make_dir
from utilvolc import volcat
from utilhysplit.runhandler import ProcessList
from utilhysplit.plotutils import map_util
import utilhysplit.evaluation.web_ensemble_plots as wep
from utilvolc import make_data_insertion as mdi
from utilhysplit.evaluation import ensemble_tools

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


def get_summary_file_df(fdir, verbose=False, hours=48, edate = datetime.datetime.now()):
    """
    fdir : str : location of summary files
    Returns
    sumdf : pandas dataframe with information from all the summary files.

    prints an error if one of the files cannot be read
    """

    # 2023 09 May (amc) use the time stamp in the filename rather than scandir

    vlist = []
    if not hours:
        logger.warning("get_summary_file_df input hours None, setting hours to 100")
        hours = 100
    strf = "s%Y%j_%H"
    sdate = edate - datetime.timedelta(hours=hours)
    s_str = sdate.strftime(strf)
    done=False
    flist = []
    while not done:
          flist.extend(glob.glob(fdir + '/*{}*.json'.format(s_str)))
          sdate += datetime.timedelta(hours=1)
          s_str = sdate.strftime(strf)
          if sdate > edate:  done=True
    for fln in flist:
       try:
           sfn = SummaryFile(fln, fdir)
       except:
           if verbose:
             print("Not summary file", fln)
             continue
       try:
           vlist.append(sfn.open_dataframe())
       except Exception as eee:
           print('ERROR READING', fln,  eee)
    sumdf = pd.concat(vlist)
    return sumdf


class VFile:
    def __init__(self, fname, fdir="./"):
        self.fname = fname
        self.fdir = fdir
        self.dtfmt = "s%Y%j_%H%M%S"
        self.vhash = {"FileName": fname}
        #self.dii = -2
        #self.djj = -1
        temp = self.parse(fname)
        self.exists = os.path.isfile(os.path.join(self.fdir, self.fname))

    def parse(self, fname):
        """
        Parses information from filename
        """
        fname = fname.replace("Full_Disk", "FullDisk")
        fname = fname.replace("FULL_DISK", "FullDisk")
        temp = fname.split("_")
        dstr = "{}_{}".format(temp[self.dii], temp[self.djj].replace(".json", ""))
        try:
            self.date = datetime.datetime.strptime(dstr, self.dtfmt)
        except:
            print("date not in correct format {} {}".format(dstr, fname))
            print(temp)
            self.date = datetime.datetime.now()
        self.vhash["date"] = self.date
        return temp


class SummaryFile(VFile):
    def __init__(self, fname, fdir="./"):
        self.dii = -2
        self.djj = -1
        super().__init__(fname, fdir)

    def open_dataframe(self):
        """
        Returns dataframe with information from
        1. summary file filename
        2. summary file (json format) contents
        2. event log file filename(s) event log filenames and URLS
           are contained in the summary file

        1. summary date : date in filename of summary file
        1. summary file : filename of summary file
        1. satellite : indicated in filename of summary file
        VOLCANO_NAME
        VOLCANO_GVP_ID
        VOLCANO_LAT
        VOLCANO_LON
        VOLCANO_COUNTRY
        VOLCANO_REGION
        VOLCANO_SUBREGION
        VAAC_REGION
        VOLCANO_TIMEFRAME
        VOLCANO_ELEVATION
        EVENT_TYPE
        LOG_URL

        """
        jsonf = open_json(os.path.join(self.fdir, self.fname))
        dataf = jsonf["VOLCANOES"]
        datra = []  # list of dictionaries
        if isinstance(dataf, dict):
            dataf = [dataf]
        if isinstance(dataf, list):
            for volc in dataf:
                if "EVENTS" in volc.keys():
                    # if there is only one event, then it will be dictionary
                    aaa = volc["EVENTS"]
                    if isinstance(aaa, dict):
                        volc.update(aaa)
                        volc.pop("EVENTS")
                        volc.update(self.add_event_file_info(aaa["LOG"]))
                        volc.update(self.add_file_name_atts())
                        datra.append(volc)
                    # if there is more than one event, it will be a list
                    elif isinstance(aaa, list):
                        for event in aaa:
                            volc2 = volc.copy()
                            volc2.update(event)
                            volc2.pop("EVENTS")
                            volc2.update(self.add_event_file_info(event["LOG"]))
                            volc2.update(self.add_file_name_atts())
                            datra.append(volc2)
            data = pd.DataFrame.from_dict(datra)
        # if isinstance(dataf,dict):
        #    data = pd.DataFrame.from_dict([dataf])
        # elif isinstance(dataf,str):
        #    data = dataf
        columns = data.columns
        newc = [x.lower() for x in columns]
        data.columns = newc

        return data

    def add_file_name_atts(self):
        aaa = {"summary date": self.date}
        aaa["summary file"] = self.fname
        aaa["satellite"] = self.sat
        return aaa

    def add_event_file_info(self, ename):
        efile = EventFile(ename)
        ehash = {}
        ehash["event date"] = efile.date
        ehash["event instrument"] = efile.instrument
        ehash["event vid"] = efile.vid
        ehash["event gid"] = efile.gid
        return ehash

    def parse(self, fname):
        """
        Parses information from filename
        """
        temp = super().parse(fname)
        self.sat = temp[1]

    def open(self):
        fname = os.path.join(self.fdir, self.fname)
        df = open_dataframe(fname, varname="VOLCANOES")
        self.df = df


def flist2eventdf(flist,inphash):
    """
    create a dataframe from list of volcat filenames
    inphash contains additional information should be in following format
    VOLCANO_NAME : str
    """
    # example usage in Raikoke2023Volcat
    vlist = []
    for fle in flist:
        try:
            temp = volcat.VolcatName(fle)
        except:
            continue
        vlist.append(temp.vhash)
    vframe = pd.DataFrame.from_dict(vlist)
    cols = vframe.columns
    # rename columns which need to be renamed.
    cols = [x if x != 'satellite platform' else 'sensor_name' for x in cols]
    cols = [x if x != 'idate' else 'observation_date' for x in cols]
    cols = [x if x != 'filename' else 'event_file' for x in cols]
    cols = [x if x != 'WMO satellite id' else 'SENSOR_WMO_ID_INITIAL' for x in cols]
    cols = [x if x != 'feature id' else 'FEATURE_ID' for x in cols]
    vframe.columns = cols
    checklist = ['VOLCANO_NAME', 'VOLCANO_LAT', 'VOLCANO_LON']
    for key in checklist:
        if key in inphash.keys():
           vframe[key.lower()] = inphash[key] 
        if key.lower() in inphash.keys():
           vframe[key.lower()] = inphash[key] 

    return vframe

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


class EventFile(VFile):
    """
    methods
         __init__
         parse
         open

    attributes

    df : pandas DataFrame
         colums are
         sensor_name
         SENSOR_MODE
         SENSOR_WMO_ID
         START_COVERAGE_TIME
         OBSERVATION_TIME
         OBSERVATION_EPOCH
         FEATURE_ID
         event_file
         EVENT_URL

    attrs : dicionary
         keys are
         PROCESSING_SYSTEM
         VOLCANO_NAME
         VOLCANO_GVP_ID
         VOLCANO_LAT
         VOLCANO_LON
         VAAC_REGION
         EVENT_TYPE
         EVENT_STATUS
         FEATURE_ID_INITIAL
         FEATURE_ID_LATEST
    """

    def __init__(self, fname, fdir="./"):
        self.dii = -5
        self.djj = -4
        self.df = pd.DataFrame()
        self.attrs = {}
        super().__init__(fname, fdir)

    def parse(self, fname):
        """
        Parses information from filename
        """
        temp = super().parse(fname)
        self.vid = temp[5]
        self.instrument = temp[1]
        self.gid = temp[-1].replace(".json", "")

    def open(self):
        fname = os.path.join(self.fdir, self.fname)
        if not self.exists:
            logger.warning("file does not exist {}".format(fname))
            return False
        jsonf = open_json(os.path.join(self.fdir, self.fname))

        # some files have 'files' and some have 'FILES'
        keys = jsonf.keys()
        fstr = [x for x in keys if 'files' in x.lower()]
        if fstr:
           fdict = jsonf[fstr[0]]
        else:
           logger.warning('correct key FILES or files  not found.')
           return False

        if isinstance(fdict, dict):
            fdict = [fdict]
        df = pd.DataFrame(fdict)

        # make all the columns lower case
        newc = [x.lower() for x in df.columns]
        df.columns = newc
        if 'observation_time' not in df.columns:
           logger.warning('observation_time not found')
           return False

        # format of the data is different in different files.
        dtfmt = "%Y-%m-%dT%H:%M:%SZ"
        dtfmt2 = "%Y-%m-%dT%H:%M:%S.0Z"
        ttest = df['observation_time'].values[0]
        if ttest[-3:] == '.0Z': dtfmt = dtfmt2
 
        dftemp2 = df.apply(
            lambda row: datetime.datetime.strptime(row["observation_time"], dtfmt),
            axis=1,
        )
        df.insert(1, "observation_date", dftemp2)
        self.df = df

        jsonf.pop(fstr[0])
        self.attrs = jsonf
        return True

    def check_feature_id(self):
        dtemp = self.df
        sns.set()
        for fid in dtemp["feature_id"].unique():
            dtemp2 = dtemp[dtemp["feature_id"] == fid]
            plt.plot(dtemp2["observation_date"], dtemp2["feature_id"], ".")
        fig = plt.gcf()
        fig.autofmt_xdate()
        ax = plt.gca()
        ax.set_ylabel("feature id")


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

        # initialize attributes.

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
        indf = self.logdf[self.logdf['status']==True]
        self.ehash = create_event(indf, self.dirhash, vlist=self.greenlist)
        for volc in self.ehash.keys():
            print("download {} netcdf files".format(volc))
            logfile = '/hysplit-users/alicec/{}_FAILED.txt'.format(volc)
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


def generate_report(vmin=None, vmax=None, **kwargs):
    import matplotlib.pyplot as plt

    # get_files()
    if "VOLCAT_DIR" in kwargs.keys():
        data_dir = kwargs["VOLCAT_DIR"]
    else:
        data_dir = "/pub/ECMWF/JPSS/VOLCAT/Files/"

    vnames = os.listdir(data_dir)
    print(vnames)
    vnames = vnames[vmin:vmax]
    print(len(vnames))
    for iii, volc in enumerate(vnames):
        fig = plt.figure(figsize=[10, 2])
        try:
            events = list_times(os.path.join(data_dir, volc))
        except Exception as eee:
            print(eee)
            print("warning in generate report for directory ", volc)
            continue
        if 'Event_Dates' in events.keys():
            plt.plot(events["Event_Dates"], events["Num_Files"], "ko")
            plt.title(volc)
            fig.autofmt_xdate()
            plt.show()
        else:
           print(volc, events.keys())

def make_sumdf(flist):
    return -1


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


def get_log_files(sumdf, verbose=False, inp={"VOLCAT_LOGFILES": "./"},**kwargs):
    """
    INPUTS
    sumdf  : pandas DataFrame : 
    verbose : boolean
    inp    : dictionary. VOLCAT_LOGFILES 
    kwargs : dictionary where key is a column name and value is list of values to keep for that column.

    Return:
       pandas dataframe with additional status column indicating whether log file was downloaded

    """
    df2 = sumdf.copy()
    ccc = sumdf.columns
    rval = []
    for kwargkey in kwargs.keys():
        if kwargkey in ccc:
            if isinstance(kwargs[kwargkey], (list, np.ndarray)):
                df2 = df2[df2[kwargkey].isin(kwargs[kwargkey])]
            else:
                df2 = df2[df2[kwargkey] == kwargs[kwargkey]]
    eurl = df2["log_url"].unique()
    if 'VOLCAT_LOGFILES' not in inp.keys():
        logger.warning('directory for volcat logfiles not inpt')
    for iii, logfile in enumerate(eurl):
        rv = get_log(logfile, verbose=verbose, VOLCAT_LOGFILES=inp["VOLCAT_LOGFILES"])
        #if iii>10: break
        rval.extend(rv)
    failed = []
    succeeded = []
    for rv in rval:
        temp = df2[df2['log_url']==logfile]
        vname = temp['volcano_name'].unique()
        if not rv[1]:
           failed.append((logfile, vname))
        else:
           succeeded.append((logfile, vname))
    for fail in failed:
        print('FAILED\n')
        print('   {} {}'.format(fail[0], fail[1]))
    else:
        if verbose:
           print('\nSUCCEEDED\n')
           for success in succeeded:
              print('   {} {}'.format(success[0], success[1]))
    temp = pd.DataFrame.from_records(rval,columns=['log','status'])
    df2 = pd.merge(df2,temp,on='log')
    #return temp
    return df2


def get_files(inp={"JPSS_DIR": "/pub/jpsss_upload"}, vaac=None, verbose=False):
    """
    Use various functions to get all available netcdf files from json event log files

    Sorts through the ftp directory to find json event summary files

    Loops through them, and downloads the corresponding json event log files

    Keeps track of the json event summary files that have already been downloaded

    Lists event log files, finds the correspondind event file urls for download

    Downloads the event netcdf files

    Keeps track of the downloaded netcdf files

    Parameters
    ----------
    vaac : string
          vaac region can specify region. If None, pulls all available files
    verbose : boolean
    inp : dictionary
         JPSS_DIR :        string
         VOLCAT_DIR :      string
         VOLCAT_LOGFILES : string
    """
    if "JPSS_DIR" in inp.keys():
        jdir = inp["JPSS_DIR"]
    else:
        jdir = "/pub/jpsss_upload/"
    if "VOLCAT_DIR" in inp.keys():
        ddir = inp["VOLCAT_DIR"]
    else:
        ddir = "/pub/ECMWF/JPSS/VOLCAT/Files/"
    if "VOLCAT_LOGFILES" in inp.keys():
        logdir = inp["VOLCAT_LOGFILES"]
    else:
        # this is location of the event log files
        # which are pulled according to the summary files which are in the jpsss upload.
        # these contain the netcdf file names and tell us which netcdfs to upload.
        # json_log.txt is just a list of all
        logdir = "/pub/ECMWF/JPSS/VOLCAT/LogFiles/"

    # Delete files from jpsss_uploads folder that is older than 7 days
    # Files are only available for 7 days on the ftp
    # Dont have permissions to delete files from jpsss_upload!
    # delete_old(jdir, verbose=verbose)

    # Finds json files added to ftp folder
    status = check_dirs(jdir, logdir, verbose=False)
    if np.all(status):
        added = new_json(jdir, logdir)
        added = sorted(added)
        i = 0
        for afiles in added:
            # while i < len(added):
            data = open_dataframe(os.path.join(jdir, afiles), varname="VOLCANOES")
            log_url = get_log_list(data)
            # Downloads json event log files
            get_log(log_url, verbose=verbose, VOLCAT_LOGFILES=logdir)
        # Logs event summary json files
        # writes names of files in jpsss directory which we have already read and pulled the
        # relevant files from.
        record_change(
            ddir=jdir,
            logdir=logdir,
            logfile="json_log.txt",
            suffix=".json",
            verbose=verbose,
        )

    status = check_dirs(logdir, ddir, verbose=True)
    # Delete files from json event log folder that are older than 7 days
    # Netcdf files are only available for 7 days on the ftp
    if np.all(status):
        delete_old(logdir, days=7, verbose=verbose)
        # TODO - delete old files in jpsss folder when permissions are set correctly.
        # delete_old(jdir, days=7,verbose=verbose)
        # Opens json event files
        # Finds event file urls for download
        # Downloads netcdf files
        # Creates list of downloaded netcdf files for reference
        log_list = sorted(list(f for f in os.listdir(logdir) if f.endswith(".json")))
        for log in log_list:
            if verbose:
                print(logdir + log)
            num_downloaded = get_nc(
                logdir + log, vaac=vaac, mkdir=True, verbose=verbose, VOLCAT_DIR=ddir
            )
        # TO DO:
        # Could create a function that moves already downloaded netcdf files to new location
        # Some sort of filing system if desired
    return check_dirs(logdir, ddir, jdir, verbose=verbose)


def new_json(jdir, logdir, logfile="json_log.txt"):
    """
    Get list of json files pushed to our system.
    Inputs
    jdir: directory containing json event summary files (string)
    logdir: directory of json_log file (string)
    logfile: name of log file (string)
    Outputs:
    sum_list: list of summary json files (list)
    """
    original, current, added, removed = determine_change(jdir, logdir, logfile, ".json")
    return added


def open_json(fname):
    """Opens json file.
    Input: full filename (with directory) (string)
    Output: dictionary"""
    f = open(fname)
    jsonf = json.load(f)
    return jsonf


def open_dataframe(fname, varname=None):
    """Opens json file, and converts to pandas dataframe
    Inputs:
    fname: full filename (with directory) (string)
    varname: VOLCANOES if opening event summary json file (strin
             Not needed for event log json files.
             FILES if looking for event netcdf files from event log files (string)
             VAAC_REGION is looking for the vaac region of the events
             varname is case insenstive so can use FILES or files.
    Output: pandas dataframe (or string depending on varname)"""
    jsonf = open_json(fname)
    if varname:
        vstr = [x for x in jsonf.keys() if x.lower() == varname]
        vstr = vstr[0]
        dataf = jsonf[vstr]
    else:
        dataf = jsonf
    if type(dataf) == list:
        data = pd.DataFrame.from_dict(dataf)
    if type(dataf) == dict:
        data = pd.DataFrame.from_dict([dataf])
    if type(dataf) == str:
        data = dataf
    return data


def delete_old(directory, days=7, verbose=False):
    """Determines the age of the files in the specified folder.
    Deletes files older than 7 days, since this is the length of time
    the files exist on the wisconsin ftp site.
    CURRENTLY NOT CREATING A LOG
    Could modify to adjust time for deletion (longer or shower than 7 days)
    Inputs:
    directory: directory of files to determine age (string)
    Outputs:
    string: number of files deleted from directory and total size of files (string)
    """
    import time
    import shutil

    # import
    now = time.time()  # current time
    deletetime = now - (days * 86400)  # delete time
    deletefiles = []  # creating list of files to delete
    for files in os.listdir(directory):
        files = os.path.join(directory, files)  # Joining path and filename
        if os.stat(files).st_mtime < deletetime:
            if os.path.isfile(files):
                deletefiles.append(files)  # Creating list of files to delete
    if verbose:
        print("Files to be deleted: " + str(deletefiles))
    size = 0.0
    count = 0
    for count, deletef in enumerate(deletefiles):
        size += os.path.getsize(deletef) / (124 * 124)
        os.remove(deletef)
    if verbose:
        return print(
            "Deleted "
            + str(count)
            + " files, totalling "
            + str(round(size, 2))
            + " MB."
        )
    else:
        return None


def get_log_list(data):
    ## 5/25/2022 now this is done somewhere else.
    """Pulls the log url from the pandas dataframe
    Inputs:
    data: pandas dataframe
    Outputs:
    logurl: list of urls for the event log files
    """
    log_url = []
    if "EVENTS" in data.keys():
        events = data["EVENTS"]
    else:
        logger.warning("no EVENTS found in dictionary")
        return log_url
    for eve in events:
        if isinstance(eve, dict):
            if "LOG_URL" in eve.keys():
                log_url.append(eve["LOG_URL"])
            else:
                logger.warning("no LOG_URL found in dictionary")
        elif isinstance(eve, list):
            for subeve in eve:
                if "LOG_URL" in subeve.keys():
                    log_url.append(subeve["LOG_URL"])
                else:
                    logger.warning("no LOG_URL found in dictionary")
    return log_url


def get_log(log_url, verbose=False, **kwargs):
    """
     Downloads desired json event log files from ftp.
     downloads files only if newer version is available using -N option on wget.
     Log files are downloaded to specified location

    Inputs:
     log_url: list of strings
              list of urls to the log files

     verbose : boolean

     VOLCAT_LOGFILES : string
              location to download logfiles to.
     Returns:
     None
    """
    # Log_dir should be changed to something more generic (/pub/volcat_logs/ ?)
    rval = []
    if "VOLCAT_LOGFILES" in kwargs.keys():
        log_dir = kwargs["VOLCAT_LOGFILES"]
    else:
        log_dir = "/pub/ECMWF/JPSS/VOLCAT/LogFiles/"
    if isinstance(log_url, str):
        log_url = [log_url]
    for gurl in log_url:
        fname = gurl.split("/")[-1]
        # wget -r -l1 -A.nc  : retrieve all nc files
        # wget -N: timestamping - retrieve files only if newer than local
        # wget -P: designates location for file download
        os.system("wget -N -P " + log_dir + " " + gurl)
        if os.path.isfile(os.path.join(log_dir, fname)):
            #logger.info("File {} downloaded to {}".format(gurl, log_dir))
            rval.append((fname,True))
        else:
            logger.warning("FAILED to download File {} ".format(gurl, log_dir))
            rval.append((fname,False))
        # if verbose:
        #    logger.info("File {} downloaded to {}".format(gurl, log_dir))
        #    print("File {} downloaded to {}".format(gurl, log_dir))
    return rval


def open_log(logdir, logfile=None):
    """Opens the event file download log file
    Inputs:
    logdir: Directory location for data log file
    logfile: name of log file (string)
    Outputs:
    original: list of files already downloaded to our server (list)
    """
    import json

    if os.path.isfile(os.path.join(logdir, logfile)):
        with open(os.path.join(logdir, logfile), "r") as f:
            original = json.loads(f.read())
        return original
    else:
        return []


def check_file(fname, directory, suffix=".nc", verbose=False):
    """Checks if file in fname exists on our servers
    Inputs:
    fname: full path filename of file (string)
    directory: directory of data file list
    suffix: file suffix (string)
    outputs:
    Boolean: True, False
    """
    # original = open_log(directory)
    original = list(f for f in os.listdir(directory) if f.endswith(suffix))
    s = fname.rfind("/")
    current = fname[s + 1 :]
    if current in original:
        # if verbose:
        #    print('File '+current+' already downloaded')
        return False
    else:
        return True


def check_dirs(*args, verbose=False):
    status = []
    for direc in args:
        if not os.path.isdir(direc):
            status.append(False)
            if verbose:
                print("{} NOT FOUND".format(direc))
            logger.warning("{} directory not found".format(direc))
        else:
            status.append(True)
    return status


def determine_change(ddir, logdir, logfile, suffix):
    """Determines which files were original, which are current, which were added, which were removed

    Inputs:
    ddir:    data directory (string)
    logdir:  location of event file download log file (string)
    logfile: name of log file (string)
    suffix:  file suffix for list criteria (string)

    Returns
    original : list  : files in log when opened
    current  : list  : all files in the jpsss directory
    added    : listi : files in jpss directory but not in log
    removed  : list  : files in log but not in jpss directory.
    """
    status = check_dirs(ddir, logdir, verbose=True)
    if not np.all(status):
        return None, None, None, None
    # Files downloaded during previous check
    original = open_log(logdir, logfile=logfile)
    # Includes files just downloaded (if any)
    current = list(fi for fi in os.listdir(ddir) if fi.endswith(suffix))
    # Determining what was added and what was removed
    added = [fs for fs in current if not fs in original]
    removed = [fs for fs in original if not fs in current]
    return original, current, added, removed


def record_change(ddir=None, logdir=None, logfile=None, suffix=".nc", verbose=False):
    """Records file changes in data directory
    Inputs:
    ddir:    data directory (string)
    logdir:  location of event file download log file (string)
    logfile: name of log file (string)
    suffix:  file suffix for list criteria (string)
    Outputs:
    """
    original, current, added, removed = determine_change(ddir, logdir, logfile, suffix)

    if added:
        h = 0
        while h < len(added):
            original.append("".join(added[h]))
            h += 1
        if verbose:
            print("Added " + str(len(added)) + " files")
    if removed:
        g = 0
        while g < len(removed):
            original.remove("".join(removed[g]))
            g += 1
        if verbose:
            print("Removed " + str(len(removed)) + " files")
    if added or removed:
        with open(logdir + "tmp_file2.txt", "w") as fis:
            fis.write(json.dumps(original))
        Helper.move(
            os.path.join(logdir, "tmp_file2.txt"), os.path.join(logdir, logfile)
        )
        # os.system('mv '+logdir+'tmp_file2.txt '+logdir+logfile)
        # os.chmod(logdir+logfile, 0o666)
        if verbose:
            print("Updates recorded to file!\n")
        return None
    else:
        if verbose:
            print("No updates to " + ddir + " folder\n")
        return None


def record_missing(mlist, mdir, mfile="missing_files.txt", verbose=False):
    """Records files that are not downloaded.
    Inputs:
    mlist: list of missing files (list)
    mdir: directory to write missing file (string)
    mfile: missing file full name (string)
    Outputs:
    text file with list of missing files
    """
    # os.chmod(mdir+mfile, 0o666)
    if os.path.exists(mdir + mfile):
        txtfile = open(mdir + mfile, "a")
    else:
        txtfile = open(mdir + mfile, "w")
    for element in mlist:
        txtfile.write(element + "\n")
    txtfile.close()
    # os.chmod(mdir+mfile, 0o666)
    if verbose:
        print("Missing files added to {}".format(os.path.join(mdir, mfile)))
    return None


def fix_volc_name(volcname):
    """Fixes the volcano name if a comma, or space appear in the name"""
    if "," in volcname:
        s = volcname.find(",")
        tmp = volcname[:s]
        tmp2 = volcname[s + 2 :]
        volcname = tmp2 + "_" + tmp
    if " " in volcname:
        volcname = volcname.replace(" ", "_")
    return volcname


def make_volcdir(data_dir, fname, verbose=False):
    """Finds volcano name from json event log file.
    If name has ',' or spaces, the name is modified.
    Example: Tigre, Isla el --> Isla_el_Tigre
    Checks if a folder by that name already exists. If it does not, the folder is generated.
    Inputs:
    data_dir: directory where data are located
    fname: name of json event log file
    Outputs:
    volcname: string
    New directory is created if it didn't exist
    """
    volcname = fix_volc_name(volc)
    make_dir(data_dir, newdir=volcname, verbose=verbose)
    return volcname


def get_nc(fname, vaac=None, mkdir=True, verbose=False, **kwargs):
    """
    Finds and downloads netcdf files in json event log files from ftp.
    Netcdf event files are download to specified location

    Inputs:
    fname : string
            filename of json event log file
    vaac : string
           vaac region for file downloads Default: None, all files are downloaded
    mkdir: boolean
           make directory of volcano name, download files to that directory
    verbose: boolean
    kwargs : dictionary
             'VOLCAT_DIR' : string : netcdf files will be downloaded to subfolder
                                     with volcano name under this directory.

    Returns:
    num_downloaded : number of files downloaded
    """
    num_missing = 0
    num_downloaded = 0
    if "VOLCAT_DIR" in kwargs.keys():
        data_dir = kwargs["VOLCAT_DIR"]
    else:
        data_dir = "/pub/ECMWF/JPSS/VOLCAT/Files/"


    volcname = fix_volc_name(open_dataframe(fname)["volcano_name"].values[0])
    make_dir(data_dir, newdir=volcname, verbose=verbose)
    # volcname = make_volcdir(data_dir, volc, verbose=verbose)
    if vaac is not None:
        # checking that file is for specified vaac region
        vaac2 = open_dataframe(fname, varname="VAAC_REGION")
        if vaac != vaac2:
            # if the regions do not agree, then files are not downloaded
            if verbose:
                print("Files in " + fname + " are not for " + vaac + " VAAC region")
            return None
    dfiles = open_dataframe(fname, varname="FILES")
    dfile_list = dfiles["event_url"].values
    missing = []
    # Checking for type - if only one file in json event log, then event_url will be string
    # Need a list type
    if type(dfile_list) == str:
        dfile_list = [dfile_list]
    i = 0
    while i < len(dfile_list):
        # Check if file exists or has already been downloaded
        # If it has not, the download file from event_url

        file_download = check_file(dfile_list[i], data_dir + volcname, verbose=verbose)
        if file_download:
            # Might want to add a check for complete download of files
            # Need to figure out a good way to do this
            # os.system('wget -a '+data_dir+'data_logfile.txt --rejected-log=' +data_dir+'nodata_logfile.txt -P'+data_dir+' '+dfile_list[i])
            os.system("wget -P" + data_dir + volcname + "/ " + dfile_list[i])
            s = dfile_list[i].rfind("/")
            dfile = dfile_list[i][s + 1 :]
            if os.path.isfile(data_dir + volcname + "/" + dfile):
                if verbose:
                    print("File " + dfile + " downloaded to " + data_dir + volcname)
                num_downloaded += 1
            else:
                missing.append(dfile_list[i])
                num_missing += 1
                # if verbose:
                #    print('File '+dfile+' NOT DOWNLOADED!')
                #    print('From json file: '+fname)
        i += 1
    # record_change(ddir=data_dir, logdir=data_dir, logfile='data_logfile.txt')
    if len(missing) > 0:
        record_missing(missing, data_dir, mfile="missing_netcdfs.txt")
        if verbose:
            print(
                "File downloads complete. {} Missing files located in missing_netcdfs.txt".format(
                    num_missing
                )
            )
    if verbose:
        print("File downloads complete. {} files downloaded".format(num_downloaded))
    return num_downloaded


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


def num_files(data_dir, verbose=False):
    """Lists the subdirectories within the given directory
    and the number of files within the subdirectories.
    Inputs:
    data_dir: directory path of parent string()
    outputs:
    dirs_num: pandas DataFrame with directories and number of volcat files within
    """
    import pandas as pd

    vdirlist = list_dirs(data_dir)
    volcname = []
    numfiles = []
    for fdir in vdirlist:
        tdir = data_dir + "{}/".format(fdir)
        fnames = glob.glob(tdir + "*.nc")
        volcname.append(fdir)
        numfiles.append(len(fnames))
        if verbose:
            print(fdir, len(fnames))
    vdf = pd.DataFrame(volcname, columns=["Volcano Name"])
    numdf = pd.DataFrame(numfiles, columns=["Num Files"])
    dirs_num = pd.concat([vdf, numdf], axis=1)
    return dirs_num


def get_latlon(data_dir):
    """Read csv file containing volcano name, latitude, longitude
    Inputs:
    data_dir: directory where Volcanoes.csv is located
    Outputs:
    volcdf: pandas dataframe of volcanos, latitude, longitude
    """
    import pandas as pd

    volcdf = pd.read_csv(data_dir + "Volcanoes.csv", sep=",")
    return volcdf


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


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
from utilhysplit.runhandler import ProcessList
from utilhysplit.plotutils import map_util
import utilhysplit.evaluation.web_ensemble_plots as wep
from utilvolc import make_data_insertion as mdi
from utilhysplit.evaluation import ensemble_tools

logger = logging.getLogger(__name__)

"""
 Classes and functions for working with VOLCAT alerts and log files

 Classes
   VFile (base class)
     EventFile
     SummaryFile

 Functions
   check_dirs
   check_file
   delete_old
   determine_change
   get_nc
   get_summary_file_df 
   num_files
   open_json
   open_dataframe
"""




def get_summary_file_df(fdir, verbose=False, hours=48, edate = datetime.datetime.now()):
    """
    fdir : str : location of summary files
    Returns
    sumdf : pandas dataframe with information from all the summary files.

    prints an error if one of the files cannot be read
    """

    # 2023 09 May (amc) use the time stamp in the filename rather than scandir
    if not isinstance(edate,datetime.datetime):
       raise Exception("input edate wrong type {}: should be \
                        datetime.datetime".format(type(edate)))


    vlist = []
    if not hours:
        logger.warning("get_summary_file_df input hours None, setting hours to 100")
        hours = 100
    strf = "s%Y%j_%H"
    sdate = edate - datetime.timedelta(hours=hours)
    s_str = sdate.strftime(strf)
    done=False
    flist = []
    iii=0
    while not done:
 
          flist.extend(glob.glob(fdir + '/*{}*.json'.format(s_str)))
          sdate += datetime.timedelta(hours=1)
          s_str = sdate.strftime(strf)
          if sdate > edate:  done=True
          iii+=1
          if iii>10: done=True
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
    if len(vlist)>0:
        sumdf = pd.concat(vlist)
    return sumdf


class VFile:
    """
    Base class for the files
    """
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
    """
    Volcat summary files
    """
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


class EventFile(VFile):
    """
    Volcat event files. These are json files which are pulled from ftp 
    and contain information about the data files. data file classes are in volcat.py

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
     list of tuples (filename, boolean)
     The boolean is True if the file was successfully downloaded and false if it was not.
    """
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

def num_files(data_dir, verbose=False):
    """Lists the subdirectories within the given directory
    and the number of files within the subdirectories.
    Inputs:
    data_dir: directory path of parent string()
    outputs:
    dirs_num: pandas DataFrame with directories and number of volcat files within

    -------------------------------------------------------------------------------
    Example usage:

    numfiles('/topdir/VOLCAT/Files/')

    Example output:

        Volcano Name    Num Files
    0   Abu            17
    1   Adatarayama    1
    2   Agua           8
    3   Agung          110
    4   Aira           104

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






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


from utilvolc.volcat_event import create_event


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
 2023 Jul 26 AMC  moved Events class to volcat_events.py  



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
        gfile = os.path.join(self.dirhash["VOLCAT_DIR"], "green_list.txt")
        with open(gfile, "r") as fid:
            greenlist = fid.readlines()
        greenlist = [x.strip() for x in greenlist]
        return greenlist

    def write_parallax_corrected(self, gridspace=None, elist=None):

        failed_list = []
        if not isinstance(elist,list):
           elist = self.ehash.keys()

        print("WRITING FOR THESE", elist)
        for volc in elist:
            print("writing parallax corrected for {}".format(volc))
            try:
                self.ehash[volc].write_parallax_corrected(gridspace=gridspace)
            except Exception as eee:
                print("FAILED to create pc corrected")
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


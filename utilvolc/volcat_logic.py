import json
import pandas as pd


def workflow():
    # check for event summary files and read
    # decide which file list json files to pull
    # these files change so will need to be re-pulled (or maybe new files are written?).
    # pull file list json files and read
    # decide which event files to pull.
    # check if file already exists.
    # os.isfile(filename)
    # check if file meets time resolution requirement (spacing every 10 minutes)
    # check if file meets other requirements as needed.
    # pull file(s)
    # link from json dictionary
    # make parallax corrected files
    # make emit-times files
    # area calculation
    # some decision logic about when and what kind of HYSPLIT runs we make.
    # are they automatic
    # are they triggered by a request from web.
    # Limited number of runs for active volcano and then link to READY website
    # where they can choose to modify inputs.
    # populate with latitude longitude what satellite retrievals are available.etc

    # run HYSPLIT data insertion
    # write CONTROL file
    # write SETUP file

    # regrid volcat to HYSPLIT grid
    # needed for inverse modeling
    # needed for evaluation
    # two possible paths
    # write regridded file for every volcat retrieval
    # write hour averaged regridded files.

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


def open_json(fname):
    """Opens json file.
    Input: full filename (with directory) (string)
    Output: dictionary"""
    f = open(fname)
    jsonf = json.load(f)
    return jsonf


def open_dataframe(fname):
    """ Opens json file, and converts to pandas dataframe
    Inputs: 
    fname: full filename (with directory) (string)
    Output: pandas dataframe"""
    jsonf = open_json(fname)
    dataf = jsonf['VOLCANOES']
    data = pd.DataFrame.from_dict(dataf)
    return data


def get_log_list(data):
    """Pulls the log url from the pandas dataframe
    Inputs:
    data: pandas dataframe
    Outputs:
    logurl: list of urls for the event log files
    """
    events = data['EVENTS']
    log_url = []
    i = 0
    while i < len(events):
        tmp = pd.DataFrame([events[i]])
        log_url.append(tmp['LOG_URL'].values[0])
        i += 1
    return log_url


def get_log(log_url, verbose=False):
    """ Downloads desired json event log files from ftp.
    Inputs:
    log_url: list of urls to the log files
    Outputs:
    Log files are downloaded to specified location
    """
    import os
    # Log_dir should be changed to something more generic (/pub/volcat_logs/ ?)
    log_dir = '/hysplit-users/allisonr/Hysplit_Tools/utilhysplit/utilvolc/data/'
    i = 0
    while i < len(log_url):
        os.system('wget -P '+log_dir+' '+log_url[i])
        if verbose:
            print('File '+log_url[i]+' downloaded to '+log_dir)
        i += 1
    return 0

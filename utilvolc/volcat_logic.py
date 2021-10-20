import json
import pandas as pd
import os


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


def new_json(jdir, directory, logfile='json_log.txt', verbose=False):
    """ Get list of json files pushed to our system
    Inputs
    jdir: directory containing json event summary files (string)
    directory: directory of json_log file (string)
    logfile: name of log file (string)
    Outputs:
    sum_list: list of summary json files (list)
    """
    # record_change(ddir=jdir, directory=directory, logfile=logfile, suffix='.json', verbose=verbose)
    original, current, added, removed = determine_change(jdir, directory, logfile, '.json')

    return added


def open_json(fname):
    """Opens json file.
    Input: full filename (with directory) (string)
    Output: dictionary"""
    f = open(fname)
    jsonf = json.load(f)
    return jsonf


def open_dataframe(fname, varname=None):
    """ Opens json file, and converts to pandas dataframe
    Inputs:
    fname: full filename (with directory) (string)
    varname: VOLCANOES if opening event summary json file (string)
                    Not needed for event log json files.
                    FILES if looking for event netcdf files from event log files (string)
    Output: pandas dataframe"""
    jsonf = open_json(fname)
    if varname:
        dataf = jsonf[varname]
    else:
        dataf = jsonf
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
    return print('Event log json files downloaded')


def open_log(directory, logfile=None):
    """Opens the event file download log file
    Inputs:
    directory: Directory location for data log file
    logfile: name of log file (string)
    Outputs:
    original: list of files already downloaded to our server (list)
    """
    import json
    with open(directory+logfile, 'r') as f:
        original = json.loads(f.read())
    return original


def check_file(fname_url, directory, verbose=False):
    """ Checks if file in fname_url exists on our servers
    Inputs:
    fname_url: full ftp url for file (string)
    directory: directory of data file list
    outputs:
    Boolean: True, False
    """
    original = open_log(directory)
    s = fname_url.rfind('/')
    current = fname_url[s+1:]
    if current in original:
        if verbose:
            print('File '+current+' already downloaded')
        return False
    else:
        return True


def determine_change(ddir, directory, logfile, suffix):
    """Determines which files were original, which are current, which were added, which were removed"""
    # Files downloaded during previous check
    original = open_log(directory, logfile=logfile)
    # Includes files just downloaded (if any)
    current = list(fi for fi in os.listdir(ddir) if fi.endswith(suffix))
    # Determining what was added and what was removed
    added = [fs for fs in current if not fs in original]
    removed = [fs for fs in original if not fs in current]
    return original, current, added, removed


def record_change(ddir=None, directory=None, logfile=None, suffix='.nc', verbose=False):
    """Records file changes in data directory
    Inputs:
    ddir: data directory (string)
    directory: location of event file download log file (string)
    logfile: name of log file (string)
    suffix: file suffix for list criteria (string)
    Outputs:
    """
    original, current, added, removed = determine_change(ddir, directory, logfile, suffix)

    if added:
        h = 0
        while h < len(added):
            original.append(''.join(added[h]))
            h += 1
        if verbose:
            print('Added '+str(len(added))+' files')
    if removed:
        g = 0
        while g < len(removed):
            original.remove(''.join(removed[g]))
            g += 1
        if verbose:
            print('Removed '+str(len(removed))+' files')
    if added or removed:
        with open(directory+'tmp_file2.txt', 'w') as fis:
            fis.write(json.dumps(original))
        os.system('mv '+directory+'tmp_file2.txt '+directory+logfile)
        return print('Updates recorded to file!\n')
    else:
        return print('No updates to '+ddir+' folder\n')


def get_nc(fname, verbose=False):
    """ Finds and downloads netcdf files in json event log files from ftp.
    Inputs:
    fname: filename of json event log file
    Outputs:
    Netcdf event files are download to specified location
    """
    import os
    # This should be changed, a specified data file location
    data_dir = '/pub/ECMWF/JPSS/VOLCAT/Files/'
    dfiles = open_dataframe(fname, varname='FILES')
    dfile_list = dfiles['EVENT_URL'].values
    i = 0
    while i < len(dfile_list):
        # Check if file exists or has already been downloaded
        # If it has not, the download file from event_url

        file_download = check_file(dfile_list[i], data_dir, verbose=verbose)
        if file_download:
            # os.system('wget -a '+data_dir+'data_logfile.txt --rejected-log=' +data_dir+'nodata_logfile.txt -P'+data_dir+' '+dfile_list[i])
            os.system('wget -P'+data_dir+' '+dfile_list[i])
            if verbose:
                print('File '+dfile_list[i]+' downloaded to '+data_dir)
        i += 1
    record_change(ddir=data_dir, directory=data_dir, logfile='data_logfile.txt')
    return print('File downloads complete')

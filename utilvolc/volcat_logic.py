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


def get_files(verbose=False):
    """ Use various functions to get all available netcdf files from json event log files
    Uses the different functions within volcat_logic.py
    Sorts through the ftp directory to find json event summary files
    Loops through them, and downloads the corresponding json event log files
    Keeps track of the json event summary files that have already been downloaded
    Lists event log files, finds the correspondind event file urls for download
    Downloads the event netcdf files
    Keeps track of the downloaded netcdf files"""

    jdir = '/pub/jpsss_upload/'
    ddir = '/pub/ECMWF/JPSS/VOLCAT/Files/'
    logdir = '/pub/ECMWF/JPSS/VOLCAT/LogFiles/'

    # Delete files from jpsss_uploads folder that is older than 7 days
    # Files are only available for 7 days on the ftp
    # Dont have permissions to delete files from jpsss_upload!
    # delete_old(jdir, verbose=verbose)

    # Finds json files added to ftp folder
    added = new_json(jdir, logdir, verbose=verbose)
    added = sorted(added)
    i = 0
    while i < len(added):
        data = open_dataframe(jdir+added[i], varname='VOLCANOES')
        log_url = get_log_list(data)
        # Downloads json event log files
        get_log(log_url, verbose=verbose)
        i += 1
    # Logs event summary json files
    record_change(ddir=jdir, logdir=logdir, logfile='json_log.txt', suffix='.json', verbose=verbose)

    # Delete files from json event log folder that are older than 7 days
    # Netcdf files are only available for 7 days on the ftp
    delete_old(logdir, verbose=verbose)

    # Opens json event files
    # Finds event file urls for download
    # Downloads netcdf files
    # Creates list of downloaded netcdf files for reference
    log_list = sorted(list(f for f in os.listdir(logdir) if f.endswith('.json')))
    x = 0
    while x < len(log_list):
        print(logdir+log_list[x])
        get_nc(logdir+log_list[x], mkdir=True, verbose=verbose)
        x += 1

    # TO DO:
    # Could create a function that moves already downloaded netcdf files to new location
    # Some sort of filing system if desired


def new_json(jdir, logdir, logfile='json_log.txt', verbose=False):
    """ Get list of json files pushed to our system
    Inputs
    jdir: directory containing json event summary files (string)
    logdir: directory of json_log file (string)
    logfile: name of log file (string)
    Outputs:
    sum_list: list of summary json files (list)
    """
    original, current, added, removed = determine_change(jdir, logdir, logfile, '.json')
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
    if type(dataf) == list:
        data = pd.DataFrame.from_dict(dataf)
    if type(dataf) == dict:
        data = pd.DataFrame.from_dict([dataf])
    return data


def delete_old(directory, verbose=False):
    """ Determines the age of the files in the specified folder.
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
    import glob
    # import
    days = 7
    now = time.time()  # current time
    deletetime = now - (days * 86400)  # delete time
    deletefiles = []  # creating list of files to delete
    for files in os.listdir(directory):
        files = os.path.join(directory, files)  # Joining path and filename
        if os.stat(files).st_mtime < deletetime:
            if os.path.isfile(files):
                deletefiles.append(files)  # Creating list of files to delete
    if verbose:
        print("Files to be deleted: "+str(deletefiles))
    count = 0
    size = 0.0
    mm = 0
    while mm < len(deletefiles):
        size = size + (os.path.getsize(deletefiles[mm]) / (124*124))
        os.remove(deletefiles[mm])
        count = count+1
        mm += 1
    if verbose:
        return print('Deleted '+str(count) + ' files, totalling '+str(round(size, 2))+' MB.')
    else:
        return print('Done')


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
        if type(events[i]) == dict:
            tmp = pd.DataFrame([events[i]])
            log_url.append(tmp['LOG_URL'].values[0])
        elif type(events[i]) == list:
            j = 0
            while j < len(events[i]):
                tmp = pd.DataFrame([events[i][j]])
                log_url.append(tmp['LOG_URL'].values[0])
                j += 1
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
    log_dir = '/pub/ECMWF/JPSS/VOLCAT/LogFiles/'
    i = 0
    while i < len(log_url):
        # wget -N: timestamping - retrieve files only if newer than local
        # wget -P: designates location for file download
        os.system('wget -N -P '+log_dir+' '+log_url[i])
        if verbose:
            print('File '+log_url[i]+' downloaded to '+log_dir)
        i += 1
    return print('Event log json files downloaded')


def open_log(logdir, logfile=None):
    """Opens the event file download log file
    Inputs:
    logdir: Directory location for data log file
    logfile: name of log file (string)
    Outputs:
    original: list of files already downloaded to our server (list)
    """
    import json
    with open(logdir+logfile, 'r') as f:
        original = json.loads(f.read())
    return original


def check_file(fname_url, directory, suffix='.nc', verbose=False):
    """ Checks if file in fname_url exists on our servers
    Inputs:
    fname_url: full ftp url for file (string)
    directory: directory of data file list
    suffix: file suffix (string)
    outputs:
    Boolean: True, False
    """
    # original = open_log(directory)
    original = list(f for f in os.listdir(directory) if f.endswith(suffix))
    s = fname_url.rfind('/')
    current = fname_url[s+1:]
    if current in original:
        if verbose:
            print('File '+current+' already downloaded')
        return False
    else:
        return True


def determine_change(ddir, logdir, logfile, suffix):
    """Determines which files were original, which are current, which were added, which were removed"""
    # Files downloaded during previous check
    original = open_log(logdir, logfile=logfile)
    # Includes files just downloaded (if any)
    current = list(fi for fi in os.listdir(ddir) if fi.endswith(suffix))
    # Determining what was added and what was removed
    added = [fs for fs in current if not fs in original]
    removed = [fs for fs in original if not fs in current]
    return original, current, added, removed


def record_change(ddir=None, logdir=None, logfile=None, suffix='.nc', verbose=False):
    """Records file changes in data directory
    Inputs:
    ddir: data directory (string)
    logdir: location of event file download log file (string)
    logfile: name of log file (string)
    suffix: file suffix for list criteria (string)
    Outputs:
    """
    original, current, added, removed = determine_change(ddir, logdir, logfile, suffix)

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
        with open(logdir+'tmp_file2.txt', 'w') as fis:
            fis.write(json.dumps(original))
        os.system('mv '+logdir+'tmp_file2.txt '+logdir+logfile)
        return print('Updates recorded to file!\n')
    else:
        return print('No updates to '+ddir+' folder\n')


def record_missing(mlist, mdir, mfile='missing_files.txt'):
    """ Records files that are not downloaded.
    Inputs:
    mlist: list of missing files (list)
    mdir: directory to write missing file (string)
    mfile: missing file full name (string)
    Outputs:
    text file with list of missing files
    """
    if os.path.exists(mdir+mfile):
        txtfile = open(mdir+mfile, 'a')
    else:
        txtfile = open(mdir+mfile, 'w')
    for element in mlist:
        txtfile.write(element + '\n')
    txtfile.close()
    return print('Missing files added to '+mdir+mfile)


def make_dir(data_dir, fname, verbose=False):
    """ Finds volcano name from json event log file.
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
    volcname = open_dataframe(fname)['VOLCANO_NAME'].values[0]
    if ',' in volcname:
        s = volcname.find(',')
        tmp = volcname[:s]
        tmp2 = volcname[s+2:]
        volcname = tmp2+'_'+tmp
    if ' ' in volcname:
        volcname = volcname.replace(' ', '_')
    if not os.path.exists(data_dir+volcname):
        os.makedirs(data_dir+volcname)
        if verbose:
            print('Directory '+data_dir+volcname+' created')
    return volcname


def get_nc(fname, mkdir=True, verbose=False):
    """ Finds and downloads netcdf files in json event log files from ftp.
    Inputs:
    fname: filename of json event log file
    mkdir: make directory of volcano name, download files to that directory (boolean)
    verbose: (boolean)
    Outputs:
    Netcdf event files are download to specified location
    """
    # This should be changed, a specified data file location
    data_dir = '/pub/ECMWF/JPSS/VOLCAT/Files/'
    volcname = make_dir(data_dir, fname, verbose=verbose)
    dfiles = open_dataframe(fname, varname='FILES')
    dfile_list = dfiles['EVENT_URL'].values
    missing = []
    # Checking for type - if only one file in json event log, then event_url will be string
    # Need a list type
    if type(dfile_list) == str:
        dfile_list = [dfile_list]
    i = 0
    while i < len(dfile_list):
        # Check if file exists or has already been downloaded
        # If it has not, the download file from event_url

        file_download = check_file(dfile_list[i], data_dir+volcname, verbose=verbose)
        if file_download:
            # Might want to add a check for complete download of files
            # Need to figure out a good way to do this
            # os.system('wget -a '+data_dir+'data_logfile.txt --rejected-log=' +data_dir+'nodata_logfile.txt -P'+data_dir+' '+dfile_list[i])
            os.system('wget -P'+data_dir+volcname+'/ '+dfile_list[i])
            s = dfile_list[i].rfind('/')
            dfile = dfile_list[i][s+1:]
            if os.path.isfile(data_dir+volcname+'/'+dfile):
                if verbose:
                    print('File '+dfile+' downloaded to '+data_dir+volcname)
            else:
                missing.append(dfile_list[i])
                if verbose:
                    print('File '+dfile+' NOT DOWNLOADED!')
                    print('From json file: '+fname)
        i += 1
    #record_change(ddir=data_dir, logdir=data_dir, logfile='data_logfile.txt')
    if len(missing) > 0:
        record_missing(missing, data_dir, mfile='missing_netcdfs.txt')
        return print('File downloads complete. Missing files located in missing_netcdfs.txt')
    else:
        return print('File downloads complete. No missing files.')


def write_emitimes(fname, emitdir, inputdict, verbose=True):
    """Writes emitimes file from the volcat netcdf file provided. Still in progress!
    Inputs:
    fname: name of netcdf file (string)
    emitdir: Directory to write emitimes files (string)
    inputdict: dictionary of inputs for writing emitimes file (number of particles, duration,layer, etc.)
    Output:
    emitimes file written to designated directory
    """
    from utilvolc import write_emitimes as we

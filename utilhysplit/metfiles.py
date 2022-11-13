# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import datetime
import logging
import sys
from os import path

"""
PRGMMR: Alice Crawford  ORG: ARL
This code written at the NOAA  Air Resources Laboratory
ABSTRACT: choosing met files for HYSPLIT control file

All dates assumed to be in UTC time.


"""
logger = logging.getLogger(__name__)

#TODO
# if using archived gdas1, does not figure out
# how many files to use correctly.
    

def gefs_suffix_list():
    base0 = 'gec'
    base = 'gep'
    suffix = []
    for num in range(0,31,1):
        if num==0:
           suffix.append(base0 + str(num).zfill(2))
        else:
           suffix.append(base + str(num).zfill(2))
    return suffix

class MetFileFinder:
   """
   METHODS
   set_forecast_directory
   find_forecast_cycle
   find_forecast
   """
   def  __init__(self, metid):
       self.metid = metid
       self.forecast_directory = '/pub/forecast'
       self.archive_directory = '/pub/archive'
       self.mstr = get_forecast_str(self.metid, self.forecast_directory)
 
   def set_forecast_directory(self, dpath):
       self.forecast_directory = dpath 
       self.mstr = get_forecast_str(self.metid, self.forecast_directory)

   def set_archives_directory(self, dpath):
       self.archive_directory = dpath 
       #self.mstr = get_archive_str(self.metid, self.forecast_directory)

   def set_ens_member(self, suffix):
       self.mstr = get_forecast_str(self.metid, self.forecast_directory)
       self.mstr += suffix
       logger.info('Setting mstr for ensemble {}'.format(self.mstr)) 

   def find_forecast_cycle(self,dstart,duration,cycle):
       """
       cycle : str '00', '06', '12', '18'
       return files in a particular forecast cycle.
       backward runs don't work quite right.
       """
       forecast_info = get_forecast_info(self.metid)
       cycle_time = forecast_info['cycle_time']
       forecast_length = forecast_info['forecast_length']
       mstr = self.mstr.replace('%H',cycle)
       # e.g.
       # if using the 12z and dstart is before 12 then must
       # get the previous day. 
       # change start date by -24 hours and
       # add cycle time + start date hour to the duration
       # forward run
       if dstart.hour < int(cycle) and duration > 0:
              dt = datetime.timedelta(hours=24)
              newdate = dstart - dt
              runtime = (duration + int(cycle) + dstart.hour)
       # backward run
       elif dstart.hour < int(cycle) and duration < 0:
              #dt = datetime.timedelta(hours=24)
              #newdate = dstart - dt
              newdate = dstart 
              runtime = duration - int(cycle)
       else:
          newdate = dstart
          runtime = duration
       #print(dstart, duration, dstart+datetime.timedelta(hours=duration))
       #print(newdate, runtime, newdate + datetime.timedelta(hours=runtime))
       mf = MetFiles(mstr,hours=24)
       files = []
       iii=0
       while not files:
          newdate = newdate - datetime.timedelta(hours=cycle_time * iii)
          files = mf.make_file_list(newdate, runtime)
          if cycle_time *iii > forecast_length: break
          iii += 1
       return files    

   def find(self, dstart, duration,hours=-1):
       #self.mstr = get_forecast_str(self.metid, self.archive_directory)
       metfiles = self.find_forecast(dstart,duration)
       if not metfiles:
           logger.info('Looking in archive for met files')
           print('Looking in archive for met files')
           metfiles = self.find_archive(dstart,duration,hours=hours)
       return metfiles

   def find_archive(self, dstart, duration,hours=-1):
       self.mstr = get_archive_str(self.metid, self.archive_directory)
       mf = MetFiles(self.mstr,hours=hours)
       mfiles = mf.get_files(dstart, duration)
       return  mfiles

   def find_forecast(self, dstart, duration):
       """
       dstart : datetime object.
       duration : integer. can be positive or negative
                  negative for backward runs.
       Returns: 
       list of files that cover dstart to dstart +/- duration.
       """
       forecast_info = get_forecast_info(self.metid)
       cycle_time = forecast_info['cycle_time']
       forecast_length = forecast_info['forecast_length']
       mstr = self.mstr
       mf = MetFiles(mstr,hours=cycle_time)
       files = []
       iii=0
       while not files:
          # start looking at start date and
          # keep going back 6 hours until find first cycle that covers
          # the time period.
          newdate = dstart - datetime.timedelta(hours=cycle_time * iii)
          files = mf.make_file_list(newdate, duration)
          # double check that start of first files is before start date.
          if files:
              metdate = datetime.datetime.strptime(files[0],mstr)
              if metdate > newdate:
                 files = []
          if cycle_time *iii > forecast_length: break
          iii += 1
       if not files:
          logger.warning('No meteorological files found {}'.format(mstr))
          return files
       files = weed_files(files,dstart,duration,self.metid,self.mstr)
       return process(files)   

def weed_files(metfiles,dstart,duration,metid,metstr):
    #TO DO. fix so files don't overlap in time if more than one.
    mhash = get_forecast_info(metid)
    metlist =[]
    for count, mfile in enumerate(metfiles): 
      logger.debug('{}'.format(mfile))
      metf = mfile.split('/')[-1]
      metdir = mfile.replace(metf,'')
      mdate = datetime.datetime.strptime(mfile, metstr)       
      edate = mdate + datetime.timedelta(hours=mhash['forecast_length']) 
      metlist.append(mfile)
      # if the end date is past the end date of the simulation then don't
      # need any more files. 
      if edate > dstart + datetime.timedelta(hours=duration):
         break 
    return metlist

def get_forecast_info(metid):
    # TO DO. set forecast_length based on metid
    mhash = {}
    mhash['cycle_time'] = 12       # default cycle time
    mhash['time_res'] = 1          # default time resolution
    mhash['forecast_length'] = 24  # default forecast length

    if 'gfs' in metid.lower():
        mhash['forecast_length'] = 84
        mhash['time_res'] = 3  #3 hour time resolution

    elif 'gefs' in metid.lower():
        mhash['cycle_time'] = 6
        mhash['forecast_length'] = 6
        mhash['time_res'] = 3  #3 hour time resolution

    elif 'nam' in metid.lower():
        mhash['forecast_length'] = 72

    elif 'gfs0p25' in metid.lower():
        mhash['forecast_length'] = 24

    elif('era5' in metid.lower()):
        mhash['forecast_length'] = 24

    elif('merra2' in metid.lower()):
        mhash['forecast_length'] = 24

    else:
        logger.error('get_forecast_info: metid not recognized')

    return mhash

def get_archive_str(metid, ARCDIR='/pub/archive'):
    if (metid.lower() == 'gfs0p25'):
        metstr = 'gfs0p25/%Y%m%d_gfs0p25'
    elif (metid.lower() == 'gfs'):
        metstr = 'gdas1/gdas1.%b%y.week'
    elif('gefs' in metid.lower()):
        metstr = 'gefs'
    elif('era5' in metid.lower()):
        metstr = 'ERA5_%Y%m%d.ARL'
    elif('merra2' in metid.lower()):
        metstr = 'MERRA2_%Y%m%d.ARL'
    else:
        logger.warning('METID not found for archive {}'.format(metid))
        sys.exit()
    return path.join(ARCDIR, metstr)

def get_forecast_str(metid, FCTDIR='/pub/forecast'):
    """Finds forecast meteorology data files
    dstart : datetime object of start date
    metdata : string (options: GEFS, GFS, GFS0p25, NAM, NAMAK, NAMHI)
    """
    if (metid.lower() == 'gfs') or (metid.lower() == 'gfs0p25'):
        met = metid.lower()
    elif('nam' in metid.lower() and 'ak' in metid.lower()):
        met = 'namsf.AK'
    elif('nam' in metid.lower() and 'hi' in metid.lower()):
        met = 'namsf.HI'
    elif('nam' in metid.lower()):
        met = 'namsf'
    elif('gefs' in metid.lower()):
        met = 'gefs'
    else:
        logger.warning('Did not recognize MET name {}.'.format(metid))
        logger.warning('Using GFS')
        met = 'gfs' 
    #metnamefinal = 'No data found'
    #        metime = dtm
    metdir = path.join(FCTDIR , '%Y%m%d/')
    metfilename = 'hysplit.t%Hz.' + met 
    if 'gfs' in metid.lower():
        metfilename += 'f'
        #metfilename = 'hysplit.' + metime.strftime('t%Hz') + '.' + met
    return metdir + metfilename

def getmetfiles(strfmt, sdate, runtime,
                altstrfmt):
    """
    INPUTS
    strfmt : string
    sate : datetime object
    runtime : integer
    altstrfmt : string (currently not used)
    RETURNS :
    list of tuples (directory, filename)
    """
    mfiles = MetFiles(strfmt)
    return mfiles.get_files(sdate, runtime)

class MetFiles:
    """
    Class for finding metfiles used in HYSPLIT CONTROL file.
    """

    def __init__(self, strfmt, hours=-1, 
                 altstrfmt=None, althours=-1, verbose=False):
        """
        INPUTS
        strfmt : string
        hours : integer
        verbose : boolean

        Attributes
        mdt : datetime timedelta object.
              how much time one met file spans in hours.
        verbose : boolean.
        strfmt : string

        """
        self.verbose = verbose
        self.strfmt = strfmt
        if hours >= 0:
            self.mdt = datetime.timedelta(hours=hours)
        else:
            self.mdt = -1
        #self.altmet = MetFiles(altstrfmt, hours=althours, verbose=verbose)         
        self.altstrfmt = altstrfmt
        self.maxfiles = 50

    def __bool__(self):
        if self.strfmt : return True
        else: return False

    def set_mdt(self, hours):
        """
        set self.mdt.
        can be used to over-ride algorithm in find_mdt method.
        """
        self.mdt = datetime.timedelta(hours=hours)

    def get_files(self, sdate, runtime):
        # MetFiles class
        """
        sdate : datetime object. start date.
        runtime : integer. hours of runtime.
        RETURNS :
        list of tuples (directory, filename)
        """
        nlist = self.make_file_list(sdate, runtime)
        return process(nlist)

    @staticmethod
    def sub_handle_hour(sdate):
        """
        INPUTS
        sdate : datetime.datetime object
        returns
        date with 0 hour.
        """
        sdate = sdate.replace(tzinfo=None)
        year = sdate.year
        month = sdate.month
        day = sdate.day
        #hour = sdate.hour
        testdate = datetime.datetime(year, month, day)
        return testdate


    def handle_hour(self, sdate):
        """
        Handles met files which span less than one day.
        Handles finding nearest forecast hour as well.
        """
        # testdate is sdate with hour=0. beginning of day.
        sdate = sdate.replace(tzinfo=None)
        testdate = self.sub_handle_hour(sdate)
        done = False
        imax = 100
        iii = 0
        while not done:
            # add self.mdt to the testdate
            newdate = testdate + self.mdt
            # if not greater than the input date
            # then keep adding.
            if newdate <= sdate:
                testdate = newdate
            # else return this date.
            else:
                done = True
            iii += 1
            if iii > imax:
                done = True
        return testdate

    def find_mdt(self, testdate):
        """
        # finds time spacing between met files by
        # seeing which spacing produces a new file name.
        # this can be problematic if there are missing files.
        # testdate = datetime.datetime(2010,1,1)
        """
        if "%H" in self.strfmt:
            mdtlist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
            testdate = self.sub_handle_hour(testdate)
        else:
            mdtlist = [1, 24, 24 * 7, 24 * 31, 24 * 356]

        file1 = testdate.strftime(self.strfmt)
        done = False
        iii = 0
        while not done:
            dttt = datetime.timedelta(hours=mdtlist[iii])
            dtt2 = testdate + dttt
            file2 = dtt2.strftime(self.strfmt)
            if file2 != file1 and path.isfile(file2):
                done = True
            iii += 1
            if iii >= len(mdtlist):
                done = True
        ##AMC ADDED TEMP
        #dttt=datetime.timedelta(hours=24)
        return dttt

    def make_file_list(self, sdate, runtime):
        """
        INPUTS
        sdate : datetime.datetime ojbect
        runtime : int (hours of runtime)
        """
        nlist = []
        sdate = sdate.replace(tzinfo=None)
        if not isinstance(self.mdt, datetime.timedelta):
            self.mdt = self.find_mdt(sdate)
        # handle backwards runs. by switching sdate and edate
        if runtime < 0:
            runtime = abs(runtime)
            end_date = sdate
            sdate = end_date - datetime.timedelta(hours=runtime)
        else:
            end_date = sdate + datetime.timedelta(hours=runtime)
        done = False
        # self.verbose=True
        if "%H" in self.strfmt:
            sdate = self.handle_hour(sdate)
        edate = sdate
        #if self.verbose:
        #   logger.inro("GETMET", sdate, edate, end_date, runtime, self.mdt)
        zzz = 0
        while not done:
            if "week" in self.strfmt:
                self.mdt = datetime.timedelta(hours=7 * 24)
                temp = parse_week(self.strfmt, edate)
            else:
                temp = edate.strftime(self.strfmt)

            # this is beginning of forecast in the file.
            mdate = datetime.datetime.strptime(temp,self.strfmt)
            # end time of this particular file.
            medate = mdate + self.mdt
            #print('file', temp, self.strfmt)
            #print('begin time of file', mdate)
            #print('end time of file', medate)
            #print('-------------')

            #temp = temp.lower()
            #print('MDT', self.mdt) 
            # also need to increment the edate to get next possible file name 
            edate = edate + self.mdt
            #if not path.isfile(temp):
            #    temp = temp.lower()
            if not path.isfile(temp):
                logger.info("WARNING " +  temp + " meteorological file does not exist")
                print("WARNING " +  temp + " meteorological file does not exist")
                #pass
                #logger.debug("WARNING " +  temp + " meteorological file does not exist")
                #temp = self.altmet.makefilelist(edate, self.altmet.mdt)
                #print("REPLACE with", temp)
            #if temp != "None":
            elif temp not in nlist:
               nlist.append(temp)
               zzz += 1
            if medate > end_date:
                done = True
            if zzz > self.maxfiles:
                done = True
                print('warning: maximum number of met files reached {}'.format(self.maxfiles))
        return nlist

def process(nlist):
    """
    # convert full path to
    # list of directories and list of filenames
    # and then zips the lists to return list of tuples.
    """
    mfiles = []
    mdirlist = []
    for temp in nlist:
        siii = [x for x, char in enumerate(temp) if char == "/"]
        siii = siii[-1]
        fname = temp[siii + 1 :]
        mdir = temp[0 : siii + 1]
        mfiles.append(fname)
        mdirlist.append(mdir)
    return list(zip(mdirlist, mfiles))

def parse_week(strfmt, edate):
    """
    helps handle weekly files.
    """
    # used if week is in the strfmt (mostly for gdas1)
    temp = edate.strftime(strfmt)
    day = int(edate.strftime("%d"))
    # week 1 is
    if day < 8:
        temp = temp.replace("week", "w1")
    # week 2 is 8 through 14
    elif day < 15:
        temp = temp.replace("week", "w2")
    elif day < 22:
        temp = temp.replace("week", "w3")
    # week 4 is 22 through  28
    elif day < 29:
        temp = temp.replace("week", "w4")
    else:
        temp = temp.replace("week", "w5")
    return temp

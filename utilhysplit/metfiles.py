# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import datetime
from os import path
import pytz

"""
PRGMMR: Alice Crawford  ORG: ARL
This code written at the NOAA  Air Resources Laboratory
ABSTRACT: choosing met files for HYSPLIT control file
"""

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

   def set_ens_member(self, suffix):
       self.mstr += suffix

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
       print(dstart, duration, dstart+datetime.timedelta(hours=duration))
       print(newdate, runtime, newdate + datetime.timedelta(hours=runtime))
       mf = MetFiles(mstr,hours=24)
       files = []
       iii=0
       while not files:
          newdate = newdate - datetime.timedelta(hours=cycle_time * iii)
          files = mf.make_file_list(newdate, runtime)
          if cycle_time *iii > forecast_length: break
          iii += 1
       return files    

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
          if cycle_time *iii > forecast_length: break
          iii += 1

       files = weed_files(files,dstart,duration,self.metid,self.mstr)
       return process(files)   

def weed_files(metfiles,dstart,duration,metid,metstr):
    mhash = get_forecast_info(metid)
    metlist =[]
    #metstr = get_forecast_str(metid)
    done=False
    iii=0
    while not done:
          metf = metfiles[iii].split('/')[-1]
          metdir = metfiles[iii].replace(metf,'')
          #metstr = get_forecast_str(metid,metdir)
          mdate = datetime.datetime.strptime(metfiles[iii], metstr)       
          edate = mdate + datetime.timedelta(hours=mhash['forecast_length']) 
          metlist.append(metfiles[iii])
          if edate > dstart + datetime.timedelta(hours=duration):
             done = True  
          iii+=1
    return metlist




def get_forecast_info(metid):
    # TO DO. set forecast_length based on metid
    mhash = {}
    mhash['cycle_time'] = 6
    if 'gfs' in metid.lower():
        mhash['forecast_length'] = 84
        mhash['time_res'] = 3  #3 hour time resolution
    if 'gefs' in metid.lower():
        mhash['forecast_length'] = 84
    if 'nam' in metid.lower():
        mhash['forecast_length'] = 72
    return mhash

def get_forecast_str(metid, FCTDIR='/pub/forecast'):
    """Finds forecast meteorology data files
    dstart : datetime object of start date
    metdata : string (options: GEFS, GFS, GFS0p25, NAM, NAMAK, NAMHI)
    """
    if (metid.lower() == 'gfs') or (metid.lower() == 'gfs0p25'):
        met = metid.lower()
    elif(metid.lower() == 'nam'):
        met = 'namsf'
    elif(metid.lower() == 'namak'):
        met = 'namsf.AK'
    elif(metid.lower() == 'namhi'):
        met = 'namsf.HI'
    elif(metid.lower() == 'gefs'):
        met = 'gefs'
    #metnamefinal = 'No data found'
    #        metime = dtm
    metdir = FCTDIR + '/%Y%m%d/'
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
        if self.verbose:
            print("GETMET", sdate, edate, end_date, runtime, self.mdt)
        zzz = 0
        while not done:
            if "week" in self.strfmt:
                self.mdt = datetime.timedelta(hours=7 * 24)
                temp = parse_week(self.strfmt, edate)
            else:
                temp = edate.strftime(self.strfmt)
            edate = edate + self.mdt
            #if not path.isfile(temp):
            #    temp = temp.lower()
            if not path.isfile(temp):
                print("WARNING", temp, " meteorological file does not exist")
                #temp = self.altmet.makefilelist(edate, self.altmet.mdt)
                #print("REPLACE with", temp)
            #if temp != "None":
            elif temp not in nlist:
               nlist.append(temp)
            if edate > end_date:
                done = True
            if zzz > self.maxfiles:
                done = True
            zzz += 1
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


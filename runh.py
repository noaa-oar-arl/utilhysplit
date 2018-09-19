# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import numpy as np
import datetime
import time
from os import path, chdir
from subprocess import call
import pickle as pickle
import sys
from .process import is_process_running
from .process import ProcessList
from .hcontrol import HycsControl
from .hcontrol import Species
from .hcontrol import ConcGrid
from .hcontrol import NameList
from .hcontrol import writelanduse

"""
NAME: ashfall_base.py
UID: r100
PRGMMR: Alice Crawford  ORG: ARL  
This code written at the NOAA  Air Resources Laboratory
ABSTRACT: This code manages HYSPLIT runs for Hanford resuspension project.
CTYPE: source code

List of classes and functions in this file:

class RunParams - defines critical parameters for each run number.

-------------------------------------------------------------------------------
functions to create directory tree based on date and to get directory from date.
-------------------------------------------------------------------------------
function get_vlist
function date2dir
function dirtree (generator function)

-------------------------------------------------------------------------------
functions and classes to create source locations.
-------------------------------------------------------------------------------
function make_sourceid - creates string from lat lon point. 
class Source  - Helper class which stores information about a source.
function source_generator (generator function) - returns a Source object
function get_landuse_sources - returns list of tuples (lat, lon) of sources which have land
                               use of 2,3 or 8 (dryland cropland, irrigated cropland, shrubland)

-------------------------------------------------------------------------------
functions and classes to create CONTROL and SETUP files
-------------------------------------------------------------------------------
function add_psizes - defines ash particle sizes and densities to be used.
function make_setup - creates SETUP.CFG needed as input to  HYSPLIT.
function getmetfiles - returns directory and name of WRF ARL files which are needed for each run.
function add_pterm_grid - add an extra concentration grid for terminating particles.
class Station - helper class which contains information about location for concentration 
function stations_info - returns a Station object.

-------------------------------------------------------------------------------
Main functions
-------------------------------------------------------------------------------

function mult_run - main function for making multiple HYSPLIT runs.

"""
def source_generator(df, area=1, altitude=10):
    """df is a pandas dataframe with index being release date and column headers showing
       location of release and values are emission rates"""
    locs = df.columns.values
    for index, row in df.iterrows(): 
        for vals in zip(locs, row):
            rate = vals[1]
            center = vals[0]
            date = index
            yield  Source(center, altitude, area, rate, date)


class Station(object):  
    """Helper class to store information about a measurement station or location where we want to model concentrations"""
    def __init__(self, name='name', sample_start="00 00 00 00", stime=1, location=(0,0), cname='9999',
                 cspan=[10,10], cdiff=[0.1,0.1]):
        """name - name of measurement location or station
           sample_start - time that concentration sampling starts
           stime - sampling time
           location - tuple of (latitude, longitude)
           cname - name for cdump output file. """
        self.name = name
        self.sample_start = sample_start
        self.stime = stime
        self.location = location
        self.cspan = cspan
        self.cdiff = cdiff
        if cname =='9999':             
           self.cname = name
        else:
           self.cname = cname


def getdegrees(degrees, minutes, seconds):
    return degrees + minutes/60.0 + seconds/3600.00


def stations_info(name, cdiff, cspan):
    """function which will return a station object""" 
    if name.lower() == 'nd':
       return Station(name = 'nd',
                      sample_start = "00 00 00 00",
                      stime = 1,
                      location = (47.0,-101.0),
                      cdiff=cdiff,
                      cspan=cspan)


class RunParams(object):
    """RunParams is a class which defines parameters for input into HYSPLIT.
    """
    def __init__(self, num, met_type='wrf27',  topdirpath='./'):
       """num : integer - defines the run number associated with the run parameters
          met_type : string
          topdirpath : string
          Sets parameters according to the value of num.
          Creates a subdirectory of topdirpath/run(num)/ if it does not already exist. 
       """
       self.topdirpath =  topdirpath        
       self.met_type = met_type
       #adds subdirectory runN under the topdirpath.  
       self.topdirpath += 'run' +   str(num) + '/' 
       #creates directory if it does not already exist
       if not path.isdir(self.topdirpath):
            callstr = 'mkdir -p ' +  self.topdirpath  
            call(callstr, shell=True)        
   
       ##First set defaults for all the parameters.  
       ###Set where the concentration outputs are. 
       self.vert = True                                    #if True makes vertical line source up to 10m.

       ###Inputs for source_generator function. 
       ##Antelope valley stack height is about 183 Meters.
       ## Coyote Station is 152 meters.
       self.altitude = [183, 200]             #altitude to realease computational particles to.

       ##antelop valley stack height is 182.9 meters.

       ##values for CONTROL file 
       self.runtime = 1*24                       #time to run HYSPLIT for (hours).
       self.ztop = 10000                         #height of model top. 
       self.psizelist= [0]  #particle sizes in microns.
       self.emission_hrs = 1                     #how long emission of each particle size lasts
       self.nbptyp = '1'                         #number of sub-bins on the particle sizes.
       self.clevels = [100]                      #vertical levels of concentration grid. (m)
       self.rate  = 1                            #emission rate multiplier.
       ###Parameters for SETUP.CFG file
       self.initd = '0'         #Puffs
       self.numpar = '10'    #number of particles/ puffs to release per hour
       self.delt = '0'          #time step is automatically determined by HYSPLIT.
       self.maxpar = 100000      #maximum number of particles to emit.

       ##Change values depending on run number input. So each run number is
       ##associated with these inputs.
       if num == 1 or num==2:
           self.clevels = [20,50,100,200,250]   #use 4 vertical level.
           self.nbptyp  =  '1'              #use 1 bin for particle sizes.       
           self.numpar = '10000'            #use 100,000 particles
           self.maxpar = '10000'        
           cdiff = [0.05 , 0.05]                        #concentration grid resolution (degrees)
           cspan = [10 ,10]                           #concentration grid span.      (degrees)
           self.stations = [stations_info('nd', cdiff=cdiff, cspan=cspan)]         ##used to define concentration grid in the CONTROL file 
   
       #convert psize list of floats to list of strings.
       self.psizestr = list(map(str, self.psizelist))




##------------------------------------------------------------------------------------##
##------------------------------------------------------------------------------------##
#####get_vlist date2dir and dirtree are functions to create directory tree based on date and
#####to get directory from date.
##------------------------------------------------------------------------------------##
##------------------------------------------------------------------------------------##

def get_vlist(sdate, dhour):
       """This function is called by the dirtree generator to help generate a directory
       tree based on date.
       Input is sdate (datetime object) and dhour (integer specifying how many directories per hour in the tree).
       If dhour = -1 then directory tree only goes down to month.
       """
       if dhour < 24:
            list1 = [sdate.year, sdate.month, sdate.day, sdate.hour]
            list2 = ['y', 'm', 'd', 'h']
            list3 = [4,2,2,2]
            vlist = list(zip(list1, list2, list3))
       elif dhour >= 24:
            list1 = [sdate.year, sdate.month, sdate.day]
            list2 = ['y', 'm', 'd']
            list3 = [4,2,2]
            vlist = list(zip(list1, list2, list3))
       elif dhour == -1:
            list1 = [sdate.year, sdate.month]
            list2 = ['y', 'm']
            list3 = [4,2]
            vlist = list(zip(list1, list2, list3))
       return vlist

def date2dir(topdirpath, sdate, dhour=1, chkdir=False):
    """given a topdirpath(string), run_num (integer), sdate (datetime) 
       and dhour(int) returns directory"""
    finaldir = topdirpath
    if not path.isdir(finaldir) and chkdir:
        print('creating directory: ' , finaldir)
        callstr = 'mkdir -p ' +   finaldir
        call(callstr, shell=True)      
    for val in get_vlist(sdate, dhour):
        finaldir  += val[1] + str(val[0]).zfill(val[2]) + '/'
    if not path.isdir(finaldir) and chkdir:
        callstr = 'mkdir -p ' +   finaldir
        call(callstr, shell=True)      
    return finaldir


def dirtree(topdirpath,  sdate, edate,
            chkdir=True, verbose = True, maxiter=1e6, dhour = 1):
    """A generator function which generates the directories for the HYSPLIT output.
       A generator function (which uses yield instead of return) can be iterated over.
       If chkdir=True then will create day level directory if it does not exist already.
       Directory is topdirpath / y(year) / m(month) / d(day) / h(hour)/
       Starts creating directories at input sdate and  creates new directory for every date spaced
       by hours specified in dhour until reach edate or maxiter.
       The topdirpath should already exist
    """
    dt = datetime.timedelta(seconds = 60*60*dhour)
    done = False
    zzz = 0
    while not done:
        finaldir = topdirpath
        for val in get_vlist(sdate, dhour):
            finaldir  += val[1] + str(val[0]).zfill(val[2]) + '/'
            if not path.isdir(finaldir) and chkdir:
                callstr = 'mkdir -p ' +   finaldir
                call(callstr, shell=True)      
        zzz+=1
        sdate += dt
        if sdate > edate or zzz> maxiter:
           done = True
        yield finaldir    

##------------------------------------------------------------------------------------##
##------------------------------------------------------------------------------------##
###The following functions and classes define the sources used  
## make_sourceid, Source, source_generator 
##source_generator can easily be modified to produce different sources.
##------------------------------------------------------------------------------------##
##------------------------------------------------------------------------------------##


def make_sourceid(lat, lon):
        """Returns string for identifying lat lon point"""
        latstr = str(abs(int(round(lat*100))))
        lonstr = str(abs(int(round(lon*100))))
        return latstr + '_' + lonstr


class Source(object):
    """Helper class. source generator returns this type of object """

    def __init__(self, center, altitude, area=1, rate=1, date=None):
        ##A source line in the HYSPLIT CONTROL file would be 
        ##latitude, longitude, altitude, rate, area.
        ##For most applications the rate is the mass per hour released and it is multplied by 
        ##the mass per hour released which is specified in the definition of each pollutant species.
 
        self.center = center              #tuple (latitutde, longitude)
        self.area = area                  
        self.rate = rate
        self.altitude = altitude
        self.nid = make_sourceid(center[0], center[1])
        self.sdate = date
    
    def __str__(self):
        rstr = 'Source --- \n'
        rstr += self.sdate.strftime("%m/%d/%Y\n")
        rstr += 'Rate: ' + str(self.rate) + '\n'
        try:
            rstr += 'Altitude: ' + str(self.altitude) + '\n'
        except:
            rstr += 'Altitude: ' + str(self.altitude[0]) + ' to ' + str(self.altitude[1]) +  '\n'
        rstr += 'ID: ' + self.nid + '\n'
        return rstr

#def reverse_source_generator(rate=1, mult_area=1, altitude=50):
#    ##makes the stations the sources. Used for reverse runs.
#    names = ['heim', 'hval', 'hvol', 'gren']
#    for nm in names:
#        station = stations_info(nm)
#        point = station.location
#        yield Source(point, [point], altitude, mult_area, rate)
   


def add_psizes(control, emission_hrs, psizelist, pshape=-1, verbose=True):
        """Defines particle properties to be used.
        density is related to particle size with smaller particles generally being more dense.
        density chosen here informed by JGR: atmospheres 10.1002/2015JD023609 Figure 3.
        control is a HycsControl object.
        emission_hrs is list of hours to emit the pollutant
        psize list is a list of the particle sizes in microns
        pshape is a list of particle shapes. 
        if verbose True then will print out information. """
        prholist = []
        pnames = []
        for sz in psizelist:
            pnames.append('P' + str(int(sz*10)).zfill(3))
            if sz == 0:
               prholist.append(0)
            elif sz <= 10:
               prholist.append(2.500)
            elif  sz < 30:
               prholist.append(2.200)
            elif  sz < 100:
               prholist.append(1.750)
            else:
               prholist.append(1.250)

        for zzz in range(0,len(pnames)): 
            particle=Species(pnames[zzz], psize=psizelist[zzz],  rate=1, 
                                     duration=emission_hrs, density=prholist[zzz], 
                                     shape=pshape, wetdep=1)
            if verbose:
               print('Particle' , zzz) 
               print(particle)
               print('-----------------------------')
            control.add_species(particle)



def make_setup(rhash, setupname='SETUP.CFG', pardumpname='PARDUMP', wdir='./' ):
        """writes the SETUP file which is needed to run HYSPLIT.
           input rhash iu a RunParams object.
        """
        namelist={}
        namelist['delt'] =  rhash.delt           #fix time step.
        namelist['initd'] = rhash.initd          #particle or puff.
        namelist['maxdim'] = '1'           
        namelist['kmix0'] = '250'                #default value is 250. controls minimum mixing depth.
        namelist['kblt'] =  '2'                  #Use Kantha Clayson for vertical mixing. 
        namelist['kbls'] =  '1'
        namelist['numpar'] =  rhash.numpar       #number of particles/puffs to release per hour.
        namelist['maxpar'] = rhash.maxpar        #maximum number of particles/puffs to simulate

        ##The pardump file can be used to restart HYSPLIT if it crashes. Although the runs
        ##are short enough in this application that rather than restart a run, we would just re-run it.
        namelist['ndump'] = '168'                #create pardump file after this many hours.
        namelist['ncycl'] = '168'                #add to pardump file every ncycl hours.
        namelist['poutf'] = '"' + pardumpname + '"'     #name of pardump file

        ##The termination grid is set in the function add_pterm_grid
        #namelist['ptrm'] = '0'                   #terminate particle outside of concentration grids
        namelist['nbptyp'] = rhash.nbptyp        #number of sub bins for each particle size.
  
        ## THIS IS NO LONGER USED.
        #namelist['p10f'] = rhash.p10f            #determines mass flux relationship in p10f to use
        #namelist['ichem'] = rhash.ichem          #resuspension run.

        nl = NameList(setupname, working_directory=wdir)
        nl.add_n(namelist)
        nl.write()


def roundtime(dto):
    """rounds input datetime to day at 00 H"""
    #return datetime.datetime(dto.year, dto.month, dto.day, 0, 0)
    return datetime.datetime(dto.year, dto.month, dto.day, 0, 0)
    

def getmetfiles(sdate, runtime, verbose=False, warn_file = 'MetFileWarning.txt', met_type='wrf27',  
                mdir='./', altdir='./'):
    '''Input start date and run time and returns list of tuples - (meteorological file directory, meteorological files)
       to use for HYSPLIT run.
       INPUTS:
       sdate : start date (datetime object)
       runtime : in hours
       verbose : print out extra messages if True
       met_type : always WRF for this project
       mdir : directory where met files should be accessed from
       altdir : directory where met files may be stored. If they do not exist in mdir, then the routine will look in
                the altdir and copy them to the mdir.

       OUTPUT:
       List of tuples with (met dir, met filename)    
       If the met files cannot be found in the mdir or the altdir then
    '''
    verbose=True
    mdirbase = mdir.strip()
    if mdirbase[-1] != '/': mdirbase += '/'
    dt =datetime.timedelta(days=1)
    mfiles=[]
    mdirlist = []
    sdate = sdate.replace(tzinfo=None)
    if runtime < 0:
       runtime = abs(runtime)
       end_date = sdate
       sdate = (end_date - datetime.timedelta(hours=runtime))
    else:
       end_date = (sdate +  datetime.timedelta(hours=runtime))
    edate = sdate
    notdone=True
    if verbose: print('GETMET', sdate, edate, end_date, runtime)
    zzz=0
    while notdone:
        yr = edate.strftime("%Y")
        if met_type=='ERA5':
            fmt =      'ERA5_%Y%m.ARL'  
            temp =  edate.strftime(fmt)
            edate = edate + datetime.timedelta(hours=24)
        elif met_type=='wrf27':
            fmt = 'wrfout_d01_%Y%m%d.ARL'
            temp =  edate.strftime(fmt)
            mdir = mdirbase +  yr + '/'
            edate = edate + datetime.timedelta(hours=24)
        else: temp='none'
        if not path.isfile(mdir + temp):
            print('WARNING', mdir + temp,  ' meteorological file does not exist')  
            with open(warn_file, 'a') as fid:
                 fid.write('WARNING ' + mdir + temp + ' meteorological file does not exist\n') 
        else:
           mfiles.append(temp) 
           mdirlist.append(mdir)
           if verbose: print('Adding',  temp, ' meteorological file')  
        if verbose: print(edate, end_date)
        if edate > end_date:
            notdone=False 
        if zzz>12:
            notdone=False
            print('WARNING: more than 12 met files')
        zzz+=1
    #return mdir, mfiles
    return list(zip(mdirlist, mfiles))


def add_pterm_grid(control,iii): 
            """adds concentration grid for removing particles to the control file object.
               center is (64.1,-19.5) and span is 5 degrees so particles will be followed if they
               stay within the 61.6 to 66.6 degree latitude and -26.5 to -12.5 degree longitude"""
            cgrid = ConcGrid("junkname", levels=[1000], centerlat=64.1, centerlon=-19.5,
            latdiff=0.2, londiff=0.2, latspan=5, lonspan=14, outdir="./", outfile='cdump.all.' + iii,
            sampletype=0, interval=(1,0))
            control.add_cgrid(cgrid)




def sphu2relh(sphu, temp, pres ):
    """converts specific humidity, sphu, to relative humidity.
       esat is saturated vapor pressure, pres is pressure, temp is temperature (Kelvin) """
    esat = sat_vap_density(temp)
    
    #Pc=220.640  #hPa
    #Tc=647.096  #kelvin
    #C1 = -7.859
    #C2 = 1.844
    #C3 = -11.7866
    #C4 = 2.6807
    #C5 = -15.9618
    #C6 = 1.80123
    #theta = 1 - temp/Tc
    #M = Tc/T * (C1*theta + C2 * theta**1.5 + C3 * theta **3  + C4*theta**3.5 + C5*theta**4+C6*theta**7.5)
    #Pws = Pc * np.exp(M)
    return sphu*(pres) / (0.622*esat)

def sat_vap_density(t, type=1):
    """calculates saturation vapour pressure used to convert sphu to relh"""
    esat = np.exp(21.4-5351.0/t)
    return esat


class MultRun(object):
    def __init__(self, run_num, topdirpath='./'):
        from monet.models.hysplit import HYSPLIT
        self.run_num = run_num 
        self.rhash = RunParams(run_num) 
        self.model = HYSPLIT()              
        self.topdirpath = topdirpath
     
 
    def add_source_generator(self, dframe):
        self.sgenerator = source_generator(dframe, altitude=self.rhash.altitude)
        
    def add_sources(self):
        filera=[]
        zzz=0
        for source in self.sgenerator:
            sdate = source.sdate  #get the date from the source
            newdir = date2dir(self.topdirpath, sdate, chkdir=True) #now create the directory from the date.
            iii = source.nid
            filera.append(newdir + 'cdump' + '.' + iii)
            zzz+=1
            if zzz> 2: break
        print(filera)   
        self.model.open_files(filera) 


def concmerge(hdir, files2merge):
    """
    hdir: string
    files2merge: list of strings
    """
    mfile = 'mfiles.txt'
    ofile = 'merged.bin'
    with open(mfile, 'w') as fid:
         for fl in files2merge:
             fid.write(fl + '\n')
    callstr =  hdir + '/conmerge -i' + mfile + ' -o' + ofile 
    #subprocess.call(callstr, shell=True)
    return ofile


def create_plume(run_num, start_date, end_date,  dframe, 
             verbose=True, topdirpath = './',  
             hysplitdir = './'):
    rhash = RunParams(run_num, met_type = met_type, topdirpath=topdirpath)
    verbose=True 
    tdirpath = rhash.topdirpath
    zzz = 0 
    sgenerator = source_generator(dframe, altitude=rhash.altitude)
    pdir = None
    cnames = []
    for source in sgenerator:
        if verbose: print(source)
    #for newdir in dirtree(rhash.topdirpath,  start_date, end_date, dhour=dhour):  
        sdate = source.sdate  #get the date from the source
        newdir = date2dir(tdirpath, sdate, chkdir=False) #now create the directory from the date.
        if newdir != pdir: zzz=0
        #sdate = datetime.datetime.strptime(newdir, dfmt)       #get the date from the directory tree.
        iii = source.nid
        cnames.append('cdump.' +   iii)
        zzz+=1
        pdir = newdir
        if verbose: print(newdir, iii)
    ofile = concmerge(newdir, cnames) 
    return ofile 
 
def mult_run(run_num, start_date, end_date,  dframe, runh=0, process_max = 24, met_type='wrf27', logname='junklog',
             verbose=True, topdirpath = './', metdirpath='./', 
             hysplitdir = './', altdirpath='./',  adderr=False):
    """ Starts and keeps track of multiple HYSPLIT runs. 
        The first HYSPLIT run will be started at start_date. CONTROL, SETUP and output files from the run 
        are placed in a directory tree Yyear/month/day/hour/.
        The sources are determined by the source_generator function.
    

        run_num    : integer 
                     passed to RunParams to determine input values into HYSPLIT.
                     the runs created are also placed under the directory run(run_num).
        start_date : datetime.
                     First date to start a HYSPLIT for. 
        end_date   : datetime 
                     Last date to create a HYSPLIT run for.
        dframe     : dataframe 
                     to use as input into the source_generator.
        runh       : integer 
                     1 run check_ustar only to create USTAR directory tree with data from xtrct_stn.
                     2 do a dry run and create directory tree and CONTROL and SETUP files but do not run HYSPLIT.
                     0 is the default and will run check_ustar as well as HYSPLIT.
        process_max : number of hysc_std processes to keep running at the same time.
        met_type    : type of meteorological data to use.
        logname     : basename to use for the the logfiles
                      logname.endlog.txt writes time at which HYSPLIT runs finish
                      logname.sumlog.txt writes time at which HYSPLIT runs in a new directory start.
                      logname.current.txt write PID of HYSPLIT runs which are currently running
                      logname.stderr.txt writes stdout and stderr info if any HYSPLIT runs do not run properly.
        topdirpath  : path to write CONTROL, SETUP files and HYSPLIT output to. 
        metdirpath  : path to location of meteorological files.
        hysplitdir  : directory where the hysplit executable is 
        altdirpath  : path to where meteorological files are stored. If files are not found in metdirpath then getmetfiles will look in
                      altdirpath and copy  the files to metdirpath.
        Writes CONTROL and SETUP files for HYSPLIT. 
 
        Each run_num gets its own directory
    """ 

    ##RunParams class sets inputs based on the run_num. 
    rhash = RunParams(run_num, met_type = met_type, topdirpath=topdirpath)
    verbose=True 
    tdirpath = rhash.topdirpath
    endlogname     = tdirpath + logname + ".endlog.txt"        #when process finished running
    startlogname   = tdirpath + logname + ".sumlog.txt"        #when start running in a new directory.
    currentlogname = tdirpath + logname + ".current.txt"       #pid currently running
    stderrlog       = tdirpath + logname + ".stderr.txt"       #writes stdout and stderr for any runs which didn't finish properly.

    hysplit_version = 'hycs_std'                               #use standard version of hysplit
     
    pidtracking = ProcessList(hysplit_version)                 #set up object for keeping track of multiple HYSPLIT runs.

    landusedir = hysplitdir.replace('exec', 'bdyfiles')
    #sfile = 'sources.txt'
    #with open(sfile, 'w') as fid:
    #    fid.write('Sources ' + datetime.datetime.now().strftime("%Y%m%d %H:%M \n"))
    #    fid.write('Run number ' + str(run_num) + '\n') 
    #    for source in source_generator(altitude=rhash.altitude, rate=rhash.rate):
    #        strlat = "{0:.3f}".format(source.latlonlist[0][0])
    #        strlon = "{0:.3f}".format(source.latlonlist[0][1])
    #        fid.write(strlat + ' ' + strlon  + '\n')

    if runh == 0:
       dryrun = False
    elif runh == 1:
       dryrun = True

    dfmt = rhash.topdirpath  + 'y%Y/m%m/d%d/h%H/'
    dhour=1
    dfmt = rhash.topdirpath  + 'y%Y/m%m/d%d/h%H/'
    print(start_date , end_date) 
    
    zzz = 0 
    sgenerator = source_generator(dframe, altitude=rhash.altitude)
    pdir = None
    for source in sgenerator:
        if verbose: print(source)
    #for newdir in dirtree(rhash.topdirpath,  start_date, end_date, dhour=dhour):  
        sdate = source.sdate  #get the date from the source
        newdir = date2dir(rhash.topdirpath, sdate, chkdir=True) #now create the directory from the date.
        if newdir != pdir: zzz=0
        #sdate = datetime.datetime.strptime(newdir, dfmt)       #get the date from the directory tree.

        if verbose: print('------------------')
        if verbose: print('Directory', newdir)
        if verbose: print('DATE' , sdate)
  
        #Get the directory and filename for the meteorological data files that will be used for runs 
        #this directory, newdir.
        metfiles = getmetfiles(sdate, int(rhash.runtime), met_type=met_type, mdir=metdirpath, altdir=altdirpath) 


        with open(startlogname, 'a') as fid:
             fid.write(datetime.datetime.now().strftime("%Y/%m/%d %H:%M  ") + \
                       'Working on: ' + sdate.strftime("%Y/%m/%d %H:%M ") +  '\n')
        #Loop through each source location.         

        #for source in source_generator(altitude=rhash.altitude, rate=rhash.rate):
        #for source in sgenerator:
        cont=True
        if cont:
            iii = source.nid
            tempdir = newdir +  'temp' + str(zzz).zfill(4) + '/'

            #create the control file object.            
            control = HycsControl(fname='CONTROL.' + iii)
            #add information to the control file object. 
            control.add_duration(rhash.runtime) 
            control.add_sdate(sdate)
            control.add_ztop(rhash.ztop)
            control.add_vmotion(0)
            for mfile in metfiles:
                control.add_metfile(mfile[0], mfile[1])
            ##NOTE that the rate in the particle section is set to 1 because the rate is set when the location is specified.
            add_psizes(control, rhash.emission_hrs, rhash.psizelist, verbose=False)

            cfile_list = []
            ##define and add concentration grids to the control file object.
            for stn in rhash.stations:
                #cname2 = 'cdump.' + stn.cname + '.' + ustar_str + '.' +  iii
                ##no longer put ustar in cdump name
                cname2 = 'cdump.' +   iii
                cgrid = ConcGrid("junkname", levels=rhash.clevels, centerlat=stn.location[0], 
                                         centerlon=stn.location[1],
                                         latdiff=stn.cdiff[0], londiff=stn.cdiff[1], 
                                         latspan=stn.cspan[0], lonspan=stn.cspan[1], 
                                         outdir="./", outfile=cname2,
                                         sampletype=0, interval=(stn.stime,0), sample_start=stn.sample_start)
                control.add_cgrid(cgrid)
                cfile_list.append(cname2)
            ##add source lines to the control file object. 
            sloc = source.center
            ##The hourly release rate is set here, on the lines for location.
            if rhash.vert:  #if want to make a line source from ground to top altitude add a source at 0 elevation.
               control.add_location(latlon = (sloc[0], sloc[1]), alt= source.altitude[0], rate = source.rate, area = source.area)
               control.add_location(latlon = (sloc[0], sloc[1]), alt= source.altitude[1], rate = source.rate, area = source.area)
            else:
               control.add_location(latlon = (sloc[0], sloc[1]), alt= source.altitude, rate = source.rate, area = source.area)
            if dryrun:   #If no HYSPLIT run is going to be done then do not create temp directory.
               tempdir = newdir

            if not path.isdir(newdir):
                callstr = 'mkdir -p ' +  newdir
                if verbose: print(callstr)
                call(callstr, shell=True)      

            if not path.isdir(tempdir):
                callstr = 'mkdir -p ' +  tempdir
                if verbose: print(callstr)
                call(callstr, shell=True)      

            chdir(newdir)
           
            #write the CONTROL and SETUP file
            make_setup(rhash, 'SETUP.' + iii, pardumpname='PARDUMP', wdir=tempdir)
            if verbose: print('working dir' , newdir)
            if verbose: print('fname' , control.fname)
            control.rename(control.fname, working_directory=tempdir)
            control.write(annotate=False)
            writelanduse(landusedir, outdir=tempdir)

            rundescrip =  newdir + ' : ' + source.nid

            pid_done = False  #Loop so only max_process hycs_std runs are running simultaneously.
                              #Will run this waiting loop until number of hycs_std runs drops below max.
            min2sleep = 0.1   #How long to sleep in each cycle of the loop.
            while not pid_done:
                pidtracking.checkprocs(verbose=False, moveoutput=True)
                num_running = is_process_running(hysplit_version)
                if verbose: print("Number of hycs_std running is " , str(num_running))
                if num_running >=  process_max:
                   time.sleep(min2sleep * 60)
                else:
                   pid_done = True
                                 
            pid_done=False 
            callstr = [hysplitdir +  hysplit_version,  iii]
            if not dryrun:
                 if verbose: print('RUNNING IN :' , newdir)
                 if verbose: print('temporary directory :' , tempdir, zzz)
                 if verbose: print(callstr)
                 npid = pidtracking.startnew(callstr, wdir=tempdir, descrip = rundescrip)
                 fhash = pidtracking.updatelist(writelog=True, logname=currentlogname)      #update the pid tracking
                 with open(endlogname, "a") as fid:
                      for pkey in list(fhash.keys()):
                          fid.write(fhash[pkey] + '\n')
                 
            ##wait 30 seconds between starting each HYSPLIT run.
            #min2sleep=0.5
            #time.sleep(min2sleep * 60)
            pdir = newdir
            zzz+=1

    ###Loop to wait for last processes to finish.
    pid_done=False
    with open(endlogname, "a") as fid:
        fid.write('Out of main loop. Waiting for final runs to complete\n')
    if verbose: print('Out of main loop. Waiting for final runs to complete')
    while not pid_done and  not dryrun:
        numprocs = pidtracking.checkprocs(verbose=False)
        #print 'numprocs', numprocs
        if numprocs > 0:
           time.sleep(10)
           fhash = pidtracking.updatelist(writelog=True, logname=currentlogname)      #update the pid tracking
           with open(endlogname, "a") as fid:
                 for pkey in list(fhash.keys()):
                     fid.write(fhash[pkey] + '\n')
        else:
            fhash = pidtracking.updatelist(writelog=True, logname=currentlogname)      #update the pid tracking
            with open(endlogname,  "a") as fid:
                 for pkey in list(fhash.keys()):
                     fid.write(fhash[pkey] + '\n')
            pid_done=True

    #If any of the runs did not close with a 0 then the stdout was saved in an array and is printed out here.
    if pidtracking.err == []:
       print('Completed. No errors returned.')
    else:
       if not adderr:
           print('WARNING: Completed and errors returned from some runs. Check stderr log.')
       with open(stderrlog, 'a') as fid:
           for message in pidtracking.err:
                fid.write('descrip  '+  str(message[3]) + '\n')
                fid.write('pid  '+  str(message[2]) + '\n')
                fid.write('poll returncode  '+  str(message[4]) + '\n')
                fid.write('stdout '+  message[0] + '\n')
                fid.write('stderr '+  message[1] + '\n')
                fid.write('-----------------------------\n') 
                fid.write('-----------------------------\n') 

    for newdir in dirtree(rhash.topdirpath,  start_date, end_date, dhour=dhour):  
        #clean up any temporary cdump files that were left over
        chdir(newdir)
        rmtemp = True  
        ##remove cdump files which are created from the large concentration grid which is there to remove particles.
        if rmtemp:
           callstr = 'rm -f cdump.temp.*'
           call(callstr, shell=True)        
           

                
    return 1        





# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import numpy as np
import datetime
import os
from subprocess import call
from os import path
#from ashfall_base_iceland import RunParams
#from ashfall_base_iceland import stations_info
#from ashfall_base_iceland import dirtree
#from ashfall_base_iceland import source_generator
#from ashfall_base_iceland import reverse_source_generator
#from ashfall_base_iceland import check_ustar   
#from ashfall_base_iceland import make_sourceid  
import sys
import os
import cPickle as pickle
import pandas as pd
#from ashfall_base_iceland import RunParams
from arlhysplit.runh import date2dir

 
"""
NAME: iceland_datem.py
UID: r203
PYTHON 2.7 
PRGRMMR: Alice Crawford ORG: NOAA ARL
ABSTRACT: Calls c2datem to extract information from HYSPLIT output files. Adds information to the c2datem output.
CTYPE: source code

FUNCTIONS
writedatem_sh : writes a shell script to run c2datem on the cdump files (HYPSLIT output binary concentration files).
                The shell script will also have lines to concatenate all the output into one text file and add extra information
                to the end of each line (infor about particle size and vertical concentration level which is originally in the file name).

writedatem : writes a datem file which tells c2datem which positions and times to extract concentrations from the cdump file for.

"""


def writedatem_sh(cfile_list, mult='1e20', mdl='./',  ofile_list=None,
                  psizes=[1],  zlra=[1],  concat=True, concatfile='model.txt',
                  add2ra=['none']):
    """Writes a .sh file to run the c2datem program on the HYSPLIT output files  
       Writes separate file for each particle size. and vertical concentration level.
       psizes : a list of particle sizes to write (particle size indices 1..N)
       zlra   : a list of vertical concentration levels to write (indice 1..N)
       concat : concatenate all the files into one file with name specified by concatfile (model.txt is default).
                information about the particle size and vertical concentration level is added to the end of each line. 
       concatfile : name of file with all the data
       add2ra : extra information to add to each line. The extra information is added using sed.
    """
    with open('datem.sh', "w") as fid:
         fid.write('#!/bin/sh \n')
         fid.write('rm -f ' + concatfile + '\n')   #removes any existing model.txt file
         fid.write('MDL='+mdl+'\n')
         fid.write('mult='+str(mult)+'\n')
         if ofile_list == None:
             ofile_list=[]
             for cfile in cfile_list:
                 ofile_list.append(cfile.replace('cdump','model'))
         outfile_list = []
         for zl in zlra: 
             fid.write('zl='+str(zl) + '\n')
             iii=0
             for cfile in cfile_list:
                 for psz1 in psizes:
                    try:
                       psz = abs(int(psz1))
                    except:
                       psz1 = 1
                    else:  #if try does not raise an exception then this code is executed.
                        if zl==-1:
                           zstr = 'zn1'
                        else:
                           zstr = '.z' + str(zl)
                        outfile = ofile_list[iii]  + '.p' + str(int(psz)) + zstr + '.txt'
                        #print iii, outfile
                        outfile_list.append(outfile)
                    ##Following block writes line in shell script to run c2datem
                    ##the -h0 specifies to not write header lines.
                    #print('WRITING', cfile)
                    fid.write('$MDL/c2datem -n -h0 -i'+cfile.strip() + ' -mdatemfile.txt -o' + outfile + '  -c$mult -z$zl')
                    fid.write(' -p' + str(int(psz)))                                   #pollutant index select for multiple species
                    fid.write('\n')
                    if concat: 
                        temp = cfile.split('.')
                        if add2ra[0] != 'none':
                            a2l = ' ' + add2ra[iii]  + ' ' + str(psz) + ' ' + str(zl)
                            fid.write("sed 's/$/" + a2l + "/' " + outfile + " >> " + concatfile  + '\n')  #add info to end of line and add to file.
                 iii+=1
         fid.write("if [ ! -s model.txt ]\n") 
         fid.write("then\n")
         fid.write("rm -f model.txt\n")
         fid.write("fi\n")
    
    return outfile_list 

def frame2datem(dfile, df,  header_str='Header', writeover=True,\
                 cnames = ['date', 'duration', 'lat', 'lon', 'obs', 'vals', 'sid', 'altitude']):
    """converts a pandas dataframe with columns names by cnames (date, duration, lat, lon, obs, vals, sid, altitude)
       to a text file in datem format.
       date should be a datetime object.
       duration should be a string format HHMM (TODO- make this more flexible?)
       lat - latitude, float
       lon - longitude, float
       obs - value of observation, float
       vals - modeled value, float
       sid  - station id, int or string
       altitude - float """
    iii=0
    if writeover:
        with open(dfile, "w") as fid:
            fid.write(header_str + ' (obs then model) ' + '\n')
            fid.write('year mn dy shr dur(hhmm) LAT LON  ug/m2 ug/m2 site_id  height \n' )
    with open(dfile, "a") as fid:
          for index, row in df.iterrows():
              fid.write(row[cnames[0]].strftime('%Y %m %d %H%M') + ' ')
              fid.write(row[cnames[1]] + ' ' )
              fid.write("%8.3f  %8.3f" % (row[cnames[2]], row[cnames[3]]))
              fid.write("%8.4f  %8.4f " % (row[cnames[4]], row[cnames[5]] ))
              if isinstance(row[cnames[6]], int):
                 fid.write("%12i" % (row[cnames[6]]))
              elif isinstance(row[cnames[6]], str):
                 fid.write("%12s  " % (row[cnames[6]]))
              fid.write("%7.2f \n" % (row[cnames[7]]))
              #fid.write(str(row[cnames[7]]) + '\n')           

def writedatem(dfile, stationlocs, sample_start, sample_end, stime, height=' 10'):
    """writes a station datem file which has times for each station location.
       This file is used by c2datem to determine what concentration values to pull from the cdump files.
       stationlocs is a list of (lat,lon) tuples. 
 
       If the -z option in c2datem is set to -1 then the height 
       indicates which level will be used. It is the actual height in meters, not
       the index of the height lev
el. 

       outputs 1 in the measurement concentration and sourceid columns.
    """
    iii=0
    with open(dfile, "w") as fid:
        fid.write('DOE ASHFALL PROJECT\n')
        fid.write('year mn dy shr dur(hhmm) LAT LON g/m2  site_id  height \n' )
        for iii in range(0, len(stationlocs)):
            sdate = sample_start
            while sdate < sample_end:
                fid.write(sdate.strftime('%Y %m %d %H%M') + ' ')
                fid.write(str(stime[iii]) + '00 ' )
                fid.write("{:0.3f}".format(stationlocs[iii][0]) + ' ' + "{:0.3f}".format(stationlocs[iii][1]) + ' ')
                fid.write('1 ')
                fid.write('1 ')
                #fid.write(height + '\n')
                fid.write(height +  '\n')
                sdate  += datetime.timedelta(hours=stime[iii])


##-------------------------------------------------------------------------------------------##
##-------------------------------------------------------------------------------------------##
##-------------------------------------------------------------------------------------------##

#def datem(tdir, 



   





"""
NAME: ashfall_pandas.py
UID: r202
PRGMMR: Alice Crawford  ORG: ARL  
PYTHON 2.7
This code written at the NOAA  Air Resources Laboratory

ABSTRACT: The functions in this file collect data from c2datem output and store them in a pandas dataframe object.
CYTPE: source code

------------------------------------------------------------
The following functions are found in this file
mk_datem_pkl :  creates daily pickle files for c2datem output (which is found in each hour directory).
panda_daily  :  pickles a panda dataframe that contains datem data for 1 day (calls panda_conc for time frame of one day).
panda_conc   :  returns a panda Dataframe with datem data for time period sdate to edate (calls read_datem_file for time period of interest).
read_datem_file : reads datem file and returns a dataframe with all the information in the datem file.

"""

##------------------------------------------------------------------------------------##


def mk_datem_pkl(run_num, sdate, edate,  topdirpath, verbose=False, logfile='log.txt'):
    """Creates daily pickle file for c2datem output in each hour directory.
    Loop to call panda_daily so daily pickle files are created for time range sdate to edate.
    run_num    :  the run number to create the pickle file for.
    topdirpath : the directory path where the run_num runs are located
    sdate to edate is the daterange to create the pickle files for
    verbose :  print out information 
    logfile :  name of logfile which will be palced in the topdirpath. 
    """
    dt = datetime.timedelta(hours=24)
    sdate = datetime.datetime(sdate.year, sdate.month, sdate.day, 0,0)
    edate = datetime.datetime(edate.year, edate.month, edate.day, 0,0)
    done = False
    while not done:
        panda_daily(sdate, run_num=run_num, topdirpath=topdirpath, verbose=verbose, logfile=logfile)
        sdate += dt
        if sdate > edate:
           done = True


def panda_daily(sdate, run_num=2, verbose=1, topdirpath='./', pkl_name='conc_daily.pkl', logfile='log.txt'):
    '''creates a binary pickle file with datem output for one day (24 hour period).
    sdate      : day to create the pickle file for
    run_num    : the run number to create the pickle file for.
    topdirpath : the directory path where the run_num runs are located
    verbose    : print out information 
    pkl_name   : name of binary file to create
    logfile    : name of logfile which will be palced in the topdirpath. 
    ''' 
    verbose=False
    dt = datetime.timedelta(hours=24)
    dftot = panda_conc(sdate, sdate + dt, run_num=run_num, verbose=verbose, topdirpath=topdirpath, logfile=logfile)
    #vpi = np.where(dftot>0)
    #print('VPI')
    #print(vpi)
    outdir = date2dir(topdirpath,  sdate, dhour=24)  
    pickle.dump(dftot, open(outdir + 'conc_daily.pkl', "wb"))
    if verbose:
       print outdir
       print dftot
    print 'panda pkl done ' , sdate, outdir

def read_datem_file(fname, zlevs, pdict,sdate, dummy=False, verbose=False, \
    colra=['date','meas_lat', 'meas_lon', 'vals','sourceid','stationid', 'level','thickness','psize','sourcedate',  \
           ]):
    """ Reads a datem file and returns a dataframe with colums described by colra
       colra : 1st must be date, second must be lat, third must be longitude 
       fname : base name of datem file to read and get information from
       zlevs : 
       pdict :
       sdate :
       dummy : if dummy= True creates and returns a dataframe with correct columns with 0's. 
    """
    tdate = []
    #mvals = []
    #levra = []
    #thickra = []
    #psize = []
    #sdatera = []
    mlat = []
    mlon = []
    #sourceid=[]
    #stationid=[]
    #colra=['date','vals','sourceid','stationid', 'level','thickness','psize','sourcedate',  \
    #       'meas_lat', 'meas_lon']
    vhash = {}
    nhash = {}
    qhash = {}
    iii=7
    jjj=0
    for val in colra:
        vhash[val].append([])
        nhash[val] = iii
        qhash[iii] = val
        iii += 1
        jjj += 1
    ##the date takes up  0,1,2,3
    ##lat is 5th
    ##lon is 6th
    if not dummy:
        with open(fname, 'r') as fid:
             for line in fid:
                 temp = line.split() 
                 vhash[colra[1]].append(float(temp[5]))  #lat of measurement
                 vhash[colra[2]].append(float(temp[6]))  #lon of measurement
                 hh = int(temp[3][0:2])
                 try:
                     mm = int(temp[3][2:4])
                 except:
                     print 'ERR in read', temp[3] 
                     sys.exit()
                 ##if value is -1 that means no valid info on the cdump grid. Meas point may be off grid.
                 #if float(temp[7]) != -1:
                 vhash[colra[0]].append(datetime.datetime(int(temp[0]), int(temp[1]), int(temp[2]), hh, mm))
                 for val in vra:
                     iii = nhash[val]
                     vhash[val].append(float(temp[iii]))                     

        if vhash[colra[0]] != []:
            vra = []
            for col in colra:
                vra.append(vhash[col])
            tempzip = zip(*vra) 
            if not tempzip:
               print 'Did not read ', fname, ' file correctly. read_datem_file function exiting'
               sys.exit() 
            df = pd.DataFrame(tempzip, columns=colra)
            if verbose: print 'read_datem_file returning df', df
            return df
        else: dummy = True
    elif dummy:
        tdate = [0]
        mvals = [0]
        ustar = [0]
        sourceid = [0]
        levra = [0]
        thickra = [0]
        psize = [0]
        sdatera = [sdate]
        u10m = [0]
        v10m = [0]
        vwnd1 = [0]
        uwnd1 = [0]
        vwnd2 = [0]
        uwnd2 = [0]
        pres  =[0]
        prss1  =[0]
        prss2  =[0]
        relh  =[0]
        mlat = [0]
        mlon = [0]
        df = pd.DataFrame(zip(tdate, mvals, ustar, sourceid, levra, thickra, psize, sdatera, u10m, v10m, uwnd1, vwnd1, \
                              uwnd2, vwnd2, pres, prss1, prss2, relh, mlat, mlon),  \
                              columns=colra)
        print 'USTAR ZZZZZZZZZ' , df.ustar.unique()
        return df



def panda_conc(sdate, edate,  run_num=2, verbose=0, 
               topdirpath = './', logfile='log.txt'):
    """ returns a pandas Dataframe with colums 
      date (dates in datem output file)
      vals (concentration from datem output file)
      sourceid (ascii latitude_longitude of release location)
      ustar (friction velocity at release location and time)
      sourcedate (date of release) 
      psize (particle size in microns)
      level (top height of leve in meters)
      thickness (thickness of level in meters) 
      topdirpath is directory where the USTAR directory and the run(number) directory are found.

      sdate   : first day to return data for
      edate   : last day to return data for
      run_num    : the run number to create the pickle file for.
      topdirpath : the directory path where the run_num runs are located
      logfile :  name of logfile which will be palced in the topdirpath. 


      calls read_datem_file  function
      """
    #rhash = RunParams(run_num, met_type='WRF', topdirpath=topdirpath)
    maxiter=1e4
    zlevs = [50,100]
    iii=1                             #create dictionary showing particle sizes in microns
    pdict={}
    for psize in ['1']:
        pdict[str(iii)] = psize
        iii+=1
    #topdirpath = rhash.topdirpath
    done = False
    prev_dir = ''
    nnn=0
    odate = sdate
    while not done:
        if sdate.hour ==0 and verbose>=1:
           print 'working on' , sdate
        newdir = date2dir(topdirpath,  sdate, dhour=1) 
        
        outdir = date2dir(topdirpath,  sdate, dhour=24)  #directory where panda concentration pkl file will be written.
        os.chdir(newdir)
        if os.path.isfile('model.txt'):
            if verbose>=1: print 'model.txt exists in' , newdir
            #print 'model.txt exists in' , newdir
            df = read_datem_file('model.txt', zlevs, pdict,sdate)
            if verbose>=1: print df
            df = df[df.vals != 0]                     #remove zero values
            if nnn==0:
               dftot = df.copy()
            else:
               dftot = pd.merge(dftot, df, how='outer') 
            nnn+=1
        else:
            with open(topdirpath + logfile, 'a') as fid:
                fid.write('no model.txt file in ' +  newdir + '\n')


        sdate = sdate + datetime.timedelta(hours=1)
        if sdate >= edate or nnn>maxiter:
           done = True
    if nnn > 0:
        return dftot
    else:
        with open(topdirpath + logfile, 'a') as fid:
            fid.write('MESSAGE: no emissions for ' +  odate.strftime('%Y %m %d %H ') + '\n')
        print 'MESSAGE: no emissions for ' ,  odate.strftime('%Y %m %d %H ') , '\n'
        return read_datem_file('model.txt', zlevs, pdict,sdate, dummy=True)


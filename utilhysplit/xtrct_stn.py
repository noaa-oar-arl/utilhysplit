#!/n-home/alicec/anaconda/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import datetime
import string
import subprocess
from os import path

import numpy as np


"""
PRGMMR: Alice Crawford  ORG: ARL  
PYTHON 2.7
This code written at the NOAA  Air Resources Laboratory
UID: r102
CTYPE: source code
ABSTRACT: manages the xtrct_stn program and outputs.

CLASSES

Xtrct_Stn_File   represents a file produced by xtrct_stn (fortran utility program distributed with HYSPLIT) 

Provides methods to
1. write a shell script that will run xtrct_stn and call the shell script.
2. Read output file produced by xtrct_stn
"""
          

class Xtrct_Stn_File(object):
    """ python wrapper for fortran utility code xtrct_stn which is distributed with HYSPLIT. 
        readfile - method to read output file from xtrct_stn
        mksh - method to write shell script to run xtrct_stn and then run the shell script."""
       
    def __init__(self, fname, dir='./' , valra=[('U10M','01'),('V10M','01')], century=2000, verbose=True):
        """fname : name of file output by xtrct_stn
           valra : list of values that are in fname
           century : fname lists year by last two digits only. century is needed to process date.
           """
        self.fname = fname
        self.dir = dir
        self.valra = valra
        self.century = century
        self.valhash = {}
        self.dates = []
        for val in valra:
            self.valhash[val] = []



    def make_dummies(self, data_ra=[-999]):
        """instead of running xtrct_stn, write a dummy file like the one xtrct_stn would write.
           This is used for testing functions which use the xtrct_stn output."""
        sdate = datetime.datetime()
        dt = datetime.timedelta(hour=1)
        iii=1
        with open(self.fname) as fid:
             fid.write(str(iii) + sdate.strftime(" %y %m %d %h"))
             iii+=1 


    def readfile(self,  verbose=False):
        """Reads file and returns True if the file exists.
           returns False if file is not found"""
        if path.isfile(self.dir + self.fname):
            with open(self.dir + self.fname, "r") as fid:
                 try:
                    vdate1 = self.line2date(fid.readline())
                 except:
                    return False
                 try:
                    vdate2 = self.line2date(fid.readline())
                 except:
                    return False
            self.dt = vdate2 - vdate1
            if verbose:
               print self.dt
            with open(self.dir + self.fname, "r") as fid:
               for line in fid:
                   temp = line.split()
                   self.dates.append( datetime.datetime(self.century + int(temp[1]),
                                          int(temp[2]), int(temp[3]) , int(temp[4])))
                   jjj=5
                   for val in self.valra:
                       self.valhash[val].append(float(temp[jjj])) 
                       jjj+=1
            return True 
        else:
            return False


    def line2date(self, line):
        """get date from a line in the xtrct_stn output and return datetime object"""
        temp = line.strip().split()
        year = int(temp[1]) + self.century
        month = int(temp[2])
        day =   int(temp[3])
        hour =  int(temp[4])
        vdate = datetime.datetime(year, month, day , hour)
        return vdate

    def mksh(self, coord, xname, mdir, metname, interpolate=False,  multra=[], runsh=False, verbose=False,
             hysplitdir = '-99'):
        """writes shell script which will run xtrct_stn
           coord is lat lon in string format e.g. '34 -119'
           mdir is the directory of the meteorological file.
           metname is the name of the meteorological file."""
       
        valra = zip(*self.valra)[0]
        levra = zip(*self.valra)[1]
 
        if interpolate:
           interp = '1 \n'
        else:
           interp = '0 \n'
        if hysplitdir == '-99':
           hysplitdir = ''
        elif hysplitdir[-1] != '/':
           hysplitdir += '/'
        with open(xname, 'w') as fid:
             if verbose: print 'writing to ', xname
             fid.write(hysplitdir + 'xtrct_stn -i << EOF \n')
             fid.write(mdir + '\n')
             fid.write(metname + '\n')
             fid.write(str(len(valra)) + '\n')
             if len(levra) != len(valra):
                levra = ['01'] * len(valra)
             if multra==[]:
                ##this block will pick the amount to multiply the value by. 
                for zval in valra:
                    if zval in  ['U10M', 'V10M', 'T02M', 'PRSS','PBLH','SHGT','PRECIP','LHTF','DSWF','LHTF']:
                        multra.append('1')
                    elif zval in ['USTR']:  
                        multra.append('100')
                    elif zval in ['SPHU']:  
                        multra.append('1000')
                    elif zval in ['TPP1', 'TPPT', 'TPP6']:
                        multra.append('10000')
                    else:
                        multra.append('1')
                    #print zval, multra
             outra = zip(valra, levra, multra)
             for val in outra:
                 fid.write(val[0].strip() + ' ' + val[1].strip() + ' ' + val[2].strip() + '\n')
             fid.write(coord + '\n')
             fid.write(interp)
             fid.write(self.fname + '\n')
             fid.write('1\n')          #record number to start with
             fid.write('99999\n')      #record number to end with
             fid.write('EOF')
        if runsh:
           callstr = "chmod u+x " + './' + xname
           if verbose: print 'CALLING', callstr
           subprocess.call(callstr, shell=True) 

           callstr = './' + xname
           if verbose: print callstr
           subprocess.call(callstr, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return callstr            

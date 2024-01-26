import os
import sys
from utilvolc import runhelper


def lookup_anchorage(key):
    dhash = {}
    dhash['311360'] = 'Shishaldin'
    dhash['311070'] = 'Gareloi'
    dhash['311120'] = 'Great_Sitkin'
    

    if key in dhash.keys():
       return dhash[key]
    else:
       return key


if __name__ == "__main__":
   helper = runhelper.Helper
   helper.remove('index.html')
   ftpname = 'ftp.ssec.wisc.edu/pub/volcat/daily_so2_composites/'
   os.system("wget -P" + './'  + " " + ftpname)
   vidlist = []
   dlist = []
   with open('index.html','r') as fid:
        for line in fid.readlines():
            if 'VOLCAT' in line:
                aaa = line.split('"')
                fname = [x for x in aaa if 'VOLCAT' in x]
                fname = fname[0]
                os.system("wget -P " + " /pub/ECMWF/JPSS/VOLCAT/SO2Files/  {}/{}".format(ftpname,fname))

import os
import sys
from utilvolc import runhelper


def lookup_anchorage(key):
    dhash = {}
    dhash['311360'] = 'Shishaldin'

    if key in dhash.keys():
       return dhash[key]
    else:
       return key


if __name__ == "__main__":
   vaac='anchorage'
   #vid = '300260'
   #vname = 'Klyuchevskoy'

   helper = runhelper.Helper
   helper.remove('{}.1'.format(vaac))
   helper.remove('index.html')
   ftpname = 'ftp.ssec.wisc.edu/pub/volcat/events/{}'.format(vaac)
   os.system("wget -P" + './'  + " " + ftpname)
   vidlist = []
   dlist = []
   with open('{}.1'.format(vaac),'r') as fid:
        for line in fid.readlines():
            a = line.split('"')
            if len(a) > 8:
               vid = a[7][0:-1]
               if vid[0] == '3':
                   vidlist.append(vid)
 
   print(vidlist) 
   for vid in vidlist:
       vname = lookup_anchorage(vid)
       helper.remove('index.html')
       ftpname = 'ftp.ssec.wisc.edu/pub/volcat/events/{}/{}/netcdf/'.format(vaac,vid)
       os.system("wget -P" + './'  + " " + ftpname)
       print('-----------', ftpname)

       fpath = '/pub/ECMWF/JPSS/VOLCAT/Files/{}/'.format(vid)
       runhelper.make_dir(fpath,newdir=None,verbose=True)

       print('VID ', vid) 
       with open('index.html','r') as fid:
            print('OPEN INDEX.HTML')
            for line in fid.readlines():
                if 'VOLCAT' in line:
                    a = line.split('"')
                    print( a[7])
                    fpath = '/pub/ECMWF/JPSS/VOLCAT/Files/{}/{}'.format(vname,a[7])
                    if not os.path.isfile(fpath):
                        dlist.append(fpath)
                        os.system("wget -P " + " /pub/ECMWF/JPSS/VOLCAT/Files/{}/  {}{}".format(vname, ftpname,a[7]))
                    else:
                        print('not getting {}'.format(fpath))
   print('downloaded the following files')
   for d in dlist: print(d)

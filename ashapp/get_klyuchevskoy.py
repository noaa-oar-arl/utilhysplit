import os
import sys
from utilvolc import runhelper

if __name__ == "__main__":

   vaac='anchorage'
   vaac='tokyo'
   vid = '300260'
   vname = 'Klyuchevskoy'

   helper = runhelper.Helper
   helper.remove('index.html')

   ftpname = 'ftp.ssec.wisc.edu/pub/volcat/events/{}/{}/netcdf/'.format(vaac,vid)
   os.system("wget -P" + './'  + " " + ftpname)

   with open('index.html','r') as fid:
        for line in fid.readlines():
            if 'VOLCAT' in line:
                a = line.split('"')
                #print( a[7])
                fpath = '/pub/ECMWF/JPSS/VOLCAT/Files/{}/{}'.format(vname,a[7])
                if not os.path.isfile(fpath):
                    print('getting {}'.format(fpath))
                    os.system("wget -P " + " /pub/ECMWF/JPSS/VOLCAT/Files/{}/  {}{}".format(vname, ftpname,a[7]))
                    #sys.exit()
                else:
                    print('not getting {}'.format(fpath))

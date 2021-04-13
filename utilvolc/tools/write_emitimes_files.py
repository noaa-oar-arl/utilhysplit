# make_netcdf_area.py
# Makes area netcdf files - necessary for emitimes file creation
from utilvolc import write_emitimes as we
from datetime import datetime, timedelta
from glob import glob

wdir = '/hysplit-users/allisonr/Soufriere/Data_Insertion/'
vdir = '/pub/ECMWF/JPSS/VOLCAT/Vincent/Ash/'
vid = None
correct_parallax = False

#vfiles = '*'+vid+'*.nc'
vfiles = '*some_vars.nc*'
volclist1 = glob(vdir+vfiles)
volclist = volclist1[:]        # Can choose subset of volcat file list if desired
x = 0
while x < len(volclist):
    if vid in volclist[x]:
        start = volclist[x].find("_s") + len("_s")
        end = volclist[x].find("_v")
        date = volclist[x][start:end]
        dates = datetime.strptime(date, '%Y%j_%H%M%S')
        tmp = we.InsertVolcat(wdir, vdir, dates, vid=vid)
        done = tmp.write_emit(correct_parallax=correct_parallax)
        print(done)
    x += 1

# write_emitimes_files_L1L2nc.py
from utilvolc import write_emitimes as we
from datetime import datetime, timedelta
from glob import glob

wdir = '/hysplit-users/allisonr/Soufriere/Data_Insertion/'
vdir = '/pub/ECMWF/JPSS/VOLCAT/Vincent/Ash/'
vid = None
correct_parallax = False
decode_times = True

# vfiles = '*'+vid+'*.nc'
vfiles = '*2021104*some_vars.nc*'
volclist1 = glob(vdir+vfiles)
volclist = volclist1[:]        # Can choose subset of volcat file list if desired
x = 0
while x < len(volclist):
    start = volclist[x].find("Disk.") + len("Disk.")
    end = volclist[x].find(".some")
    date = volclist[x][start:end]
    dates = datetime.strptime(date, '%Y%j_%H%M%S')
    tmp = we.InsertVolcat(wdir, vdir, dates, vid=vid)
    done = tmp.write_emit(correct_parallax=correct_parallax, decode_times=decode_times)
    print(done)
    x += 1

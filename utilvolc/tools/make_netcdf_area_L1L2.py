# make_netcdf_area_L1L2.py
# Makes area netcdf files - necessary for emitimes file creation
from utilvolc import write_emitimes as we
from datetime import datetime, timedelta
from glob import glob

vdir = '/pub/ECMWF/JPSS/VOLCAT/Vincent/Ash/'
vid = None
correct_parallax = False
decode_times = True
write = True

vfiles = '*.nc'
volclist1 = glob(vdir+vfiles)
volclist = volclist1[:]        # Can choose subset of volcat file list if desired
dates = []
x = 0
while x < len(volclist):
    start = volclist[x].find("Disk.") + len("Disk.")
    end = volclist[x].find(".some")
    date = volclist[x][start:end]
    dates.append(datetime.strptime(date, '%Y%j_%H%M%S'))
    x += 1
print('Starting: ')
print(datetime.now())
j = len(volclist)
y = 0
while y <= len(dates):
    print(dates[y])
    tmp = we.InsertVolcat('./', vdir, dates[y], vid=vid)
    area = tmp.get_area(correct_parallax=correct_parallax, decode_times=decode_times, write=write)
    print(str(y)+' files of '+str(j)+' completed!')
    print(datetime.now())
    y += 1

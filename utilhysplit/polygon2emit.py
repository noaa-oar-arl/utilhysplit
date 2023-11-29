import numpy as np
import datetime
from utilhysplit import emitimes


   



def points2emit(ptlist,height,time,rate,area=0,heat=0,species=1,dt=5,filename='EMITIMES.txt'):
    """
    ptlist = list of tuples with (longitude, latitude) coordinates
    height = altitude in meters
    dt = duration of release in minutes
    """
    print('here')
    # make cycle_time at least one hour
    cycle_time = np.max([dt/60.0,1])
    efile = emitimes.EmiTimes(filename=filename)
    efile.add_cycle(time,cycle_time)

    # if record is added and no cycle found for it then print warning.
    duration = emitimes.minutes2duration(dt)
    ht = height
    rate = rate
    area = area
    heat = heat
    species= species
    for ppp in ptlist:
        lon = ppp[0]
        lat = ppp[1]
        if np.abs(lat>90.0):
           print('Warning latitude, longitude order may be mixed up')
        efile.add_record(time,duration,lat,lon,ht,rate,area,heat,species)
    efile.write_new(filename=filename)
    return efile


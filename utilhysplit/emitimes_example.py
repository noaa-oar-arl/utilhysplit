import numpy as np
import datetime
import emitimes

cycle_time=12
d1 = datetime.datetime(2022,1,1,0)
dt = datetime.timedelta(hours=cycle_time)
d2 = d1 + dt
d3 = d2 + dt

efile = emitimes.EmiTimes(filename='EMITIMES.txt')

# add emission cycles

efile.add_cycle(d1,cycle_time)
efile.add_cycle(d2,cycle_time)
efile.add_cycle(d3,cycle_time)

# add records
duration = "0010"
lat = 45
lon = 160
ht = 500
rate = 0.01
area = 0
heat = 0
species = 1

dtt = datetime.timedelta(hours=1)
d0 = d1
for iii in np.arange(0,3*cycle_time):
    efile.add_record(d0,duration,lat,lon,ht,rate,area,heat,species)
    d0 = d0 + dtt


efile.write_new(filename='EMIT.txt')





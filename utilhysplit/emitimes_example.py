import numpy as np
import datetime
import emitimes

def example2():
    cycle_time=12
    d1 = datetime.datetime(2022,1,1,0)
    dt = datetime.timedelta(hours=cycle_time)
    d2 = d1 + dt
    d3 = d2 + dt
    d4 = d1 + dt


    efile = emitimes.EmiTimes(filename='EMITIMES.txt')
    efile.add_cycle(d1,cycle_time)
    efile.add_cycle(d2,cycle_time)
    efile.add_cycle(d3,cycle_time)

    # emission cycles must be added in order.
    # if you try to add a cycle 
    # that starts before last cycle end time, you will get a message
    # and that cycle will not be added..
    efile.add_cycle(d4,cycle_time)

    # if record is added and no cycle found for it then print warning.
    d1 = datetime.datetime(2021,1,1,0)
    duration = "0010"
    lat = 45
    lon = 160
    ht = 500
    rate = 0.01
    area = 0
    heat = 0
    species=1
    efile.add_record(d1,duration,lat,lon,ht,rate,area,heat,species)

def example1():

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
    # records are automatically added to the correct emission cycle
    duration = "0010"
    lat = 45
    lon = 160
    ht = 500
    rate = 0.01
    area = 0
    heat = 0
    species = 1

    # here each record is for a different time.
    dtt = datetime.timedelta(hours=1)
    d0 = d1
    for iii in np.arange(0,3*cycle_time):
        efile.add_record(d0,duration,lat,lon,ht,rate,area,heat,species)
        d0 = d0 + dtt

    # this part will generate more lines in only one of the cycles
    # to 
    for iii in np.arange(0,cycle_time):
        ht = ht + 1
        efile.add_record(d2,duration,lat,lon,ht,rate,area,heat,species)
        d0 = d0 + dtt

    # write the file
    efile.write_new(filename='EMIT.txt')

    # Produces this output
    """
    YYYY MM DD HH DURATION(hhhh) #RECORDS #spnum P001 
    YYYY MM DD HH MM DURATION(hhmm) LAT LON HGT(m) RATE(/h) AREA(m2) HEAT(w)  
    2022 01 01 00  0012 12
    2022 01 01 00 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 01 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 02 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 03 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 04 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 05 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 06 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 07 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 08 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 09 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 10 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 11 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 12  0012 12
    2022 01 01 12 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 13 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 14 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 15 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 16 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 17 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 18 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 19 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 20 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 21 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 22 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 01 23 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 02 00  0012 12
    2022 01 02 00 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 02 01 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 02 02 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 02 03 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 02 04 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 02 05 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 02 06 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 02 07 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 02 08 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 02 09 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 02 10 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    2022 01 02 11 00 0010 45.0000 160.0000 500 1.00e-02 0.00e+00 0.00e+00 
    """

example2()

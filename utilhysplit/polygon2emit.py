import numpy as np
import datetime
import shapely.geometry as sgeo
from utilhysplit import emitimes
from utilhysplit import geotools



def polygon2emit(polygon, height, thickness, time, dy=0.1,filename='EMITIMES.txt',concentration=2e-3):
    """
    concentration should be given in g/m3.
    default is 2e-3 g/m3 (2 mg/m3)
    """
    if isinstance(polygon,list):
       polygon = sgeo.Polygon(polygon)
    lats = [[x[1] for x in polygon.exterior.coords]]
    lats = np.array(lats)
    meanlat = lats.mean()
    ratio = np.cos(meanlat*np.pi/180.0)
    dx = dy/ratio
    deg2meters = 111.1e3
    dym = dy*deg2meters
    dxm = dx*deg2meters*ratio
    area =dxm*dym
    points = geotools.polygon2points(polygon,dx=dx,dy=dy,res=0)
    points = [(p.x,p.y) for p in points]
    dt = 5          # 5 minute release
                    # approx 1000m thickness
    mass = concentration * area*thickness  # gives mass in grams
    rate = mass*60/dt # grams/hour
    species=1
    heat=0
    heightlist = [height-thickness, height]
    efile = points2emit(points,heightlist,time,rate,area,heat,species,dt,filename)
    return efile


def points2emit(ptlist,height,time,rate,area=0,heat=0,species=1,dt=5,filename='EMITIMES.txt'):
    """
    ptlist = list of tuples with (longitude, latitude) coordinates
    height = altitude in meters
    dt = duration of release in minutes
    """
    # make cycle_time at least one hour
    cycle_time = np.max([dt/60.0,1])
    efile = emitimes.EmiTimes(filename=filename)
    efile.add_cycle(time,cycle_time)

    if isinstance(height,(float,int)):
       height = [height]


    # if record is added and no cycle found for it then print warning.
    duration = emitimes.minutes2duration(dt)
    rate = rate
    area = area
    heat = heat
    species= species
    for ppp in ptlist:
        lon = ppp[0]
        lat = ppp[1]
        if np.abs(lat>90.0):
           print('Warning latitude, longitude order may be mixed up')
        for ht in height:
            efile.add_record(time,duration,lat,lon,ht,rate,area,heat,species)
    efile.write_new(filename=filename)
    return efile


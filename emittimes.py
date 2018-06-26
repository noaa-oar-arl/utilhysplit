#!/opt/Tools/anaconda3/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
from math import *
import sys 
from scipy.io import netcdf
from pylab import *
import numpy as np
import numpy.ma as ma
#import matplotlib.pyplot as plt
#import string
#import datetime
#from matplotlib.path import Path
#import shapely.geometry as sgeo
#from mpl_toolkits.basemap import Basemap
#from hysplit import *
import datetime
#from mpl_toolkits.basemap import Basemap

##Note. the place holder lat lon values in the CONTROL file must be within the meteorological grid-space.
##otherwise HYSPLIT will throw up an error and stop.

##For a line source, in the second line the emissions must be zero. Emission for the source defined in first line only:

"""class EmitTimes 
   function circpts
   function sat2emit"""



def circpts(x0, y0, r, dtheta=10):
    theta = np.array(list(range(0,360,dtheta))) * np.pi/180
    xp = r * np.cos(theta) + x0
    yp = r * np.sin(theta) + y0
    #plt.plot(xp,yp,'b.')
    #plt.show()
    #print 'ptx' , x0 , xp
    #print 'pty' , y0 , str(yp) 
    return np.array(xp), np.array(yp) 



class EmitLine(object):

   def __init__(self, date, duration, lat, lon, height, rate, area=0, heat=0):
       self.date = date
       self.duration = duration
       self.lat = lat
       self.lon = lon
       self.heigt=height
       self.rate = rate
       self.area = area
       self.heat = heat

   def __str__(self):
       returnstr = self.date.strftime("%Y %m %d %H %M ")
       returnstr += self.duration + ' '
       returnstr += str(self.lat) + ' '
       returnstr += str(self.lon) + ' '
       reuturnstr += str(self.height) + ' '
       reuturnstr += str(self.rate) + ' '
       reuturnstr += str(self.area) + ' '
       reuturnstr += str(self.heat) + ' \n'
       return returnstr

class EmitTimes(object):

   def __init__(self, filename='EMITIMES.txt'):
       self.filename=filename
       self.recordra=[]
       self.nrecs = 0

   def parse_header(self, header):
       temp = header.split()
       year = int(temp[0])
       month = int(temp[1])
       day = int(temp[2])
       hour = int(temp[3])
       #minute = int(temp[4])
       dhour = int(temp[4])
       nrecs = int(temp[5])
       self.sdate = datetime.datetime(year, month, day, hour)
       self.dt = datetime.timedelta(hours=dhour) 
       return nrecs

   def write_new(self, filename, duration='0024'):
       datestr = self.sdate.strftime('%Y %m %d %H ')
       with open(filename, 'w') as fid:
            fid.write(self.header_str())
            fid.write(datestr + ' ' + duration + ' ' + str(self.nrecs) + '\n')
             
            for record in self.recordra:
                fid.write(str(record)) 

   def parse_record(self, record):
       temp = record.split()
       year = int(temp[0])
       month = int(temp[1])
       day = int(temp[2])
       hour = int(temp[3])
       dhour = int(temp[5][0:2])
       dmin = int(temp[5][-2:])
       duration = temp[5]
       sdate = datetime.datetime(year, month, day, hour)
       lat = float(temp[6])
       lon = float(temp[7])
       ht = float(temp[8])
       rate = float(temp[9])
       try:
          area = float(temp[10])
       except:
          area=0
       try:
          heat = float(temp[11])
       except:
          heat=0
       return EmitLine(sdate, duration, lat, lon, ht, rate, area, heat) 


   def read_file(self):
       recordra=[]
       with open(self.filename, 'r') as fid:
            fid.readline()
            fid.readline()
            header = fid.readline()
            nrecs =  self.parse_header(header)
            for line in range(0,nrecs):
                temp = fid.readline()
                recordra.append(self.parse_record(temp))
       self.recordra.extend(recordra)
       self.nrecs = nrecs         

   def filter_records(self,llcrnr, urcrnr):
       iii=0
       rrr=[]
       for record in self.recordra:
           print('lat', record.lat, llcrnr[1], urcrnr[1])
           print('lon', record.lon, llcrnr[0], urcrnr[0])

           if record.lat < llcrnr[1] or record.lat > urcrnr[1]:
              rrr.append(iii)
              print('lat out')
           elif record.lon< llcrnr[0] or record.lon > urcrnr[0]:
              rrr.append(iii)
              print('lon out')
           iii+=1
       print(rrr)
       for iii in sorted(rrr, reverse=True):
           self.recordra.pop(iii)
           self.nrecs -= 1

   def header_str(self):
      returnval =  'YYYY MM DD HH    DURATION(hhhh) #RECORDS \n'
      returnval += 'YYYY MM DD HH MM DURATION(hhmm) LAT LON HGT(m) RATE(/h) AREA(m2) HEAT(w)  \n'
      return returnval 


   def write_header(self, sdate, duration='0100'):
       self.sdate = sdate
       self.datestr = sdate.strftime('%Y %m %d %H %M')
       self.duration = duration
       with open(self.filename, 'w') as emitfile:
            emitfile.writelines(self.header_str())

   def ash(self):
       pollnum = 4
       pollpercents=[0.008,0.068,0.254,0.670]
       return pollnum, pollpercents

   def cyl_source(self, lat, lon, radius=50000, dr=10000, pollnum=1, pollpercents=[1], dh=[0,10000], rate=1):
       """for default operation ash setting use pollnum=-1"""
       latlist = np.array([lat])
       lonlist = np.array([lon])
       emitfile = open(self.filename, 'a') 
      
       if pollnum==-1:
          pollnum, pollpercents=self.ash()
 
       numcircs = int(radius / dr)
       if radius%dr != 0:
          print(('Warning. Radius not evenly divisible by dr. Radius to ' , numcircs*dr))
       if numcircs == 0:
          print('Warning. dr larger than radius') 
          numcircs = 1
          dr = radius
       #map = Basemap(projection="sinu", lon_0=lon, resolution="c")
       #x0,y0 = list(map(lon,lat))
       #x1,y1 = list(map(lon+1,lat+1))
       dt = 1000
       #t0,t1 = list(map(x0+dt, y0+dt, inverse=True))
       dtheta = 40
       for n in range(1,numcircs+1):
           dth = dtheta / n   #keep arc length between points the same.
           if dth <=5:
              dth = 5
           xp, yp = circpts(x0,y0, dr*n, dtheta = dth)
           #lonp , latp = list(map(xp, yp, inverse=True))
           latlist= np.concatenate((latlist, latp))
           lonlist = np.concatenate((lonlist, lonp))
           #print 'CONCAT' , latlist
       heat = '0'
       duration = self.duration
      
       nrecords = 2 * pollnum * len(latlist)
       nnnrecords = 2*len(latlist)   #This is number of source lines that need to be written to CONTROL file.
       dtstr = self.sdate.strftime("%Y %m %d %H")
       emitfile.write(dtstr + '  9999  ' + str(nrecords) + '\n')

       ##Lat lon position (outer loop)
         ##bottom height, top height (next loop)
           ##particle (inner loop)
                       
       for ilat in range(len(latlist)):
           latstr = '%.3f' % latlist[ilat]
           lonstr = '%.3f' % lonlist[ilat]
           for iht in range(0,2):
               htstr = str(dh[iht]) 
               for ipolls in range(pollnum):
                   if iht==0:
                      rate_str = '%0.3e' % (pollpercents[ipolls] * rate / nrecords * 2 * pollnum)
                   elif iht==1:
                      rate_str = '0'
                   emitfile.write(self.datestr +  ' '  + duration + '  ' + latstr + ' ' + 
                          lonstr + ' ' + htstr + ' ' + rate_str + ' ' + 
                          '0   ' +  heat + '\n')
       #map2 = Basemap(projection="tmerc", lon_0=lonlist[0], lat_0= latlist[0], llcrnrlat=np.min(latlist)-1, urcrnrlat=np.max(latlist)+1,
       #                                llcrnrlon=np.min(lonlist)-1, urcrnrlon=np.max(lonlist)+1)
       #xp, yp = map2(lonlist, latlist)
       #map2.drawmeridians(np.arange(0,360,1), labels=[0,0,0,1])
       #map2.drawparallels(np.arange(0,360,1), labels=[1,0,0,0])
       #map2.drawcoastlines()
       #map2.plot(xp,yp, '-bo')
       #plt.show()
       emitfile.close()
       return nnnrecords


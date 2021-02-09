#cylsource.py
#Writes multiple line sources in cylinder around volcanic vent
#Makes figure of line source locations
import numpy as np
import datetime
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

##Note. the place holder lat lon values in the CONTROL file must be within the meteorological grid-space.
##otherwise HYSPLIT will throw up an error and stop.

##For a line source,  emission for the source defined in first line

"""class EmitTimes 
   function circpts"""

def circpts(x0, y0, r, dtheta=10):
    theta = np.array(range(0,360,dtheta)) * np.pi/180
    xp = r * np.cos(theta) + x0
    yp = r * np.sin(theta) + y0
    return np.array(xp), np.array(yp)

class EmitTimes(object):

   def __init__(self, filename='EMITIMES.txt'):
       """write EMITIMES file header"""

       self.filename = filename
       emitfile = open(filename, 'w')
       emitfile.writelines('YYYY MM DD HH DURATION(hhhh) #RECORDS \n'
                       'YYYY MM DD HH MM DURATION(hhmm) LAT LON HGT(m) RATE(/h) AREA(m2) HEAT(w)  \n')
       emitfile.close() 

   def ash(self):
       """ Default volcanic ash 4 particle divisions"""
       pollnum = 4
       pollpercents=[0.008,0.068,0.254,0.670]
       return pollnum, pollpercents

   def cyl_source(self, sdate, lat, lon, radius=30000, dr=10000, duration = '0010', pollnum=1, pollpercents=[1], dh=[0,10000], rate=1, heat='0',image = True):
       """ Creating cylinder volcanic ash source
       for default operation ash setting use pollnum=-1 """
       latlist = np.array([lat])
       lonlist = np.array([lon])
       emitfile = open(self.filename, 'a') 
       print('open file', self.filename) 
       if pollnum == -1:
          pollnum, pollpercents=self.ash()
 
       numcircs = int(radius / dr)
       print('Num circles', numcircs)
       if radius%dr != 0:
          print( 'Warning. Radius not evenly divisible by dr. Radius to ' , numcircs*dr)
       if numcircs == 0:
          print( 'Warning. dr larger than radius' )
          numcircs = 1
          dr = radius
       map = Basemap(projection="sinu", lon_0=lon, resolution="c")
       x0, y0 = map(lon,lat)
       x1, y1 = map(lon+1,lat+1)
       dt = 1000
       t0,t1 = map(x0+dt, y0+dt, inverse=True)
       dtheta = 40
       for n in range(1,numcircs+1):
           dth = int(dtheta / n)   #keep arc length between points the same.
           if dth <=5:
              dth = 5
           xp, yp = circpts(x0,y0, dr*n, dtheta = dth)
           lonp, latp = map(xp, yp, inverse=True)
           latlist = np.concatenate((latlist, latp))
           lonlist = np.concatenate((lonlist, lonp))
           
       nrecords = 2 * pollnum * len(latlist)
       dtstr = sdate.strftime("%Y %m %d %H")
       emitfile.write(dtstr+'  '+str(duration)+' '+ str(nrecords)+'\n')
       
       datestr = sdate.strftime('%Y %m %d %H %M')
                       
       for ilat in range(len(latlist)):
           latstr = '%.3f' % latlist[ilat]
           lonstr = '%.3f' % lonlist[ilat]
           for iht in range(0,2):
               htstr = str(dh[iht]) 
               for ipolls in range(pollnum):
                   if iht==0:
                      # divide the amount of mass 
                      rate_str = '%0.3e' % (pollpercents[ipolls] * rate / nrecords * 2 * pollnum)
                   elif iht==1:
                      rate_str = '0.0'
                   emitfile.write(datestr+  ' '  +duration+ '  ' +latstr+ ' ' +lonstr+ ' ' +htstr+' '+rate_str+' 0.0 ' + heat+'\n')
       if image == True:
            map2 = Basemap(projection="tmerc", lon_0=lonlist[0], lat_0= latlist[0], llcrnrlat=np.min(latlist)-1, urcrnrlat=np.max(latlist)+1,llcrnrlon=np.min(lonlist)-1, urcrnrlon=np.max(lonlist)+1)
            xp, yp = map2(lonlist, latlist)
            map2.drawmeridians(np.arange(0,360,1), labels=[0,0,0,1])
            map2.drawparallels(np.arange(0,360,1), labels=[1,0,0,0])
            map2.drawcoastlines()
            map2.plot(xp,yp, '-bo')
            plt.show()
       emitfile.close()

#Allison updates (making code able to do umbrella plumes)
   def calc_cyl(self, lat, lon, radius = 30000, dr = 10000):
       """ Calculate the source points in concentric circles from volcanic source
       Inputs: latitude, longitude, radius (m), dr (radius should be evenly divisible by dr)
       Returns list of latitudes and longitudes for point sources """ 
       latlist = np.array([lat])
       lonlist = np.array([lon])
       #Determine number of circles around volcanic vent
       numcircs = int(radius / dr)
       print('Num circles', numcircs)
       if radius%dr != 0:
           print( 'Warning. Radius not evenly divisible by dr. Radius to ' , numcircs*dr)
       if numcircs == 0:
           print( 'Warning. dr larger than radius' )
           numcircs = 1
           dr = radius
       #Calculating lats and lons in circle
       map = Basemap(projection="sinu", lon_0=lon, resolution="c")
       x0, y0 = map(lon,lat)
       dt = 1000
       dtheta = 40
       for n in range(1,numcircs+1):
           dth = int(dtheta / n)   #keep arc length between points the same.
           if dth <=5:
               dth = 5
           xp, yp = circpts(x0, y0, dr*n, dtheta = dth)
           lonp, latp = map(xp, yp, inverse=True)
           latlist = np.concatenate((latlist, latp))
           lonlist = np.concatenate((lonlist, lonp))
       return latlist, lonlist

   def calc_nrecs(self, latlist, pollnum = 1, umbrella = 1):
       """ Calculates the number of records. 
       If umbrella == 2, column divided into 2, emissions at upper half only
       If umbrella == 3, column divided into 3, emissions at upper third only
       nrecs = 2 * umbrella * pollnum * len(latlist)
       Inputs: latlist, pollnum, umbrella
       For default volcanic ash (4 particle sizes) use pollnum = -1"""

       if pollnum == -1:
           pollnum, pollpercents=self.ash()
       nrecords = (umbrella + 1) * pollnum * len(latlist)
       return nrecords

   def write_data(self, sdate, latlist, lonlist, nrecords, pollnum, duration = '0010', pollpercents = [1], height = [0,10000], rate = 1):
       """ Opens the emittimes file and writes data in emitimes file. 
       Inputs: latlist, lonlist (use calc_cyl), nrecords, duration (string: '0010' = 10 min), 
       pollnum (use calc_nrecs), pollpercents, rate, 
       height (str list that should account for umbrella heights)"""

       emitfile = open(self.filename, 'a') 
       print('open file', self.filename) 
       #Writing first line of data
       dtstr = sdate.strftime('%Y %m %d %H')
       emitfile.write(dtstr+'  '+str(duration)+' '+ str(nrecords)+'\n')
       datestr = sdate.strftime('%Y %m %d %H %M')
       #Writing records
       for ilat in range(len(latlist)):
           latstr = '%.3f' % latlist[ilat]
           lonstr = '%.3f' % lonlist[ilat]
           h = 0
           tt = len(height)
           if tt < 2:
               print('Warning: need at least two heights in list')
               print('Defaulting to height = [0,10000]')
               height = [0,10000]
               tt = len(height)
           while h < len(height):
               hstr = str(height[h])
               for ip in range(pollnum):
                   if h == tt-2:
                       # divide the amount of mass 
                       rate_str = '%0.3e' % (pollpercents[ip] * rate / nrecords * 2 * pollnum)
                   else:
                       rate_str = '0.0'
                   emitfile.write(datestr+' '+str(duration)+' '+latstr+' '+lonstr+' '+hstr+' '+rate_str+' 0.0 0.0 \n')
               h += 1
       emitfile.close()                      
           
def example(sdate, duration, lat, lon, altitude, filename, radius, dr):
    efile = EmitTimes(filename = filename)
    efile.cyl_source(sdate, lat, lon, radius = radius, dr = dr, pollnum = 1, dh = altitude)


### EXAMPLE HOW TO CALL
#dt = datetime.datetime(2019,6,21,18)
#duration = 0100
#lat = 48.292
#lon = 153.385
#altitude = [551,12000]

#example(dt, duration, lat, lon, altitude, filename, radius, dr)









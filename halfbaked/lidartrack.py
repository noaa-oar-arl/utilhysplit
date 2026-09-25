#!/n-home/alicec/anaconda/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
from math import *

import numpy as np
import shapely.geometry as sgeo
from pylab import *

##notes - fundamental sampling resolution of caliop lidar is 30 meters vertical and 333 meters horizontal.
##06/08/15 modified to read pardump binary file (rather than ascii file created by par2asc. 
##Also  modified to use pandas dataframes rather than the particle class.

#classes: LidarTrack


def getproj(coords):
    """returns a list which contains the distance of each point from the beginning of the line.
       coords should be array of (lon,lat) points"""
    #map = Basemap(projection="sinu", lon_0=0,  resolution='c')
    #x , y = map(zip(*coords)[0], zip(*coords)[1])
    #xy = zip(x,y)
    d1 = [xy[0]]
    pr = [0]
    cntr = 0
    l1old = 0
    length = 0
    for pt in xy[1:]:
        d1.append(pt) 
        l1 = sgeo.linestring.LineString(coordinates = d1)
        length += l1.length
        pr.append(length)
        #pr.append(length/1000.0)
        if cntr % 1000 ==0:
          print 'working on point ' , cntr , ' out of ' , len(xy[1:]) , 'distance ' , length/1000.0 , ' km'  , l1.length / 1000.0
        d1 = [pt]
        cntr +=1

    ##This method is slower than above  
    #ln = sgeo.linestring.LineString(coordinates=xy)
    #for pt in xy[1:]:
    #    if cntr % 1000 ==0:
    #       print 'working on point ' , cntr , ' out of ' , len(xy[1:])
    #    ln = sgeo.linestring.LineString(coordinates=xy)
    #    pnt = sgeo.Point(pt[0],pt[1])
    #    pr.append(ln.project(sgeo.Point(pnt)))
    #    cntr +=1

    return pr   

def find_dist(coords, pnt, bdist=0):
    """returns distance from beginning of line to that point. Coords and pnt are in lat lon.
       distance returneds is in km.""" 
    #map = Basemap(projection="sinu", lon_0=0,  resolution='c')
    #x , y = map(zip(*coords)[0], zip(*coords)[1])
    #xp , yp = map(pnt[0], pnt[1])
    #xy = zip(x,y)   
 
    ln = sgeo.linestring.LineString(coordinates = xy)
    pnt = sgeo.Point(xp, yp) 
    dist = ln.project(sgeo.Point(pnt))

    if bdist ==1:
        ##This section computes the distance of the point to the nearest point on the line
        pnt1 = ln.interpolate(dist)   
        ln2 = sgeo.linestring.LineString([pnt, pnt1])
        bffr_dist = ln2.length
        return dist/1000.0 , bffr_dist/1000.0
    else:
        return dist / 1000.0


class LidarTrack(object):
    """A class which can represent the path of a space based lidar.
       the line attribute is a shapely sgeo linestring object of the lat-lon coordinates of the path.
       the line2 attribute is  a sgeo linestring object defined by the lat-lon coordinates with a sinusoidal projection done by basemap.
       the line_plus attribute is a shapely object which is the line plus a buffer. 
       testpoint is a method which tests if a lat lon point is within the line_plus.
       dist2pnt is a method which returns the distance from a sgeo point to the line2 object.
    """
 
    def __init__(self, coords, bffr=0.05):
        #self.drange  =  drange       
        self.coords = coords               #list of lat/lon coordinates
        self.line = sgeo.linestring.LineString(coordinates = coords)     #line given by lat-lon
        #print 'BUFFER VALUE' , bffr
        map = Basemap(projection="sinu", lon_0=0,  resolution='c')
        x , y = map(zip(*coords)[0], zip(*coords)[1])
        self.map = map
        self.line2 = sgeo.linestring.LineString(coordinates = zip(x,y))  #line give by projected coordinated
        self.line_plus = self.line.buffer(bffr)

    def range(self):
        lon  = self.line.xy[0]
        lat  = self.line.xy[1]
        rlon = [min(lon), max(lon)]
        rlat = [min(lat), max(lat)]
        return rlon, rlat


    def plot(self, map=''):
        lonmin = -90
        lonmax = 0
        latmin = -80
        latmax = 30 
        #if map == '':
        #   self.map = Basemap(projection='cea', lon_0=-72, lat_0=041, llcrnrlon=lonmin, urcrnrlon=lonmax , llcrnrlat=latmin, urcrnrlat=latmax)
        #   self.map.drawlsmask(land_color="#E0E0E0", ocean_color='white', lakes=True)
        #else:
        #   self.map = map  
        #x , y  = self.map(self.line.xy[0].tolist(), self.line.xy[1].tolist())
        #mkr = '-mo' 
        #self.map.plot(x, y, mkr, markersize=1)
        #map.drawlsmask(land_color="#E0E0E0", ocean_color='white', lakes=True)

    def testpoint(self, lon, lat):
        pnt = sgeo.Point(lon, lat) 
        test = pnt.intersects(self.line_plus)
        return test

    def dist2pnt(self, pnt): 
        """
        pnt : shapely.geometry  Point object
        """
        x , y = self.map(pnt.x, pnt.y)
        dist = self.line2.project(sgeo.Point(x,y))
        return dist 

#!/n-home/alicec/anaconda/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import datetime
import string
import sys
from math import *

import matplotlib as mpl
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
#import ashgeo
#from pyvolcat.modis import ModisNC
#from pyhysplit.hysplit import Pardump
#from mytools import sortxyz
import pandas as pd
import shapely.geometry as sgeo
from matplotlib.path import Path
#from mpl_toolkits.basemap import Basemap
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
from scipy.io import netcdf

from utilhysplit import ashgeo

##notes - fundamental sampling resolution of caliop lidar is 30 meters vertical and 333 meters horizontal.
##06/08/15 modified to read pardump binary file (rather than ascii file created by par2asc. 
##Also  modified to use pandas dataframes rather than the particle class.

class InvGroup:
    def __init__(self,key,sorti,htbnds,latbnds,lonbnds):
        self.key = key
        self.sorti = sorti
        self.htbnds = htbnds
        self.latbnds = latbnds
        self.lonbnds = lonbnds


#def plot_geom(chull,clr='-r'):
    #for geom in chull:
        #x,y = ashgeo.plotpoly(geom)
        #plt.plot(x,y,clr)

def plot_pars(pdump,mass2,vloc=None):
    
    alpha = 10
    ch, ep = obs2concavehull(mass2,thresh=0.1,alpha=alpha)
    
    val = 'lon'
    valy = 'lat'
    d1 = pd.to_datetime(mass2.time.values)
    
    print('plotting all particle locations {}'.format(d1.strftime("%Y %m/%d %H:%M UTC")))
    begin = pdump[pdump['date']==d1]
    plt.plot(begin[val],begin[valy],'r.',MarkerSize=0.1)
    if vloc: plt.plot(vloc[0],vloc[1],'b^')
    plt.show()

    # second plot plots only particles that are inside the concave hull.
    sorti = match2obs(pdump,mass2,thresh=0.1,alpha=alpha)
    plt.contourf(mass2.longitude, mass2.latitude, mass2.values,cmap='Blues')
    abegin = begin[begin['sorti'].isin(sorti)]
    plt.plot(abegin[val],abegin[valy],'k.',MarkerSize=1)
    # works if ch is a multi-polygon object.
    #try:
    #    plot_geom(ch)
    #except:
    #    x,y = ashgeo.plotpoly(ch)
    #    plt.plot(x,y,'-r')
    #if vloc: plt.plot(vloc[0],vloc[1],'r^')
    #plt.show()
    return sorti

def obs2concavehull(obs,thresh=0.1,alpha=1):
    """
    INPUTS
    obs : xarray
    RETURN
    chull : shapely.geometry.polygon.Polygon
    edge_points : list of numpy arrays.
    """
    longitude = obs.longitude.values.flatten()
    latitude = obs.latitude.values.flatten()
    mval = obs.values.flatten()
    vpi = np.where(mval > thresh)
    #mpts = ashgeo.make_multi(longitude[vpi],latitude[vpi])
    #chull, edge_points = ashgeo.concave_hull(mpts,alpha=alpha)
    return chull, edge_points

def roundtime(time):
    time = pd.to_datetime(time)
    hour = time.hour
    if time.minute > 30: hour+=1
    return datetime.datetime(time.year, time.month, time.day, hour, 0)


def match2obs(df,obs,thresh=0.2, alpha=1):
    # returns sorti values of particles which are within the concave
    # hull defined by the observations (mass)
    gdf = match2obs_sub(df,obs,thresh=0.2,alpha=alpha)
    matchedi = gdf[gdf['test']].sorti.values
    return matchedi


def method1_histogram(dfnew, log=True):
    binx = len(dfnew.time_elapsed.unique())-1
    #biny = np.ceil(len(dfnew.ht.unique())/100.0)
    biny = np.ceil(np.max(dfnew.ht)/100.0) 
    if log: cb = plt.hist2d(dfnew.time_elapsed, dfnew.ht, bins=[binx,biny],norm=mpl.colors.LogNorm())
    else: cb = plt.hist2d(dfnew.time_elapsed, dfnew.ht, bins=[binx,biny])
    plt.colorbar(cb[3])
    

    



class ParMatched:

    def __init__(self,pardump, obs):
        # obs should be one time period.
 
        self.df = pardump
        self.obs = obs
        self.get_matchedi()

    def get_matchedi(self):
        matchedi = match2obs(self.df,self.obs,thresh=0.1)
        self.matchedi = matchedi

    def method1(self):
        # dfnew is a pandas dataframe in which
        # the first location of the matched particles is stored.
        
        sortused = []
        sortnew = []
        datelist = []
        dflist = []
        for date in self.df.date.unique():
            new = self.df[self.df['date'] == date]
            # only keep points which were matched
            anew = new[new['sorti'].isin(self.matchedi)]
            # take out points which were matched in previous dates
            anew = anew[~anew['sorti'].isin(sortused)]
            # add sorti values to sort used list.
            sortused.extend(anew['sorti'])
            dflist.append(anew)  
            datelist.append(date)
        dfnew = pd.concat(dflist,axis=0)
        def changedate(x):
            dt = x - datelist[0]
            return dt.days*24.0 + dt.seconds/3600.0
        dfnew['time_elapsed'] = dfnew['date'].apply(changedate) 
        self.dfmatched = dfnew
        return dfnew

    def plot_matched_dist(self):
        return -1

class ParInv:

    def __init__(self,df,d0):
        self.df = df
        self.d0 = d0
        self.grouphash = {}

    def add_group(self,groupkey, htbnds, latbnds, lonbnds):
        # find particles at d0 which are within the bounds.
        # return list of sorti for particles within those bounds.
        temp = self.df[self.df['date'] == self.d0]
        if htbnds:
            temp = temp[temp['ht'] <= htbnds[1]]
            temp = temp[temp['ht'] >= htbnds[0]]
        if latbnds:
            temp = temp[temp['lat'] <= latbnds[1]]
            temp = temp[temp['lat'] >= latnds[0]]
        if lonbnds:
            temp = temp[temp['lon'] <= lonbnds[1]]
            temp = temp[temp['lon'] >= lonbnds[0]]
        sorti = temp.sorti.unique()
        self.grouphash[groupkey] = InvGroup(groupkey,sorti,htbnds,latbnds,lonbnds)
        return sorti

    def return_group(self,groupkey,date):
        # pandas data frame with particles within the group at time date.
        sorti = self.grouphash[groupkey].sorti 
        temp = self.df[self.df['date'] == date]
        temp = temp[temp.sorti.isin(sorti)]
        return temp
         

def get_slice(df, t1, t2, val='ht', comparison=['ge','le']):
    """
    returns dataframe with points between t1 and t2 for
    column with name val.
    """
    if comparison[0] == 'ge':
        df2 = df[df[val]>=t1]
    else:
        df2 = df[df[val]>t1]
    if comparison[1] == 'le':
        df2 = df2[df2[val]<=t2]
    else:
        df2 = df2[df2[val]<t2]
    return df2


class ParGroup(object):
    """group of particle positions at the same time.
       setlidar method finds particles in the group which are close to a lidar track.
       plotlidar method plots particle position as function of altitude and distance along
                 lidar track.
    """

    def __init__(self,time, parframe):
        self.time = time
        self.df = parframe                    #A dataframe with particle positions / information
        self.dtfmt = "%Y%m%d%H%M"
        self.showfmt = "%Y %m %d %H:%M UTC"
        temp1 = parframe.min(axis=0)
        temp2 = parframe.max(axis=0)
        self.rlat = [temp1.lat , temp2.lat]
        self.rlon = [temp1.lon , temp2.lon]
        self.rht =  [temp1.ht , temp2.ht]
        #self.rlat = [par.lat,par.lat]
        #self.rlon = [par.lon,par.lon]
        #self.rht = [par.ht, par.ht] 
        self.lidar = []

    def hthist(self, fignum=2, bins=100, ttl=''):
        # histogram of heights in the dataframe.
        figh = plt.figure(fignum)
        ax = figh.add_subplot(1,1,1)
        plt.hist(self.df['ht'] , bins=bins, normed=1)
        plt.savefig(ttl + 'hist.jpg')


    def setlidar(self, lidartrack, prnt=0):
        """input is a LidarTrack object"""
        
        # projection
        pmap = Basemap(projection="sinu", lon_0=0, resolution='c')
        fmap = lambda x: sgeo.Point(pmap(x.x, x.y))
        fintersects= lambda x: x.intersects(lidartrack.line_plus)

        fdist =    lambda x: lidartrack.line2.project(x)         #input should be point in projected coordinates
        fnearpnt = lambda x: lidartrack.line2.interpolate(x)     #input should be distance from beginning of line.
        fldist =   lambda x , y: sgeo.linestring.LineString([x,y]).length 

        self.df['pnt'] = map(sgeo.Point, zip(self.df['lon'], self.df['lat']))  #lat lon coordinates
        self.df['pnt2'] = self.df['pnt'].apply(fmap)                             #projected coordinates
        self.df['blidar'] = self.df['pnt'].apply(fintersects)                    #true or false

        #print 'finding lidar group ' , len(r['pnt'].tolist())

        lidar_group = self.df.groupby(['blidar'], axis=0)                          
        lidar = lidar_group.get_group(True)
        #print 'found lidar group ' , len(lidar['pnt'].tolist())
        lidar['ldist'] = lidar['pnt2'].apply(fdist)                 #distance from the beginning of the line.
        lidar['lpnt']  = lidar['ldist'].apply(fnearpnt)              #closest point on the lidar line.
        lidar['ld']    = lidar['pnt2'].combine(lidar['lpnt'], func=fldist)  #this is the distance from pt and line

        ##could probably drop the point objects now.
        #self.df.drop(['pnt', 'pnt2'] , inplace=True, axis=1)

        self.lidar = lidar
        if prnt ==1:
           fid = open('par_lidar.txt','w')
           indx=0
           for pnt in lidar['pnt']:
               fid.write(str(pnt.x) + ' ' + str(pnt.y) + ' ' + str(lidar['ht'][indx]) + ' ' +  str(lidar['sorti'][indx]) + '\n') 
           #              + str(lidar['ldist'][indx]) + '\n') 
               indx += 1
           fid.close() 
        self.df.drop(['pnt', 'pnt2'] , inplace=True, axis=1)
        return lidar

    def plotlidar(self, fignum=10, xlocs=[], lidartrack=[], ttl='', ax=[], sv=1, prtcls=pd.DataFrame(), hrange=[]):
        fs=42
        #rcParams['figure.figsize']=9,10
        rcParams['figure.figsize']=9,9
        rcParams['axes.labelsize']=fs 
        rcParams['xtick.labelsize']=fs
        rcParams['ytick.labelsize']=fs
        fig = plt.figure(fignum)
        ax = fig.add_subplot(1,1,1)
        plot(self.lidar['ldist'],self.lidar['ht']/1000.0, '.', color='#CC33FF', markersize=4)
        if not prtcls.empty:
           plot(prtcls['ldist'],prtcls['ht']/1000.0, '+', color='#000066', markersize=10)

        if xlocs != []:
           dxlocs = []
           xlabels= []
           i = 0
           for loc in  xlocs:
               pntloc = sgeo.Point(loc[0],loc[1])
               dxlocs.append(lidartrack.dist2pnt(pntloc))    
               #print 'dx' , dxlocs   
               #dxlocs.append(find_dist(lidartrack.coords, loc))
               if i%2 == 0:
                  xlabels.append(str(loc[1]) + '\n' + str(loc[0]))
               else:
                  xlabels.append('')
               i+=1
               #print 'labels' , xlabels
           ax.set_xticks(dxlocs)
           ax.set_xticklabels(xlabels, fontsize=fs)     
        ax.set_ylabel('height (km)')
        ax.set_xlabel('Lat / Lon')
        ax.grid(b=True, which='both', color='0.6', linestyle='-', linewidth=2)
        if hrange != []:
           ax.set_ylim(hrange[0], hrange[1])
        #plt.title('HYSPLIT' + self.time.strftime(self.showfmt))
        ttl.replace(' ', '')
        #plt.subplots.adjust_left(0.15)
        plt.tight_layout()
        plt.savefig(ttl+'.jpg')
        return ax

    def add_par(self, parframe):
        self.df = pd.concat([self.df, parframe], axis=0)
        #update min and max lat and lon
        temp1 = parframe.min(axis=0)
        temp2 = parframe.max(axis=0)
        self.rlat = [temp1.lat , temp2.lat]
        self.rlon = [temp1.lon , temp2.lon]
        self.rht =  [temp1.ht , temp2.ht]
        #print 'RANGE of ParGroup' , self.rlat , self.rlon


    def __str__(self):
        """needs to be modified to use dataframe"""
        str1 = 'TIME ' + self.time.strftime(self.dtfmt) + "\n"
        for par in self.df:
            str1 += str(par) + "\n"
        return str1

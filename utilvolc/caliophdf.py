#!/n-home/alicec/anaconda/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
from math import *
import sys 
from pylab import *
import numpy as np
import numpy.ma as ma
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import imshow
from matplotlib.ticker import MultipleLocator
import datetime
from pyhdf.SD import SD, SDC
from pyhdf.HDF import *
from pyhdf.VS import *
import pytz
from utils import cmap
#from pyhysplit.par2ascplot import *
from mytools import roundtime
#from kernely import kernel2d
from scipy import ndimage
from mytools import crop_contour_plot
from lidartrack import LidarTrack
import pandas as pd

#classes LidarFile (lidar_file)
#        LidarFileList (lidar_file_list)
#        CalipsoL1


class LidarFile(object):
   """name of file containing calipso data and date range of the file.
   data in the file is a calipso_L1 object. 
   Can also have a partimes object associated with info from a pardump file associated."""
      
   def __init__(self,fname, drange, dir='./'):
       self.fname = fname
       
       self.drange = drange
       self.dir = dir
       self.data =[]
       self.pardump =[]
       self.lidartrack =[]

   def __str__(self):
       returnval = self.dir + '\n' +  self.fname + '    ' + self.drange[0].strftime("%Y-%m-%d T%H:%M") + ' '
       returnval += self.drange[1].strftime("%Y-%m-%d T%H:%M") + '  \n '
       return returnval

   def cropdata(self,xrange=[29.09,60.23], yrange=[0,25]):
       """Only crops height / altitude data.
          If crop latitudes then run into problems making sure the particles are plotted correctly from saved file.
       """
       if self.data ==[]:
          self.get_data()
       trange = [0,0]
       print 'CROPPINT lidar data. SHAPE original' , self.data.mab532.shape
       ##crops by latitude values.  commented out due to problems making compatible with saved particle positions.
       #vp = np.where(np.logical_and(self.data.lat > xrange[0], self.data.lat < xrange[1]))
       #self.data.lat = self.data.lat[vp]
       #self.data.lon = self.data.lon[vp]
       #self.data.time = self.data.time[vp]
       #self.data.mab532 = self.data.mab532[vp[0][0]:vp[0][-1],]
       #self.data.mab532 = self.data.mab532[vp,:]
       #shp = self.data.mab532.shape
       #self.data.mab532 = self.data.mab532.reshape(shp[1],shp[2])

       vp = np.where(np.logical_and(self.data.altitudes > yrange[0], self.data.altitudes < yrange[1]))
       self.data.altitudes = self.data.altitudes[vp]  
 
       #self.data.mab532 = self.data.mab532[:,vp[0][0]:vp[0][-1]]
       self.data.mab532 = self.data.mab532[:,vp]
       shp = self.data.mab532.shape
       self.data.mab532 = self.data.mab532.reshape(shp[0],shp[2])
       print 'SHAPE final' , self.data.mab532.shape , 'lat range' , np.min(self.data.lat), np.max(self.data.lat)


   def get_data(self):
       self.data = CalipsoL1(self.dir + self.fname)

   def get_lidartrack(self, bffr=0.01):
       print 'LIDAR bffr', bffr
       self.lidartrack = LidarTrack(self.data.get_coords(), self.drange, bffr=bffr)

   def  pshow(self):
       if self.data ==[]:
          self.get_data()
       self.data.plot532mesh(xaxis='time')

   def compare(self, pardump, show=0, ttl='', outpt='jpg', bffr=0.01):
       if self.data ==[]:
          self.get_data()
       if self.pardump ==[]:
          self.load_pardump(pardump)
       if self.lidartrack ==[]:
          self.get_lidartrack(bffr=bffr)
       if ttl=='':
          ttl=self.fname.replace('hdf','')
          if '.netcdf' in pardump or '.nc' in pardump:
              ttl += '.modis'
          else:
              temp = pardump[0].split('/')
              temp = temp[-1]
              ttl += temp.replace('txt','')
       prtcls = self.pardump.setlidar(self.lidartrack)
       print 'got particles. plotting figure'
       ax, xlocs = self.data.plot532mesh(prtcls=prtcls, ttl=ttl, outpt=outpt)
       ax, xlocs = self.data.plot532mesh(fignum=11 ,ttl=ttl, outpt=outpt)
       if not show:
          plt.close("all")

   def check_file(self, pardump, ldr=1, hrange=[], range=[], type='stere', bffr=0.01):
       """input a pardump file or netcdf (modis retrieval) 
          to compare to. creates plot which plots particle position and lidar track position.
          purpose is to see if lidar track co-incides with HYSPLIT output.
          if ldr=0 then will just plot particle position and no lidar file.
          hrange and range can be used to set the height range and the lat lon range of the plot.
          If they are [] then the ranges will be automatically set.
          range should be a list [lonmin, lonmax, latmin, latmax]
          The type of projection can be set to cyl (default) or stere."""
       if self.data ==[]:
          self.get_data()
       if self.pardump ==[]:
          self.load_pardump(pardump)
       if self.lidartrack ==[] and ldr==1:
          self.get_lidartrack(bffr=bffr)


      
       ##decides whether file is a netcdf (which means MODIS retrieval) or
       ##a HYSPLIT pardump file.
       if '.netcdf' in pardump[0] or '.nc' in pardump[0]:
           typstr='modis.'
           temp = pardump[0].split('/')
           temp = temp[-1]
           typstr = temp.replace('netcdf','') + 'modis'
       elif 'pardump' in pardump[0].lower():
           temp = pardump[0].split('/')
           temp = temp[-1]
           typstr = temp.replace('txt','')
       else:
           typstr = ''
       ##ttl is the title of the plot.
       ttl=self.fname.replace('hdf','') +  typstr
       if ldr !=1: 
          ttl = typstr
          #ttl=ttl.replace('CAL_LID_L1-ValStage1-V3-01.','parfile') 
          #ttl=ttl.replace('ZN_Subset','') 
          #ttl=ttl.replace('ZD_Subset','') 

       print 'typstr' , typstr , ttl, pardump[0] , 'pardump' , pardump
       self.pardump.range()
       if ldr==1:
          self.pardump.plot(drange=self.pardump.drange, lidartrack=self.lidartrack, ttl=ttl, range=range, hrange=hrange, type=type)
       else:
          self.pardump.plot(drange=self.pardump.drange, lidartrack=[], ttl=ttl, range=range, hrange=hrange, type=type)
       plt.close("all")


   def load_pardump(self, pardump):
       if type(pardump) == str:
          pardump = [pardump] 
       rd1 = roundtime(self.drange[0])
       rd2 = roundtime(self.drange[1])
       eruption = ParTimes()
       print pardump , type(pardump)
       if '.netcdf' in pardump[0] or '.nc' in pardump[0]:
           eruption.readmodis(pardump[0], rd1)
       elif '.txt' in pardump[0]:
           eruption.readpardump(pardump[0], drange = [rd1, rd2])
       else:
           for pd in pardump:
               print 'reading ' , pd
               eruption.readpardump_binary(pd, drange = [rd1, rd2])
       self.pardump =  eruption 


class LidarFileList:
   """reads a text file with list of calipso data files and determines date range from the filename"""   
 
   def __init__(self, fname='', dir='./'):
        """can initialize with name of single hdf file or with a text file with list of hdf files"""
        self.fn_list = []                               #list of lidar_file objects
        if fname != '':
           if fname[-4:] == '.hdf':
              self.add_hdf(fname, dir=dir)
           else:
              self.get_files(fname, dir=dir)

   def __str__(self):
         returnval = 'Number of files: ' + str(len(self.fn_list)) + '\n'
         for fn in self.fn_list:
             returnval += str(fn)
         return returnval

   def check_files(self, pardump, ldr=1, range=[], hrange=[], bffr=0.01):
        for fn in self.fn_list:
            fn.check_file(pardump, ldr=ldr, range=range, hrange=hrange, bffr=bffr)

   def compare(self, pardump, show=0, ttl='', bffr=0.01):
        for fn in self.fn_list:
            fn.compare(pardump, ttl=ttl, bffr=bffr)


   def add_hdf(self, hdf_name, drange=[], dir='./'):
     
        if drange==[]:                        #determine date range from file-name
            temp = hdf_name.split('.')
            tempdate = temp[1].replace('ZN_Subset','')
            tempdate = tempdate.replace('ZD_Subset','')
            min2 = float(tempdate[-2:])
            tempdate = tempdate[:-3]
            min1 = float(tempdate[-2:])
            if min2 > min1:
               dt = datetime.timedelta(seconds= (min2-min1) * 60)
            else:      
               dt = datetime.timedelta(seconds= min2 * 60)
            dt1 = datetime.datetime.strptime(tempdate, "%Y-%m-%dT%H-%M")
            dt2 = dt1 + dt
            #dt1 = dt2 
            drange = [dt1, dt2]
        self.fn_list.append(LidarFile(hdf_name , drange, dir=dir))
   
   def pshow(self):
        for fn in self.fn_list:
            fn.pshow()     

   def get_files(self, fname, dir='./'):
        """input text file with list of hdf calipso data file names"""
        fid = open(dir + fname, 'r')
        test = True
        while test:
          wl = fid.readline()
          if not wl:
             break
          self.add_hdf(wl.strip(), dir=dir)


def calipso_time2dt(timestamp):
    """Convert float which is seconds since Jan 1, 00:00 , 1993 to date."""
    d1 = datetime.datetime(1993, 1, 1, 0 ,0)
    dt = datetime.timedelta(seconds=timestamp)
    return d1 + dt

class CalipsoL2(object):
      """reads data from a calipso L1 file"""
      def __init__(self, fname, verbose=0, datatype='532'):
            self.fname = fname
            self.missing = -9999  
            self.get_vdata(verbose=1)
            self.mab532 = self.get_sd(verbose, datatype='532')
            #self.mabcr = self.get_sd(verbose, datatype='color ratio')

      def get_vdata(self, verbose=0):
            """retrieves altitude data from metadata"""
            fid = HDF(self.fname)
            vs = fid.vstart()
            vd = vs.attach('metadata', write=0)
            alt = vd.field('Lidar_Data_Altitudes')
            vd.setfields('Lidar_Data_Altitudes')
            altdata = vd.read(nRec=1)
            if verbose:
                print vd.attrinfo()
                print vd.fieldinfo()
                print 'INFO' , alt.attrinfo()
                print type(altdata) , len(altdata) , len(altdata[0][0])
                print altdata[0][0][0] , altdata[0][0][-1]
            self.altitudes = np.array(altdata[0][0][:])
            vd.detach()
            fid.close()
  
      def get_sd(self, verbose=1, datatype='532'):
            """retrieves data from SD (scientific data) set."""
            dset = SD(self.fname,SDC.READ)
            a = dset.datasets()
            self.attributes  = dset.attributes()
            verbose = 1
            if verbose:
               print 'ATTRIBUTES'
               print self.attributes 
               print a
               for key in self.attributes.keys():
                    print key 
                    #temp = dset.select(key)
                    #print temp.info()
            print 'Done'
            self.lat = dset.select('Latitude')[:,:]
            self.lon = dset.select('Longitude')[:,:]
            #self.lat = self.lat.reshape(self.lat.shape[0])
            #self.lon = self.lon.reshape(self.lon.shape[0])
            print 'latitude' , self.lat.shape 
            print self.lat
            print 'longitude', self.lon.shape
              
            self.time = dset.select('Profile_Time')[:,:]
            #ztime = dset.select('Profile_UTC_Time')[:,:]
            if datatype == '532':
               sd_data = dset.select('Particulate_Depolarization_Ratio_Profile_532')
            elif datatype == '1064':
               sd_data = dset.select('Attenuated_Backscatter_1064')
            elif datatype == 'depolarization':
               sd_data = dset.select('Depoloarization_Gain_Ratio_532')
            elif datatype == 'color ratio':
               sd_data = dset.select('Particulate_Color_Ratio')
               #sd_data = dset.select('Integrated_Attenuated_Total_Color_Ratio')
               #sd_data = dset.select('Integrated_Volume_Depolarization_Ratio')
               #sd_data = dset.select('Integrated_Attenuated_Backscatter_1064')
            ma_sd_data = ma.masked_equal(sd_data[:,:] , -9999)
            #self.mab532 = ma.masked_equal(b532[:,:], -9999)
            #self.b532= b532[:,:]
            if verbose:
                print '\n \n DIR scientific dataset ' , datatype, dir(sd_data)
                print sd_data.attributes() , '\n'
                #print 'Number of time locations' , len(time)
            print ma_sd_data.shape
            plt.imshow(sd_data[:,:])
            plt.show()
            return ma_sd_data

        
class Caliop(object):
      """reads data from a caliop hdf data file.  This is the base class."""

      def __init__(self, fname, verbose=0):
            self.fname = fname
            self.missing = -9999  
            self.get_vdata(verbose)

      def describe(self):
          self.get_vdata(verbose=1)
          self.get_sd(verbose=1, datatype='')
 
      def get_distance(self):
            self.dist = np.array(proj(self.get_coords()))
 
      def get_coords(self, verbose=0):
          coords = zip(self.lon.tolist(), self.lat.tolist())
          print 'COORDS' , coords[0]
          return coords

      def get_vdata(self, verbose=0):
            """retrieves altitude data from metadata"""
            fid = HDF(self.fname)
            vs = fid.vstart()
            vd = vs.attach('metadata', write=0)
            alt = vd.field('Lidar_Data_Altitudes')
            vd.setfields('Lidar_Data_Altitudes')
            altdata = vd.read(nRec=1)
            if verbose:
                print 'MetaData '
                print '-------------------------------'
                print vd.attrinfo()
                print vd.fieldinfo()
                print 'INFO' , alt.attrinfo()
                print type(altdata) , len(altdata) , len(altdata[0][0])
                print altdata[0][0][0] , altdata[0][0][-1]
                print '-------------------------------'
     
            self.altitudes = np.array(altdata[0][0][:])
            vd.detach()
            fid.close()
        
      def get_sd(self, verbose=1, datatype=''):
            """retrieves data from SD (scientific data) set."""
            dset = SD(self.fname,SDC.READ)
            a = dset.datasets()
            self.attributes  = dset.attributes()
            verbose = 1
            if verbose:
               print 'Scientific datasets ATTRIBUTES'
               print '-------------------------------'
               print self.attributes 
               print a
               for key in self.attributes.keys():
                    print key 
                    #temp = dset.select(key)
                    #print temp.info()
               print '-------------------------------'
            #self.lat = dset.select('Latitude')[:,:]
            #self.lon = dset.select('Longitude')[:,:]
            #self.lat = self.lat.reshape(self.lat.shape[0])
            #self.lon = self.lon.reshape(self.lon.shape[0])
            #print self.lat.shape
            #print self.lon.shape
            #self.time = dset.select('Profile_Time')[:,:]
            #ztime = dset.select('Profile_UTC_Time')[:,:]
            if datatype != '':
                sd_data = dset.select(datatype)
                ma_sd_data = ma.masked_equal(sd_data[:,:] , -9999)
                if verbose:
                    print '\n \n DIR scientific dataset ' , datatype, dir(sd_data)
                    print sd_data.attributes() , '\n'
                    #print 'Number of time locations' , len(time)
                return ma_sd_data
            else:
                return 0

      def smoothdata(self, data, sigma=2, trunc=2):
                smoothdata = ndimage.filters.gaussian_filter(data.T, sigma, truncate=trunc, mode='nearest')                
                return smoothdata


      def plot532hist(self, verbose=0, fignum=1):
                plt.figure(fignum)
                plt.hist(self.mab532.compressed())
  

      def plot532mesh(self, verbose=0, fignum=10, xaxis='distance', prtcls=pd.DataFrame(), 
                      ttl='', cbar=0, outpt='jpg', clrs='grey'):
                """Plots 532 nm backscatter"""
                datatype='532'
                forpaper=False
                if datatype=='532':
                   sdata = self.mab532
                   clabel = 'Total Attenuated \n Backscatter  532 nm (km$^{-1}$ sr$^{-1}$)'
                elif datatype == 'color ratio':
                   sdata = self.mabcr
                   clabel = 'Color Ratio'
                    
                cmapdir = '/n-home/alicec/python/cmap/'
                if forpaper:
                    fs = 28
                    rcParams['figure.figsize'] = 11 , 6
                    rcParams['axes.labelsize'] =  fs
                    rcParams['xtick.labelsize'] = fs
                    rcParams['ytick.labelsize'] = fs
                    rcParams['xtick.major.width'] = 2
                    rcParams['xtick.major.size'] = 15
                    rcParams['ytick.major.width'] = 2
                    rcParams['ytick.major.size'] = 15
                    rcParams['ytick.minor.width'] = 2
                    rcParams['ytick.minor.size'] = 10
                    ticksp = 3000
                else:
                    fs = 8
                    rcParams['figure.figsize'] =  4, 2
                    rcParams['axes.labelsize'] =  fs
                    rcParams['xtick.labelsize'] = fs
                    rcParams['ytick.labelsize'] = fs
                    ticksp = 2000
                    #print mab532.shape
                #shp1 = self.mab532.shape[0]
                #shp2 = self.mab532.shape[1]
                shp1 = sdata.shape[0]
                shp2 = sdata.shape[1]
                if xaxis == 'distance' and self.dist==[]:
                   self.get_distance()
                       
                #This block calculates xtick locations and labels 
                lbl_loc = []
                lbl_nm = [] 
                icntr = 0
                for tm in self.time:
                    #dt.append(calipso_time2dt(tm[0]))
                    if icntr % ticksp == 0:
                       latlonstr = "{:.2f}".format(self.lat[icntr]) + '\n' + "{:.2f}".format(self.lon[icntr])
                       dt = calipso_time2dt(tm[0])
                       #lbl_nm.append(latlonstr + '\n' + dt[icntr].strftime("%m%d %H:%M"))            
                       lbl_nm.append(latlonstr + '\n' + dt.strftime("%d %H:%M"))        
                       if xaxis == 'time':    
                          lbl_loc.append(tm)       
                       else:
                          lbl_loc.append(self.dist[icntr]) 
                    icntr+=1 

                y = self.altitudes.reshape(shp2)

                if xaxis == 'time':
                   x = self.time.reshape(shp1)
                else:
                   x = self.dist.reshape(shp1)
                fig = plt.figure(fignum)
                ax = fig.add_subplot(1,1,1)
                if clrs == 'browse':
                    if datatype == '532':
                       cmp = cmap(cmapdir + 'calipso-backscatter.cmap')
                    elif datatype == 'color_ratio':
                       cmp = cmap(cmapdir + 'calipso-cratio.cmap')
                else:
                    cmp = cmap(cmapdir + 'calipso-backscatter-amc.cmap')
                cm =  mpl.colors.ListedColormap(cmp['colors']/255.0)
                cm.set_under(cmp['under']/255.0)
                cm.set_over(cmp['over']/255.0)
                cm.set_bad(cmp['bad']/255.0)
                norm = mpl.colors.BoundaryNorm(cmp['bounds'], cm.N)
                #cp = ax.pcolormesh(x, y, self.smoothdata(self.mab532), cmap=cm, norm = norm, vmin=-4)
                cp = ax.pcolormesh(x, y, self.smoothdata(sdata), cmap=cm, norm = norm, vmin=-4)
                if not prtcls.empty:
                    print 'printing particles'
                  # ptdist = []
                  # ptht = []
                    ##magenta points
                    #ax.plot(prtcls['ldist'],prtcls['ht']/1000.0, '.', color='#CC33FF', markersize=4)
                    ##bright green crosses
                    #ax.plot(prtcls['ldist'],prtcls['ht']/1000.0, '+', color='#99FF33', markersize=4)
                    ##bright green crosses
                    #ax.plot(prtcls['ldist'],prtcls['ht']/1000.0, '.', color='#FF0000', markersize=4, markeredgewidth=2, alpha=0.1, rasterized=True)
                    ax.plot(prtcls['ldist'],prtcls['ht']/1000.0, '.', color='#FF0000', markersize=4 )
                ax.set_ylabel('Height (km)', fontsize=fs)
                ax.set_xlabel('Lat / Lon / Day HH:MM UTC', fontsize=fs)
                ax.set_xticks(lbl_loc)
                ax.set_yticks([5,15])
                minorLocator=MultipleLocator(5)
                ax.yaxis.set_minor_locator(minorLocator)
                ax.set_xticklabels(lbl_nm, fontsize=fs)
                plt.tight_layout()
                if cbar:
                    cbar = plt.colorbar(cp)
                    cbar.set_label(clabel , fontsize=fs)
                    cbar.ax.tick_params(labelsize=fs)
                print 'saving figure'
                plt.subplots_adjust(left=0.15)
                plt.savefig(ttl + 'lidar.' + str(fignum).zfill(2) + '.' + outpt.strip(), dpi=400) 
                #return ax , lbl_loc



class CaliopL1(Caliop):
      """reads data from a calipso L1 file"""

      def __init__(self, fname, verbose=0, datatype='532'):
            super(CaliopL1, self).__init__(fname, verbose=verbose)
            self.mab532 = self.get_sd(verbose, datatype='532')
            self.dist=[]


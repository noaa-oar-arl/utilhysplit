#plot_volcat.py
#Hard coded plotting of volcat data from Reventador eruption
#15min output files
import datetime
import os
from os import walk
import os.path as path
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import seaborn as sns
import numpy as np
import numpy.ma as ma
import pandas as pd
from utilvolc import volcat
#import utilvolc.reventador_volcat as volcat


def make_test_file(tdir, daterange=None,verbose=True):
    vidA = 'v300250'  #Bezy data
    vnlistA = volcat.find_volcat(tdir, daterange=daterange, vid=vidA)
    for key in vnlistA.keys():
        vname = vnlistB[key].fname
        dset = volcat.open_dataset(os.path.join(tdir,vname),pc_correct=False)
        return dset 


def plot_volcat(tdir, daterange=None, verbose=True):
    vidB = 'v290260'  # other data
    vidA = 'v300250'  #Bezy data
    vnlistA = volcat.find_volcat(tdir, daterange=daterange, vid=vidA)
    vnlistB = volcat.find_volcat(tdir, daterange=daterange, vid=vidB)
    fignum = 1
    for key in vnlistA.keys():
        fig = plt.figure(fignum)
        vname = vnlistA[key].fname
        dset = volcat.open_dataset(os.path.join(tdir,vname),pc_correct=False)
        mass = volcat.get_mass(dset)
        mass.isel(time=0).plot.pcolormesh(x='longitude',y='latitude',cmap='Blues')
        if key in vnlistB.keys():
            vname = vnlistB[key].fname
            dset = volcat.open_dataset(os.path.join(tdir,vname),pc_correct=False)
            mass = volcat.get_mass(dset)
            mass.isel(time=0).plot.pcolormesh(x='longitude',y='latitude',cmap='Reds')
        plt.title(key.strftime("%Y %m/%d %H:%M"))
        fignum += 1
        #return mass
        #mass.isel(time=0).plot.pcolormesh(x='logintude',y='latitude')
        #plt.show()
         
#def open_dataset(fname):
#    dset = xr.open_dataset(fname,mask_and_scale=False, decode_times=False)
    # use parallax corrected values of latitude and longitude
#    dset = volcat._get_latlon(dset,'pc_latitude','pc_longitude')
#    dset = volcat._get_time(dset)
#    return dset 

class ReventadorVolcat:
    # for easy retrieval of data for eruption on 2019 Feb 25

    def __init__(self,tdir):     
        # location of data
        self.directory=tdir
        # hours available
        self.hr=['15','16','17','18','19','20','21','22']
        # minutes available
        self.minute=['00','15','30','45']
        self.basename='geocatL2.GOES-16.Full_Disk.2019056.'

    def create_name(self, hr, mm):
        fname =  self.basename  +hr + mm +'30.hdf'
        if path.isfile(path.join(self.directory, fname)):
           return path.join(self.directory, fname)
        else:
           print('could not find ', path.join(self.directory, fname))
           return fname

    def get_dset(self, hr, mn):
        """
        hr : str
        mn : str
        input hour and minute of dataset to get.
        Returns xarray with data for that time period.

        """
        print(self.create_name(hr,mn))
        dset = volcat.open_dataset(self.create_name(hr, mn))
        return dset

    def generate_dsets(self, hlist=None, mlist=None):
        if not hlist: hlist = self.hr
        if not mlist: mlist = self.minute
        for hr in hlist:
            for mm in mlist:
                try:
                  dset = self.get_dset(hr, mm)
                except:
                  print('FILE cannot be opened ', self.create_name(hr,mm))
                  continue
                yield dset




def animate_revenatdor():
    fignum=1
    rev = ReventadorVolcat()
    for dset in rev.generate_dsets():
        fig = plt.figure(fignum)
        try:
            mass = volcat.get_mass(dset)
        except:
            continue
        makemap(mass,fignum)
        fignum += 1


def makemap(mass,fignum):
        sns.set_style('white')
       
        #mass.isel(time=0).plot.pcolormesh(x='longitude',y='latitude')
        img_proj=ccrs.PlateCarree()
        m=plt.axes(projection=img_proj)
        m.add_feature(cfeat.LAND)
        m.add_feature(cfeat.COASTLINE)
        m.add_feature(cfeat.BORDERS)
        #Plots volcat data
        plt.pcolormesh(mass.longitude, mass.latitude, mass.isel(time=0).values, 
                       vmin=0.1,vmax=5,
                       cmap='autumn', transform=ccrs.PlateCarree())
        cb2 = plt.colorbar()
        gl = m.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                         linewidth=2,color='gray',alpha=0.5)
        gl.xlabels_top=False
        gl.ylabels_right=False
        gl.xlocator = mticker.FixedLocator([-78,-77.5,-77])
        gl.ylocator = mticker.FixedLocator([-0.8,-0.4,0,0.4])
        gl.xformatter=LONGITUDE_FORMATTER
        gl.yformatter=LATITUDE_FORMATTER
        gl.xlabel_style = {'size':12, 'color':'gray'}
        gl.ylabel_style = {'size':12, 'color':'gray'}
        cb2.set_label('g/m$^2$',fontsize=12)
        cb2.ax.tick_params(labelsize=10)
        ax = plt.gca()
        ax.set_aspect(1.0, adjustable="box")
        
        #Plots Reventador Volcano symbol
        plt.scatter(-77.6558,-0.0775,c='black',marker='^')
        #Zoomed in region
        plt.xlim(-78.2, -76.8)
        plt.ylim(-1.,0.5)
        #ttl = 'Ash Top Height (km) \n from VOLCAT at ' 
        #ttl =  d1.strftime("%Y/%m/%d %H:%M")
        d1 = pd.to_datetime(mass.time.values[0])
        ttl = d1.strftime("%Y %b %d %H:%M UTC")
        #ttl = str(mass.time.values[0])
        plt.title(ttl,size=10)
        plt.tight_layout()
        figname = 'rev.{:02d}.png'.format(fignum)
        print('saving figure ' + figname)
        plt.savefig(figname)
        #plt.savefig('Images/VOLCAT_AshTopHeight_'+hr[iii]+mm[jjj]+'30.png')
        #plt.show()
        #plt.close()   
 


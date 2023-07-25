import datetime

import cartopy
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from matplotlib.colors import BoundaryNorm
from monetio.models import hysplit


class LabelData:
    def __init__(self,time,descrip,units,source='',tag=''):
        timefmt = "%m/%d/%Y %H:%M UTC"
        self.model = 'NOAA HYSPLIT'
        try:
            self.time = time.strftime(timefmt)
            temp = time
        except:
            temp = pd.to_datetime(str(time))
            self.time = temp.strftime(timefmt)
        temp2 = temp + datetime.timedelta(hours=1)
        self.time = '{} to {}'.format(self.time, temp2.strftime(timefmt))
        self.size=15
        self.descrip = descrip
        self.units = 'units: {}'.format(units)
        self.source = 'source: {}'.format(source)
        self.tag = tag

def label_ax(ax,label,transform):
    # add labels to plot meant for labels.
    gl = ax.gridlines(crs=transform, draw_labels=False,
                      linewidth=0, color='gray', alpha=0)
    gl.top_labels = False
    gl.bottom_labels = False
    gl.right_labels = False
    gl.left_labels = False
    size = label.size
    va = 'bottom'
    ha = 'left'
    rot='horizontal'
    xloc=0.05
    yspc = 0.1
    ax.text(xloc,0.9,label.model,va=va,ha=ha,rotation=rot,
            rotation_mode='anchor', transform=ax.transAxes,size=size)
    ax.text(xloc,0.9-yspc,label.time,va=va,ha=ha,rotation=rot,
            rotation_mode='anchor', transform=ax.transAxes,size=size)
    ax.text(xloc,0.9-2*yspc,label.descrip,va=va,ha=ha,rotation=rot,
            rotation_mode='anchor', transform=ax.transAxes,size=size)
    ax.text(xloc,0.9-3*yspc,label.units,va=va,ha=ha,rotation=rot,
            rotation_mode='anchor', transform=ax.transAxes,size=size)
    ax.text(xloc,0.9-4*yspc,label.source,va=va,ha=ha,rotation=rot,
            rotation_mode='anchor', transform=ax.transAxes,size=size)
    ax.text(xloc,0.9-5*yspc,label.tag,va=va,ha=ha,rotation=rot,
            rotation_mode='anchor', transform=ax.transAxes,size=size)
     

def get_transform():
    transform = cartopy.crs.PlateCarree()
    return transform


def set_levels(mass):
    levels = [0.2,2,5,10]
    if np.max(mass) < 1:
       levels = [0.01,0.05,0.1,0.5,1]
    elif np.max(mass) < 10:
       levels = lev1
    elif np.max(mass) < 50:
       levels.extend([20,30,40,50])
    elif np.max(mass) < 100:
       levels.extend([20,30,40,50,60,70,80,90,100])
    elif np.max(mass) < 1000:
       levels.extend([50,100,250,500,750,1000])
    else:
       levels = [0.01,0.05,0.1,0.5,1]
    return levels

def massload_plot(revash,enslist,name='None',vlist=None):
    setup_plot()
    mass = massload_ensemble_mean(revash,  enslist)
    transform = get_transform()
    iii=0
    for time in mass.time.values:
        source = 'Longitude: {:0.2f} Latitude {:0.2f}'.format(vlist[0],vlist[1])
        label = LabelData(time,
                          'Column mass loading (ensemble mean)',
                          'g/m$^2$',source)
        
        mass2 = mass.sel(time=time).isel(source=0)
        x = mass2.longitude
        y = mass2.latitude
        #z = mass2.where(mass2!=0)
       
        figname = name.replace('zzz',"{}".format(iii))
        z = mass2/1000.0  #convert to g/m2
        levels = set_levels(z)
        sub_massload_plot(x,y,z,transform,levels,label,figname,vlist)      
        iii+=1
 
def sub_massload_plot(x,y,z,transform,levels,labeldata,name='None',
                      vlist=None):
    nrow = 2
    ncol = 1 #second row is for text information.
    fig,axra = plt.subplots(nrows=nrow,ncols=ncol,figsize=(10,10),
                             constrained_layout=False,
                             subplot_kw={'projection':transform})
    ax = axra[0]
    dax = axra[1]
    cmap = plt.get_cmap('tab20b')
    norm = BoundaryNorm(levels,ncolors=cmap.N,clip=False)
    cb2 = ax.contourf(x,y,z,levels=levels,cmap=cmap,transform=transform)
    #cb2 = ax.pcolormesh(x,y,z,cmap=cmap,transform=transform,norm=norm)
    if vlist:
        ax.plot(vlist[0],vlist[1],'r^')
    format_plot(ax,transform)
    label_ax(dax,labeldata,transform)
    cb = plt.colorbar(cb2)
    cb.set_label('g/m$^2$')
    if name != 'None':
       plt.savefig(name)
 
def ATLtimeloop(revash,enslist,thresh,level,vlist=None,
                name='None',norm=True,clevels=[1,10,20,40,60,80,95],
                title='HYSPLIT Applied Threshold Levels',
                adjust=0):
    """
    creates plot for each time period.
    vlist : list of volcano coordinates to plot [longitude,latitude]
    name : str
           name for saving figures. if include zzz in name will replace
           with number for each time period.

    adjust: if >0 then will adjust threshold downwards if max concentration
            below the threshold which would result in empty plots.
            New threshold will be maxval / adjust.
    """
    iii=0
    if adjust > 0:
       thresh = check_thresh(np.max(revash),thresh,adjust)

    for time in revash.time.values:
        revash2 = revash.sel(time=time)
        source = 'Longitude: {:0.2f} Latitude {:0.2f}'.format(vlist[0],vlist[1])
        label = LabelData(time,
                          'Probability of exceedance',
                          '%',source)
        rtot = ATL(revash2,enslist,thresh,level,norm=norm)
        #figname = name.replace('zzz',"{:02d}".format(iii))
        figname = name.replace('zzz',"{}".format(iii))
        if norm: rtot = rtot * 100
        title2 = '{}\n{}'.format(title, label.time)
        plotATL(rtot,vlist,name=figname,levels=clevels,thresh=thresh,
                title=title2)
        iii+=1

def check_thresh(cmax,thresh,adjust):
    """
    set threshold at half the maximum value.
    """
    thresh2 = thresh
    if cmax < thresh:
       thresh2 = float(cmax / float(adjust))
    print(cmax, adjust, thresh2)
    return thresh2


def set_ATL_text(nrow):
    if nrow >= 6:
       yplace = 0.8
       xplace = 0.15
       size=10
    else:
       yplace = 0.90
       xplace = 0.15
       size=20
      
    return xplace, yplace, size

def plotATL(rtot, vlist,name='None',levels=[1,5,10,20],
            thresh=0.2,title='HYSPLIT'):
    """
    creates plot for one time period.
    creates subplot for each vertical level.
    """
    setup_plot()
    x = rtot.longitude
    y = rtot.latitude
    temp = rtot.max(dim='z')
    zvals = rtot.z.values
    nz = len(zvals)
    if nz>1: 
       ncol = 2  
       nrow = int(np.ceil(nz/ ncol))
    else:
       nrow = 1
       ncol = 1
    z = temp.where(temp!=0)
    transform = get_transform()
    fig,axarr = plt.subplots(nrows=nrow,ncols=ncol,figsize=(10,10),
                             constrained_layout=False,
                             subplot_kw={'projection':transform})
    if nrow > 1: 
        axlist = axarr.flatten()
    else:
        axlist = [axarr]
    iii=0
    for ax in axlist:
        z = rtot.isel(z=iii)
        z = z.where(z!=0)
        label = meter2FL(rtot.z.values[iii])
        cb = ATLsubplot(ax,x,y,z,transform,label,vlist,levels=levels,nrow=nrow)
        iii+=1
    cb2 = plt.colorbar(cb) 
    #cb2.set_label('mg/m$^3$')
    cb2.set_label('Percent')
    fig.suptitle(title,size=10) 
    if name != 'None':
       plt.savefig(name)

def meter2FL(meters):
    return 'FL{:2.0f}'.format(meters/30.48)

def setup_plot():
    mpl.rcParams.update(mpl.rcParamsDefault)
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.use('pdf')
    plt.style.use('seaborn-poster') 
    #plt.style.use('fivethirtyeight') 

def format_plot(ax,transform):
    #ax.add_feature(cartopy.feature.LAND) 
    ax.add_feature(cartopy.feature.BORDERS) 
    ax.coastlines('50m')
    # This allows latitude and longitude axis
    # to have different scales. 
    ax.set_aspect('auto', adjustable=None)
    # this will increase data limit to keep aspect ratio 1
    #ax.set_aspect(1, adjustable='datalim')
    # this will adjust axxes to  keep aspect ratio 1
    # when this is used, text is often mis-placed.
    #ax.set_aspect(1, adjustable='box')
    gl = ax.gridlines(crs=transform, draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.bottom_labels = True
    gl.right_labels = False
    gl.left_labels = True
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'gray'}
    gl.ylabel_style = {'size': 10, 'color': 'gray'}


def ATLsubplot(ax, x,y,z,transform,label='',vlist=None,
               levels=[1,5,10,15,20],name='None',nrow=6):
    """
    ax : axes
    x :
    y :
    z :
    transform :
    label : str
    vlist : [longiude, latitude]
    levels : list of floats/ints, levels for contouring.
    name : str save figure to this name. 
    """
    xplace, yplace, size = set_ATL_text(nrow)
    cmap = plt.get_cmap('viridis')
    norm = BoundaryNorm(levels,ncolors=cmap.N,clip=False)
    cb2 = ax.pcolormesh(x,y,z,cmap=cmap,transform=transform,norm=norm)
    ax.plot(vlist[0],vlist[1],'r^')
    format_plot(ax,transform)
    ax.text(xplace,yplace,label,va='bottom',ha='center',rotation='horizontal',
            rotation_mode='anchor', transform=ax.transAxes,size=size,
            backgroundcolor='white')
    #ax.plot(vlist[0],vlist[1],'r^')
    if name != 'None':
       plt.savefig(name)
    return cb2 


def massload_ensemble_mean(revash,  enslist):
    r2 = revash.sel(ens=enslist)
    mass = hysplit.hysp_massload(r2)
    massmean = mass.mean(dim='ens')
    return massmean

def ATL(revash, enslist, thresh=0.2, level=1,norm=False):
     """
     Returns array with number of ensemble members above
     given threshold at each location.
     norm : boolean
            if True return percentage of members.
            if False return number of members.
     """
     #import matplotlib.pyplot as plt
     #sns.set_style('whitegrid')
     source=0
     if isinstance(level, int):
        level=[level]
     iii=0
     for lev in level:
         #r2 = revash.sel(time=time)
         r2 = revash.sel(ens=enslist)
         r2 = r2.sel(z=lev)
         if 'source' in r2.coords:
             r2 = r2.isel(source=source)
         # place zeros where it is below threshold
         r2 = r2.where(r2>=thresh)
         r2 = r2.fillna(0)
         # place onces where it is above threshold
         r2 = r2.where(r2<thresh)
         r2 = r2.fillna(1)
         # ensemble members were above threshold at each location. 
         r2 = r2.sum(dim=['ens'])
         if iii==0: 
            rtot = r2
            rtot.expand_dims('z')
         else:
            r2.expand_dims('z')
            rtot = xr.concat([rtot,r2],'z')
         iii+=1 
     # This gives you maximum value that were above concentration 
     # at each location.
     if norm:
         nmembers = len(enslist)
         rtot = rtot / nmembers
     return rtot

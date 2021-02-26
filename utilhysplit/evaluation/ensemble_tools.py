from matplotlib.colors import BoundaryNorm
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import numpy.ma as ma
import monetio.models.hysplit as hysplit
import datetime
import seaborn as sns
import os
import matplotlib.pyplot as plt
from utilhysplit.evaluation.statmain import cdf

"""
FUNCTIONS
The functions in this file are for
manipulating concentration xarray objects which have dimesions of lat,lon,z,time,ens,source

maxconc   : returns maximum concentration along a dimension(x,y, or z)
topheight : returns height of level
ATL       : applied threshold level
ens_mean  : creates array with ensemble mean values.

PLOTTING
plotmaxconc : creates plot of maximum concentrations along a dimension (x,y,z)
getclrslevs : 
plotATL :
plot_ens_mean :
ashcmap :
draw_map :

For reading writing netcdf file
makenc :
readnc :

Use Examples
maketestra :
exampleATL
"""

def listvals(dra):
    """
    returns 1d list of values in the data-array
    """
    vals = dra.values
    vshape = np.array(vals.shape)
    return vals.reshape(np.prod(vshape))

def ens_cdf(dra, timelist=None, source=0, threshold=0):
    """
    produces plots of cumulative distribution functions.
    dra : xarray DataArray produced by combine_dataset function.
    timelist : list of times in the time coordinate to produce plots for
    source : int : index of source to use.
    threshold : float : produce CDF for values > threshold.
    """
    # select the source of interest
    from utilhysplit.evaluation.statmain import cdf
    fig = plt.figure(1)
    ax = fig.add_subplot(1,1,1)
    dra = dra.isel(source=source)
    enslist = dra.ens.values
    clrs = ['r','y','g','c','b','k']
    cdflist = []
    iii=0
    for ens in enslist:
        #print('working on ', ens)
        subdra= dra.sel(ens=ens)
        mass = hysplit.hysp_massload(subdra)
        if not isinstance(timelist,list) and not isinstance(timelist,np.ndarray):
           timelist = [mass.time.values[0]]
        for tm in timelist:
            #print('working on ', tm)
            tmass = mass.sel(time=tm)
            # create cdf from values above threshold
            sdata, y = cdf([x for x in listvals(tmass) if x > threshold])
            # plot as step function.
            ax.step(sdata, y, '-'+clrs[iii])
        if iii > len(clrs)-1: iii=0
        iii+=1
        ax.set_xscale('log')
    return -1

def draw_map(fignum, fs=20):
    proj = ccrs.PlateCarree()
    fig = plt.figure(fignum, figsize=(16,9))
        
    ax = fig.add_subplot(1,1,1)
    ax = plt.axes(projection=proj)
    gl = ax.gridlines(draw_labels=True, linewidth = 0.2, color="gray")
    gl.labels_right=False
    gl.labels_top=False
    #ocean = cfeature.NaturalEarthFeature(category='physical',name='ocean',
    #        scale='50m', facecolor = cfeature.COLORE['water'])

    #ax.add_feature(states.edgecolor='gray')
    ax.add_feature(cfeature.LAND)
    return ax

def exampleATL(revash):
    """
    output of maketestra can be used.
    """
    #if not revash:
    #    revash = maketestra()
    enslist = revash.ens.values
    d1 = datetime.datetime(2019,2,25,20)
    level = revash.z.values
    # remove deposition layer if it exists.
    if 0 in level:
       level = np.delete(level, np.where(level=0))
    # mult factor should convert unit mass into mg. 
    mult = 1e20
    fignum=1
    fig = plt.figure(fignum)
    ax = draw_map(fignum)
    plotATL(ax, revash*mult, time=d1, enslist=enslist, level=level )

def maketestra():
    flist=[]
    dname='/pub/ECMWF/JPSS/reventador/ens4/'
    basefname='cdump.004a.A'
    elist = [0,1,2,3,4,5,6,7,8,9]
    source = 'S1'
    for eee in elist:
        fname = basefname + 'e' + str(eee)
        flist.append((os.path.join(dname,fname), source, '4ERA5e' + str(eee)))
    #revash = monetio.models.hysplit.combine_datatset(flist, drange=None)
    revash = hysplit.combine_dataset(flist)
    return revash


def makenc(dset, ncname):
    """
    time, x, y, z have attributes which cause trouble
    when trying to write a netcdf file.
    """
    attlist = ['area','proj4_srs']
    for att in attlist:
        try:
            del dset.time.attrs[att]
        except:
            pass
        try:
            del dset.x.attrs[att]
        except:
            pass
        try:
            del dset.y.attrs[att]
        except:
            pass
        try:
            del dset.z.attrs[att]
        except:
            pass
    dset.to_netcdf(ncname, format='NETCDF4')
    return -1

def readnc(name):
    """
    reads in as dataset.
    convert to dataarray
    remove variable as a dimension 
    """
    dset = xr.open_dataset(name)
    dra = dset.to_array()
    dra = dra.isel(variable=0)
    return dra

def ashcmap(t=1):
     """
     possible colormaps for plotting volcanic ash.
     """
     import matplotlib
     #cmap = matplotlib.cm.cubehelix(np.linspace(0,1,20)) 
     #cmap = matplotlib.cm.RdGy(np.linspace(0,1,20)) 
     #cmap = matplotlib.colors.ListedColormap(cmap[10:,:-1])
     if t==1:
         cmap = 'RdGy_r'
     elif t==2:
         cmap=sns.diverging_palette(200,1,as_cmap=True, center='light',l=10)
     elif t==3:
         cmap=sns.diverging_palette(200,1,as_cmap=True, center='dark',l=70)
     return cmap

#_______________________________________________________________________
"""
FUNCTIONS
maxconc : returns maximum concentration along a dimension(x,y, or z)
"""

def maxconc(revash, time, enslist, level, dim='y'):
    """
    Returns maximum concetration along a dimension, x, y, or z.
    revash : xarray DataArray
    time   : datetime object
    enslist : list of integers, or integer
    level : list of integers (indices of vertical levels to choose).
            should NOT include the deposition layer if dim='y' or dim='x'.
    """
    source=0
    #L returns maximum concentration at any level.
    r2 = revash.sel(time=time)
    if 'ens' in revash.dims:
        r2 = r2.isel(ens=enslist)
    if 'source' in revash.dims:
        r2 = r2.isel(source=source)
    r2 = r2.isel(z=level)
    rtot = r2.max(dim=dim)
    # find maximum among all ensemble members as well.
    if 'ens' in revash.dims:
        if len(enslist) > 1:
            rtot = rtot.max(dim='ens')
        else:
            rtot = rtot.isel(ens=0)
    return rtot

def topheight(revash, time,ens, level,thresh=0.01,tp=0):
     source=0
     if isinstance(level, int):
        level=[level]
     iii=0
     for lev in level:
         lev_value = revash.z.values[lev]
         r2 = revash.sel(time=time)
         r2 = r2.isel(ens=ens)
         r2 = r2.isel(source=source, z=lev)
         # place zeros where it is below threshold
         r2 = r2.where(r2>=thresh)
         r2 = r2.fillna(0)
         # place ones where it is above threshold
         r2 = r2.where(r2<thresh)
         r2 = r2.fillna(lev_value)
         if iii==0: 
            rtot = r2
            rtot.expand_dims('z')
         else:
            r2.expand_dims('z')
            rtot = xr.concat([rtot,r2],'z')
         rht = rtot.max(dim='z')
     return  rht


def ATL(revash, enslist=None, sourcelist=None, thresh=0.2, level=None, norm=False):
    """
     revash: xr dataarray produced by combine_dataset or by hysp_massload
     enslist : list of values to use fo 'ens' coordinate
     sourcelist : list of values to use for 'source' coordinate
     level : integer or list of integers of vertical levels to use.
     

     Returns array with number of ensemble members above
     given threshold at each location.
     norm : boolean
            if True return percentage of members.
            if False return number of members.
     """
    rev = revash.copy()
    # if a massloading dataset with no z coordinate is input.
    # then expand dimension to z so can process.

    if 'z' not in rev.coords:
        rev = rev.expand_dims('z')

    if enslist:
        rev = rev.isel(ens=enslist)

    if sourcelist:
        rev = rev.isel(source=sourcelist)

    if isinstance(level, int):
        level = [level]
    if not level:
       mlen = len(rev.z.values) 
       level = np.arange(0,mlen) 
    # stack ensemble and source dimension into one.
    rev = rev.stack(ensemble=("ens","source"))

    # not sure if this loop over the levels is needed.
    for iii, lev in enumerate(level):
        #rev2 = rev.sel(ens=enslist)
        rev2 = rev.isel(z=lev)
        if "source" in rev2.coords:
            rev2 = rev2.isel(source=source)
        # place zeros where it is below threshold
        rev2 = rev2.where(rev2 >= thresh)
        rev2 = rev2.fillna(0)
        # place onces where it is above threshold
        rev2 = rev2.where(rev2 < thresh)
        rev2 = rev2.fillna(1)
        # ensemble members were above threshold at each location.
        rev2 = rev2.sum(dim=["ensemble"])
        if iii == 0:
            rtot = rev2
            rtot.expand_dims("z")
        else:
            rev2.expand_dims("z")
            rtot = xr.concat([rtot, rev2], "z")
    # This gives you maximum value that were above concentration
    # at each location.
    if norm:
        nmembers = len(enslist)
        rtot = rtot / nmembers
    return rtot

def ens_mean(revash,time,enslist,level=1):
    """
    revash : xarray
    time : datetime object : date in time coordinate
    enslist: list of integers.
    level : int : index of level in z coordinate.
    """
    source=0
    sns.set()
    m2 = revash.sel(time=time)
    m2 = m2.isel(ens=enslist,source=source)
    m2 = m2.mean(dim=['ens'])
    if level > 0:
        m2 = m2.isel(z=level)
    #plt.plot(-77.656, -0.077, 'k^')     
    #cmap = sns.choose_colorbrewer_palette('colorblind', as_cmap=True)
    #cb2 = plt.contourf(m2.longitude, m2.latitude, m2, cmap=cmap) 
    #plot_rev_data(time,tp='mass')
    #plt.colorbar(cb2) 
    #plt.tight_layout()
    return m2

 


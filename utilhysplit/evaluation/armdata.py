import os
import datetime
import numpy as np
import matplotlib.pyplot as plt 
import xarray as xr

def get_variance(dset):
    wvar = dset.w_variance
    # remove first two levels 
    # only use levels about 80 m from instrument.
    wvar = wvar.sel(height=wvar.height[2:133])
    wvar = wvar ** 0.5
    return wvar

def plot_variance(dset, levels=None, coarsen=(6,1), ax=None):
    if not ax: ax = plt.gca()
    from matplotlib.colors import BoundaryNorm
    wvar = get_variance(dset)
    wvar = wvar.fillna(0)
    wvar = wvar.coarsen(time=coarsen[0], height=coarsen[1], boundary='trim').mean()
    
    #wvar.plot.pcolormesh(x='time', y='height')
    lev, norm, cmap = get_vmix_colors()
    if not levels: levels = lev
    cb = ax.pcolormesh(wvar.time, wvar.height/1000.0, wvar.T,
                        cmap = cmap, norm=norm)
    #plt.colorbar(cb)
    plt.tight_layout()
    #plt.show() 
    return cb

def get_vmix_colors(levels=None):
    from matplotlib.colors import BoundaryNorm
    if not isinstance(levels,list):
        levels = [0,0.002,0.004,0.008]
        levels.extend([0.04,0.06,0.08,0.1])
        levels.extend([0.6,0.8,1.0,1.2,1.4,1.6,1.8,2,2.5,3])
    #cmap = plt.get_cmap('bone')
    #cmap = plt.get_cmap('copper')
    cmap = plt.get_cmap('viridis')
    #cmap = plt.get_cmap('cividis')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    return levels, norm, cmap
 

class ARMdata:

    def __init__(self):
        self.stations = ['C1','E32','E41']
        # zzz to be replaced by a station name.
        self.basename = 'sgpdlprofwstats4news'
        self.tdir = '/pub/Scratch/fantinen/DATA_ARM/wvar/%Y/'
        self.fmt = self.basename + 'zzz' + '.c1.%Y%m%d.000000.nc'
        self.fmt2 = self.basename + 'zzz' + '.c1.%Y%m%d.000500.nc'

    def get_locations(self,sampledate):
        drange = [sampledate, sampledate]
        lochash = {}
        for stn in self.stations:
            dset = self.get_data(stn, drange)
            #return dset
            lochash[stn] = (float(dset.lat.values), float(dset.lon.values))
        return lochash     

    def get_data(self,stn,drange):
        flist = self.make_fname_list(stn,drange)
        dset = xr.open_mfdataset(flist)
        dset.load()
        return dset 

    def make_fname_list(self,stn, drange):
        dt = datetime.timedelta(hours=24)         
        flist = []
        date = drange[0]
        iii=0
        while date <= drange[1]:
          flist.append(self.make_fname(stn,date))
          date += dt
          iii+=1
          if iii> 1000: break
        # remove files which don't exist
        flist = [x for x in flist if x != None]
        return flist
 
    def make_fname(self, stn, date):
        fmt = self.tdir + self.fmt.replace('zzz',stn)
        fname = datetime.datetime.strftime(date, fmt)
        print(fmt, fname)
        if os.path.isfile(fname):
            return fname 
        else:
            fmt = self.tdir + self.fmt2.replace('zzz',stn)
            fname = datetime.datetime.strftime(date, fmt)
            if os.path.isfile(fname):
                return fname 
            else:
                return None

    @staticmethod
    def open_dataset(fname):
        return  xr.open_dataset(fname)
       
  

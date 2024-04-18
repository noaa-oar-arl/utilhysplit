import datetime
import logging
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr

from utilhysplit.plotutils.map_util import get_transform, format_plot, PlotParams
from ashapp import utils
import utilhysplit.evaluation.web_ensemble_plots as wep


logger = logging.getLogger(__name__)
utils.setup_logger(level=logging.WARNING)

class ConcentrationMontage:

    def __init__(self, cdumpra, 
                 columns=3,cmap='viridis',
                 ):

        self._source = None
        self._ens = None
        self.cdump = cdumpra
        self._unit = 'mg m$^{-3}$'
        self.columns = columns
        self.cmap = cmap
        self.get_sloc_fromatts()
        self.central_longitude = 0    
        sns.set_style('white')

    @property
    def cdump(self):
        return self._cdump

    @cdump.setter
    def cdump(self,cset):
        if isinstance(cset,xr.core.dataset.Dataset):
           varnames = list(cset.data_vars.keys())
           if len(varnames)>1: 
              logger.warning('Multiple variables in dataset')
              logger.warning(varnames)
           self._cdump = cset[varnames[0]]
        elif isinstance(cset,xr.core.dataarray.DataArray):
           self._cdump = cset
        else:
           logger.warning('input not of correct type {}'.format(type(cset)))
       
        if 'source' in self._cdump.dims:
            slist = list(self._cdump.source.values)
            slist = [str(x) for x in slist]
            rstr = ','.join(slist)
            self._source = slist[0]
            logger.warning('Sources found as a dimension {}. choosing the first one'.format(rstr))
            self._cdump = self._cdump.isel(source=0) 
        if 'ens' in self._cdump.dims:
            slist = list(self._cdump.ens.values)
            slist = [str(x) for x in slist]
            rstr = ','.join(slist)
            self._ens = slist[0]
            logger.warning('Ensemble found as a dimension {}. choosing the first one'.format(rstr))
            self._cdump = self._cdump.isel(ens=0) 

    @property
    def sloc(self):
        return self._sloc

    @sloc.setter
    def sloc(self,sloc):
        if isinstance(sloc,tuple):
           self._sloc = sloc
        elif isinstance(sloc, [list,np.ndarray]):
           self._sloc = tuple(sloc)
        else:
           self._sloc=None   

    @property
    def cmap(self):
        return self._cmap

    @cmap.setter
    def cmap(self,cmapstr):
        self._cmap = plt.get_cmap(cmapstr)

    @property
    def times(self):
        return self._cdump.time.values

    @property
    def levels(self):
        return self._cdump.z.values

    @property
    def numpages(self):
        return len(self._cdump.time.values) 

    @property
    def rows(self):
        nnn = len(self.levels)
        return int(np.ceil(nnn/self.columns))

    @property
    def columns(self):
        return self._columns

    @columns.setter
    def columns(self,ccc):
        self._columns = int(ccc)


    @property
    def figsize(self):
        iii = 4 * self.columns
        jjj = 3 * self.rows
        return(iii,jjj)

    @property   
    def unit(self):
        return self._unit

    @property
    def central_longitude(self):
        return self._central_longitude
    
    @central_longitude.setter
    def central_longitude(self,central_longitude):
        if central_longitude <= 180 and central_longitude >= -180:
            self._central_longitude = central_longitude
        else:
            logger.warning('central longitude out of range -180 to 180')
            self._central_longitude=0 

    @property
    def transform(self):
        return get_transform(central_longitude=self.central_longitude)

    @property
    def data_transform(self):
        # data always has central_longitude=0
        return get_transform(central_longitude=0)


    def setup_figure(self,fignum=1):
        # add an extra row for text.
        fig,axra = plt.subplots(figsize=self.figsize,nrows=self.rows+1,ncols=self.columns,
                                constrained_layout=False, subplot_kw={'projection':self.transform})
        return fig, axra

    def makepages(self,pagenums=None,savename='concplot.png'):
        temp = savename.split('.')
        base = temp[0]
        ptype = temp[1]
        if isinstance(pagenums,list):
           plist = pagenums
        else:
           plist = self.times
        for iii, pnum in enumerate(plist):
            savename = '{:s}_{:02d}.{:s}'.format(base,iii,ptype)
            self.savepage(pnum,savename)

    def savepage(self,pagenum,savename):
        fig = self.plotpage(pagenum)
        fig.savefig(savename)

    def plotpage(self,pagenum):
        if isinstance(pagenum,int):
            ctemp = self.cdump.isel(time=pagenum)
        elif isinstance(pagenum,(datetime.datetime,np.datetime64)):
            ctemp = self.cdump.sel(time=pagenum)
          

        # levels and boundaries stay the same on the page.
        levels = wep.set_levels(ctemp.values)

        pm = PlotParams(ctemp,levels[0])
        self.central_longitude = pm.central_longitude
        xmax = pm.xmax
        xmin = pm.xmin
        ymax = pm.ymax
        ymin = pm.ymin
        yticks = pm.yticks
        xticks = pm.xticks

        norm = BoundaryNorm(levels,ncolors=self.cmap.N,clip=False)
        fig, axlist = self.setup_figure(fignum=pagenum)
        axtext = []
        for iii, ax in enumerate(axlist.flatten()):
            if iii >= len(self.levels): 
               ax.axis('off')
               axtext.append(ax)
               continue
            ctemp2 = ctemp.isel(z=iii)          
            xxx = ctemp2.longitude.values
            yyy = ctemp2.latitude.values
            zzz = ctemp2.values
            zzz = np.where(zzz>levels[0],zzz,np.nan)
            cb = ax.pcolormesh(xxx,yyy,zzz,transform=self.data_transform,norm=norm,cmap=self.cmap)   
            ax.set_xlim(xmin,xmax)
            ax.set_ylim(ymin,ymax)

            (row,column) = self.get_row_column(iii)
            if column==0: ylabel=True
            else: ylabel=False
            if row+1==self.rows: xlabel=True
            else: xlabel=False

            xplace,yplace,size = self.set_text(row)
            label = wep.meter2FL(ctemp2.z.values)
            ax.text(
                xplace,
                yplace,
                label,
                va="bottom",
                ha="center",
                rotation="horizontal",
                rotation_mode="anchor",
                transform=ax.transAxes,
                size=size,
                backgroundcolor="white",
            )

            if isinstance(self.sloc, tuple):
               ax.plot(self.sloc[0],self.sloc[1],'r^',markersize=10)

            format_plot(ax,self.data_transform,fsz=15,xticks=xticks,yticks=yticks,xlabel=xlabel,ylabel=ylabel)



        cb2 = fig.colorbar(cb,ax=axlist[:,self.columns-1])
        cb2.ax.tick_params(labelsize=size)
        cb2.set_label(self.unit,fontsize=size)
        time = pd.to_datetime(ctemp.time.values)
        figlabel = attrs2label(self.cdump.attrs, time, self._source, self._ens)
        axt = axtext[0]
        axt.text(0,0.0,figlabel,
                va="bottom",
                ha="left",
                rotation="horizontal",
                rotation_mode="anchor",
                transform=axt.transAxes,
                backgroundcolor="white",
                fontsize=size)
        #for iii, axi in enumerate(axtext):
        #    axi.text(0,0,str(iii),va='bottom',ha='left',rotation='horizontal',
        #             transform = axi.transAxes,backgroundcolor='white',fontsize=size) 
        return fig


    def set_text(self,nrow):
        yplace = 0.95
        xplace = 0.80
        size = 15
        return xplace, yplace, size


    def get_row_column(self,iii):
        column = iii%self.columns
        row = int(np.floor(iii/self.columns))
        return (row,column)


    def get_sloc_fromatts(self):
        attrs = self.cdump.attrs
        attkeys = list(attrs.keys())
        check=True
        if 'latitude' in attkeys:
           lat = attrs['latitude']
        else:
           check=False
        if 'longitude' in attkeys:
           lon = attrs['longitude']
        else:
           check=False
        if check: self.sloc = (lon,lat) 


def attrs2label(atthash,time,source=None,ens=None):
    label = 'HYSPLIT\n'
    label += time.strftime("%Y %m/%d %H:%M UTC ")
    label += ' 1h average\n'
    keys = list(atthash.keys())
    for xxx in ['VolcanoName','Meteorological Model ID']:
        if xxx in keys:
            label += '{}: {}\n'.format(xxx, atthash[xxx])
    for xxx in ['jobid','jobname']:
        if xxx in keys:
            label += '{}: {}  '.format(xxx, atthash[xxx])
    if isinstance(source,str):
       label += '\n Source: {}'.format(source)
    if isinstance(ens,str):
       label += '  MemberID: {}'.format(ens)
    return label 




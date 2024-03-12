import datetime
import logging
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr

from utilhysplit.plotutils.map_util import get_transform, format_plot, PlotParams
import  ashapp.utils as utils
from utilhysplit.evaluation import web_ensemble_plots as wep
from utilhysplit.evaluation import ensemble_polygons, ensemble_tools

logger = logging.getLogger(__name__)
utils.setup_logger(level=logging.WARNING)

class VAAMontageATL:

    def __init__(self, cdumpra, 
                 columns=3,cmap='Reds',
                 vaathresh=0.2,probthresh=0.01,
                 ):

        self._source = None
        self._ens = None
        # assume cdump has units of mg/m3
        self.cdump = cdumpra
        # convert to applied threshold level
        #self._atl = ensemble_tools.ATL(self.cdump,thresh=vaathresh)
        self._unit = 'Percent above threshold (maximum in column)'
        self.columns = columns
        self.cmap = cmap
        self.sloc = self.get_sloc_fromatts()
        self.central_longitude = 0    
        self.vaathresh=vaathresh
        self.probthresh=probthresh
        sns.set_style('white')
        # sets figsize based on columns and rows.
        self.figsize=None
        self.control_member = 0

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
        #if 'ens' in self._cdump.dims:
        #    slist = list(self._cdump.ens.values)
        #    slist = [str(x) for x in slist]
        #    rstr = ','.join(slist)
        #    self._ens = slist[0]
        #    logger.warning('Ensemble found as a dimension {}. choosing the first one'.format(rstr))
        #    self._cdump = self._cdump.isel(ens=0) 

    @property
    def sloc(self):
        return self._sloc

    @sloc.setter
    def sloc(self,sloc):
        if isinstance(sloc,tuple):
           self._sloc = sloc
        elif isinstance(sloc, (list,np.ndarray)):
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
        nnn = len(self.times)
        return int(np.ceil(nnn/self.columns))

    @property
    def columns(self):
        return self._columns

    @columns.setter
    def columns(self,ccc):
        self._columns = int(ccc)


    @property
    def figsize(self):
        return self._figsize

    @figsize.setter
    def figsize(self, figsize=None):
        if isinstance(figsize,tuple):
            self._figsize = figsize
        else:
            iii = 4 * self.columns
            jjj = 4 * self.rows
            self._figsize = (iii,jjj)

    @property   
    def unit(self):
        return self._unit

    @property
    def central_longitude(self):
        return self._central_longitude
    
    @central_longitude.setter
    def central_longitude(self,central_longitude):
        if central_longitude <= 180 and central_longitude >= -180:
            self._central_longitude = int(central_longitude)
        else:
            logger.warning('central longitude out of range -180 to 180')
            self._central_longitude=0 
        # always make it positive 180.
        if central_longitude == -180:
           self._central_longitude = 180

    @property
    def transform(self):
        return get_transform(central_longitude=self.central_longitude)

    @property
    def data_transform(self):
        # data always has central_longitude=0
        return get_transform(central_longitude=0)


    def setup_figure(self,fignum=1):
        # add an extra row for text.
        fig,axra = plt.subplots(figsize=self.figsize,nrows=self.rows,ncols=self.columns,
                                constrained_layout=False, subplot_kw={'projection':self.transform})
        return fig, axra


    def savepage(self,savename):
        fig = self.plotpage()
        fig.savefig(savename)

    @property
    def atl(self):
        self._atl = ensemble_tools.ATL(self.cdump,thresh=self.vaathresh)
        return self._atl * 100

    @property
    def vaathresh(self):
        return self._vaathresh

    @vaathresh.setter
    def vaathresh(self,thresh):
        if isinstance(thresh,(float,int)):
           self._vaathresh= thresh
        else:
           self._vaathresh = float(thresh)


    def levelgetter(self,zvals):
        temp = np.where(zvals>0,zvals,np.nan)
        minval = np.nanmin(temp)
        print(minval)
        minx = int(np.floor(minval))
        if minx < 5: dx = 5 
        else: dx = minx
        values = np.arange(0,100,5)
        values[0] = minx
        return values
        
    def plotpage(self):
        ctemp = self.atl
        control = self.cdump.isel(ens=self.control_member)
        # levels and boundaries stay the same on the page.
        # levels = wep.set_levels(ctemp.values)
        levels = self.levelgetter(ctemp.values)

        pm = PlotParams(ctemp,levels[0])
        self.central_longitude = pm.central_longitude
        xmax = pm.xmax
        xmin = pm.xmin
        ymax = pm.ymax
        ymin = pm.ymin
        yticks = pm.yticks
        xticks = pm.xticks

        norm = BoundaryNorm(levels,ncolors=self.cmap.N,clip=False)

        fig, axlist = self.setup_figure(fignum=1)
        axtext = []
        for iii, ax in enumerate(axlist.flatten()):
            logger.info('working on {}'.format(iii))
            if iii >= len(self.times): 
               ax.axis('off')
               axtext.append(ax)
               continue
            ctemp2 = ctemp.isel(time=iii).max(dim='z')   

            xxx = ctemp2.longitude.values
            yyy = ctemp2.latitude.values
            zzz = ctemp2.values
            zzz = np.where(zzz>levels[0],zzz,np.nan)
            cb = ax.pcolormesh(xxx,yyy,zzz,transform=self.data_transform,norm=norm,cmap=self.cmap)   
            ax.set_xlim(xmin,xmax)
            ax.set_ylim(ymin,ymax)

            # Block for polygons
            cdump = self.cdump.isel(time=iii) 
            #top, bottom = ensemble_tools.topheight(ctemp.isel(time=iii),time=None,thresh=self.probthresh)
            top, bottom = ensemble_tools.topheight(control.isel(time=iii),time=None,thresh=self.vaathresh)
            top_poly = ensemble_polygons.HeightPolygons(cmap='winter')  
            skip=False
            try:
                top_poly.process(top,alpha=0.1)
            except  Exception as eee:
                logger.warning(eee)
                skip=True
            if not skip:
                tpoly = top_poly.merge(key='high')
                bottom_poly = ensemble_polygons.HeightPolygons(cmap='winter')  
                bottom_poly.process(bottom,alpha=0.1)
                bottom_poly = top_poly.merge(key='low')
                tpoly = top_poly.merge_with(bottom_poly)
                handles, labels = tpoly.plot(ax=ax,pbuffer=0.15,legend=False, linewidth=2,transform=self.data_transform)
                #ax.legend(handles,labels,fontsize=10,loc='lower left')
                lbb = list(set(labels))
                lbb = ''.join(lbb)
                print(lbb)

            (row,column) = self.get_row_column(iii)
            if column==0: ylabel=True
            else: ylabel=False
            if row+1==self.rows: xlabel=True
            else: xlabel=False

            xplace,yplace,size = self.set_text(row)
            label = pd.to_datetime(ctemp2.time.values)
            label = label.strftime('%y %m/%d %H UTC')
            label += '  ' + lbb
            try:
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
            except:
                pass

            if isinstance(self.sloc, tuple):
               ax.plot(self.sloc[0],self.sloc[1],'y^',markersize=10)

            format_plot(ax,self.data_transform,fsz=10,xticks=xticks,yticks=yticks,xlabel=xlabel,ylabel=ylabel)



        cb2 = fig.colorbar(cb,ax=axlist[:,self.columns-1])
        cb2.ax.tick_params(labelsize=size)
        cb2.set_label(self.unit,fontsize=size)
        time = pd.to_datetime(ctemp.time.values)
        figlabel = attrs2label(self.cdump.attrs, datetime.datetime.now(),thresh=self.vaathresh,probthresh=self.probthresh )
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
        xplace = 0.50
        size = 10
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
        if check: rval = (lon,lat) 
        else: rval = None
        return rval

def attrs2label(atthash,time,source=None,ens=None, thresh='Unknown', probthresh=0):
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
    label += '\n'
    #label += 'Polygon shows area where at least {}% ensemble members predict \n ash concentration above {} mg/m$^3$ \n'.format(probthresh,thresh)
    label += 'Polygon shows where ash concentration above {} mg/m$^3$ \n for control run \n'.format(thresh)
    label += 'Top and bottom levels of polygon indicated in top right'.format(thresh)
    return label 




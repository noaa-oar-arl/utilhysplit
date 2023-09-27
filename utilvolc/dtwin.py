import logging
import json
import pandas as pd
import numpy as np
import os
import datetime
import seaborn as sns
import matplotlib.pyplot as plt
import shapely.geometry as sgeo
from utilvolc.runhelper import Helper
from utilvolc.runhelper import list_dirs
from utilvolc.runhelper import make_dir
from utilvolc import volcat
import utilvolc.iwxxmVAA as ixa
from utilvolc import volcat_plots as vp
from utilvolc import qva_logic
from utilhysplit import geotools
from utilhysplit import volcat2object


def process(work, icollection, ax=None):
    """
    work : WorkFlow object (see volcat_logic.py)
    icollection : iwxxmCollection object (see iwxxmVAA.py)
    """

    fignum=1
    elist = eruption_list(work, icollection)
    for eruption in elist:
        print('WORKING ON {} {}'.format(eruption[0], eruption[1]))
        eve = work.ehash[eruption[1]]
        #eve.write_parallax_corrected()
        eve.get_volcat_events()
        vplot = vp.VolcatPlots(eve.events)
        vplot.make_arrays()
        #vplot.describe_plot()
        vaa = ixa.create_new_collection(icollection, eruption[0])
        #if not ax:
        fig = plt.figure(fignum)
        ax = fig.add_subplot(1,1,1)
        vplot.sub_plot_area(ax)
        ax2 = ax.twinx()
        d1,d2 = vaa.time_series(ax=ax2,vname=eruption[0])
        #d2 = datetime.datetime.now()
        #d1 = d2 - datetime.timedelta(hours=14*24)
        #ax.set_xlim(d1,d2)
        for tic in ax.get_xticklabels():
            tic.set_rotation(45)
        #ax.set_xticklabels(ax.get_xticks(),rotation=45)
        #ax2.set_xticklabels(ax.get_xticks(),rotation=45)
        plt.show() 
        fignum+=1


def eruption_list(work, icollection,hours=24*14):
    """
    work : an instance of the WorkFlow class 
    icollection : an instance of the iwxxmCollection class

    Returns
    eruption_list : tuple : (name in iwxxm file, name in volcat file)

    PURPOSE: 
    the name of the volcano can be slightly different in the volcat vs.
    iwxxm files. e.g. RUIZ vs. Ruiz, Navdo del

    """

    eruption_list = []
    vaas = icollection.get_volcano_list()
    volcat_names = work.get_all_volcanos(hours=hours)
    for vname in volcat_names:
        for vaa in vaas:
            if vaa.lower() in vname.lower():
               eruption_list.append((vaa,vname))
    return eruption_list


class Eruption:
    """
    class for integrating information from different sources together.
    Begin with VOLCAT and VAA
    """

    def __init__(self, vname=None, daterange=None):
        self.name = vname
        self.daterange = daterange
        self.ixc = ixa.iwxxmCollection()
        self.vcat = qva_logic.Events()

    def add_vaas(self,icc):
        self.ixc = icc


    def add_volcat(self,vcat):
        self.vcat = vcat        


    def volcat2polygon(self,event):
           vmass = volcat.get_mass(event,clip=True)
           temp = vmass.isel(time=0)
           dlon = temp.longitude.values[0][1] - temp.longitude.values[0][0]
           dlat = temp.latitude.values[1][0] - temp.latitude.values[0][0]
           chull, edges = volcat2object.obs2concavehull(temp)
           chull2 = chull.buffer(0.5*np.abs(dlat))
           return chull, chull2

    def plotboth2(self,dt=6,vaatime=0,pc=True):
        # for each VAA plot all volcat data that falls in dt hour prior.
        for xvaa in self.ixc.ilist:
            if vaatime==0:
                otime = xvaa.obs['date']
            elif vaatime==1:
                otime = xvaa.forecast['Forecast0']['date']
            elif vaatime==2:
                otime = xvaa.forecast['Forecast1']['date']
            elif vaatime==3:
                otime = xvaa.forecast['Forecast2']['date']
            xvaa.plot_vaa()
            if pc: vlist = self.vcat.events
            else: vlist = self.vcat.pcevents
            for events in vlist:
                etime = pd.to_datetime(events.time.values[0])
                diff = etime-otime
                diff = diff.days*24 + diff.seconds/3600
                if np.abs(diff) < dt:
                   print('vaa', otime, 'volcat', etime, 'match')    
                else:
                   continue
                chull, chull2 = self.volcat2polygon(events)
                print(type(chull))
                #if isinstance(chull, sgeo.polygon.Polygon):
                #    plt.plot(*chull.exterior.xy,'-c')
                plt.plot(*chull2.exterior.xy,'-y',linewidth=5,alpha=0.75)
            plt.show()

 
    def plotboth(self,dt=0.5,pc=True,vaatime=0):
        # for each volcat retrieval plot the VAA that is close.
        tdelt = datetime.timedelta(hours=dt)
        if pc: vlist = self.vcat.events
        else: vlist = self.vcat.pcevents
        for events in vlist:
            etime = pd.to_datetime(events.time.values[0])
            for xvaa in self.ixc.ilist:
                if vaatime==0:
                    otime = xvaa.obs['date']
                elif vaatime==1:
                    otime = xvaa.forecast['Forecast0']['date']
                elif vaatime==2:
                    otime = xvaa.forecast['Forecast1']['date']
                elif vaatime==3:
                    otime = xvaa.forecast['Forecast2']['date']
                diff = etime-otime
                diff = diff.days*24 + diff.seconds/3600
                if np.abs(diff) < dt:
                   print('vaa', otime, 'volcat', etime, 'match')    

                   vmass = volcat.get_mass(events,clip=True)
                   temp = vmass.isel(time=0)
                   dlon = temp.longitude.values[0][1] - temp.longitude.values[0][0]
                   dlat = temp.latitude.values[1][0] - temp.latitude.values[0][0]
                   cb = plt.pcolormesh(temp.longitude.values - 0.5*dlon,
                                       temp.latitude.values-0.5*dlat,
                                       temp.values)
                   chull, edges = volcat2object.obs2concavehull(temp)
                   chull2 = chull.buffer(0.5*np.abs(dlat))
                   print(type(chull))
                   if isinstance(chull, sgeo.polygon.Polygon):
                       plt.plot(*chull.exterior.xy,'-c')
                   plt.plot(*chull2.exterior.xy,'-r')
                   xvaa.plot_vaa()
                   plt.show()
                   
 




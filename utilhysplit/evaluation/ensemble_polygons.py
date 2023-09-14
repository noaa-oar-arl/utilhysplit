import matplotlib.pyplot as plt
import numpy as np
import utilhysplit.evaluation.web_ensemble_plots as wep
from utilhysplit import geotools
from utilhysplit.plotutils import colormaker

def FL_colors(cmap='viridis'):
    levlist = range(50,650,50)
    levlist = [int(x*100/3.28084) for x in levlist]
    nclr = len(levlist)
    cm = colormaker.ColorMaker(cmap,nclr,ctype='rgb')
    zzz = zip(levlist,cm())
    return dict(zzz) 



def get_hull(z,thresh1=0.1,thresh2=1000,alpha=10):
    lat = z.longitude.values.flatten()
    lon = z.latitude.values.flatten()
    zzz = z.values.flatten()
    tlist = list(zip(lat,lon,zzz))
    tlist = [x for x in tlist if ~np.isnan(x[2])]
    tlist = [x for x in tlist if x[2]>=thresh1]
    tlist = [x for x in tlist if x[2]<=thresh2]
    lon = [x[1] for x in tlist]
    lat = [x[0] for x in tlist]
    numpts = len(lon)
    mpts = geotools.make_multi(lon,lat)
    if numpts >= 4: 
        ch, ep = geotools.concave_hull(mpts,alpha=alpha)
    else:
        ch = mpts.convex_hull
        ep = None

    return ch, ep



class HeightPolygons:

    def __init__(self):
        self.colorhash = FL_colors(cmap='cividis')
        self.hts = self.colorhash.keys()
        self.ch_hash = {}
  
    def process(self,dset,alpha=10):
        hlevs = wep.set_height_levels(dset.values)
        hlevs.sort()
        if hlevs[0]==0: hlevs=hlevs[1:]
        self.hlevs = hlevs
        for hhh in self.hts:
            diff = [np.abs(x-hhh) for x in hlevs]
            if np.min(diff)>10: continue
            thresh1 = hhh-10
            thresh2 = hhh+10
            try:
                ch, ep = get_hull(dset,thresh1,thresh2,alpha)
            except Exception as eee:
                print(eee)
                print('error for level', wep.meterev2FL(hhh), hhh,  dset)
                continue
            # add a buffer around the shape.
            # usually 0.25 degrees because that is resolution of output.
            self.ch_hash[hhh] = ch

    def plot(self,vloc,pbuffer,ax,legend=True):
        for hhh in self.ch_hash.keys():
            if pbuffer>0: ch = self.ch_hash[hhh].buffer(pbuffer)
            else: ch = self.ch_hash[hhh]
            for x,y in geotools.plotpoly(ch):
                ax.plot(y,x,linewidth=5,color=self.colorhash[hhh],label=wep.meterev2FL(hhh))
        ax.plot(vloc[1],vloc[0],'c^',markersize=20)
        handles, labels = ax.get_legend_handles_labels()
        if legend:
            ax.legend(handles,labels)
        return handles, labels
        



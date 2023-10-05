import matplotlib.pyplot as plt
import numpy as np
import utilhysplit.evaluation.web_ensemble_plots as wep
from utilhysplit import geotools
from utilhysplit.plotutils import colormaker
import shapely.geometry as sgeo
from shapely.ops import unary_union
#from sklearn.cluster import DBSCAN

def FL_colors(levlist=None,cmap='viridis'):
    levlist = range(50,650,50)
    levlist = [int(x*100/3.28084) for x in levlist]
    nclr = len(levlist)
    cm = colormaker.ColorMaker(cmap,nclr,ctype='rgb')
    zzz = zip(levlist,cm())
    return dict(zzz) 

def other_colors(levlist,cmap='viridis'):
    nclr = len(levlist)
    cm = colormaker.ColorMaker(cmap,nclr,ctype='rgb')
    zzz = zip(levlist,cm())
    return dict(zzz) 


def get_hull(z,thresh1=0.1,thresh2=1000,alpha=10):

    lat = z.longitude.values.flatten()
    lon = z.latitude.values.flatten()
    zzz = z.values.flatten()
    tlist = list(zip(lat,lon,zzz))

    # get lat lon values for values above thresh1 and below thresh2 and non nan.
    tlist = [x for x in tlist if ~np.isnan(x[2])]
    tlist = [x for x in tlist if x[2]>=thresh1]
    tlist = [x for x in tlist if x[2]<=thresh2]
    lon = [x[1] for x in tlist]
    lat = [x[0] for x in tlist]

    # create the polygons
    numpts = len(lon)
    mpts = geotools.make_multi(lon,lat)
    if numpts >= 4: 
        ch, ep = geotools.concave_hull(mpts,alpha=alpha)
    else:
        ch = mpts.convex_hull
        ep = None

    return ch, ep


class HeightPolygonSet:

    def __init__(self):
        self._top = HeightPolgyon(cmap='cividis')
        self._bottom = HeightPolygon(cmap='Blues')

        self._merged = HeightPolygon(cmap='Reds')

    @property
    def top(self):
        return self._top

    @property
    def bottom(self):
        return self._bottom

    @top.setter
    def top(self,dset):
        self._top.process(dset)

    @bottom.setter
    def bottom(self,dset):
        self._bottom.process(dset)


    def merge(self):
        newtop = self._top.merge()
        newbottom = self._bottom.merge()


class HeightPolygons:

    def __init__(self, levlist=None,
                 cmap='cividis',
                 colorfunction=FL_colors):
        self.colorhash = colorfunction(levlist,cmap=cmap)
        self.hts = self.colorhash.keys()
        self.ch_hash = {}
        self.fs = 20

 
    def empty(self):
        if not self.ch_hash: return True
        else: return False            
 
    def process(self,dset,alpha=10):
        """
        dset : xarray DataArray
        """
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


    def cluster(self,dist=2):
        keylist = list(self.ch_hash.keys())
        polylist = list(self.ch_hash.values())
        clist = []
        checklist = list(zip(keylist,polylist))
        for poly in checklist:
            cont=True
            for ccc in clist: 
                for c in ccc:
                    if c[0]==poly[0]: cont=False
            if cont:
                cluster = [x for x in checklist if poly[1].distance(x[1])<dist]
                clist.append(cluster)
        return clist

        #centroids = [p.centroid for p in polygons]
        #coords = [(c.x,c.y) for c in centroids]
        #dscan = DBSCAN(eps=eps_distance, min_samples=2)
        #clusters = dscan.fit(coords)           
         

    def merge(self,key='high'):
        """
        Returns a HeightPolygons object with all the polygons merged.
        """
        keylist = list(self.ch_hash.keys())
        polygons = self.ch_hash.values()
        newpoly = unary_union(polygons)
        #newpoly = newpoly.convex_hull
        #newpoly = newpoly.envelope
        low = float(np.min(keylist))
        hi = float(np.max(keylist))
        cfunc = FL_colors
        if key=='high': newkey = hi
        elif key=='low': newkey = low
        else:
            newkey = '{}_{}'.format(wep.meterev2FL(low),wep.meterev2FL(hi))
            cfunc = other_colors
        newhash = {newkey:newpoly}
        new = HeightPolygons(levlist=[newkey],colorfunction=cfunc)
        new.ch_hash = newhash
        return new 

    #def merge(self):

    def merge_with(self,other):
        self_polygons = list(self.ch_hash.values())
        other_polygons = list(other.ch_hash.values())
    
        keylist = list(self.ch_hash.keys())
        keylist.extend(list(other.ch_hash.keys()))
          
        low = float(np.min(keylist))
        hi = float(np.max(keylist))

        self_polygons.extend(other_polygons)
        
        newpoly = unary_union(self_polygons)
        newkey = '{}_{}'.format(wep.meterev2FL(low),wep.meterev2FL(hi))
        cfunc = other_colors
        newhash = {newkey:newpoly}
        new = HeightPolygons(levlist=[newkey],colorfunction=cfunc)
        new.ch_hash = newhash
       
        return new 


    def plot(self,vloc,pbuffer,ax,legend=True,linewidth=5):
        for hhh in self.ch_hash.keys():
            if isinstance(hhh,(float,int)): label = wep.meterev2FL(hhh)
            else: label = hhh
            if pbuffer>0: ch = self.ch_hash[hhh].buffer(pbuffer)
            else: ch = self.ch_hash[hhh]
            for x,y in geotools.plotpoly(ch):
                ax.plot(y,x,linewidth=linewidth,color=self.colorhash[hhh],label=label)
        ax.plot(vloc[1],vloc[0],'c^',markersize=20)
        handles, labels = ax.get_legend_handles_labels()
        if legend:
            ax.legend(handles,labels,fontsize=self.fs)
        return handles, labels
        





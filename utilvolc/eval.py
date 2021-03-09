import xarray as xr
import matplotlib.pyplot as plt
import cartopy
import numpy as np

def getDI(x):
    if not x: return False
    if 'DataInsertion' in x: return True
    return False


class UpdateProbs:

    def __init__(self):
        self.probE1 = 0.999
        self.probE2 = 0.001

        self.prior = xr.DataArray()
        self.poly = xr.DataArray()
        self.post = xr.DataArray()

        self.prior_minlev = self.set_prior_minlev(0.001)

    def set_prior_minlev(self,minlev):
        # prior may have 0 values which cannot
        # be updated. Set to low value instead.
        self.prior_minlev = minlev
        

    def set_confidence(self,clevel):
        #probE1 : probability that analyst detects ash above
        #         threshold give that there is ash above threshold
        #probE2 : probability that analyst detects ash above
        #         threshold give that there is NOT ash above threshold


        # perfect will over-ride the prior probability.
 
        # otherwise we assume the analyst is pretty good at detecting ash
        # 90% of the time there is ash, the analyst will find it.

        # The confidence is in how confident the analyst is that the ash
        # is above threshold. If they are confidenct we still give them a
        # 25 % chance that they got it wrong.
        # medium confidence gives them a 50% chance and low confidence
        # a 75% chance. 

        if clevel == 'perfect':
            probE1=0.999
            probE2=0.0001
        elif clevel == 'high':
            probE1=0.90
            probE2=0.25
        elif clevel == 'medium':
            probE1=0.90
            probE2=0.50
        elif clevel == 'low':
            probE1=0.90
            probE2=0.75

        self.probE1 = probE1
        self.probE2 = probE2
  
    def makePE1(self, poly):
        temp = poly * self.probE1
        # creates nan where it is 0.
        temp2 = temp.where(temp!=0)  
        # puts 1's in place where polygon is NOT.
        temp2 = temp2.fillna(1)
        return temp2

    def makePE2(self, poly):
        temp = poly * self.probE2
        # creates nan where it is 0.
        temp2 = temp.where(temp!=0)  
        # puts 1's in place where polygon is NOT.
        temp2 = temp2.fillna(1)
        return temp2

    def process_prior(self, prior):
        """
        replace 0 values with minlev values
        """
        temp = prior.where(prior!=0)
        temp = temp.fillna(self.prior_minlev)
        return temp

    def update_prob(self, prior, poly):
        # poly is just 1's and 0's.
        # This is given by how confident the forecaster is
       
        # poly2 has points with P(A) that is prob
        # that event A will happen under these conditions. 
        prior = self.process_prior(prior,minlev=0.001) 
        probE1 = self.makePE1(poly)
        probE2 = self.makePE2(poly)
        
        val1 = prior * probE1  
        val2 = (1-prior)*probE2
        newprob = val1 / (val1 + val2) 
        return newprob

def make_ensra(xdset):
    dset = xdset.filter_by_attrs(runtype=getDI)
    ensra = dset.to_array(dim='ens')
    return ensra

def get_ATL(xdset,thresh=0.2):
    dset = xdset.filter_by_attrs(runtype=getDI)
    ensra = dset.to_array(dim='ens')
    time = ensra.time.values[0]  
    enslist = ensra.ens.values
    levels= ensra.z.values
    return ATL(ensra, time, enslist,thresh,levels)


def plotATL(rtot, vlist):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    x = rtot.longitude
    y = rtot.latitude
    temp = rtot.max(dim='z')
    zvals = rtot.z.values
    nz = len(zvals)
    if nz>1: 
       nrow= 2  
       ncol = int(np.ceil(nz/ nrow))
    else:
       nrow = 1
       ncol = 1
 
    z = temp.where(temp!=0)
    transform = cartopy.crs.PlateCarree()
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
        cb = ATLsubplot(ax,x,y,z,transform,label,vlist)
        iii+=1
    plt.colorbar(cb) 

def meter2FL(meters):
    return 'FL{:2.0f}'.format(meters/30.48)

def ATLsubplot(ax, x,y,z,transform,label='',vlist=None,
               levels=[1,5,10,15,20]):
    from matplotlib.colors import BoundaryNorm
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
      
    cmap = plt.get_cmap('viridis')
    norm = BoundaryNorm(levels,ncolors=cmap.N,clip=False)
    cb2 = ax.pcolormesh(x,y,z,cmap=cmap,transform=transform,norm=norm)
    ax.plot(vlist[0],vlist[1],'r^')
    ax.add_feature(cartopy.feature.OCEAN) 
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
    #ax.text(0.1,0.1,label,transform=transform) 
    ax.text(0.5,-0.15,label,va='bottom',ha='center',rotation='horizontal',
            rotation_mode='anchor', transform=ax.transAxes,size=15)
    #ax.plot(vlist[0],vlist[1],'r^')
    return cb2 

def ATL(revash, time, enslist, thresh=0.2, level=1):
     """
     Returns array with number of ensemble members above
     given threshold at each location.
     """
     #import matplotlib.pyplot as plt
     #sns.set_style('whitegrid')
     source=0
     if isinstance(level, int):
        level=[level]
     iii=0
     for lev in level:
         r2 = revash.sel(time=time)
         r2 = r2.sel(ens=enslist)
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
     return rtot


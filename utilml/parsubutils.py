import datetime
import os
import sys

#import matplotlib.pyplot as plt
#import textwrap
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
import seaborn as sns
import xarray as xr
from sklearn.cluster import KMeans
from sklearn.mixture import BayesianGaussianMixture as BGM
from sklearn.mixture import GaussianMixture as GMM

#import reventador as reventador
import monet
import utilhysplit.pardump as pardump
from Tools import volcMER


"""

classes


functions

"""


def compare_pdump(pdumpdf, time, fig=None, plot=False):
    """
    pdumpdf : pandas dataframe.
    time : datetime object
    """
    import matplotlib.pyplot as plt
    #lat, lon, vmass, vht = reventador.get_rev_data(time)
    vpi = np.where(vmass>0)
    latgood = lat[vpi]
    longood = lon[vpi]
    buff = 0.005
    latmin = np.min(latgood)
    latmax = np.max(latgood)
    lonmin = np.min(longood)
    lonmax = np.max(longood)
    try:
        df = pdumpdf[pdumpdf['lat']<=latmax+buff]
        df = df[df['lat']>=latmin-buff]
        df = df[df['lon']<=lonmax+buff]
        df = df[df['lon']>=lonmin-buff]
    except:
        df = pdumpdf.copy() 
    if plot:
        plt.title(time.strftime("%Y %m %d %H:%Mz"))
        ax = fig.add_subplot(3,1,1)
        ax2 = fig.add_subplot(3,1,2)
        ax3 = fig.add_subplot(3,1,3)
        ax.pcolormesh(lon, lat, vmass)
        ax2.plot(df['lon'], df['lat'], 'k.',MarkerSize=1)
        ax2.contourf(lon, lat, vmass)
        ax3.plot(df['lon'],df['ht'], 'k.') 
    return df 

def get_thickness(df, t1, t2):
    df2 = df[df['ht']>t1]
    df2 = df2[df2['ht']<=t2]
    return df2


class ParConc:

    def __init__(self,pdict, sp=[1], mult=1,nnn=None,
                 dd=0.1, dh=0.1, bnds=None, method='gmm'):
        self.pdict = pdict
        self.sp = sp
        self.mult = mult
        self.nnn = nnn
        self.dd = dd
        self.dh = dh
        self.bnds = bnds
        self.method = method
        self.concra = xr.DataArray(None)   
        self.mfitlist = []  
        self.pdat = False
        self.pdatdf = pd.DataFrame()
        self.msl = False
        self.maxnnn= nnn
        self.time_average = 0
        self.sourcename = 'unknown'

    def add_source_description(self, sname):
        self.sourcename = sname

    def save(self, outname):
        atthash = {}
        atthash['dd'] = self.dd
        atthash['dh'] = self.dh
        atthash['nG'] = self.nG
        atthash['method'] = self.method
        atthash['mult'] = self.mult
        atthash['species'] = str.join(',',map(str,self.sp))
        atthash['Time_Average (minutes)'] = self.time_average
        atthash['description'] = self.sourcename
        self.concra = self.concra.assign_attributes(atthash)
        self.concra.to_netcdf(outname, format='NETCDF4')

    def load(self, fname):
        dset = xr.open_dataset(fname)
        dra = dset.to_array()
        dra = dra.isel(variable=0)
        return dra

    def set_nnn(self,nnn):
        self.nnn = nnn
        self.maxnnn = nnn 

    def add_pardat_file(self,df):
        self.pdat= True
        self.msl = True 
        self.pdatdf = df

    def set_msl(self):
        if self.pdat: self.msl=True
        else: print('WARNING: cannot set msl. no pdat file')
 
    def make_bnds(self, dfa):
         bnds = {}
         bnds['latmin'] = np.min(dfa['lat'].values)
         bnds['latmax'] = np.max(dfa['lat'].values)
         bnds['lonmin'] = np.min(dfa['lon'].values)
         bnds['lonmax'] = np.max(dfa['lon'].values)
         return bnds 

    def key2datetime(self,time):
        strfmt = "%Y%m%d%H%M"
        return datetime.datetime.strptime(time, strfmt)

    def timeave(self, tmave=60, stime=None, etime=None):
        # NOT DONE
        self.time_average=tmave
        iii=0
        jjj=0
        templist=[]
        atime = stime
        macc=0 # minutes accumulated.
        if not etime:
            tlist = list(self.pdict.keys())
            etime = self.key2datetime(tlist[-1]) 
        self.mfitlist=[] 
        for time in self.pdict.keys():
            print(stime,etime, self.key2datetime(time))
            if iii==0 and not stime: 
               stime = self.key2datetime(time) 
               macc = 0
               atime = stime
            else:
               ptime = self.key2datetime(time) 
               # if before start time then go to nex time
               if ptime < stime: 
                  print('wait')
                  continue
               # if after end time then exit loop
               if ptime > etime: 
                  print('stop')
                  break
               macc = (ptime-atime).seconds/60.0
               #if ((ptime-stime).seconds/60) % delta !=0:
               #   continue               
               if jjj==0: 
                  dfnew = self.pdict[time]
               else:
                  dfnew = pd.concat([dfnew, self.pdict[time]])
               if macc < tmave:
                  print('macc', macc)
                  jjj += 1
                  continue
               else:
                  print('end macc', macc)
                  tw = jjj
                  macc=0
                  jjj=0
                  atime = ptime
            print('Lenght of file', len(dfnew))
            print('time weighting', tw)
            nnn = self.checkn(dfnew)
            mfit, conc = par2conc(dfnew, 
                     self.key2datetime(time),
                     self.sp,
                     self.mult / float(tw),
                     self.method,
                     self.dd, 
                     self.dh, 
                     nnn,  
                     self. bnds)
            if iii==0:
              print('creating concra')
              concra = conc
            else:
              print('appending concra')
              concra = xr.concat([concra, conc],dim='time')
            self.mfitlist.append(mfit)
            iii+=1
            if macc > tmave: 
               print('WARNING: macc too large')
               break 
            macc=0
            jjj=0
        self.concra = concra     
        return  concra  
 

    def checkn(self, df):
        ncheck = len(df)
        if ncheck/ self.nnn > 10:
           rval = self.nnn
        elif ncheck < 50:
           rval = 1
        else:
           rval = ncheck / 50.0
        rval = np.min([rval, self.maxnnn])
        return rval 
        
    def timeloop(self, delta=60, stime=None, etime=None):
        """
        stime : datetime :  first time to extract
        delta : int : extract time every delta minutes
        """
        iii=0
        self.mfitlist=[]
        if not etime:
            tlist = list(self.pdict.keys())
            etime = self.key2datetime(tlist[-1]) 
        for time in self.pdict.keys():
            if iii==0 and not stime: 
               stime = self.key2datetime(time) 
            else:
               ptime = self.key2datetime(time) 
               if ptime < stime: continue
               if ptime > etime: break
               if ((ptime-stime).seconds/60) % delta !=0:
                  continue               
            if self.pdat:
               print('Merging pdat', time, len(self.pdict[time]))
               dfnew = merge_pdat(self.pdict[time], self.pdatdf)
            else:
               dfnew = self.pdict[time]
            print('Lenght of file', len(dfnew))
            nnn = self.checkn(dfnew)
            mfit, conc = par2conc(dfnew, 
                     self.key2datetime(time),
                     self.sp,
                     self.mult,
                     self.method,
                     self.dd, 
                     self.dh, 
                     nnn,  
                     self. bnds)
            if iii==0:
              concra = conc
            else:
              concra = xr.concat([concra, conc],dim='time')
            self.mfitlist.append(mfit)
            iii+=1
        self.concra = concra     
        return  concra  
 

    def check_ones(self):
        for mfit in self.mfitlist:
            for key in mfit.keys():
                print(key, mfit[key].one)        
                print(key, mfit[key].htunits)

    def mass_loading(self, splist=None, thresh=1e-10):
        iii=0
        for time in self.concra.time.values:
            if iii>0: break
            tempra = self.concra.sel(time=time)
            tempra = threshold(tempra)
            if 'sp' in tempra.coords:
                tempra = tempra.sum(dim='sp')
            tempra = tempra.sum(dim='z')
            tempra = tempra * self.dh* 1000
            if iii>0: break
            cb = plt.pcolormesh(tempra.x, tempra.y, tempra)
            plt.colorbar(cb)
            iii+= 1

def threshold(cra, tval=3, tp='log',fillna=True):
    """
    Apply a threshold to a DataArray
    if tp=='log' then tval indicates how many orders of magnitude
    shoudl be retained.
    """
    
    if tp == 'log':
       maxlog = np.max(np.log10(cra)) 
       minlog = np.min(np.log10(cra)) 
       thresh = 10**(maxlog - tval)
    else:
       thresh = tval
    cra = cra.where(cra>=thresh)
    if fillna: cra = cra.fillna(0)
    return cra

def merge_pdat(pardf, datdf):
    """
    pardf is a pandas DataFrame representing pardump file
    datdf is a pandas DataFrame representing a stilt particle.dat file
    """
    #SOMETHING WRONG HERE.
    if 'agl' in pardf.columns.values:
        print('WARNING: switching inputs for merge_pdat')
    #    temp = datdf.copy()
    #    datdf = pardf.copy()
    #    pardf = temp 
    ddf = datdf.copy()
    t1 = pardf.date.unique()
    #print(t1)
    ddf = ddf[ddf['date'].isin(t1)]
    print('lngth ddf', len(pardf), len(ddf), len(datdf))
    if 'lat' in ddf.columns.values: ddf=ddf.drop(['lat'],axis=1)
    if 'lon' in ddf.columns.values: ddf=ddf.drop(['lon'],axis=1)
    #datdf.drop(['lat','lon'],axis=1,inplace=True)
    #merge_cols = ['date','sorti','lat','lon']
    merge_cols = ['date','sorti']
    # merge left means rows in datdf that do not match pardf are left out.
    dfnew = pd.merge(pardf,
                     ddf,
                     how='left',
                     left_on= merge_cols,
                     right_on= merge_cols,
                     )
    print('dfnew length', len(dfnew))
    return dfnew
                   
def par2fit(pdumpdf,  
         mult=1,
         method='gmm',
         nnn=None,
         msl=True,
         wcp=1e3):  #used for Bayesian Gaussian Mixture

    """
    INPUT
    pdumpdf : pandas DataFrame
              pandas DataFrame which has particle locations and mass.
              locations are indicated by 'lon', 'lat', 'ht' columns.
              mass is indicated by 'pmass' column.
              'ht' should be in meters and is converted to km for fitting.
              'lon', 'lat' are in degrees.
        
    method : string
             options are 'gmm', 'bgm', 'kde'
    
    mult : float
           MassFit object has a mass attribute which is computed here by
           summing over the mass on all the particles and multiplying by this
           number.

    nnn  : if method 'gmm'  or 'gbm' number of clusters to use.
           if None will try to figure out for gmm.
           if method 'kde' it is the bandwidth

    wcp  : for bgm

    returns an instance of the MassFit class.
    """

    htmult = 1/1000.0 #convert height to km
    bic=True 
    mfithash = {}
    iii=0
    df2 = pdumpdf.copy()
    #df2 = df2[df2['poll']==species]
    mass = df2['pmass'].sum() * mult
    lon = df2['lon'].values
    lat = df2['lat'].values
    hval = 'ht'
    if 'agl' in df2.columns.values:
        if msl:
            df2['msl'] = df2.apply(lambda row:
                         row['agl'] + row['grdht'], axis=1)
            hval = 'msl'
        else: 
            hval = 'agl'
    ht = df2[hval].values * htmult
    xra = get_xra(lon,lat,ht)
    if method=='gmm':
        if not nnn:
            nnn1,nnn2 = find_n(xra, plot=True)
            print('bic: ', nnn1,' aic: ',nnn2)
            if bic: nnn=nnn1
            else: nnn=nnn2
        if nnn>=len(xra): nnn= np.floor(len(xra)/2.0)
        if nnn==0: nnn=1
        gmm = get_gmm(n_clusters=nnn)
    elif method=='bgm':
        gmm = get_bgm(n_clusters=nnn, wcp=wcp)
    elif method=='kde':
        if not nnn:
           print('par2fit error: bandwidth, (nnn) must be specified for KDE method')
           sys.exit() 
        gmm = get_kde(bandwidth=nnn) 
    else:
        print('Not valid method ', method)
        sys.exit()
    print('par2fit: MassFit', len(xra))
    mfit = MassFit(gmm, xra, mass)
    return mfit



def par2conc(pdumpdf,  
         time=None,
         sp=[1], 
         mult=1,
         method='gmm',
         dd=0.01,
         dh=0.1,
         nnn=None,
         bnds=None,
         thk=None,
         msl=True,
         wcp=1e3):
    """
    Three dimensional. Can use multiple species.
    pdumpdf : pandas dataframe representing pardump file.
    mult : multiplicative factor for particle mass.
    method : kde or gmm (gaussian mixture model)
    nnn  : if method 'gmm' number of clusters to use.
           if None will try to figure out.
           if method 'kde' it is the bandwidth
    dd   : horizontal resolution of output in degrees.
    dh   : vertical resolution of output in km.
    OUTPUT:
    mfithash : dictionary of MassFit objects.
    conctotal : xarray DataArray
    """
    buf=0.2
    htmult = 1/1000.0 #convert height to km
    bnds=None
    bic=True 
    if thk:
        dfa = get_thickness(pdumpdf, thk[0], thk[1])
        thickness = thk[0] - thk[1]
    else:
        dfa = pdumpdf.copy()
    mfithash = {}
    iii=0
    if not bnds:
        bnds = {}
        bnds['latmin'] = np.min(dfa['lat'].values)
        bnds['latmax'] = np.max(dfa['lat'].values)
        bnds['lonmin'] = np.min(dfa['lon'].values)
        bnds['lonmax'] = np.max(dfa['lon'].values)
        bnds['htmin'] = np.min(dfa['ht'].values*htmult)
        bnds['htmax'] = np.max(dfa['ht'].values*htmult)
    for species in sp:
        df2 = dfa[dfa['poll']==species]
        mfit = par2fit(df2,method,dd,dh,nnn,msl,wcp)
        dra = mfit.get_conc(dd, dh, buf, time=time, bnds=bnds, mass=mass)
        #conc = massload / float(thickness)
        if iii==0:
           conctotal = dra
        else:
           conctotal = xr.concat([conctotal,dra],dim='sp')
        iii+=1
        mfithash[(time,iii)] = mfit
    return mfithash, conctotal

def makeplot(lon, lat, conc, levels=None):
    from matplotlib.colors import BoundaryNorm
    if not levels: levels = np.arange(0,5,0.1)
    cmap = plt.get_cmap('PiYG') 
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    cb = plt.pcolormesh(lon, lat, conc, 
                        cmap=cmap, norm=norm)
    plt.colorbar(cb) 

def cdump_plot(cdump, t2, d2, ax=None, title='None'):
    cdump = cdump.sel(z=t2, time=d2)   
    ax = ax or plt.gca()
    levels = list(np.arange(0,5,0.1))
    makeplot(cdump.longitude, cdump.latitude, cdump, levels)
    levels = [0.01,0.5,1,1.5,2,3,4]
    #plt.contour(cdump.longitude, cdump.latitude, cdump, levels=levels)


def draw_ellipse(position, covariance, ax=None, **kwargs):
    """Draw an ellipse with a given position and covariance"""
    from matplotlib.patches import Ellipse
    ax = ax or plt.gca()
    
    # Convert covariance to principal axes
    if covariance.shape == (2, 2):
        U, s, Vt = np.linalg.svd(covariance)
        angle = np.degrees(np.arctan2(U[1, 0], U[0, 0]))
        width, height = 2 * np.sqrt(s)
    else:
        angle = 0
        width, height = 2 * np.sqrt(covariance)
    
    # Draw the Ellipse
    for nsig in range(1, 4):
        ax.add_patch(Ellipse(position, nsig * width, nsig * height,
                             angle, **kwargs))

def get_kde(bandwidth, kernel='gaussian'):
    from sklearn.neighbors import KernelDensity
    kde = KernelDensity(bandwidth=bandwidth, kernel=kernel)
    return kde

def make_bnds(dfa):
     bnds = {}
     bnds['latmin'] = np.min(dfa['lat'].values)
     bnds['latmax'] = np.max(dfa['lat'].values)
     bnds['lonmin'] = np.min(dfa['lon'].values)
     bnds['lonmax'] = np.max(dfa['lon'].values)
     return bnds 

class MassFit():

    def __init__(self,gmm, xra, mass=1):
        """
        gmm : GuassianModelMixture object 
              OR KED or BGM 
             
        xra : list of (longitude, latitude, height) points
        mass: float : mass represented by latitude longitude points
        """
        self.gmm = gmm
        self.xra = xra  
        self.fit=True
        print('MassFit: fitting')
        try:
            self.gfit = gmm.fit(xra)
        except:
            print('MassFit: fitting failed')
            self.fit=False
        print('MassFit: done fitting')
        self.mass = mass 
        self.htunits = self.get_ht_units()

    def get_ht_units(self):
        if self.xra.shape[1] == 3:
           maxht = np.max(self.xra[:,2])
           if maxht > 1000: 
              return 'm'
           else:
              return 'km' 
        else:
           return 'none'

    def err(self):
        """
        not currently useful.
        This give probability that point belongs to cluster
        to which it was assigned.
        """
        resp = self.gfit.predict_proba(self.xra)

    def scatter(self,ax=None, labels=True, dim='ht'):
        """
        create scatter plot
        """
        ax = ax or plt.gca()
        xra = self.xra
        z = self.gfit.predict(xra)
        resp = self.gfit.predict_proba(self.xra)

        if dim=='ht': 
           xxx = xra[:,0]
           yyy = xra[:,1]
        elif dim=='lat': 
           xxx = xra[:,0]
           yyy = xra[:,2]
        elif dim=='lon': 
           xxx = xra[:,1]
           yyy = xra[:,2]

        if labels: 
            ax.scatter(xxx, yyy, c=z,s=1)  
        else:
            ax.scatter(xxx, yyy, c=resp) 
        ax.axis('equal')
        return z


    def plot_gaussians(self, ax=None, dim='ht'):
        """
        plot gaussians as eillipses.
        Needs to be modified for 3d.
        """
        ax = ax or plt.gca()
        gfit = self.gfit
        wfactor = 0.5 /self.gfit.weights_.max()
        for pos, covar, www in zip(gfit.means_, gfit.covariances_, gfit.weights_): 
            draw_ellipse(pos, covar, ax, alpha=www * wfactor) 

    #def get_conc(self, dd, buf=0.2):
    #    if self.thickness:
    #       lonra, latra, massload = self.get_massload(dd, buf)
    #       conc = massload / self.thickness
    #    return lonra, latra, conc


    def get_ht_ra(self, htmin, htmax, dh):
        """

        """
        # first convert to km.
        if self.htunits == 'm':
           htmin = htmin/1000.0
           htmax = htmax/1000.0
           dh = dh / 1000.0
        #htmin = np.floor(htmin)
        #htmax = np.ceil(htmax)
        if htmin <= 0.01:
           htmin = -0.02
        htra = np.arange(htmin, htmax, dh)
        # then convert back to m.
         

        if self.htunits == 'm': htra = htra * 1000.0
        return htra

    def get_lra(self, lmin, lmax, dd):
        """
        helps make sure that multiple grids can be added
        in the xarray.
        TO DO : needs to be generalized to all values of dd.
                dh should be added as well.
        """
        if dd in [0.1,0.05,0.02]:
           qqq=10
        elif dd in [0.01,0.005,0.002]:
           qqq=100
        else:
           qqq=1
        lmin = (np.floor(lmin * qqq)) /qqq
        lmax = (np.ceil(lmax * qqq)) /qqq
        lra = np.arange(lmin, lmax, dd)
        return lra

         
    def totalgrid(self, dd,dh,buf):
        """
        returns lat, lon, ht arrays based
        on data in the xra.
        """
        lon = self.xra[:,0]
        lat = self.xra[:,1]
        ht = self.xra[:,2]
        latmin = np.min(lat)-buf
        latmax = np.max(lat)+buf
        lonmin = np.min(lon)-buf
        lonmax = np.max(lon)+buf
        htmin = np.min(ht)-buf
        htmax = np.max(ht)+buf
        latra = self.get_lra(latmin, latmax, dd)
        lonra = self.get_lra(lonmin, lonmax, dd)
        htra =  self.get_ht_ra(htmin, htmax, dh)
        return latra, lonra, htra

    def partialgrid(self,lat,lon,ht,dd,dh,buf):
        """
        returns lat, lon, ht arrays.
        dd is resolution in horizontal directions
        dh is resolution in vertical direction.
        buf[0] is 1/2 span in horizontal directions.
        buf[1] is 1/2 span in vertical directions.
        Will allow height array to go below 0.

        """
        bufx = buf[0]
        if bufx < dd: bufx = dd
        bufh = buf[1]
        if bufh < dd: bufh = dh
        latmin = lat-bufx
        latmax = lat+dd+bufx
        lonmin = lon-bufx
        lonmax = lon+dd+bufx
        htmin = ht - bufh
        htmax = ht + dh + bufh
        latra = self.get_lra(latmin, latmax, dd)
        lonra = self.get_lra(lonmin, lonmax, dd)
        htra =  self.get_ht_ra(htmin, htmax, dh)
        #latra = np.array([latmin,lat,latmax])
        #lonra = np.array([lonmin,lon,lonmax])
        #htra = np.array([htmin,ht,htmax]) 
        return latra, lonra, htra 
 
    def get_conc(self,
                 dd,
                 dh, 
                 buf=0.2, 
                 bnds=None, 
                 time=None, 
                 mass=None,
                 lat=None,
                 lon=None,
                 ht=None,
                 verbose=False):
        if not mass: 
           mass = self.mass
        if not lat:
           latra, lonra, htra = self.totalgrid(dd,dh,buf) 
        else:
           latra, lonra, htra = self.partialgrid(lat,lon,ht, dd,dh,buf)
        x,y,z = np.meshgrid(lonra, latra, htra)
        xra2 = np.array([x.ravel(), y.ravel(),z.ravel()]).T 
        score = self.gfit.score_samples(xra2)
        one = np.exp(score) * dd**2 * dh
        if verbose: print('ONE is ', one.sum()) 
        self.one = one.sum()
        deg2meter = 111e3
        volume = dh * (dd*deg2meter)**2
        if self.htunits == 'km': volume = volume * 1000.0
        self.conc = np.exp(score)*dd**2 * dh * mass / volume
        corder = [(latra,'y'), (lonra,'x'), (htra,'z')]
        rshape = []
        coords=[]
        dims=[]
        for ccc, dim in corder:
            rshape.append(ccc.shape[0])
            coords.append(ccc)
            dims.append(dim)
        conc2 = self.conc.reshape(rshape)
        dra = xr.DataArray(conc2,
                           coords=coords,
                           dims=dims)
        if time: 
           dra['time'] = time
           dra = self.dra.expand_dims(dim='time') 
        self.dra = monet.monet_accessor._dataset_to_monet(dra,lat_name='y',lon_name='x')
        return self.dra         

    def plotconc(self,dra):
        dra2 = dra.sum(dim='z')
        plt.pcolormesh(dra2)
        
    def get_massload(self, dd, buf=0.2, bnds=None):
        """
        Return mass loading on grid.
        INPUTS
        dd : float : grid spacing
        buf : float : how much to add onto sides of grid.
        bnds : dict : gives box with lat and lon bounds.
        OUTPUTS
        lonra : 2d array
        latra : 2d array
        massload : 2d array
        """
        if bnds:
            latmin = bnds['latmin']
            latmax = bnds['latmax']
            lonmin = bnds['lonmin']
            lonmax = bnds['lonmax']
        else:
            lon = self.xra[:,0]
            lat = self.xra[:,1]
            latmin = np.min(lat)-buf
            latmax = np.max(lat)+buf
            lonmin = np.min(lon)-buf
            lonmax = np.max(lon)+buf
        latra = np.arange(latmin, latmax, dd)
        lonra = np.arange(lonmin, lonmax, dd)
        x,y = np.meshgrid(lonra, latra)
        xra2 = np.array([x.ravel(), y.ravel()]).T 
        score = self.gfit.score_samples(xra2)
        score = score.reshape(len(latra), len(lonra))
        one = np.exp(score) * dd**2
        print('ONE is ', one.sum()) 
        deg2meter = 111e3
        area = (dd*deg2meter)**2
        massload = np.exp(score)*dd**2 * self.mass /area 
        return lonra, latra, massload

def use_gmm(gmm, xra, mass=1, label=True, ax=None):
     ax = ax or plt.gca()
     gfit = gmm.fit(xra)
     labels = gfit.predict(xra)
     if label: 
        ax.scatter(xra[:,0], xra[:,1], c=labels,s=1)  
     else:
        ax.scatter(xra[:,0], xra[:,1]) 
     ax.axis('equal')
     wfactor = 0.5 /gmm.weights_.max()
     #for pos, covar, www in zip(gfit.means_, gfit.covariances_, gfit.weights_): 
     #    draw_ellipse(pos, covar, alpha=www * wfactor) 
     #Z1 = gfit.score_samples(xra)
     #print(Z1.shape, xra.shape)
     #ax.contour(xra[:,0], xra[:,1], Z1)  
     
     lon = xra[:,0]
     lat = xra[:,1]

     latmin = np.min(lat)
     latmax = np.max(lat) 
     lonmin = np.min(lon)
     lonmax = np.max(lon) 
     dd = 0.05
     buf = 0.2 
     latra = np.arange(latmin-buf, latmax+buf, dd)
     lonra = np.arange(lonmin-buf, lonmax+buf, dd)
     #xra2 = np.array(list(zip(lonra,latra)))
     x,y = np.meshgrid(lonra, latra)
     xra2 = np.array([x.ravel(), y.ravel()]).T 
     score = gfit.score_samples(xra2)
     score = score.reshape(len(latra), len(lonra))
     print(score.shape, latra.shape, lonra.shape)
     ## score is the probability density so the following should
     ## sum to 1.
     one = np.exp(score) * dd**2
     print('ONE is ', one.sum()) 
     mass2 = np.exp(score) * dd**2 * mass
     print('mass is ', mass2.sum())
     #clevs = [-10,-5,-1, 0, 1,5,10,50,100]
     deg2meter = 111e3
     area = (dd*deg2meter)**2
     massload = np.exp(score)* dd**2 * mass /area 
     massload = np.exp(score)* dd**2 * mass / area
     #cb = ax.contour(lonra, latra, np.exp(score)* dd**2 * mass / area)
     cb = ax.contour(lonra, latra, massload)
     plt.colorbar(cb)
     return massload, score, gfit



def plot_gmm(lonra, latra, massload):
     clevs = [0.1,1,5,10,20,50]
     if clevs:
         cb = ax.contour(lonra, latra, np.exp(score)* dd**2 * mass / area,
                      levels=clevs)  
     else:
         cb = ax.contour(lonra, latra, np.exp(score)* dd**2 * mass / area,
                        )  
     plt.colorbar(cb)

def find_n(xra, n_max=0, step=1, plot=True):
    fig = plt.figure(1)
    ax = fig.add_subplot(1,1,1)
    min_par_num = 50 
    parnum = len(xra)
    n_max = int(parnum / min_par_num)
    print('find_n function NMAX', n_max, len(xra))
    if n_max==0: n_max = 21
    if n_max>30: n_max=40
    print('find_n function NMAX', n_max, len(xra))
    n_components = np.arange(1, n_max, step)
    models = [GMM(n, covariance_type='full', random_state=0).fit(xra)
              for n in n_components]
    aic = [m.aic(xra) for m in models]
    bic = [m.bic(xra) for m in models]
    if plot:
        ax.plot(n_components, [m.bic(xra) for m in models], label='BIC')
        ax.plot(n_components, [m.aic(xra) for m in models], label='AIC')
        ax.legend(loc='best')
        #ax.xlabel('n_components')
    return  np.argmin(bic), np.argmin(aic)

def get_xra(lon,lat, ht=None):
    if isinstance(ht, np.ndarray):
        xra = np.array(list(zip(lon,lat,ht)))
    else:    
        xra = np.array(list(zip(lon,lat)))

    return xra


def get_bgm(n_clusters=10, wcp=1.0e3, tol=None):
    """
    # wcp: higher number puts more mass in the center and will lead to more
    # more components being active
    """
    covartype = 'full'
    #init_params = 'random' #alternative is kmeans
    init_params = 'kmeans' #alternative is kmeans
    #wcpt = 'dirichlet_process'       #wcp shouldbe (float, float)
    wcpt = 'dirichlet_distribution' #wcp should be float 
    warm_start = False
    verbose = False
    n_init = 1
    max_iter=500
    reg_covar = 1e-5  #default is 1e-6
    if not tol: tol=1e-3
    gmm = BGM(n_components=n_clusters, 
              covariance_type=covartype,
              n_init = n_init,
              weight_concentration_prior_type = wcpt,
              init_params = init_params,
              max_iter=max_iter,
              verbose=verbose,
              reg_covar = reg_covar,
              tol=tol)
    return gmm


def get_gmm(n_clusters=0):
    # possibilities
    # full, tied, diag, spherical
    covartype = 'full'
    # number of initializations to perform. defaults to 1.
    # best results are kept
    n_init = 1
    #
    warm_start = False
    verbose = False
    gmm = GMM(n_components=n_clusters, 
              random_state=42,
              covariance_type=covartype,
              n_init = n_init,
              warm_start = False,
              verbose=verbose)
    return gmm

def cluster_pars(xra, n_clusters=0):
    use_gmm=True
    use_kmeans=False
    parnum = len(xra)
    if n_clusters ==0:
        min_par_num = 50 
        n_clusters = int(parnum / min_par_num)
    #xra = np.array(list(zip(lon,lat)))
    if n_clusters < 1:
       n_clusters = 1
    if n_clusters > 50:
       n_clusters = 50
    if use_kmeans:
        kpredict = KMeans(n_clusters=n_clusters, random_state=0).fit_predict(xra)
        plt.scatter(xra[:,0], xra[:,1], c=kpredict)  
    if use_gmm:
        # possibilities
        # full, tied, diag, spherical
        covartype = 'full'
        # number of initializations to perform. defaults to 1.
        # best results are kept
        n_init = 1
        #
        warm_start = False
        verbose = False
        gmm = GMM(n_components=n_clusters, 
                  random_state=42,
                  covariance_type=covartype,
                  n_init = n_init,
                  warm_start = False,
                  verbose=verbose)
        #score, gmm = plot_gmm(gmm, xra)
        #gmm = GMM(n_components=n_clusters).fit(ra)
        #glabels = gmm.predict(ra)
        #plt.scatter(ra[:,0], ra[:,1], c=glabels)  
    return gmm 

def subdivide_box(bnds, nnn):
    jra = np.arange(0,nnn+1)
    dlen = (bnds[2] - bnds[0]) / float(nnn)
    latra = bnds[0] + jra * dlen

    jra = np.arange(0,nnn+1)
    dlen = (bnds[3] - bnds[1]) / float(nnn)
    lonra = bnds[1] + jra * dlen

def temp():
    rC = RevPars('PARDUMP.C')    
    rC.read_pardump()
    pdict = rC.pdict
    key1 = "201902252000"
    dfc = pdict[key1]
    mpts = par2conc(dfc, obs, 1, 9000)

class VolcPar:
    def __init__(self, fdir='./', fname='PARDUMP.A'):
        fname = fname
        self.tname = os.path.join(fdir, fname)
        self.strfmt = "%Y%m%d%H%M"
        self.start = datetime.datetime.now()
        self.delta = datetime.timedelta(hours=1)
        self.delt=5
        self.ymax=9
        self.pdict = {}

    def read_pardump(self, century=2000):
        pd = pardump.Pardump(fname=self.tname) 
        self.pdict = pd.read(century=century)

    def key2time(self):
        datelist=[]
        for key in self.pdict.keys:
            datelist.append(datetime.datetime.strptime(key, self.strfmt))
        #self.datetlist = datelist
        return datelist

    def getbytime(self,time):
        # returns a dataframe for that time.
        dstr = time.strftime(self.strfmt)
        return self.pdict[dstr]

    def findsource(self, sorti):
        import pandas as pd
        done=False
        iii=0
        while not done:
           d1 = self.start + iii*self.delta
           print('find source', d1)
           df1 = self.getbytime(d1)
           print('Ages', df1['age'].unique())
           # keep only the new particles.
           df1 = df1[df1['age']==self.delt]         
           
           # if no particles released then
           if df1.empty and iii!=0: 
              done=True
              print('empty', iii)
              continue
           df1 = df1[df1['sorti'].isin(sorti)]
           if iii==0:
              dfsource = df1
           else:
              dfsource = pd.concat([dfsource, df1])
           iii+=1
        return dfsource

    def plotsource(self, dfhash):
        x = []
        y = []
        for key in dfhash.keys():
            sorti = dfhash[key]['sorti']            
            dfsource = self.findsource(sorti)
            x.extend(dfsource['date'])
            y.extend(dfsource['ht'])
        x2 = self.time2int(x)
        # put time into minutes since start
        x2 = np.array(x2)/60.0     
        # put height into km
        y = np.array(y)/1000.0
        xbins = np.arange(0,200,5)
        ybins = np.arange(0,4*9)/4.0
        cb = plt.hist2d(x2,y, bins=[xbins,ybins])
        plt.colorbar(cb[3])
        #sns.heatmap(x,y)
        return x2, y

    def time2int(self, timelist):
        newlist = []
        tmin = np.min(timelist)
        for ttt in timelist:
            val = (ttt-tmin)
            newlist.append(val.seconds)
        return newlist 


def average_mfitlist(mfitlist,masslist,dd=None,dh=None,buf=None,lat=None, lon=None, ht=None):
    iii=0
    for mfit in mfitlist:
        conc = mfit.get_conc(dd=dd, 
                           dh=dh, 
                           buf=buf,
                           mass=masslist[iii],
                           lat=lat,
                           lon=lon,
                           ht=ht)
        if iii==0:
           concra = conc
        else:
           concra = xr.concat([concra, conc],dim='time')
        iii+=1
    try:
        concra = concra.mean(dim='time')
    except:
        pass
    return concra

def key2datetime(time):
    strfmt = "%Y%m%d%H%M"
    return datetime.datetime.strptime(time, strfmt)

def datetime2key(time):
    strfmt = "%Y%m%d%H%M"
    return time.strftime(strfmt)

def subset_pdict(pdict,stime,tmave=60,verbose=False):
    plist = combine_pdict(pdict,stime,tmave,makelist=True, verbose=verbose)
    return plist 

def combine_pdict(pdict, stime, tmave=60, makelist=False, verbose=False):
    """
    pdict : dictionary of dataframes
    stime : datetime
    tmave : integer (minutes)
    outputs

    """
    iii=0
    jjj=0
    templist=[]
    etime = stime + datetime.timedelta(minutes=tmave)
    skey = datetime2key(stime)
    ekey = datetime2key(etime)
    print('combine_pdict ENDkey', ekey, stime, etime, tmave)
    keylist = list(pdict.keys())
    keylist.sort()
    if skey in pdict.keys():
        firsti = keylist.index(skey)
        dfnew = pdict[keylist[firsti]]
    else:
        print('WARNING: key does not exist', skey)
        print('WARNING: Using ', keylist[0])
        firsti = 0
        dfnew = pdict[keylist[firsti]]
    if ekey in pdict.keys():
        lasti = keylist.index(ekey)
    else:
        print('WARNING: key does not exist', ekey)
        print('WARNING: Using ', keylist[-1])
        lasti=len(keylist)
    #else:
    #    print('WARNING: key does not exist')
    #    return 0, pd.DataFrame()
        #continue   
    jjj=1
    for iii in range(firsti+1, lasti):
        key = keylist[iii]
        if verbose: print(key)
        if key in pdict.keys():
            dfnew = pd.concat([dfnew, pdict[key]])
        else:
            print('WARNING: key does not exist')
        templist.append(pdict[key])
        jjj+=1
    if makelist:
       return templist
    else:
       return jjj, dfnew

class KasPar(VolcPar): 
    def __init__(self, fname='PARDUMP.A', fdir=None):
        if not fdir:
            fdir = '/pub/Scratch/alicec/KASATOCHI/2020/cylens/'
        fname = fname
        self.tname = os.path.join(fdir, fname)
        self.strfmt = "%Y%m%d%H%M"
        self.start = datetime.datetime(2008,8,10,16,00)
        self.delta = datetime.timedelta(hours=1)
        self.delt=5
        self.ymax=9
        self.pdict = {}

def height_correction(zaprime, zter, msl=True):
    zmdl = 25000
    za = zaprime * (zmdl-zter)/float(zmdl)
    if msl: za += zter
    return za

class RevPar(VolcPar):
    def __init__(self, fname='PARDUMP.A', fdir=None):
        if not fdir:
            fdir = '/pub/ECMWF/JPSS/reventador/pens/'
        fname = fname
        self.tname = os.path.join(fdir, fname)
        self.strfmt = "%Y%m%d%H%M"
        self.start = datetime.datetime(2019,2,25,16,30)
        self.delta = datetime.timedelta(seconds=5*60)
        self.delt=5

    def get_data(self, time):
        lat, lon, dmass, dht = get_rev_data(time)

def process_under(under):
    lastz = under.z.values[-1]
    under = under.sum(dim='z')
    under = under.assign_coords(z=lastz)
    under = under.expand_dims('z') 
    return under


def shift_underground(dra):
    """
    Takes all mass that is underground and puts it in the first level.
    """
    import math
    iii=0
    height=1e-12
    for val in dra.z.values:
        if height < val: break
        if math.isclose(height,val): break
        iii+=1 
    under = dra.isel(z=np.arange(0,iii+1))
    under = process_under(under)
    above = dra.isel(z=np.arange(iii+1,len(dra.z.values),1))
    new = xr.concat([under,above],dim='z')    
    return new


def reflect_underground(dra):
    """
    reflects mass that is in first 3 levels underground
    onto first three levels above ground.
    """
    import math
    iii=0
    height=0
    for val in dra.z.values:
        if height < val: break
        if math.isclose(height,val): break
        iii+=1 
    print('INDEX', iii)
    under1 = dra.isel(z=[iii-1,iii])
    print('under1', under1)
    under1 = process_under(under1)
    under = under1
    jjj=iii+1
    if iii-2 > 0:
       under2 = dra.isel(z=[iii-2,iii+1])
       print('under2', under2)
       under2 = process_under(under2)
       jjj=iii+2
       under = xr.concat([under,under2],dim='z')
    if iii-3 > 0:
       under3 = dra.isel(z=[iii-3,iii+2])
       print('under3', under3)
       under3 = process_under(under3)
       jjj=iii+3
       under = xr.concat([under,under3],dim='z')


    #under = under1
    print('under', under)
    above = dra.isel(z=np.arange(jjj,len(dra.z.values),1))
    print('above',above.z.values)

    #lastz = under.z.values[-1]
    #under = under.sum(dim='z')
    
    #under = under.assign_coords(z=lastz)
    #under = under.expand_dims('z') 
    
    new = xr.concat([under,above],dim='z')    
 
    return new
    

def xslice(dra, longitude ):
    latitude = dra.latitude.values[0][0]
    xi, yi = dra.monet.nearest_ij(lat=latitude, lon=longitude) 
    return dra.isel(x=xi)

def yslice(dra, latitude ):
    longitude = dra.longitude.values[0][0]
    xi, yi = dra.monet.nearest_ij(lat=latitude, lon=longitude) 
    return dra.isel(y=yi)

def zslice(dra, height ):
    import math
    iii=0
    for val in dra.z.values:
        if height < val: break
        if math.isclose(height,val): break
        iii+=1 
        print(iii, val, height)
    return dra.isel(z=iii)

from sklearn.mixture import BayesianGaussianMixture as BGM
from sklearn.mixture import GaussianMixture as GMM
from sklearn.cluster import KMeans
from scipy.stats import multivariate_normal
import monet
import matplotlib.pyplot as plt
import monetio.models.pardump as pardump
from Tools import volcMER
import sys
import os
import seaborn as sns
import datetime
import pandas as pd
import xarray as xr
#import matplotlib.pyplot as plt
#import textwrap
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import numpy as np
import numpy.ma as ma
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
#import reventador as reventador


"""
Functions and classes to convert pardump file
output to concentrations using mixture models or 
kernel density estimation.
the KDE functionality is very basic.


classes
-------------
MassFit
Par2Conc

functions:
-------------
par2fit : returns a MassFit object.
draw_ellipse
get_kde
get_gmm
get_bgm
get_xra
find_n

average_mfitlist
threshold
make_bnds

use_gmm

plot_gmm
par2conc
get_thickness
merge_pdat
makeplot
cdump_plot

"""


def get_thickness(df, t1, t2):
    df2 = df[df['ht'] > t1]
    df2 = df2[df2['ht'] <= t2]
    return df2


def threshold(cra, tval=3, tp='log', fillna=True):
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
    cra = cra.where(cra >= thresh)
    if fillna:
        cra = cra.fillna(0)
    return cra


def merge_pdat(pardf, datdf):
    """
    pardf is a pandas DataFrame representing pardump file
    datdf is a pandas DataFrame representing a stilt particle.dat file
    """
    # SOMETHING WRONG HERE.
    if 'agl' in pardf.columns.values:
        print('WARNING: switching inputs for merge_pdat')
    #    temp = datdf.copy()
    #    datdf = pardf.copy()
    #    pardf = temp
    ddf = datdf.copy()
    t1 = pardf.date.unique()
    # print(t1)
    ddf = ddf[ddf['date'].isin(t1)]
    print('lngth ddf', len(pardf), len(ddf), len(datdf))
    if 'lat' in ddf.columns.values:
        ddf = ddf.drop(['lat'], axis=1)
    if 'lon' in ddf.columns.values:
        ddf = ddf.drop(['lon'], axis=1)
    # datdf.drop(['lat','lon'],axis=1,inplace=True)
    #merge_cols = ['date','sorti','lat','lon']
    merge_cols = ['date', 'sorti']
    # merge left means rows in datdf that do not match pardf are left out.
    dfnew = pd.merge(pardf,
                     ddf,
                     how='left',
                     left_on=merge_cols,
                     right_on=merge_cols,
                     )
    print('dfnew length', len(dfnew))
    return dfnew


def par2fit(pdumpdf,
            mult=1,
            method='gmm',
            pfit=None,
            nnn=None,
            msl=True,
            wcp=1e3):  # used for Bayesian Gaussian Mixture
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
             p_bgm use a previous bgm fit (input with pfit).
             The function will make a copy of the pfit and set
             warm_start=True.


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

    htmult = 1/1000.0  # convert height to km
    bic = True
    mfithash = {}
    iii = 0
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
    xra = get_xra(lon, lat, ht)
    if method == 'gmm':
        if not nnn:
            nnn1, nnn2 = find_n(xra, plot=True)
            print('bic: ', nnn1, ' aic: ', nnn2)
            if bic:
                nnn = nnn1
            else:
                nnn = nnn2
        if nnn >= len(xra):
            nnn = np.floor(len(xra)/2.0)
        if nnn == 0:
            nnn = 1
        gmm = get_gmm(n_clusters=nnn)
    elif method == 'bgm':
        gmm = get_bgm(n_clusters=nnn, wcp=wcp)
    elif method == 'kde':
        if not nnn:
            print('par2fit error: bandwidth, (nnn) must be specified for KDE method')
            sys.exit()
        gmm = get_kde(bandwidth=nnn)
    elif method == 'p_bgm':
        gmm = copy_fit(pfit, method='bgm')
        gmm.warm_start = True
    elif method == 'p_gmm':
        gmm = copy_fit(pfit, method='gmm')
        gmm.warm_start = True
    else:
        print('Not valid method ', method)
        sys.exit()
    print('par2fit: MassFit', len(xra))
    mfit = MassFit(gmm, xra, mass)
    return mfit


class ParArgs:

    def __init__(self,
                 stime,
                 tmave,
                 splist=None,
                 sorti=None,
                 htmin=None,
                 htmax=None):
        self.stime = stime
        self.tmave = tmave
        self.splist = splist
        self.sorti = sorti
        self.htmin = htmin
        self.htmax = htmax


def fixlondf(df, colname='lon', neg=True):
    if not neg:
        df[colname] = df.apply(lambda row: fixlon(row[colname]), axis=1)
    else:
        df[colname] = df.apply(lambda row: fixlon(row[colname])-360, axis=1)
    return df


#def fixlonra(ra):
#    newlon = []
#    for lon in ra:
#        newloeturn newlon


def fixlon(x):
    if x < 0:
        return 360+x
    else:
        return x


class Par2Conc:

    def __init__(self, df):
        self.df = df      # pandas DataFrame
        self.fitlist = []  # collection of MassFit objects.
        self.dra = None   # array with concentrations.

    def choose_time(self,
                    time_average=None):
        """ 
        time_average minutes
        """
        self.tmave = time_average

    def subsetdf(self,
                 stime,
                 tmave,
                 splist=None,
                 sorti=None,
                 htmin=None,
                 htmax=None):
        #jjj, dfnew = combine_pdict(self.pdict, pd.to_datetime(date), tmave)
        pardf = self.df.copy()
        d1 = stime
        d2 = d1 + datetime.timedelta(minutes=tmave)
        pdn = pardf[pardf.date < d2]
        pdn = pdn[pdn.date >= d1]
        if splist:
            pdn = pdn[pdn.poll.isin(splist)]
        if sorti:
            pdn = pdn[pdn.sorti.isin(sorti)]
        if htmax:
            pdn = pdn[pdn.ht <= htmax]
        if htmin:
            pdn = pdn[pdn.ht >= htmin]
        # make sure longitudes aren't split.
        pdn = fixlondf(pdn)
        return pdn

    def timeloop(self, timestep, parargs):
        iii == 0
        mfitlist = []
        for time in timelist:
            df = subsetdf(time, tmave, parargs)
            if iii == 0:
                mfit = par2fit(df, method='bgm', nnn=50)
            else:
                mfit = par2fit(df, method='p_bgm', pfit=pfit, nnn=50)
            pfit = mfit.gfit
            mfitlist.append(mfit)

    def addfit(self,
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

        return -1


def fit_timeloop(pardf, nnn, maxht=None, mlist=None, method='gmm',
                 warm_start=True,
                 ):
    """
    pardf : dataframe with particle positions.
    Returns:
    submlist : list of MassFit objects.

    creates a separate fit for each time period in the pardf file.

    """
    jjj = 0
    submlist = []
    pmethod = 'p_' + method
    #masslist = []
    # fit each unique date in the period.
    print('in sub')
    dlist = pardf.date.unique()
    dlist.sort()
    for ndate in dlist:
        pdn = pardf[pardf.date == ndate]
        print('sub', pdn)
        if maxht:
            pdn = pdn[pdn['ht'] < maxht]
        if not mlist:
            if jjj == 0:
                mfit = par2fit(pdn, nnn=nnn, method=method)
                print('sub first', mfit)
            else:
                mfit = par2fit(pdn, nnn=nnn, method=pmethod, pfit=pfit)
                print('sub mfit with previous', mfit)
        else:
            mfit = mlist[jjj]
        if not mfit.fit:
            continue
        pfit = mfit.gfit
        # masslist.append(pdn['pmass'].sum())
        submlist.append(mfit)
        if warm_start:
            jjj += 1
    return submlist


def makeplot(lon, lat, conc, levels=None):
    from matplotlib.colors import BoundaryNorm
    if not levels:
        levels = np.arange(0, 5, 0.1)
    cmap = plt.get_cmap('PiYG')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    cb = plt.pcolormesh(lon, lat, conc,
                        cmap=cmap, norm=norm)
    plt.colorbar(cb)


def cdump_plot(cdump, t2, d2, ax=None, title='None'):
    cdump = cdump.sel(z=t2, time=d2)
    ax = ax or plt.gca()
    levels = list(np.arange(0, 5, 0.1))
    makeplot(cdump.longitude, cdump.latitude, cdump, levels)
    levels = [0.01, 0.5, 1, 1.5, 2, 3, 4]
    #plt.contour(cdump.longitude, cdump.latitude, cdump, levels=levels)


def draw_ellipse(position, covariance, ax=None, **kwargs):
    """Draw an ellipse with a given position and covariance"""
    from matplotlib.patches import Ellipse
    ax = ax or plt.gca()
    print(position, covariance)
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
        print('HERE', nsig, width, height, angle)
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

    """
    For gmm
    If we can keep giving it different xra values then we can 
    """

    def __init__(self, gmm, xra, mass=1):
        """
        gmm : GuassianModelMixture object 
              OR KED or BGM 

        xra : list of (longitude, latitude, height) points
        mass: float : mass represented by latitude longitude points
        """
        self.gmm = gmm
        self.xra = xra
        self.fit = True
        print('MassFit: fitting')
        try:
            self.gfit = gmm.fit(xra)
        except:
            print('MassFit: fitting failed')
            self.fit = False
        print('MassFit: done fitting')
        self.mass = mass
        self.htunits = self.get_ht_units()

    def get_ht_units(self):
        if self.xra.shape[1] == 3:
            maxht = np.max(self.xra[:, 2])
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

    def scatter(self, ax=None, labels=True, dim='ht'):
        """
        create scatter plot
        """
        ax = ax or plt.gca()
        xra = self.xra
        z = self.gfit.predict(xra)
        resp = self.gfit.predict_proba(self.xra)

        if dim == 'ht':
            xxx = xra[:, 0]
            yyy = xra[:, 1]
        elif dim == 'lat':
            xxx = xra[:, 0]
            yyy = xra[:, 2]
        elif dim == 'lon':
            xxx = xra[:, 1]
            yyy = xra[:, 2]

        if labels:
            ax.scatter(xxx, yyy, c=z, s=1)
        else:
            ax.scatter(xxx, yyy, c=resp)
        ax.axis('equal')
        return z

    def plot_means(self, dim='ht'):
        for pos, covar, www in self.generate_gaussians(dim):
            plt.scatter(pos[0], pos[1], s=100*www)
            print(www)

    def quickconc(self, mult):
        # quick conc is SLOW.
        # Create empty xarray with  latitude - longitude grid
        df = pd.DataFrame()

        dlist = []
        latlist = []
        lonlist = []
        zlist = []

        buflat = 1
        buflon = 1
        bufht = 1
        ddd = 0.25 / 2.0
        fdd = ddd * 2.0



        # get positions of centers of Gaussians
        zlist = self.sort_fits()
        means = list(zip(*zlist))[0]
        temp = list(zip(*means))
        lats = temp[1]
        lons = temp[0]
        hts = temp[2]
        latr = [np.min(lats), np.max(lats)]
        lonr = [np.min(lons), np.max(lons)]
        htr = [np.min(hts), np.max(hts)]

        # find lower left corner.
        latmin = np.floor(latr[0]) - buflat
        lonmin = np.floor(lonr[0]) - buflon
        zmin = np.floor(htr[0]) - bufht

        latmax = np.ceil(latr[1]) + buflat
        lonmax = np.ceil(lonr[1]) + buflon
        zmax = np.ceil(htr[1]) + bufht

        lat = latmin
        lon = lonmin
        ht = zmin

        done = False
        iii = 0
        while not done:
            if iii % 100 == 0:
                print(iii, lon, lat, ht)
            conc = self.calcconc(ddd, ddd, pos=[lon, lat, ht], check=False,
                                 verbose=False)
            dlist.append(conc*mult)
            iii += 1
            lat += ddd * 2
            if lat > latmax:
                lat = latmin
                lon += ddd * 2
                if lon > lonmax:
                    lon = lonmin
                    if iii > 5000000:
                        done = True
                        print('stopping', lat, lon, ht)
                        zmax = ht
                    ht += ddd * 2
                    if ht > zmax:
                        done = True
        latlist = np.arange(latmin, latmax + fdd, fdd)
        lonlist = np.arange(lonmin, lonmax + fdd, fdd)
        zlist = np.arange(zmin, zmax + fdd, fdd)

        dra = np.array(dlist)
        dra = dra.reshape(len(zlist),len(lonlist),len(latlist))
        dset = xr.DataArray(dra, coords=[zlist,lonlist,latlist],
                            dims=['z','longitude','latitude',])

        return dset

    def get_sublist(self, pos):
        zlist = self.sort_fits()
        cutoff = 1
        for mean, cov, www in zlist:
            ddd=0
            for iii in [0,1,2]:
                ddd += (mean[iii]-pos[iii])**2 
                 
    def getcenter(self):
        zlist = self.sort_fits()
        mean = zlist[0][0]
        lat = mean[1] 
        lon = mean[0]
        ht = mean[2]  
        return (lon,lat,ht)

    def getcenters(self):
        # get positions of centers of Gaussians
        zlist = self.sort_fits()
        means = list(zip(*zlist))[0]
        temp = list(zip(*means))
        lats = temp[1]
        lons = temp[0]      
        hts = temp[2]
        latr = [np.min(lats),np.max(lats)]
        lonr = [np.min(lons),np.max(lons)]
        htr = [np.min(hts),np.max(hts)]
      
    def make_pos_list(self, dd, dh):
        zlist = self.sort_fits()
        plist = []
        means = list(zip(*zlist))[0]
        dblat = 1
        dblon = 1
        dbz = 1               
        

 

    def calcconc(self, dd, dh, pos, check=True, verbose=True):
        """
        calculates concentration in volume centerd on (lon,lat,ht) position
        width, lenght, height defined by dd and dh.
        """
        # calculate probablitilies with CDF's.
        iii=0
        totalprob = 0
        lat = pos[1]
        lon = pos[0]
        ht = pos[2]

        #if verbose: print('Position', lat,lon,ht)
        if check:
            # retrieve values at sampling grid locations.
            score = self.gfit.score_samples([[lon, lat, ht]])
            # this should approximate the "area" under the curve.
            # which for small volumes should be very close
            # to what you get from using the CDF method.
            prob = np.exp(score) * (2*dd)**2 * 2*dh

        for mean, cov, www in self.generate_gaussians(dim=None):
            # The cdf method computes the cumulative distribution function.
            # at various points. For a 1d Gaussian - two points
            # cdf at point 1 gives area under the curve from -infinity to that
            # point.
            # cdf at point 2 gives area under curve from -infinity side to that  point.
            # For a 2d Gaussian - need 4 points.
            # For 3d Gaussian - need 8 points
            # and finds "area" by subtracting cdf at 
            volume = getvolume(mean,cov,lon,lat,ht,dd,dh, verbose) 
            #if verbose: print('prob, volume, weight', volume * www, volume, www)
            #if verbose: print('\n--------------')
            totalprob += volume * www
            #iii+=1

        deg2meter = 111e3
        vol = 2*dh*1000 * (2*dd*deg2meter)**2
        conc = totalprob * self.mass / vol

        if check:
            if dd <= 0.05:
                vv = 0.005
            if dd <= 0.1:
                vv = 0.005
            elif dd <= 0.5:
                vv = 0.01
            elif dd <= 1:
                vv = 0.01
            else:
                vv = 0.01
            conc2 = self.get_conc(vv, vv, buf=[dd, dh], lat=lat, lon=lon, ht=ht)
            print('get_conc function', conc2.mean(), self.one)
            print('conc this function', conc, totalprob)
            print('conc 1 pt estimation ', prob * self.mass / vol, prob)
        return conc

    def sort_fits(self):
        gfit = self.gfit
        zlist = zip(gfit.means_, gfit.covariances_, gfit.weights_)
        zzz = list(zlist)
        # sort according to weight, biggest to smallest.
        zlist = sorted(zzz, key=lambda x: x[2], reverse=True)
        return zlist

    def generate_gaussians(self, dim='ht'):
        """
        """
        gfit = self.gfit
        wfactor = 0.5 / self.gfit.weights_.max()
        if dim == 'ht':
            c1 = 0
            c2 = 1
        if dim == 'lon':
            c1 = 1
            c2 = 2
        if dim == 'lat':
            c1 = 0
            c2 = 2
        zlist = self.sort_fits()
        for pos, covar, www in zlist:
            if dim:
                position = np.array([pos[c1], pos[c2]])
                one = covar[c1][c1]
                two = covar[c1][c2]
                three = covar[c2][c1]
                four = covar[c2][c2]
                covariance = np.array([[one, two], [three, four]])
            else:
                position = pos
                covariance = covar
            yield position, covariance, www

    def plot_gaussians(self, ax=None, dim='ht'):
        """
        plot gaussians as ellipses.
        """
        ax = ax or plt.gca()
        gfit = self.gfit
        wfactor = 0.5 / self.gfit.weights_.max()
        if dim == 'ht':
            c1 = 0
            c2 = 1
        if dim == 'lon':
            c1 = 1
            c2 = 2
        if dim == 'lat':
            c1 = 0
            c2 = 2
        for pos, covar, www in zip(gfit.means_, gfit.covariances_, gfit.weights_):
            position = np.array([pos[c1], pos[c2]])
            one = covar[c1][c1]
            two = covar[c1][c2]
            three = covar[c2][c1]
            four = covar[c2][c2]
            covariance = np.array([[one, two], [three, four]])
            draw_ellipse(position, covariance, ax, alpha=www * wfactor)

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

        if self.htunits == 'm':
            htra = htra * 1000.0
        return htra

    def get_lra(self, lmin, lmax, dd):
        """
        helps make sure that multiple grids can be added
        in the xarray.
        TO DO : needs to be generalized to all values of dd.
                dh should be added as well.
        """
        if dd in [0.1, 0.05, 0.02]:
            qqq = 10
        elif dd in [0.01, 0.005, 0.002]:
            qqq = 100
        else:
            qqq = 1
        lmin = (np.floor(lmin * qqq)) / qqq
        lmax = (np.ceil(lmax * qqq)) / qqq
        lra = np.arange(lmin, lmax, dd)
        return lra

    def totalgrid(self, dd, dh, buf):
        """
        returns lat, lon, ht arrays based
        on data in the xra.
        """
        lon = self.xra[:, 0]
        lat = self.xra[:, 1]
        ht = self.xra[:, 2]
        latmin = np.min(lat)-buf
        latmax = np.max(lat)+buf
        lonmin = np.min(lon)-buf
        lonmax = np.max(lon)+buf
        htmin = np.min(ht)-buf
        htmax = np.max(ht)+buf
        latra = self.get_lra(latmin, latmax, dd)
        lonra = self.get_lra(lonmin, lonmax, dd)
        htra = self.get_ht_ra(htmin, htmax, dh)
        return latra, lonra, htra

    def partialgrid(self, lat, lon, ht, dd, dh, buf):
        """
        returns lat, lon, ht arrays.
        dd is resolution in horizontal directions
        dh is resolution in vertical direction.
        buf[0] is 1/2 span in horizontal directions.
        buf[1] is 1/2 span in vertical directions.
        Will allow height array to go below 0.

        """
        bufx = buf[0]
        if bufx < dd:
            bufx = dd
        bufh = buf[1]
        if bufh < dd:
            bufh = dh
        latmin = lat-bufx
        latmax = lat+dd+bufx
        lonmin = lon-bufx
        lonmax = lon+dd+bufx
        htmin = ht - bufh
        htmax = ht + dh + bufh
        latra = self.get_lra(latmin, latmax, dd)
        lonra = self.get_lra(lonmin, lonmax, dd)
        htra = self.get_ht_ra(htmin, htmax, dh)
        #latra = np.array([latmin,lat,latmax])
        #lonra = np.array([lonmin,lon,lonmax])
        #htra = np.array([htmin,ht,htmax])
        return latra, lonra, htra

    def get_conc(self,
                 dd,
                 dh,
                 buf=0.2,
                 time=None,
                 mass=None,
                 lat=None,
                 lon=None,
                 ht=None,
                 verbose=False):

        if not mass:
            mass = self.mass

        # if lat,lon,ht not specified then find conc over whole area.
        if not lat:
            latra, lonra, htra = self.totalgrid(dd, dh, buf)
        else:
            latra, lonra, htra = self.partialgrid(lat, lon, ht, dd, dh, buf)

        # create the sampling grid
        x, y, z = np.meshgrid(lonra, latra, htra)
        xra2 = np.array([x.ravel(), y.ravel(), z.ravel()]).T

        # retrieve values at sampling grid locations.
        score = self.gfit.score_samples(xra2)

        # multipy by volume to get probability in that volume.
        prob = np.exp(score) * dd**2 * dh

        # sum over whole volume should be one
        self.one = prob.sum()
        if verbose:
            print('ONE is ', self.one)

        # multiply probability by total mass
        # to get the mass in that volume.
        # Divide by volume (in m^3) to get
        # concentration.
        deg2meter = 111e3
        volume = dh * (dd*deg2meter)**2
        if self.htunits == 'km':
            volume = volume * 1000.0

        #self.conc = np.exp(score)*dd**2 * dh * mass / volume
        self.conc = prob * mass / volume

        # reshape the array.
        corder = [(latra, 'y'), (lonra, 'x'), (htra, 'z')]
        rshape = []
        coords = []
        dims = []
        for ccc, dim in corder:
            rshape.append(ccc.shape[0])
            coords.append(ccc)
            dims.append(dim)
        conc2 = self.conc.reshape(rshape)

        # create the xarray
        dra = xr.DataArray(conc2,
                           coords=coords,
                           dims=dims)
        # add time dimensions
        if time:
            dra['time'] = time
            dra = self.dra.expand_dims(dim='time')

        # modify xarray.
        self.dra = monet.monet_accessor._dataset_to_monet(dra, lat_name='y', lon_name='x')

        return self.dra

    def plotconc(self, dra):
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
            lon = self.xra[:, 0]
            lat = self.xra[:, 1]
            latmin = np.min(lat)-buf
            latmax = np.max(lat)+buf
            lonmin = np.min(lon)-buf
            lonmax = np.max(lon)+buf
        latra = np.arange(latmin, latmax, dd)
        lonra = np.arange(lonmin, lonmax, dd)
        x, y = np.meshgrid(lonra, latra)
        xra2 = np.array([x.ravel(), y.ravel()]).T
        score = self.gfit.score_samples(xra2)
        score = score.reshape(len(latra), len(lonra))
        one = np.exp(score) * dd**2
        print('ONE is ', one.sum())
        deg2meter = 111e3
        area = (dd*deg2meter)**2
        massload = np.exp(score)*dd**2 * self.mass / area
        return lonra, latra, massload


def use_gmm(gmm, xra, mass=1, label=True, ax=None):
    ax = ax or plt.gca()
    gfit = gmm.fit(xra)
    labels = gfit.predict(xra)
    if label:
        ax.scatter(xra[:, 0], xra[:, 1], c=labels, s=1)
    else:
        ax.scatter(xra[:, 0], xra[:, 1])
    ax.axis('equal')
    wfactor = 0.5 / gmm.weights_.max()
    # for pos, covar, www in zip(gfit.means_, gfit.covariances_, gfit.weights_):
    #    draw_ellipse(pos, covar, alpha=www * wfactor)
    #Z1 = gfit.score_samples(xra)
    #print(Z1.shape, xra.shape)
    #ax.contour(xra[:,0], xra[:,1], Z1)

    lon = xra[:, 0]
    lat = xra[:, 1]

    latmin = np.min(lat)
    latmax = np.max(lat)
    lonmin = np.min(lon)
    lonmax = np.max(lon)
    dd = 0.05
    buf = 0.2
    latra = np.arange(latmin-buf, latmax+buf, dd)
    lonra = np.arange(lonmin-buf, lonmax+buf, dd)
    #xra2 = np.array(list(zip(lonra,latra)))
    x, y = np.meshgrid(lonra, latra)
    xra2 = np.array([x.ravel(), y.ravel()]).T
    score = gfit.score_samples(xra2)
    score = score.reshape(len(latra), len(lonra))
    print(score.shape, latra.shape, lonra.shape)
    # score is the probability density so the following should
    # sum to 1.
    one = np.exp(score) * dd**2
    print('ONE is ', one.sum())
    mass2 = np.exp(score) * dd**2 * mass
    print('mass is ', mass2.sum())
    #clevs = [-10,-5,-1, 0, 1,5,10,50,100]
    deg2meter = 111e3
    area = (dd*deg2meter)**2
    massload = np.exp(score) * dd**2 * mass / area
    massload = np.exp(score) * dd**2 * mass / area
    #cb = ax.contour(lonra, latra, np.exp(score)* dd**2 * mass / area)
    cb = ax.contour(lonra, latra, massload)
    plt.colorbar(cb)
    return massload, score, gfit


def plot_gmm(lonra, latra, massload):
    clevs = [0.1, 1, 5, 10, 20, 50]
    if clevs:
        cb = ax.contour(lonra, latra, np.exp(score) * dd**2 * mass / area,
                        levels=clevs)
    else:
        cb = ax.contour(lonra, latra, np.exp(score) * dd**2 * mass / area,
                        )
    plt.colorbar(cb)


def find_n(xra, n_max=0, step=1, plot=True):
    fig = plt.figure(1)
    ax = fig.add_subplot(1, 1, 1)
    min_par_num = 50
    parnum = len(xra)
    n_max = int(parnum / min_par_num)
    print('find_n function NMAX', n_max, len(xra))
    if n_max == 0:
        n_max = 21
    if n_max > 30:
        n_max = 40
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
        # ax.xlabel('n_components')
    return np.argmin(bic), np.argmin(aic)


def get_xra(lon, lat, ht=None):
    if isinstance(ht, np.ndarray):
        xra = np.array(list(zip(lon, lat, ht)))
    else:
        xra = np.array(list(zip(lon, lat)))

    return xra


def compare_fits(fit1, fit2, method='gmm'):
    return -1


def copy_fit(bgm, method='bgm'):
    n_clusters = bgm.n_components
    covartype = bgm.covariance_type
    n_init = bgm.n_init
    max_iter = bgm.max_iter
    tol = bgm.tol
    verbose = True
    if method == 'bgm':
        wcpt = bgm.weight_concentration_prior_type
        reg_covar = bgm.reg_covar
        init_params = bgm.init_params
        tol = bgm.tol
        copy = BGM(n_components=n_clusters,
                   covariance_type=covartype,
                   n_init=n_init,
                   weight_concentration_prior_type=wcpt,
                   init_params=init_params,
                   max_iter=max_iter,
                   verbose=verbose,
                   reg_covar=reg_covar,
                   tol=tol)
        copy.weight_concentration_prior_ = bgm.weight_concentration_prior_
        copy.weight_concentration_ = bgm.weight_concentration_
        copy.mean_precision_prior = bgm.mean_precision_prior
        copy.mean_prior_ = bgm.mean_prior_
        copy.mean_precision_ = bgm.mean_precision_
        copy.covariance_prior_ = bgm.covariance_prior_
        copy.degrees_of_freedom_prior_ = bgm.degrees_of_freedom_prior_
        copy.degrees_of_freedom_ = bgm.degrees_of_freedom_
    if method == 'gmm':
        copy = GMM(n_components=n_clusters,
                   random_state=42,
                   covariance_type=covartype,
                   max_iter=max_iter,
                   n_init=n_init,
                   tol=tol,
                   verbose=verbose)
    copy.means_ = bgm.means_
    copy.covariances_ = bgm.covariances_
    copy.weights_ = bgm.weights_
    copy.precisions_ = bgm.precisions_
    copy.precisions_cholesky_ = bgm.precisions_cholesky_
    copy.converged_ = bgm.converged_
    copy.n_iter_ = bgm.n_iter_
    copy.lower_bound_ = bgm.lower_bound_
    return copy


def get_bgm(n_clusters=10, wcp=1.0e3, tol=None):
    """
    # wcp: higher number puts more mass in the center and will lead to more
    # more components being active
    """
    covartype = 'full'
    # init_params = 'random' #alternative is kmeans
    init_params = 'kmeans'  # alternative is kmeans
    # wcpt = 'dirichlet_process'       #wcp shouldbe (float, float)
    wcpt = 'dirichlet_distribution'  # wcp should be float
    warm_start = False
    verbose = False
    n_init = 1
    max_iter = 500
    reg_covar = 1e-5  # default is 1e-6
    if not tol:
        tol = 1e-3
    gmm = BGM(n_components=n_clusters,
              covariance_type=covartype,
              n_init=n_init,
              weight_concentration_prior_type=wcpt,
              init_params=init_params,
              max_iter=max_iter,
              verbose=verbose,
              reg_covar=reg_covar,
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
              n_init=n_init,
              warm_start=False,
              verbose=verbose)
    return gmm


def cluster_pars(xra, n_clusters=0):
    use_gmm = True
    use_kmeans = False
    parnum = len(xra)
    if n_clusters == 0:
        min_par_num = 50
        n_clusters = int(parnum / min_par_num)
    #xra = np.array(list(zip(lon,lat)))
    if n_clusters < 1:
        n_clusters = 1
    if n_clusters > 50:
        n_clusters = 50
    if use_kmeans:
        kpredict = KMeans(n_clusters=n_clusters, random_state=0).fit_predict(xra)
        plt.scatter(xra[:, 0], xra[:, 1], c=kpredict)
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
                  n_init=n_init,
                  warm_start=False,
                  verbose=verbose)
        #score, gmm = plot_gmm(gmm, xra)
        #gmm = GMM(n_components=n_clusters).fit(ra)
        #glabels = gmm.predict(ra)
        #plt.scatter(ra[:,0], ra[:,1], c=glabels)
    return gmm

# def subdivide_box(bnds, nnn):
#    jra = np.arange(0,nnn+1)
#    dlen = (bnds[2] - bnds[0]) / float(nnn)
#    latra = bnds[0] + jra * dlen
#
#    jra = np.arange(0,nnn+1)
#    dlen = (bnds[3] - bnds[1]) / float(nnn)
#    lonra = bnds[1] + jra * dlen


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
        self.delt = 5
        self.ymax = 9
        self.pdict = {}

    def read_pardump(self, drange=None, century=2000):
        self.df  = pardump.open_dataset(fname=self.tname, drange=drange,
                   century=century)
        #pd = pardump.Pardump(fname=self.tname)
        #self.df = pd.read(century=century)

    # def key2time(self):
    #    datelist=[]
    #    for key in self.pdict.keys:
    #        datelist.append(datetime.datetime.strptime(key, self.strfmt))
    #    #self.datetlist = datelist
    #    return datelist

    # def getbytime(self,time):
    #    # returns a dataframe for that time.
    #    dstr = time.strftime(self.strfmt)
    #    return self.pdict[dstr]

    def findsource(self, sorti):
        done = False
        iii = 0
        while not done:
            d1 = self.start + iii*self.delta
            print('find source', d1)
            df1 = self.getbytime(d1)
            print('Ages', df1['age'].unique())
            # keep only the new particles.
            df1 = df1[df1['age'] == self.delt]

            # if no particles released then
            if df1.empty and iii != 0:
                done = True
                print('empty', iii)
                continue
            df1 = df1[df1['sorti'].isin(sorti)]
            if iii == 0:
                dfsource = df1
            else:
                dfsource = pd.concat([dfsource, df1])
            iii += 1
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
        xbins = np.arange(0, 200, 5)
        ybins = np.arange(0, 4*9)/4.0
        cb = plt.hist2d(x2, y, bins=[xbins, ybins])
        plt.colorbar(cb[3])
        # sns.heatmap(x,y)
        return x2, y

    def time2int(self, timelist):
        newlist = []
        tmin = np.min(timelist)
        for ttt in timelist:
            val = (ttt-tmin)
            newlist.append(val.seconds)
        return newlist


def average_mfitlist(mfitlist, dd=None, dh=None, buf=None,
                     lat=None,
                     lon=None,
                     ht=None):
    """
    returns xarray DataArray
    """
    concra = combine_mfitlist(mfitlist, dd, dh, buf, lat, lon, ht)
    concra = concra.mean(dim='time')
    return concra


def combine_mfitlist(mfitlist, dd=None, dh=None, buf=None,
                     lat=None,
                     lon=None,
                     ht=None,
                     ):
    """
    mfitlist : list of MassFit objects.

    finds concentrations from each fit and combines into
    one xarray along dimension called 'time'. Although in
    some cases that dimension may represent something other than time.
    e.g. ensemble member number.

    returns xarray DataArray
    """
    iii = 0
    concra = xr.DataArray(None)
    for mfit in mfitlist:
        conc = mfit.get_conc(dd=dd,
                             dh=dh,
                             buf=buf,
                             lat=lat,
                             lon=lon,
                             ht=ht)
        if iii == 0:
            concra = conc
        else:
            concra = xr.concat([concra, conc], dim='time')
        iii += 1
    if 'time' not in concra.dims:
        concra = concra.expand_dims('time')
    return concra


def height_correction(zaprime, zter, msl=True):
    zmdl = 25000
    za = zaprime * (zmdl-zter)/float(zmdl)
    if msl:
        za += zter
    return za


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
    iii = 0
    height = 1e-12
    for val in dra.z.values:
        if height < val:
            break
        if math.isclose(height, val):
            break
        iii += 1
    under = dra.isel(z=np.arange(0, iii+1))
    under = process_under(under)
    above = dra.isel(z=np.arange(iii+1, len(dra.z.values), 1))
    new = xr.concat([under, above], dim='z')
    return new


def reflect_underground(dra):
    """
    reflects mass that is in first 3 levels underground
    onto first three levels above ground.
    """
    import math
    iii = 0
    height = 0
    for val in dra.z.values:
        if height < val:
            break
        if math.isclose(height, val):
            break
        iii += 1
    print('INDEX', iii)
    under1 = dra.isel(z=[iii-1, iii])
    print('under1', under1)
    under1 = process_under(under1)
    under = under1
    jjj = iii+1
    if iii-2 > 0:
        under2 = dra.isel(z=[iii-2, iii+1])
        print('under2', under2)
        under2 = process_under(under2)
        jjj = iii+2
        under = xr.concat([under, under2], dim='z')
    if iii-3 > 0:
        under3 = dra.isel(z=[iii-3, iii+2])
        print('under3', under3)
        under3 = process_under(under3)
        jjj = iii+3
        under = xr.concat([under, under3], dim='z')

    #under = under1
    print('under', under)
    above = dra.isel(z=np.arange(jjj, len(dra.z.values), 1))
    print('above', above.z.values)

    #lastz = under.z.values[-1]
    #under = under.sum(dim='z')

    #under = under.assign_coords(z=lastz)
    #under = under.expand_dims('z')

    new = xr.concat([under, above], dim='z')

    return new


class GridClass:

    def __init__(lat, lon, ht, dlat, dlon, dz, nlat, nlon, nz):

        self.llcrnr_lat = lat
        self.llcrnr_lon = lon
        self.llcrnr_ht = ht

        self.dlon = dlon
        self.dlat = dlat
        self.dz = dz

        self.nlon = nlon
        self.nlat = nlat
        self.nz = nz

    def makegrid(self):
        """
        """
        # checked HYSPLIT code. the grid points
        # do represent center of the sampling area.
        slat = self.llcrnr_lat * self.dlat
        slon = self.llcrnr_lon * self.dlon
        sht = self.llcrnr_ht * self.dz
        lat = np.arange(slat, slat + self.nlat * self.dlat, self.dlat)
        lon = np.arange(slon, slon + self.nlon * self.dlon, self.dlon)
        ht = np.arange(sht, sht + self.nz * dz, dz)
        #lonlist = [lon[x - 1] for x in xindx]
        #latlist = [lat[x - 1] for x in yindx]
        mgrid = np.meshgrid(lon, lat, ht)
        return mgrid


def getvolume(mean, cov, x, y, z, dh, dz, verbose=False):
    """
    mean and cov of gaussian
    rv is multivariate_normal object
    x, y, z : float. center position.
    dh : float. half width/ length
    dz : float. half height.
    Returns 
    volume under portion of 3d Gaussian.

    """
    rv = multivariate_normal(mean, cov)

    a1 = rv.cdf([x+dh, y+dh, z+dz])
    a5 = rv.cdf([x-dh, y+dh, z-dz])
    a8 = rv.cdf([x-dh, y+dh, z+dz])
    a4 = rv.cdf([x+dh, y+dh, z-dz])

    a6 = rv.cdf([x+dh, y-dh, z+dz])
    a2 = rv.cdf([x-dh, y-dh, z-dz])
    a3 = rv.cdf([x+dh, y-dh, z-dz])
    a7 = rv.cdf([x-dh, y-dh, z+dz])

    v1 = (a1+a5) - (a8+a4)
    v2 = (a6+a2) - (a3+a7)

    if verbose:
        print(a1, a5)
    if verbose:
        print(a8, a4)
    if verbose:
        print(a6, a2)
    if verbose:
        print(a7, a3)
    volume = v1 - v2
    if verbose:
        print('volume, V1, V2', volume, '=', v1, '-',  v2)
    return volume

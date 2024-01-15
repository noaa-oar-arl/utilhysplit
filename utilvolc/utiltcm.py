import logging
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
from utilvolc import plottcm

logger = logging.getLogger(__name__)

"""
Classes
    ParametersIn  - for Parameters_in.dat file
    InvEstimatedEmissions - for estimated emissions out.dat output from inverse
    InvOut2dat  - for the out2.dat output from inverse.
Functions
    readoutdat
    make_emissions
"""

# 2023 Dec 04 (amc) changed InverseOutDat to InvEstimatedEmissions
# 2023 Dec 04 (amc) changed self.outdat to self.emissions
# 2023 Dec 04 (amc) added make_emissions function
# 2024 Jan 13 (amc) removed InverseDat class
# here is a difference

class ParametersIn:
    """
    simple way of changing the n_ctrl, nx_ctrl and
    lbfgs_nbd values in Parmeters_in.dat file.
    """

    def __init__(self, fname):
        # read existing file.
        with open(fname, "r") as fid:
            self.lines = fid.readlines()

    def change_and_write(self, fname, nx_ctrl):
        self.change(nx_ctrl)
        self.write(fname)

    def change(self, nx_ctrl):
        # change relevant lines.
        nx_ctrl = int(nx_ctrl)
        self.lines[1] = " N_ctrl={}\n".format(nx_ctrl)
        self.lines[2] = " Nx_ctrl={}\n".format(nx_ctrl)
        self.lines[40] = " lbfgs_nbd={}*0/\n".format(nx_ctrl)

    def write(self, fname):
        with open(fname, "w") as fid:
            [fid.write(line) for line in self.lines]


def readoutdat(wdir,fname):
    """
    name : str : filename of out.dat file which has estimated release rates in same
                 order as tcm columns. Currently first column is just a dummy.

    Returns :
    df : pandas dataframe/
    """
    if os.path.isfile(os.path.join(wdir, fname)):
        df = pd.read_csv(os.path.join(wdir, fname), sep="\s+", header=None)
    else:
        print("File not found {}".format(os.path.join(wdir, fname)))
        df = pd.DataFrame()
    return df


class InvEstimatedEmissions:

    def __init__(self, fname="out.dat",sourcehash=None,columns=None):
        self.fname = fname
        # properties
        self._df = pd.DataFrame()    # dataframe which holds data from out.dat
        self.sourcehash = sourcehash # dictionary with informtion about what the columns mean.
        self.columns = columns       # column names

    @property
    def df(self):
        return self._df.copy()

    @df.setter
    def df(self,df):
        self._df = df

    @property
    def sourcehash(self):
        return self._sourcehash

    @sourcehash.setter
    def sourcehash(self,sourcehash):
        if isinstance(sourcehash,dict):
           self._sourcehash=sourcehash
        else:
           self._sourcehash = {}

    @property
    def columns(self):
        return self._columns

    @columns.setter
    def columns(self,columns):
        if isinstance(columns,(list,np.ndarray)):
           self._columns = columns
        else:
           self.columns = []

    def read(self, wdir):
        self.df = readoutdat(wdir, self.fname)

    def make_emissions(self):
        return make_emissions(self.sourcehash,self.columns,self.df)

    def emis_hist(self, units="g/h",log=True):
        """
        histogram of emissions
        """
        if log:
            if units == "g/h":
                vals = np.log10(self.df[1].values)
                nbin =1
                vals = np.where(vals>=-2,vals,-5)
        else:
            vals = self.df[1].values
            nbin=10
        nmin = int(np.floor(np.nanmin(vals)))
        nmax = int(np.ceil(np.nanmax(vals)))
        nbins = len(np.arange(nmin, nmax, nbin))
        if nbins > 100: nbins = 100
        plt.hist(vals, bins=nbins)
        ax = plt.gca()
        ax.set_xlabel("Log emission (g/h)")
        print(nmax)
        return ax

    def plot(self,log=True,thresh=0):
        edf = self.make_emissions()
        plottcm.plot_emissions(edf,log=log,thresh=thresh)

    def plot_timeseries(self,log=True,marker='o'):
        """
        plots a time series of the emissions
        """
        edf = self.make_emissions()
        plottcm.plot_emissions_timeseries(edf,log,marker)

    def plot_profile(self,log=True,marker='o'):
        """
        plots a time series of the emissions
        """
        edf = self.make_emissions()
        plottcm.plot_emissions_profile(edf,marker=marker)

    def write_emit(self,vlat,vlon,threshold=50,area=1,name='EMIT.txt',date_cutoff=None):
        from utilvolc.tcm_emit import construct_efile
        phash = {1:'PASH'}
        edf = self.make_emissions()
        time = edf['date'].values
        ht = edf['ht'].values
        mass = edf['mass'].values
        vals = list(zip(time,ht,mass))
        efile = construct_efile(vals,vlat,vlon,area=area,emis_threshold=50,name=name,
                                date_cutoff=date_cutoff,phash=phash)
        print('writing efile {}'.format(name))
        efile.write_new(name) 


class InverseOut2Dat:

    def __init__(self, fname="out2.dat"):
        self.fname = fname
        self.df = pd.DataFrame()  # dataframe which holds data from out.dat

    def read(self, wdir):
        df = readoutdat(wdir, self.fname)
        df.columns = ['index','observed','model']
        self.df = df 

    def plot_conc(self, cmap="viridis"):
        if np.all(np.isnan(self.df.model.values)):
           logger.warning('plotting failed. All values in the model column are nan')
           return 

        sns.set()
        sns.set_style("whitegrid")
        df = self.df
        # plt.plot(df['observed'],df['model'],'k.',MarkerSize=3)
        cb = plt.hist2d(
            df["observed"],
            df["model"],
            cmap=cmap,
            norm=mpl.colors.LogNorm(),
            bins=[20, 20],
        )
        cbar = plt.colorbar(cb[3])
        cbar.ax.set_ylabel("Number of Points")
        ax = plt.gca()
        ax.set_xlabel("observed")
        ax.set_ylabel("modeled")
        nval = np.max(df["observed"])
        # plot 1:1 line
        plt.plot([0, nval], [0, nval], "--b", linewidth=1)
        return ax


def make_emissions(sourcehash, tcm_columns, dfdat):
    """
    makeoutdat for InverseAshPart class.
    dfdat : pandas dataframe output by InvEstimatedEmissions class get_emis method.
    Returns
    vals : tuple (date, height, emission mass, particle size)
    """
    # matches emissions from the out.dat file with
    # the date and time of emission.
    # uses the tcm_columns array which has the key
    # and the sourehash dictionary which contains the information.
    datelist = []
    htlist = []
    valra = []
    psizera = []
    for val in zip(tcm_columns, dfdat[1]):
        shash = sourcehash[val[0]]
        datelist.append(shash["sdate"])
        htlist.append(shash["bottom"])
        valra.append(val[1])
        psizera.append(val[0][1])
    vals = list(zip(datelist, htlist, valra, psizera))
    emission_df = pd.DataFrame.from_records(vals,columns=['date','ht','mass','psize'])
    return emission_df


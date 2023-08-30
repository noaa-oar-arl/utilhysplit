import logging
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns

logger = logging.getLogger(__name__)

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


class InverseOutDat:

    def __init__(self, fname="out.dat"):
        self.fname = fname
        self.df = pd.DataFrame()  # dataframe which holds data from out.dat

    def read(self, wdir):
        self.df = readoutdat(wdir, self.fname)

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

class InverseDat:
    def __init__(self, wdir, fname="out.dat", fname2="out2.dat"):
        """
        # out.dat has estimated release rates in same order as tcm columns.
        # out2.dat has observed(2nd column) and modeled(3rd column) mass loadings.
        """
        self.wdir = wdir
        self.outdat = InverseOutDat(fname)
        self.out2dat = InverseOut2Dat(fname2)
        self.read()

    def read(self):
        """
        name : str : filename of out.dat file which has estimated release rates in same
                     order as tcm columns. Currently first column is just a dummy.

        Returns :
        df : pandas dataframe/
        """
        self.outdat.read(self.wdir)
        self.out2dat.read(self.wdir)

    def get_emis(self, name=None):
        """
        """
        return self.outdat.df 

    def emis_hist(self, units="g/h", log=True):
        ax = self.outdat.emis_hist(units=units, log=log)
        return ax

    def plot_conc(self, cmap="viridis"):
        ax =  self.out2dat.plot_conc(cmap=cmap)
        return ax


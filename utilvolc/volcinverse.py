import datetime
import logging
import os
from subprocess import call

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
import seaborn as sns
import xarray as xr

import utilhysplit.hysplit_gridutil as hgu
from monetio.models import hysplit
from utilhysplit import hcontrol
from utilhysplit.evaluation import ensemble_tools, plume_stat
from utilhysplit.metfiles import MetFileFinder
from utilhysplit.plotutils import colormaker
from utilvolc import volcat
from utilvolc.runhelper import Helper

from utilvolc.inversioninterface import PairedDataInterface
from utilvolc.inversioninterface import RunInversionInterface
from utilvolc.inversioninterface import TCMInterface
from utilvolc.volcpaired import VolcatHysplit
from utilvolc import invhelper
from utilvolc.volctcm import TCM

# def save_emis_gefs(inverse):
#    tlist = 'gec00'
#    inverse.save_emis(eii=0, savename = 'gec00')
#    for nlist in np.arange(1,30):
#        ttt = 'gep{:2d}'.format(nlist)
#        inverse.save_emis(eii=nlist

logger = logging.getLogger(__name__)


def check_type_basic(check):
    if isinstance(check, int):
        return True
    if isinstance(check, float):
        return True
    if isinstance(check, str):
        return True
    return False


# def check_type_iterable(check):


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


def read_emis_df(savename):
    """
    reads dataframe saved by the make_outdat_df method.
    index is height
    columns are date
    values are emissions.
    """
    try:
        emisdf = pd.read_csv(savename, index_col=0, header=1)
    except:
        print("trying this", savename)
        # emisdf = pd.read_csv(savename,index_col='date',header='ht',dtype={'ht':'float'})
        emisdf = pd.read_csv(savename, index_col="date", header="ht")
    return emisdf.dropna()


class InverseOutDat:
    def __init__(self, wdir, fname="out.dat", fname2="out2.dat"):
        # out.dat has estimated release rates in same order as tcm columns.
        # out2.dat has observed(2nd column) and modeled(3rd column) mass loadings.
        self.wdir = wdir
        self.subdir = wdir  # initially set the sub directory to the working directory.

        self.df = pd.DataFrame()  # dataframe which holds data from out.dat
        self.df2 = pd.DataFrame()  # dataframe which holds data from out2.dat

        self.fname = fname  # name of out.dat file
        self.fname2 = fname2  # name of out2.dat file

    def read_out(self, name):
        """
        name : str : filename of out.dat file which has estimated release rates in same
                     order as tcm columns. Currently first column is just a dummy.

        Returns :
        df : pandas dataframe/
        """

        # out.dat has estimated release rates in same order as tcm columns.
        wdir = self.wdir
        if os.path.isfile(os.path.join(wdir, name)):
            df = pd.read_csv(os.path.join(wdir, name), sep="\s+", header=None)
        else:
            print("File not found {}".format(os.path.join(wdir, name)))
            df = pd.DataFrame()
        return df

    def get_emis(self, name=None):
        """
        Reads in out.dat file
        """
        if not name:
            name = self.fname
        else:
            self.fname = name
        df = self.read_out(name)
        self.df = df
        return df

    def emis_hist(self, units="g/h"):
        if units == "g/h":
            vals = np.log10(self.df[1].values)
        nmin = int(np.floor(np.min(vals)))
        nmax = int(np.ceil(np.max(vals)))
        nbins = len(np.arange(nmin, nmax, 1))
        plt.hist(vals, bins=nbins)
        ax = plt.gca()
        ax.set_xlabel("Log emission (g/h)")

    def get_conc(self, name=None):
        """
        Reads in out2.dat file
        """
        if not name:
            name = self.fname2
        df = self.read_out(name)
        df.columns = ["index", "observed", "model"]
        self.df2 = df
        return df

    def plot_conc(self, cmap="viridis"):
        sns.set()
        sns.set_style("whitegrid")
        df = self.df2
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

    def get_vals(self):
        return -1



class InvDirectory:
    def __init__(self):
        # directories
        # use set_directory method to change
        default = "./"
        self._modeldir = default  # location of netcdf files with HYSPLIT output (cdump)
        self._obsdir = default  # location of observation data
        self._wdir = default  # location to write outputs, create new directories
        self._subdir = self.wdir  # name of subdirectory in working directory
        self._execdir = default  # directory where inversion code executable is located
        self._hysplitdir = default  # directory where HYSPLIT executables are lcated.
        self._subdirlist = [self._subdir] # list of subdirectories associated 

    @staticmethod
    def check(testdir):
        if os.path.isdir(testdir):
            return True
        else:
            return False

    @property
    def obsdir(self):
        return self._obsdir

    @obsdir.setter
    def obsdir(self, idir):
        if self.check(idir):
            self._obsdir = idir
        else:
            logger.warning("Directory does not exist {}".format(idir))

    @property
    def modeldir(self):
        return self._modeldir

    @modeldir.setter
    def modeldir(self, idir):
        if self.check(idir):
            self._modeldir = idir
        else:
            logger.warning("Directory does not exist {}".format(idir))

    @property
    def wdir(self):
        return self._wdir

    @wdir.setter
    def wdir(self, idir):
        if self.check(idir):
            self._wdir = idir
        else:
            logger.warning("Directory does not exist {}".format(idir))

    @property
    def subdir(self):
        return self._subdir

    # creates the subdirectory if it does not exist:w
    @subdir.setter
    def subdir(self, idir):
        subdir = os.path.join(self._wdir, idir)
        if self.check(subdir):
            self._subdir = subdir
        else:
            callstr = "mkdir -p" + os.path.join(self._wdir,subdir) 
            call(callstr, shell=True)
        self.add_subdir(subdir)
 
    @property
    def hysplitdir(self):
        return self._hysplitdir

    @hysplitdir.setter
    def hysplitdir(self, idir):
        if self.check(idir):
            self._hysplitdir = idir
        else:
            logger.warning("Directory does not exist {}".format(idir))

    @property
    def execdir(self):
        return self._execdir

    @execdir.setter
    def execdir(self, idir):
        if self.check(idir):
            self._execdir = idir
        else:
            logger.warning("Directory does not exist {}".format(idir))

    def __str__(self):
        """
        InverseAshEns class
        """
        rstr = "Working directory, wdir :{}".format(self.wdir)
        rstr += "\n"
        rstr += "execdir :{}".format(self.execdir)
        rstr += "\n"
        rstr += "hysplitdir :{}".format(self.hysplitdir)
        rstr += "\n"
        rstr += "subdir :{}".format(self.subdir)
        rstr += "\n"
        rstr += "modeldir :{}".format(self.modeldir)
        rstr += "\n"
        rstr += "obsdir :{}".format(self.obsdir)
        rstr += "\n"
        return rstr

    def add_subdir(self,subdir):
        if subdir not in self._subdirlist:
            self._subdirlist.append(subdir) 

    def set_directory(self, wdir, execdir, modeldir, obsdir, hysplitdir, subdir=None):
        """
        InverseAshEns class
        """
        self.modeldir = modeldir
        self.obsdir = obsdir
        self.wdir = wdir
        if not subdir:
            self.subdir = wdir
        else:
            self.subdir = subdir
        self.execdir = execdir
        self.hysplitdir = hysplitdir


class RunInversion(RunInversionInterface):
    def __init__(self, inp):

        # keys which need to be in inp
        # WORK_DIR is from the configfile and is where the cdump file is.
        # so this is where the model_dir is as well
        self._ilist = ["WORK_DIR", "INV_DIR", "HYSPLIT_DIR", "OBS_DIR"]

        self._inp = {}
        self.inp = inp

        self._directories = InvDirectory()
        self.set_directory(self.inp)

        self._paired_data = VolcatHysplit(inp, self._directories.modeldir, inp["fname"],
                                          self._directories.obsdir)
        self._tcm = TCM()

    def set_directory(self, inp):
        self._directories.set_directory(
            inp["WORK_DIR"], inp["INV_DIR"], inp["WORK_DIR"], inp["OBS_DIR"], inp["HYSPLIT_DIR"],
            None
        )

    @property
    def inp(self):
        return self._inp

    # updates does not overwrite.
    @inp.setter
    def inp(self, inp):
        if "configfile" in inp.keys():
            if "configdir" in inp.keys():
                from utilvolc.invhelper import get_inp_hash

                inp2 = get_inp_hash(inp["configdir"], inp["configfile"])
                inp.update(inp2)
        self._inp.update(inp)
        for icc in self._ilist:
            if icc not in self._inp.keys():
                logger.warning("Missing required input {}".format(icc))
                print("RunInversion: Missing required input {}".format(icc))

    @property
    def paired_data(self):
        return self._paired_data

    @property
    def tcm(self):
        return self._tcm

    @property
    def directories(self):
        return self._directories

    def copy(self):
        return -1

    def print_summary(self):
        return -1

    def get_times(self):
        kyy = ['start_date', 'durationOfSimulation']
        shash =  {key:self.inp[key] for key in kyy}
        dlist = invhelper.create_dlist(shash)
        return dlist

    def setup_paired_data(self):
        dlist = self.get_times()
        print('Times to prepare {}'.format(len(dlist)))
        for time in dlist:
            d1 = time[0].strftime("%Y %m/%d %HZ")
            d2 = time[1].strftime("%Y %m/%d %HZ")
            logger.info('Preparing paired data for {} {}'.format(d1,d2))
            print('Preparing paired data for {} {}'.format(d1,d2))
            self._paired_data.prepare_one_time(time)
                

    def setup_tcm(self,tiilist, 
                  concmult=1,
                  remove_cols=True,
                  remove_rows=False,
                  remove_sources=False,
                  remove_ncs = 0,
                  tcm_name='tcm'):
        pdata = []
        columns = self._paired_data.cdump.ens.values
        self.tcm.columns = columns

        basetag = self.inp['jobname']
        self.directories.add_subdir(basetag)
        tcmname = '{}.{}.txt'.format(self.inp['jobname'],tcmname)
        tcmname = os.path.join(self.directories.subdir,tcmname) 
        tag = invhelper.create_runtag(basetag,tiilist,remove_cols,remove_rows,remove_sources,remove_ncs)
        if not isinstance(tcm_name,str):
           print('TCM output to {}'.format(tag))
           tcm_name = tag

        # get the cdump and volcat data
        for tii in tiilist:
            model = self._paired_data.cdump_hash[tii]
            obs = self._paired_data.volcat_avg_hash[tii]
            pdata.append((model,obs))
        self._tcm.make_tcm_mult(pdata,concmult,remove_cols,remove_rows,remove_sources,remove_ncs)
        self._tcm.write(tcm_name)
        return -1

    def run_tcm(self):
        #InverseAshEns 
        """
        """
        out_name1 = "out.dat"
        out_name2 = "out2.dat"
        inp_name = "TCM_sum.csv"
        cmd = os.path.join(self.execdir, "new_lbfgsb.x")
        new_name1, new_name2 = self.make_tcm_names()
        for iii, tcm in enumerate(self.tcm_names):
            print("run_tcm tag", self.taglist[iii])
            os.chdir(self.subdir)
            print("working in ", os.getcwd())
            params = ParametersIn(os.path.join(self.execdir, "Parameters_in.dat.original"))
            params.change_and_write(
                os.path.join(self.subdir, "Parameters_in.dat"), self.n_ctrl_list[iii]
            )

            Helper.remove(inp_name)

            # remove any existing output files.
            Helper.remove(out_name1)
            Helper.remove(out_name2)

            # run

            Helper.copy(tcm, inp_name)
            print(cmd)
            # Helper.execute_with_shell(cmd)
            Helper.execute(cmd)
            # move output files to new names
            Helper.move(out_name1, new_name1[iii])
            Helper.move(out_name2, new_name2[iii])
            Helper.move("fort.188", "fort.188.{}".format(self.taglist[iii]))


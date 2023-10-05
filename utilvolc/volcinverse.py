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
from utilvolc.runhelper import make_dir

from utilvolc.inversioninterface import PairedDataInterface
from utilvolc.inversioninterface import RunInversionInterface
from utilvolc.inversioninterface import TCMInterface
from utilvolc.volcpaired import VolcatHysplit
from utilvolc import invhelper
from utilvolc.volctcm import TCM
from utilvolc import plottcm

# def save_emis_gefs(inverse):
#    tlist = 'gec00'
#    inverse.save_emis(eii=0, savename = 'gec00')
#    for nlist in np.arange(1,30):
#        ttt = 'gep{:2d}'.format(nlist)
#        inverse.save_emis(eii=nlist

logger = logging.getLogger(__name__)

"""
 Classes
     InvDirectory
     RunInversion
 Functions
     check_type_basic
     read_emis_df
"""

# 2023 Aug 21 (amc) ParametersIn moved to utiltcm.py


def check_type_basic(check):
    if isinstance(check, int):
        return True
    if isinstance(check, float):
        return True
    if isinstance(check, str):
        return True
    return False


# def check_type_iterable(check):


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
        if make_dir(self._wdir, idir):
            self._subdir = os.path.join(self._wdir,idir)
        self.add_subdir(idir)
 
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

        self._inp = {}  #dictionary with information from the config file.
        self.inp = inp

        self._directories = InvDirectory()
        self.set_directory(self.inp)

        self._paired_data = VolcatHysplit(inp, self._directories.modeldir, inp["fname"],
                                          self._directories.obsdir)
        self._tcm = TCM()

        # can be constructed in the inp setter from inp['configdir'] and inp['configfile'].
        self._sourcehash = {}

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
                print('making sourcehash')
                self._sourcehash = invhelper.get_sourcehash(inp['configdir'],inp['configfile'])

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

    #def get_times(self):
    #    kyy = ['start_date', 'durationOfSimulation']
    #    shash =  {key:self.inp[key] for key in kyy}
    #    dlist = invhelper.create_dlist(shash)
    #    return dlist

    def setup_paired_data(self, ilist=None):
        #dlist = self.get_times()

        dlist = self._paired_data.get_times()

        if isinstance(ilist,list):
           dlist = dlist[ilist[0]:ilist[1]]
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
        tag = invhelper.create_runtag(basetag,tiilist,remove_cols,remove_rows,remove_sources,remove_ncs)
        self.directories.subdir = tag
        print('subdir', self.directories.subdir) 


        tcm_name = '{}.{}.txt'.format(self.inp['jobname'],tcm_name)
        tcm_name = os.path.join(self.directories.subdir,tcm_name) 

        if not isinstance(tcm_name,str):
           print('TCM output to {}'.format(tag))
           tcm_name = tag

        # get the cdump and volcat data
        for tii in tiilist:
            model = self._paired_data.cdump_hash[tii]
            obs = self._paired_data.volcat_avg_hash[tii]
            pdata.append((model,obs))

        print('tag is ', tag)
        self._tcm.make_tcm_mult(pdata,concmult,remove_cols,remove_rows,remove_sources,remove_ncs)
        self._tcm.write(tcm_name)
        return -1

    # TODO also a funciton in volctcm.
    def make_outdat(self):
        """
        makeoutdat for InverseAshPart class.
        dfdat : pandas dataframe output by InverseOutDat class get_emis method.
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
        dfdat = self.tcm.output.get_emis()
        # if not isinstance(self.tcm_columns, (list, np.ndarray)):
        #   self.tcm_columns = self.make_tcm_columns()
        for val in zip(self.tcm.columns, dfdat[1]):
            shash = self._sourcehash[val[0]]
            datelist.append(shash["sdate"])
            htlist.append(shash["bottom"])
            valra.append(val[1])
            psizera.append(val[0][1])
        vals = list(zip(datelist, htlist, valra, psizera))
        return vals

    def plot_outdat_ts(self, log=False, fignum=1, unit="kg/s", ax=None):
        # InverseAshPart class

        # plots time series of MER. summed along column.
        # dfdat : pandas dataframe output by InverseOutDat class get_emis method.
        if not ax:
            sns.set()
            sns.set_style("whitegrid")
            fig = plt.figure(fignum, figsize=(10, 5))
            ax = fig.add_subplot(1, 1, 1)
        df = self.make_outdat_df()
        ax, ts = plottcm.plot_outdat_ts_psize_function(df,ax=ax,log=log)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels)
        return ax, ts

    # to do -also a function in volctcm.py
    def make_outdat_df(self, savename=None, part="basic"):
        # InverseAsh class
        """
        dfdat : pandas dataframe output by InverseOutDat class get_emis method.
                this is a list of tuples (source tag), value from emissions
        """
        vals = self.make_outdat()
        vals = list(zip(*vals))
        ht = vals[1]
        time = vals[0]
        # emit = np.array(vals[2])/1.0e3/3600.0
        emit = np.array(vals[2])

        # this is for particle size. Used by the InverseAshPart class.
        if len(vals) == 4:
            psize = vals[3]
            data = zip(time, ht, psize, emit)
            if part == "index":
                iii = 0
                cols = [1, 2]
            elif part == "cols":
                iii = 1
                cols = [0, 2]
            colnames = ["date", "ht", "psize", "mass"]
        # this is for only height and time
        else:
            data = zip(time, ht, emit)
            iii = 1
            cols = 0
            colnames = ["date", "ht", "mass"]
        dfout = pd.DataFrame(data)
        if part == "condense" and len(vals) == 4:
            dfout.columns = colnames
            dfout = dfout.groupby(["date", "ht"]).sum()
            dfout = dfout.reset_index()
            dfout.columns = [0, 1, 2]
            iii = 1
            cols = 0

        if part == "basic":
            dfout.columns = colnames
            return dfout
        dfout = dfout.pivot(columns=cols, index=iii)
        if savename:
            print("saving  emissions to ", savename)
            dfout.to_csv(savename)
        return dfout

    def plot_outdat_profile(self, fignum=1, unit="kg", ax=None):
        # InverseAshPart class
        dfdat = self.tcm.output.get_emis()
        if not ax:
            sns.set()
            sns.set_style("whitegrid")
            fig = plt.figure(fignum, figsize=(10, 5))
            ax = fig.add_subplot(1, 1, 1)
        df = self.make_outdat_df(dfdat, part="basic")
        ax, ts = plottcm.plot_outdat_profile_psize_function(df,ax=ax)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels)
        return ax, ts

    def run_hysplit(self):
        from ashapp.maindispersion import MainEmitTimes
        jobid = self.directories.subdir.split('/')[-1]
        inp = self.inp.copy()
        inp['durationOfSimulation']=24
        inp['WORK_DIR'] = self.directories.subdir
        arun = MainEmitTimes(inp,jobid)
        arun.doit()


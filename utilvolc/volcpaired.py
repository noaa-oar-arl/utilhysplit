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
from ashapp import level_setter

#from utilvolc.inversioninterface import PairedDataInterface

logger = logging.getLogger(__name__)

"""
Classes
    VolcatHysplit
"""

def add_qva_levels_to_dset(self,dset):
    lsetter = level_setter.LevelSetter()
    


#class VolcatHysplit(PairedDataInterface):
class VolcatHysplit():

    def __init__(
        self,
        #inp,
        #tdir,
        #fname,
        #vdir,
        #vid,
        #configdir="./",
        #configfile=None,
        verbose=False,
        ensdim="ens",
    ):

        """
        configfile : full path to configuration file.

        fname : str. full file path to netcdf file with model data.
                xarray DataArray or DataSet with cdump information
        """
        # InverseAsh class

        #self.tdir = tdir    # directory for hysplit output
        #self.fname = fname  # name of hysplit output

        self.vdir = None  # directory for volcat data
        #self.vid = vid  # volcano id.

        #self.n_ctrl = 0  # determined in the write_tcm method.
        # needed for Parameters_in.dat input into inversion algorithm

        # keep volcat arrays for different averaging times.
        self.volcat_hash = {}

        # prepare_one_time method adds data to these dictionaries.
        # the data in cdump_hash and volcat_avg_hash have been aligned.
        # and are on the same grid.
        self.volcat_avg_hash = {}
        self.volcat_ht_hash = {}
        self.cdump_mass_hash = {}

        # multiplication factor if more than 1 unit mass released.
        self.concmult = 1

        # get_cdump method loads the cdump file.
        self.cdump = None
        #self.get_cdump(tdir, fname, verbose, ensdim)

        # add_config_info creates these two dictionaries from the configuration.
        #self.sourcehash = {}
        #self.inp = {}
        #self.add_config_info(configdir, configfile)
        self._inp = {}

        # These values are for if spatial coarsening is used.
        self.coarsen = None
        self.original_volcat_avg_hash = {}
        self.original_volcat_ht_hash = {}
        self.original_cdump_mass_hash = {}

        self.set_bias_correction(slope=None, intercept=None)

        # particle size information

        # tcm_columns will be filled in later.
        #self.tcm_columns = None

    #property
    def modeloutput(self):
        return self._modeloutput

    #property
    def inp(self):
        return self._inp

    #property
    def obs(self):
        return self._obs

    def close_arrays(self):
        """ """
        # InverseAsh class
        self.cdump.close()
        for key in self.volcat_avg_hash.keys():
            self.volcat_avg_hash[key].close()
        # self.massload.close()

    def copy(self):
        """ """
        # InverseAsh class
        iacopy = InverseAsh(self.tdir, self.fname, self.vdir, self.vid, self.n_ctrl)
        iacopy.volcat_avg_hash = self.volcat_avg_hash
        iacopy.volcat_cdump_hash = self.cdump_mass_hash
        icacopy.concmult = self.concmult
        return icacopy

    def print_summary(self):
        """ """
        # InverseAsh class
        print("Observations availalbe in volcat_avg_hash")
        print(self.volcat_avg_hash.keys())
        print("times in cdump file")
        self.print_times()

    def add_cdump_dset(self,dset,ensdim='ens'):
        self.cdump = self.process_cdump(dset,ensdim)

    def add_qva_levels(self):
        qvalevs = level_setter.get_qva_levels()
        self.add_levels(qvalevs)

    def add_levels(self, levels):
        lkey = 'Level top heights (m)'
        z = self.cdump.z.values
        for zval in z:
            if zval not in levels:
               levels.append(zval)
        levels.sort()

        if lkey in self.cdump.attrs.keys():
           logger.info('current levels', self.cdump.attrs[lkey])
           logger.info('changing to ', levels)
        self.cdump = self.cdump.assign_attrs({lkey:list(levels)})        


    def add_model(self,tdir,fname,verbose=False,ensdim='ens'):
        """
        adds another cdump file along the time dimension.
        """
        cdump = self.get_cdump_sub(tdir,fname,verbose,ensdim)
        return xr.concat([self.cdump,cdump],dim='time')

    def get_cdump(self,tdir,fname,verbose=False,ensdim='ens'):
        self.cdump = self.get_cdump_sub(tdir,fname,verbose,ensdim)

    def get_cdump_sub(self, tdir, fname, verbose=False, ensdim="ens"):
        """ """
        # InverseAsh class
        # hysplit output. xarray.
        print("working on {}".format(fname))
        if isinstance(fname, str):
            cdump = xr.open_dataset(os.path.join(tdir, fname), engine="netcdf4")
            logger.info("opening {} {}".format(tdir, fname))
        else:
            cdump = fname
        if not hysplit.check_grid_continuity(cdump):
            logger.info("Grid not continuous")
        self.process_cdump(cdump,ensdim)

    def process_cdump(self,cdump,ensdim):
        if not isinstance(cdump, xr.core.dataarray.DataArray):
            # turn dataset into dataarray
            temp = list(cdump.keys())
            cdump = cdump[temp[0]]
        # get rid of source dimension
        if ensdim == "ens":
            cdump = cdump.isel(source=0)
        elif ensdim == "source":
            cdump = cdump.isel(ens=0)
        cdump, dim = ensemble_tools.preprocess(cdump)
        if dim == "source":
            cdump = cdump.rename({"ens": "metid"})
            cdump = cdump.rename({"source": "ens"})
        return cdump.fillna(0)
        #self.cdump = cdump.fillna(0)
        # print('HERE')
        # self.cdump = cdump

    def add_config_info(self, configdir, configfile):
        """ """
        # the ens dimension holds is key to what emission source was used.
        # the sourcehash is a dictionary
        # key is the ensemble number
        # values is another dictionary with
        # sdate: begin emission
        # edate: end emission
        # bottom : lower height of emission
        # top : upper height of emission.
        if configfile:
            if not os.path.isfile(os.path.join(configdir, configfile)):
                configfile = None
        if configfile:
            self.sourcehash = get_sourcehash(configdir, configfile)
            self.inp = get_inp_hash(configdir, configfile)
        else:
            self.sourcehash = {}
            self.inp = {}

    def add_inp(self, configdir, configfile):
        """ """
        # InverseAsh class
        self.inp = get_inp_hash(configdir, configfile)

    def set_concmult(self, concmult):
        """ """
        # InverseAsh class
        self.concmult = concmult

    def get_volcat(self, daterange, verbose=False):
        """ """
        # InverseAsh class
        vdir = self.vdir
        #vid = self.vid
        vid = None
        tii = self.time_index(daterange[0])
        # if tii not in self.volcat_hash.keys():
        done = True
        if done:
            das = volcat.get_volcat_list(
                vdir,
                daterange=daterange,
                vid=vid,
                decode_times=True,
                verbose=verbose,
                include_last=False,
            )
            self.volcat_hash[tii] = das
        else:
            das = self.volcat_hash[tii]
        # create one dataset with dimension of time.
        if len(das) > 0:
            # vset = xr.concat(das, dim="time")
            vset = volcat.combine_regridded(das)
        else:
            print("No volcat files found ")
            return xr.DataArray()
        return vset

    def clip(self, dummy, buf=0):
        """ """
        # InverseAsh class
        # clip ra according to where dummy has 0's.
        aaa = np.where(dummy > 0)
        a1 = np.min(aaa[0])
        a2 = np.max(aaa[0])
        b1 = np.min(aaa[1])
        b2 = np.max(aaa[1])

        if a2 + buf < dummy.y.values[-1]:
            a2 = a2 + buf

        if b2 + buf < dummy.x.values[-1]:
            b2 = b2 + buf

        a1 = a1 - buf
        b1 = b1 - buf
        if a1 < 0:
            a1 = 0
        if b1 < 0:
            b1 = 0

        return a1, a2, b1, b2

    def print_times(self):
        """ """
        # InverseAsh class
        timelist = [pd.to_datetime(x) for x in self.cdump.time.values]
        for time in timelist:
            print(time.strftime("%Y %m %d %H:%Mz"))

    def get_times(self):
        """
        tii : integer
        """
        # TODO check whether time is beginning or end of averaging period.
        timelist = [pd.to_datetime(x) for x in self.cdump.time.values]
        dt = timelist[1] - timelist[0]
        tlist = [[x,x+dt] for x in timelist] 
        return tlist

    def get_time(self, tii):
        """
        tii : integer
        """
        # InverseAsh class
        timelist = [pd.to_datetime(x) for x in self.cdump.time.values]
        return timelist[tii]

    def time_index(self, time):
        """ """
        # InverseAsh class
        timelist = [pd.to_datetime(x) for x in self.cdump.time.values]
        try:
            iii = timelist.index(time)
        except:
            # print('timelist', timelist)
            iii = None
        return iii

    def compare_plotsA(self, daterange, zii=None, tii=None, levels=None):
        # InverseAsh class
        """
        must input either daterange or tii.

        if zii is None then sum along ensemble dimension showing coverage of all HYSPLIT runs.
        For the inversion runs, the release from different heights is shown by the ens dimension.

        """

        fig = plt.figure(1, figsize=(10, 5))
        ax1 = fig.add_subplot(1, 1, 1)
        if tii is None:
            tii = self.time_index(daterange[0])
        print("tii", tii)
        cdump = self.concmult * self.cdump_mass_hash[tii]
        volcat = self.volcat_avg_hash[tii]
        if not zii:
            csum = cdump.sum(dim="ens")
        else:
            csum = cdump.isel(ens=zii)
            print(cdump.ens.values[zii])
            # print(csum.ens.values)
            # print(self.sourcehash[str(csum.ens.values)])
        # print(cdump.time)
        # print(csum.coords)
        # volcat.plot.pcolormesh(x='longitude',y='latitude',levels=levels,ax=ax1)
        # cdump.sum(dim='ens').plot.contour(x='longitude',y='latitude',ax=ax2)
        try:
            # plt.pcolormesh(csum.longitude, csum.latitude, np.log10(csum),cmap='Reds',shading='nearest')
            cbm = ax1.pcolormesh(
                csum.x, csum.y, np.log10(csum), cmap="Reds", shading="nearest"
            )
        except:
            print("FAILED max value", np.max(csum))
            print("------------------")
        # plt.pcolormesh(csum.longitude, csum.latitude, csum,cmap='Reds',shading='nearest')
        # cb= csum.plot.pcolormesh(x='longitude',y='latitude',cmap='viridis',ax=ax2)
        cb = ax1.pcolormesh(
            volcat.x,
            volcat.y,
            np.log10(volcat),
            cmap="Blues",
            shading="nearest",
            alpha=0.5,
        )
        # cb = plt.scatter(volcat.longitude, volcat.latitude, c=np.log10(volcat),s=2,cmap='Blues')
        # cb = plt.scatter(volcat.longitude, volcat.latitude, c=volcat.values,s=2,cmap='viridis',levels=levels)
        # vals = np.log10(volcat)
        # cb = plt.contour(volcat.x, volcat.y, vals.fillna(0),cmap='viridis',levels=[0,1,10,100])
        plt.colorbar(cbm)
        plt.tight_layout()
        # return csum.copy()
        return csum.copy()


    def set_bias_correction(self, slope, intercept, dfcdf=pd.DataFrame()):
        # InverseAsh class
        """ """
        self.slope = slope
        self.intercept = intercept
        self.dfcdf = dfcdf

    def manage_cdump(
        self,
        tii,
        coarsen=None,
        slope=None,
        intercept=None,
        dfcdf=pd.DataFrame(),
        cii=None,
        coarsen_max=None,
    ):
        # InverseAsh class
        """
        tii : datetime.datetime OR integer
        coarsen :  int : if > 0 then use the xarray coarsen method to coarsen the data using mean.
        coarsen_max :  int : if > 0 then use the xarray coarsen method to coarsen the data using max.

        bias correction.
        For bias correction, first will try to use the dfcdf DataFrame if input.
        If it is empty will check to see if the dfcdf attribute is not empty.
        If the dfcdf attribute is not empty then it will use that.
        If the dfcdf attribute and dfcdf input are both empty then will
        use the slope and intercept inputs. If those are both None then no bias correction is applied.

        cii   : datetime.datetime. used with dfcdf
        dfcdf : pandas DataFrame with information for bias correction.
        slope : float : if not None used for bias correction if  dfcdf is empty
        intecept : float : if not None used for bias correction if  dfcdf is empty
        """

        # make sure that tii is datime and iii is integer - index for tii.
        if isinstance(tii, int):
            iii = tii
            tii = self.get_time(iii)
        elif isinstance(tii, datetime.datetime):
            iii = self.time_index(tii)

        cdump = self.cdump_mass_hash[iii] * self.concmult

        if not cii:
            cii = tii

        # cii = tii-datetime.timedelta(hours=2)
        # if pandas dataframe available use that.
        if dfcdf.empty and not self.dfcdf.empty:
            dfcdf = self.dfcdf
        if not dfcdf.empty:
            # create xarray with slope and intercept with coordinate of ens.
            temp = dfcdf[dfcdf["time"] == cii]
            temp = temp[["slope", "intercept", "ens"]]
            temp = temp.set_index("ens")
            xrt = temp.to_xarray()
            # apply slope and intercept correction
            cdump_prev = cdump.copy()
            # print('cdump A', cdump.ens.values)
            cdump = cdump * (1 - xrt.slope)
            cdump = cdump - xrt.intercept
            f2 = xr.where(cdump > 0, cdump, 0)
            # print('xrt', xrt.ens.values)
            # print('f2', f2.ens.values)
            # print('cdump', cdump.ens.values)
            # print('prev', cdump_prev.ens.values)
            f2 = xr.where(cdump_prev < 0.0001, 0, f2)
            cdump = f2
        # otherwise use slope and intercept if those are specified.
        else:
            if not slope and self.slope:
                slope = self.slope
            cdump_prev = cdump
            if slope:
                # print('multiplying cdump by {}'.format(1-slope))
                cdump = cdump * (1 - slope)
            if not intercept and self.intercept:
                intercept = self.intercept
            if intercept:
                # print('shifting cdump by {}'.format(-1*intercept))
                cdump = cdump - intercept
                # if intercept is positive then remove negative values.
                f2 = xr.where(cdump > 0, cdump, 0)
                # if intercept is negative then do not add to values that were previously 0.
                f2 = xr.where(cdump_prev < 0.0001, 0, f2)
                cdump = f2

        # coarsen comes last.
        # find the mean value
        if coarsen:
            cdump = cdump.coarsen(x=coarsen, boundary="trim").mean()
            cdump = cdump.coarsen(y=coarsen, boundary="trim").mean()
        # use maximum value in scene
        elif coarsen_max:
            cdump = cdump.fillna(0)
            cdump = cdump.coarsen(x=coarsen_max, boundary="trim").max()
            cdump = cdump.coarsen(y=coarsen_max, boundary="trim").max()

        return cdump

    def get_pair(
        self,
        tii,
        coarsen=None,
        slope=None,
        intercept=None,
        dfcdf=pd.DataFrame(),
        cii=None,
        coarsen_max=None,
    ):
        # InverseAsh class
        """ """
        if isinstance(tii, int):
            iii = tii
        elif isinstance(tii, datetime.datetime):
            iii = self.time_index(tii)
        volcat = self.volcat_avg_hash[iii]
        cdump = self.manage_cdump(
            tii, coarsen, slope, intercept, dfcdf, cii, coarsen_max
        )
        if coarsen:
            volcat = volcat.coarsen(x=coarsen, boundary="trim").mean()
            volcat = volcat.coarsen(y=coarsen, boundary="trim").mean()
        elif coarsen_max:
            volcat = volcat.coarsen(x=coarsen_max, boundary="trim").max()
            volcat = volcat.coarsen(y=coarsen_max, boundary="trim").max()
        return volcat, cdump

    def match_volcat(self, forecast):
        # InverseAsh class
        """ """
        time = pd.to_datetime(forecast.time.values)
        tii = self.time_index(time)
        volcat = self.volcat_avg_hash[tii]
        return volcat

    def compare_forecast_dist(self, forecast, thresh=None):
        # InverseAsh class
        """ """
        from utilhysplit.evaluation import statmain

        volcat = self.match_volcat(forecast)
        evals = forecast.values
        vvals = volcat.values
        exval, eyval = statmain.nancdf(evals.flatten(), thresh)
        vxval, vyval = statmain.nancdf(vvals.flatten(), thresh)
        fig = plt.figure(1)
        ax = fig.add_subplot(1, 1, 1)
        ax.step(exval, eyval, "-r", label="HYSPLIT")
        ax.step(vxval, vyval, "-b", label="Volcat")
        ks1, ks2 = statmain.kstestnan(evals.flatten(), vvals.flatten(), thresh)
        try:
            print(
                "Kolmogorov-Smirnov Parameter {} {}".format(
                    np.max(np.abs(ks1)), np.max(np.abs(ks2))
                )
            )
        except:
            return 1
        return np.max(np.abs(ks1))

    def remove_nans(self, data, thresh=None):
        # InverseAsh class
        """ """
        # remove nans
        data2 = data[~np.isnan(data)]
        if thresh:
            data2 = data2[data2 > thresh]
            # apply threshold.
            # vpi = data2 < thresh
            # data2[vpi] = np.nan
            # data2 = data2[~np.isnan(data2)]
        return data2

    def compare(self, time,thresh=None):
        # InverseAsh class
        """ """
        tii = self.time_index(time)
        print("tii", tii)
        volcat = self.volcat_avg_hash[tii]
        forecast = self.cdump_mass_hash[tii]
        #volcat = self.match_volcat(forecast)
        evals = self.remove_nans(forecast.values.flatten(), thresh)
        vvals = self.remove_nans(volcat.values.flatten(), thresh)
        # number of pixles above threshold in each.
        forecast_parea = len(evals)
        volcat_parea = len(vvals)
        return forecast_parea, volcat_parea

    def compare_forecast_alt(
        self, forecast, cmap="Blues", ptype="pcolormesh", vloc=None
    ):
        # InverseAsh class
        """ """
        # forecast should be an xarray in mass loading format with no time dimension.
        sns.set()
        sns.set_style("whitegrid")
        fig = plt.figure(figsize=[15, 5])
        ax1 = fig.add_subplot(1, 3, 1)
        ax2 = fig.add_subplot(1, 3, 2)
        ax3 = fig.add_subplot(1, 3, 3)

        time = pd.to_datetime(forecast.time.values)
        tii = self.time_index(time)
        print("tii", tii)
        volcat = self.volcat_avg_hash[tii]
        evals = forecast.values
        vpi = evals < 0.001
        evals[vpi] = np.nan
        norm = get_norm(volcat, forecast)

        # evals = np.log10(evals)
        # vvals = np.log10(volcat.values)
        vvals = volcat.values
        vpi = vvals < 0.001
        vvals[vpi] = np.nan
        # clevels = [0.2,2,5,10]
        clevels = [0.2, 10]
        if ptype == "pcolormesh":
            cb = ax1.pcolormesh(
                volcat.longitude,
                volcat.latitude,
                vvals,
                norm=norm,
                cmap=cmap,
                shading="nearest",
            )
            cb2 = ax2.pcolormesh(
                forecast.longitude,
                forecast.latitude,
                evals,
                norm=norm,
                cmap=cmap,
                shading="nearest",
            )
            # cb2 = ax3.pcolormesh(forecast.longitude, forecast.latitude,evals,norm=norm, cmap=cmap,shading='nearest')
            ax3.contourf(
                volcat.longitude,
                volcat.latitude,
                volcat.values,
                levels=clevels,
                cmap="Reds",
            )
            ax3.contour(
                forecast.longitude,
                forecast.latitude,
                evals,
                levels=clevels,
                cmap="viridis",
            )
        plt.title(time.strftime("%Y %m/%d %H:%M UTC"))
        plt.colorbar(cb, ax=ax1)
        plt.colorbar(cb2, ax=ax2)
        ylim = ax1.get_ylim()
        ax2.set_ylim(ylim)
        ax3.set_ylim(ylim)
        xlim = ax1.get_xlim()
        ax2.set_xlim(xlim)
        ax3.set_xlim(xlim)
        if vloc:
            ax1.plot(vloc[0], vloc[1], "y^")
            ax2.plot(vloc[0], vloc[1], "y^")
            ax3.plot(vloc[0], vloc[1], "y^")
        return fig

    def compare_forecast(self, time,cmap="Blues", ptype="pcolormesh", vloc=None):
        # InverseAsh class
        """
        forecast should be an xarray in mass loading format with no time dimension.
        """
        #time = pd.to_datetime(forecast.time.values)
        tii = self.time_index(time)
        print("tii", tii)
        volcat = self.volcat_avg_hash[tii]
        forecast = self.cdump_mass_hash[tii]
        for ens in forecast.ens.values:
            fff = forecast.sel(ens=ens)  
            print(ens)
            fig = compare_forecast(volcat,fff,cmap,ptype,vloc) 
            plt.show()
        return fig

    def compare_plots(self, daterange, levels=None):
        # InverseAsh class
        """ """
        fig = plt.figure(1, figsize=(10, 5))
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)
        tii = self.time_index(daterange[0])
        cdump = self.concmult * self.cdump_mass_hash[tii]
        volcat = self.volcat_avg_hash[tii]
        volcat = xr.where(volcat>0,volcat,np.nan)
        csum = cdump.max(dim="ens")
        csum = xr.where(csum>0,csum,np.nan)
        volcat.plot.pcolormesh(x="longitude", y="latitude", levels=levels, ax=ax1)
        # cdump.sum(dim='ens').plot.contour(x='longitude',y='latitude',ax=ax2)
        # plt.pcolormesh(csum.longitude, csum.latitude, np.log10(csum),cmap='Reds')
        # plt.pcolormesh(csum.longitude, csum.latitude, csum,cmap='Reds',shading='nearest')
        cb = csum.plot.pcolormesh(x="longitude", y="latitude", cmap="viridis", ax=ax2,levels=levels)
        # cb = plt.pcolormesh(volcat.longitude, volcat.latitude, np.log10(volcat),cmap='Blues',levels=levels)
        # cb = plt.scatter(volcat.longitude, volcat.latitude, c=np.log10(volcat),s=2,cmap='Blues')
        # cb = plt.scatter(volcat.longitude, volcat.latitude, c=volcat.values,s=2,cmap='viridis',levels=levels)
        # cb = plt.contour(volcat.longitude, volcat.latitude, np.log10(volcat),cmap='Blues')
        # plt.colorbar(cb)
        plt.tight_layout()
        return ax1, ax2

    def prepare_times(self):
        dlist = self.get_times()
        for time in dlist:
            self.prepare_one_time(time)

    def add_volcat_hash(self, volcat_hash):
        """ """
        # allow to update from another instance of the class.
        self.volcat_avg_hash.update(volcat_hash)

    def add_cdump_hash(self, cdump_hash):
        # allow to update from another instance of the class.
        """ """
        self.cdump_mass_hash.update(cdump_hash)

    def change_sampling_time(self, value):
        if value not in ["start", "end"]:
            print("time description not recognized {}".format(value))
            print("value should be {} or {}".format("start", "end"))
        else:
            self.cdump.attrs.update({"time description": value})

    def get_sampling_time(self, default="start"):
        if "time description" in self.cdump.attrs.keys():
            if "start" in str.lower(self.cdump.attrs["time description"]):
                return "start"
            if "end" in str.lower(self.cdump.attrs["time description"]):
                return "end"
        else:
            return default

    def mass_load_modified(self):
        # logic for only getting mass loading in certain levels.

        # have an array of top heights.
        #
        return -1

    def coarsen(self, num=3):
        if not self.coarsen:
            self.original_cdump_hash = self.cdump_hash.copy()
            self.original_volcat_avg_hash = self.volcat_avg_hash.copy()
        self.coarsen = num
        for key in self.cdump_mass_hash.keys():
            c2 = self.original_cdump_hash[key].coarsen(x=num).mean()
            c2 = c2.coarsen(y=num).mean()
            v2 = self.original_volcat_avg_hash[key].coarsen(x=num).mean()
            v2 = v2.coarsen(y=num).mean()
            self.cdump_hash[key] = c2
            self.volcat_hash[key] = v2

    def clear_one_time(self, daterange, tii=None):
        """
        remove entry from paired dictionary
        """

        if not isinstance(tii,int):
            if self.get_sampling_time() == "start":
                # use first date in daterange list
                model_tii = 0
            elif self.get_sampling_time() == "end":
                # use second date in daterange list
                model_tii = 1
            tii = self.time_index(daterange[model_tii])
        if tii in self.cdump_mass_hash.keys():
            del self.cdump_hash[tii]
        if tii in self.volcat_hash.keys():
            del self.volcat_hash[tii]


    def prepare_one_time(self, daterange, htoptions=0, zvals=None, verbose=False):
        """
        daterange : list of datetime objects
        htoptions : boolean : capability not fully implemented.
        zvals : list of vertical levels (in meters) to use to calculate mass loading.


        Populates the folling dictionaries:
        self.cdump_hash[tii] = cdump_a
        self.volcat_avg_hash[tii] = avra
        self.volcat_ht_hash[tii] = maxvhra
        """
        # currently must coarsen all data at the same time.
        if self.coarsen:
            print(
                "Warning: Adding new data after some data has already been coarsened."
            )
        vdir = self.vdir
        #vid = self.vid
        # key for hashes is determined from times in cdump file.
        # check whether cdump time is start or end of sampling time.
        if self.get_sampling_time() == "start":
            # use first date in daterange list
            model_tii = 0
        elif self.get_sampling_time() == "end":
            # use second date in daterange list
            model_tii = 1
        # if st == "end":
        #    model_tii = 1
        # else:
        #    model_tii = 0
        tii = self.time_index(daterange[model_tii])

        if tii in self.cdump_mass_hash.keys():
           print('Time {} already prepared'.format(tii))
           return self.cdump_mass_hash[tii], self.volcat_hash[tii]

        if not isinstance(tii, int):
            print("No time found for {}".format(daterange[0]))
            return None, None
        if isinstance(zvals, np.ndarray):
            zvals = list(zvals)
        cdump_a = hysplit.hysp_massload(
            self.cdump.sel(time=daterange[model_tii]), zvals=zvals
        )
        if htoptions == 1:
            cdump_b = self.cdump.sel(time=daterange[model_tii])

        vset = self.get_volcat(daterange, verbose)
        buf = 5
        # clip the volcat array before aligning.
        try:
            vra = vset.ash_mass_loading
        except:
            return None, None

        # also get the volcat observed height.
        try:
            vhra = vset.ash_cloud_height
        except:
            return None, None
        vset.close()
        a1, a2, b1, b2 = self.clip(vra.sum(dim="time"), buf=buf)
        vra = vra.transpose("time", "y", "x")
        # vra has dimensions of time, y, x
        vra = vra[:, a1:a2, b1:b2]

        vhra = vhra.transpose("time", "y", "x")
        vhra = vhra[:, a1:a2, b1:b2]
        # clip the cdump array before aligning.
        #if "ens" in cdump_a.coords:
        if "ens" in cdump_a.dims:
            cdump_a = cdump_a.transpose("ens", "y", "x")
            if htoptions == 1:
                cdump_b = cdump_b.transpose("ens", "z,", "y", "x")
            dummy = cdump_a.sum(dim="ens")
        else:
            cdump_a = cdump_a.transpose("y", "x")
            if htoptions == 1:
                cdump_b = cdump_b.transpose("z", "y", "x")
            dummy = cdump_a

        try:
            a1, a2, b1, b2 = self.clip(dummy, buf=5)
        except:
            print("dummy cannot be clipped", dummy)
        if "ens" in cdump_a.dims:
            cdump_a = cdump_a[:, a1:a2, b1:b2]
            if htoptions == 1:
                cdump_b = cdump_b[:, :, a1:a2, b1:b2]
        else:
            cdump_a = cdump_a[a1:a2, b1:b2]
            if htoptions == 1:
                cdump_b = cdump_b[:, a1:a2, b1:b2]

        # align the grids.
        tol = 0.01
        if not hgu.compare_grids(cdump_a, vra, verbose=False, tolerance=tol):
            if not hgu.compare_grids_compat(cdump_a, vra):
                avra = vra.fillna(0).mean(dim="time")
                self.volcat_avg_hash[tii] = avra
                print("prepare_one_time: grids are not compatible")
                return False
            if verbose:
                print("grids are compatible. changing grids")
            vatt = hgu.find_grid_specs(vra)
            dlat = vatt["Latitude Spacing"]
            dlon = vatt["Longitude Spacing"]
            crnr = vatt["llcrnr longitude"]
            newatt = hgu.create_global_grid(dlat, dlon, crnr)
            cdump_a = hgu.change_grid(cdump_a, newatt)
            vhra = hgu.change_grid(vhra, newatt)
            # print(hgu.find_grid_specs(cdump_a))
            # print(vatt)
        if hgu.compare_grids(cdump_a, vra, verbose=True, tolerance=tol):
            dummy, vhra = hgu.align_grids(cdump_a, vhra, tolerance=tol)
            cdump_a, vra = hgu.align_grids(cdump_a, vra, tolerance=tol)
        else:
            print("----------")
            hgu.compare_grids(cdump_a, vra, verbose=True)
            print("prepare_one_time: grids cannot be aligned")
            print(daterange)
            # print(vset)
            print("----------")
            return cdump_a, vra

        # take the average of the volcat data.
        avra = vra.fillna(0).mean(dim="time")
        maxvhra = vhra.fillna(0).max(dim="time")

        self.cdump_mass_hash[tii] = cdump_a
        self.volcat_avg_hash[tii] = avra
        self.volcat_ht_hash[tii] = maxvhra
        return cdump_a, avra

    def compare_time_ave(self, daterange):
        """
        creates plots illustrating how time averaging affects
        volcat data.
        """
        sns.set()
        vset = self.get_volcat(daterange)
        vset = vset.ash_mass_loading
        for time in vset.time.values:
            temp = vset.sel(time=time)
            #a1, a2, b1, b2 = self.clip(temp.fillna(0))
            #temp = temp[a1:a2, b1:b2]
            temp.plot.pcolormesh(x="longitude", y="latitude")
            plt.show()
        avra = vset.fillna(0).mean(dim="time")
        avra2 = vset.mean(dim="time")
        #a1, a2, b1, b2 = self.clip(avra)
        #avra = avra[a1:a2, b1:b2]
        #avra2 = avra2[a1:a2, b1:b2]
        print("Average with nans set to 0")
        avra = xr.where(avra<=0,np.nan,avra)
        avra.plot.pcolormesh(x="longitude", y="latitude")
        plt.title('average of values with 0')      
        plt.show()
        print("Average of values above 0")
        avra2.plot.pcolormesh(x="longitude", y="latitude")
        plt.title('average of values above 0')      
        plt.show()
        diff = avra2.fillna(0) - avra
        diff.plot.pcolormesh(x="longitude", y="latitude")
        plt.title('difference')


def compare_forecast(volcat,forecast, cmap="Blues", ptype="pcolormesh", vloc=None):
    # InverseAsh class
    """
    forecast should be an xarray in mass loading format with no time dimension.
    """
    sns.set()
    sns.set_style("whitegrid")
    fig = plt.figure(figsize=[15, 5])
    ax1 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 3, 2)
    ax3 = fig.add_subplot(1, 3, 3)

    time = pd.to_datetime(forecast.time.values)
    #tii = self.time_index(time)
    #print("tii", tii)
    #volcat = self.volcat_avg_hash[tii]
    evals = forecast.values
    vpi = evals < 0.001
    evals[vpi] = np.nan
    norm = get_norm(volcat, forecast)

    # evals = np.log10(evals)
    # vvals = np.log10(volcat.values)
    vvals = volcat.values
    vpi = vvals < 0.001
    vvals[vpi] = np.nan
    # clevels = [0.2,2,5,10]
    clevels = [0.02,0.2,2,5,10,20]
    if ptype == "pcolormesh":
        cb = ax1.pcolormesh(
            volcat.longitude,
            volcat.latitude,
            vvals,
            norm=norm,
            cmap=cmap,
            shading="nearest",
        )
        cb2 = ax2.pcolormesh(
            forecast.longitude,
            forecast.latitude,
            evals,
            norm=norm,
            cmap=cmap,
            shading="nearest",
        )
        # cb2 = ax3.pcolormesh(forecast.longitude, forecast.latitude,evals,norm=norm, cmap=cmap,shading='nearest')
        ax3.contourf(
            volcat.longitude,
            volcat.latitude,
            volcat.values,
            levels=clevels,
            cmap="Reds",
        )
        ax3.contour(
            forecast.longitude,
            forecast.latitude,
            evals,
            levels=clevels,
            cmap="viridis",
        )
    plt.title(time.strftime("%Y %m/%d %H:%M UTC"))
    plt.colorbar(cb, ax=ax1)
    plt.colorbar(cb2, ax=ax2)
    ylim = ax1.get_ylim()
    ax2.set_ylim(ylim)
    ax3.set_ylim(ylim)
    xlim = ax1.get_xlim()
    ax2.set_xlim(xlim)
    ax3.set_xlim(xlim)
    if vloc:
        ax1.plot(vloc[0], vloc[1], "y^")
        ax2.plot(vloc[0], vloc[1], "y^")
        ax3.plot(vloc[0], vloc[1], "y^")
    return fig

def get_norm(model, r2, logscale=False):
    # InverseAsh class
    """ """
    if isinstance(r2, np.ndarray):
        rval = r2
    else:
        rval = r2.values

    m_max = np.nanmax(model)
    v_max = np.nanmax(rval)
    m_min = np.nanmin(model)
    v_min = np.nanmin(rval)
    p_min = np.nanmin([m_min, v_min])
    p_max = np.nanmax([m_max, v_max])
    lower_color = 0.8
    if logscale:
        norm = mpl.colors.LogNorm(vmin=lower_color * p_min, vmax=p_max)
    else:
        norm = mpl.colors.Normalize(vmin=p_min, vmax=p_max)
    return norm

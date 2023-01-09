import datetime
import os
from subprocess import call

import cartopy.crs as ccrs
import cartopy.feature as cfeat
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
from utilvolc.ash_inverse import InverseAsh, InverseAshEns
from utilvolc.runhelper import Helper


class InverseEns(InverseAshEns):
    def placeholder(self):
        return -1

    def add_invlist(self, tdirlist, vdir, vid, configdir, configfile, ensdim, verbose):
        for hruns in zip(tdirlist, self.fnamelist):
            self.invlist.append(
                InverseSO2(
                    hruns[0],
                    hruns[1],
                    vdir,
                    vid,
                    configdir,
                    configfile,
                    ensdim=ensdim,
                    verbose=verbose,
                )
            )

    def add_phash(self):
        self.invlist[0].add_phash()
        self.phash = self.invlist[0].phash

    def add_obs(self, obslist):
        for iii, hruns in enumerate(self.invlist):
            hruns.add_obs(obslist[iii])

    def make_tcm_mult(
        self, remove_cols=True, remove_rows=True, remove_sources=None, remove_ncs=0
    ):
        for hrun in self.invlist:
            tcm = hrun.make_tcm(remove_cols, remove_rows, remove_sources, remove_ncs)


class InverseSO2(InverseAsh):
    def __init__(
        self,
        tdir,
        fname,
        vdir,
        vname,
        configdir="./",
        configfile=None,
        ensdim="ens",
        verbose=False,
    ):
        """
        configfile : full path to configuration file.
        """

        self.tdir = tdir  # directory for hysplit output
        self.fname = fname  # name of hysplit output

        self.vdir = vdir  # directory for volcat data
        self.vname = vname  # volcano id.

        self.n_ctrl = 0  # determined in the write_tcm method.
        # needed for Parameters_in.dat input into inversion algorithm

        # keep volcat arrays for different averaging times.
        self.volcat_hash = {}
        # prepare_one_time method adds data to these dictionaries.
        # the data in cdump_hash and volcat_avg_hash have been aligned.
        # and are on the same grid.
        self.volcat_avg_hash = {}
        self.volcat_ht_hash = {}
        self.cdump_hash = {}

        # multiplication factor if more than 1 unit mass released.
        self.concmult = 1
        self.get_cdump(tdir, fname, verbose, ensdim)
        self.add_config_info(configdir, configfile)
        self.add_phash()

    def add_phash(self, phash=None):
        psizelist = ["pSO2", "pSO4"]
        spnum = np.arange(1, len(psizelist) + 1)
        self.phash = dict(zip(psizelist, spnum))

    def close_arrays(self):
        self.cdump.close()
        # self.massload.close()

    def copy(self):
        iacopy = InverseAsh(self.tdir, self.fname, self.vdir, self.vid, self.n_ctrl)
        iacopy.volcat_avg_hash = self.volcat_avg_hash
        iacopy.volcat_cdump_hash = self.cdump_hash
        icacopy.concmult = self.concmult
        return icacopy

    def print_summary(self):
        print("Observations availalbe in volcat_avg_hash")
        print(self.volcat_avg_hash.keys())
        print("times in cdump file")
        self.print_times()

    def get_cdump(self, tdir, fname, verbose=False, ensdim="ens"):
        # hysplit output. xarray.
        cdump = xr.open_dataset(os.path.join(tdir, fname), engine="netcdf4")
        if not hysplit.check_grid_continuity(cdump):
            print("Grid not continuous")
        if verbose:
            print("opened", tdir, fname)
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
        self.cdump = cdump.fillna(0)

    def add_inp(self, configdir, configfile):
        self.inp = get_inp_hash(configdir, configfile)

    def set_concmult(self, concmult):
        self.concmult = concmult

    def get_volcat(self, daterange, verbose=False):
        vdir = self.vdir
        vid = self.vid
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
            vset = xr.concat(das, dim="time")
        else:
            print("No volcat files found ")
            return xr.DataArray()
        return vset

    def clip(self, dummy, buf=0):
        return None

    def print_times(self):
        timelist = [pd.to_datetime(x) for x in self.cdump.time.values]
        for time in timelist:
            print(time.strftime("%Y %m %d %H:%Mz"))

    def get_time(self, tii):
        timelist = [pd.to_datetime(x) for x in self.cdump.time.values]
        return timelist[tii]

    def time_index(self, time):
        timelist = [pd.to_datetime(x) for x in self.cdump.time.values]
        try:
            iii = timelist.index(time)
        except:
            # print('timelist', timelist)
            iii = None
        return iii

    def get_massload(self, vdata, cdump):
        # compute mass loading only for times which have matching observations.
        times = vdata.time.unique()
        times = [pd.to_datetime(x) for x in times]
        times = [(x.year, x.month, x.day, x.hour) for x in times]
        times = list(set(times))
        times = [pd.datetime(x[0], x[1], x[2], x[3]) for x in times]
        cdump = cdump.sel(time=times)
        # century = int(x[0]/100.0)*100
        massload = hysplit.hysp_massload(cdump)
        return massload

    def make_tcm(
        self, remove_cols=False, remove_rows=False, remove_sources=None, remove_ncs=0
    ):

        vdata = (
            self.vdata
        )  # pandas dataframe with columns of 'mass', 'lon', 'lat', 'time'
        mass = self.get_massload(
            self.vdata, self.cdump
        )  # mass loading xarray from HYSPLIT
        tcm = []
        iii = 0
        # dictionary with key being index of cdump and value being tcm row.
        thash = {}
        for iindex, vrow in vdata.iterrows():
            if iii % 1000 == 0:
                print("working on ", iii, len(vdata), vrow)

            # get concentrations valid for observation time.
            tlist = [x for x in mass.time.values if x < vrow.time]
            tlist.sort(reverse=True)
            tii = tlist[0]
            mtemp = mass.sel(time=tii)

            # get concentrations valid at lat-lon point
            try:
                xii, yii = mtemp.monet.nearest_ij(lon=vrow.lon, lat=vrow.lat)
            # if error is returned trying changing longitude to negative or positive value.
            except:
                fixed = True
                if vrow.lon < 0:
                    lon = 360 + vrow.lon
                    try:
                        xii, yii = mtemp.monet.nearest_ij(lon=lon, lat=vrow.lat)
                    except:
                        print("ERROR longitude not within grid ", vrow.lon, lon)
                        fixed = False
                if vrow.lon > 0:
                    lon = vrow.lon - 360
                    try:
                        xii, yii = mtemp.monet.nearest_ij(lon=lon, lat=vrow.lat)
                    except:
                        print("ERROR longitude not within grid ", vrow.lon, lon)
                        fixed = False
                if not fixed:
                    continue
            trow = mtemp.isel(x=xii, y=yii)
            # append observation onto end of row.
            tcm_row = np.append(trow, vrow.mass)
            # check to see if other obs corresponded to same model data.
            if (tii, xii, yii) in thash.keys():
                prev_value = thash[(tii, xii, yii)]
                new_value = tcm_row
                # print('-------------------------------------')
                # print('SAME ', prev_value[-1], new_value[-1])
                # print('-------------------------------------')
                # replace with average value
                total_value = (prev_value[-1] + new_value[-1]) / 2.0
                total_row = np.append(trow, total_value)
                thash[(tii, xii, yii)] = total_row
            else:
                thash[(tii, xii, yii)] = tcm_row
            # if iii==0:
            #   tcm = tcm_row
            # elif iii==1:
            #   tcm = np.append([tcm],[tcm_row],axis=0)
            # else:
            #   tcm = np.append(tcm,[tcm_row],axis=0)
            iii += 1
        # convert the hash to a numpy array.
        tcm = np.array(list(thash.values()))

        self.tcm_columns = mtemp.ens.values
        self.tcm = tcm
        return tcm

    # def plot_tcm(self):
    #    cb = plt.pcolormesh(np.log10(self.tcm),cmap='tab20')
    #    plt.colorbar(cb)

    # def write_tcm(self, name):
    #    astr = ''
    #    sep = ' '
    #    hstr = '' # header string
    # print(self.tcm.shape)\
    # this is number of columns minus 1.
    # and needs to be input into Parameters.in.dat
    #    self.n_ctrl = self.tcm.shape[1]-1
    # print('N_ctrl {}'.format(self.tcm.shape[1]-1))
    # print('output file {}'.format(name))
    #    for iii, line in enumerate(self.tcm):
    #        for jjj, val in enumerate(line):
    #            if iii==0:
    # write a dummy header line.
    #               hstr += '43637.750' + sep
    #            if not np.isnan(val): astr += '{:1.5e}'.format(val)
    #            else: astr += '{:1.4e}'.format(0)
    #            astr += sep
    #        astr += '\n '
    #        if iii==0: hstr += '\n'
    #    with open(name, 'w') as fid:
    #        fid.write(hstr + astr)
    #    return hstr + astr

    def make_outdat(self, dfdat):
        """
        make_outdat for InverseAsh class.
        There is a duplicate method in the InverseAshPart class.
        dfdat : pandas dataframe output by InverseOutDat class get_emis method.
        Returns
        vals : tuple (date, height, emission mass)
        """
        # matches emissions from the out.dat file with
        # the date and time of emission.
        # uses the tcm_columns array which has the key
        # and the sourehash dictionary which contains the information.
        datelist = []
        htlist = []
        valra = []
        for val in zip(self.tcm_columns, dfdat[1]):
            shash = self.sourcehash[val[0]]
            datelist.append(shash["sdate"])
            htlist.append(shash["bottom"])
            valra.append(val[1])
        vals = list(zip(datelist, htlist, valra))
        return vals

    def make_outdat_df(self, dfdat, savename=None, part="basic"):

        # dfdat : pandas dataframe output by InverseOutDat class get_emis method.
        # this is a list of tuples (source tag), value from emissions
        vals = self.make_outdat(dfdat)
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

    def plot_outdat(
        self, vals, log=False, fignum=1, cmap="Blues", unit="kg/s", thresh=0
    ):
        """
        vals is output by make_outdat.
        """
        fig = plt.figure(fignum, figsize=(10, 5))
        vals = list(zip(*vals))
        sns.set()
        sns.set_style("whitegrid")
        # output in kg/s?/
        if unit == "kg/s":
            emit = np.array(vals[2]) / 1.0e3 / 3600.0
        elif unit == "kg/h":
            emit = np.array(vals[2]) / 1.0e3
        elif unit == "g/h":
            emit = np.array(vals[2]) / 1.0
        vpi = np.where(emit < thresh)
        emit[vpi] = 0
        ht = np.array(vals[1]) / 1e3
        if not log:
            cb = plt.scatter(vals[0], ht, c=emit, s=100, cmap=cmap, marker="s")
        else:

            cb = plt.scatter(
                vals[0], ht, c=np.log10(emit), s=100, cmap=cmap, marker="s"
            )
            # cb = plt.pcolormesh(vals[0],ht,emit,cmap=cmap)
        cbar = plt.colorbar(cb)
        cbar.ax.set_ylabel(unit)
        fig.autofmt_xdate()

    def plot_outdat_profile(
        self, dfdat, fignum=1, unit="kg", ax=None, clr="--ko", alpha=1, lw=1
    ):
        if not ax:
            sns.set()
            sns.set_style("whitegrid")
            fig = plt.figure(fignum, figsize=(10, 5))
            ax = fig.add_subplot(1, 1, 1)
        df = self.make_outdat_df(dfdat, part="cols")
        ts = df.sum(axis=1)

        sns.set()
        # the one represents one hour time period.
        yval = ts.values * 1 / 1e12
        ax.plot(yval, ts.index.values / 1000.0, clr, alpha=alpha, linewidth=lw)
        ax.set_xlabel("Tg of mass emitted", fontsize=15)
        ax.set_ylabel("Height (km)", fontsize=15)
        totalmass = yval.sum()
        # print('total {} Tg'.format(totalmass))
        return ax, df

    def plot_outdat_ts(
        self, dfdat, log=False, fignum=1, unit="kg/s", ax=None, clr="--ko"
    ):
        # plots time series of MER. summed along column.
        # dfdat : pandas dataframe output by InverseOutDat class get_emis method.
        if not ax:
            sns.set()
            sns.set_style("whitegrid")
            fig = plt.figure(fignum, figsize=(10, 5))
            ax = fig.add_subplot(1, 1, 1)
        df = self.make_outdat_df(dfdat, part="index")
        sns.set()
        ts = df.sum()
        if unit == "kg/s":
            yval = ts.values / 3.6e6
        elif unit == "g/h":
            yval = ts.values
        ax.plot([x[1] for x in ts.index.values], yval, clr)
        # fig.autofmt_xdate()
        ax.set_ylabel("MER {}".format(unit), fontsize=15)
        return ax, df

    def plot_out2dat_times(self, df2, cmap="viridis"):
        return -1

    def plot_out2dat_scatter(self, tiilist, df2, vloc, cmap="Blues"):
        if isinstance(tiilist, int):
            tiilist = [tiilist]
        sns.set()
        ppp = 0
        # tii = self.time_index(daterange[0])
        modelall = df2["model"].values
        nnn = 0
        for tii in tiilist:
            print(self.cdump.time.values[tii])
            fig = plt.figure(figsize=[10, 5])
            ax1 = fig.add_subplot(1, 2, 1)
            ax2 = fig.add_subplot(1, 2, 2)
            volcat = self.volcat_avg_hash[tii]
            shape = volcat.shape
            lon = self.lonlist[tii - 1]
            lat = self.latlist[tii - 1]
            model = modelall[nnn : nnn + len(lon)]
            volcat = self.volcat_avg_hash[tii]
            r2 = volcat.where(volcat > 0)
            cb = ax1.scatter(lon, lat, c=model, cmap=cmap, s=10, marker="o")
            cb2 = ax2.scatter(
                volcat.longitude,
                volcat.latitude,
                c=r2.values,
                s=10,
                cmap=cmap,
                marker="o",
            )
            nnn = len(lon)
            if vloc:
                ax1.plot(vloc[0], vloc[1], "y^")
                ax2.plot(vloc[0], vloc[1], "y^")
            plt.colorbar(cb, ax=ax1)
            plt.colorbar(cb2, ax=ax2)

    def plot_out2dat(self, tiilist, df2, cmap="viridis", vloc=None, ptype="pcolormesh"):
        # tiilist needs to be same order as for tcm.
        # df2 is from InverseOutDat get_conc method.
        # doesn't work when only later time periods are used!
        if isinstance(tiilist, int):
            tiilist = [tiilist]
        sns.set()
        ppp = 0
        # tii = self.time_index(daterange[0])
        for tii in tiilist:
            print(self.cdump.time.values[tii])
            fig = plt.figure(figsize=[10, 5])
            ax1 = fig.add_subplot(1, 2, 1)
            ax2 = fig.add_subplot(1, 2, 2)
            volcat = self.volcat_avg_hash[tii]
            shape = volcat.shape
            model = df2["model"].values
            temp = model[ppp : ppp + shape[0] * shape[1]]
            model = temp.reshape(shape[0], shape[1])
            vpi = np.where(model < 0.01)
            model[vpi] = np.nan
            ppp = ppp + shape[0] * shape[1]
            r2 = volcat.where(volcat > 0)
            m_max = np.nanmax(model)
            v_max = np.nanmax(r2.values)
            m_min = np.nanmin(model)
            v_min = np.nanmin(r2.values)
            p_min = np.nanmin([m_min, v_min])
            p_max = np.nanmax([m_max, v_max])
            norm = mpl.colors.Normalize(vmin=p_min, vmax=p_max)
            print(np.nanmax(model), np.nanmax(r2.values))
            if ptype == "pcolormesh":
                cb = ax1.pcolormesh(
                    volcat.longitude,
                    volcat.latitude,
                    model,
                    norm=norm,
                    cmap=cmap,
                    shading="nearest",
                )
                cb2 = ax2.pcolormesh(
                    volcat.longitude,
                    volcat.latitude,
                    r2.values,
                    norm=norm,
                    cmap=cmap,
                    shading="nearest",
                )
            # cb = ax1.scatter(volcat.longitude, volcat.latitude,c=np.log10(model),cmap=cmap,s=50,marker='s')
            else:
                # lon = self.lonlist[tii-1][0:14]
                # lat = self.latlist[tii-1][0:21]
                # xv, yv = np.meshgrid(lon,lat)
                lon = self.lonlist[tii - 1][0:14]
                lat = self.latlist[tii - 1][0:21]
                xv, yv = np.meshgrid(lon, lat)
                print("vshape", volcat.longitude.shape)
                print(xv.shape, yv.shape)
                print("model shape", model.shape)
                cb = ax1.scatter(xv, yv, c=model, cmap=cmap, s=10, marker="o")
                cb2 = ax2.scatter(
                    volcat.longitude,
                    volcat.latitude,
                    c=r2.values,
                    s=10,
                    cmap=cmap,
                    marker="o",
                )
            plt.colorbar(cb, ax=ax1)
            # cb2 = ax2.scatter(volcat.longitude, volcat.latitude, c=r2.values,s=10,cmap=cmap,marker='o')
            plt.colorbar(cb2, ax=ax2)
            if vloc:
                ax1.plot(vloc[0], vloc[1], "y^")
                ax2.plot(vloc[0], vloc[1], "y^")
            plt.show()
        return volcat

    def compare_plotsA(self, daterange, zii=None, tii=None, levels=None):
        """
        must input either daterange or tii.
        if zii is None then sum along ensemble dimension showing coverage of all HYSPLIT runs.
        For the inversion runs, the release from different heights is shown by the ens dimension.

        """

        fig = plt.figure(1, figsize=(10, 5))
        ax1 = fig.add_subplot(1, 1, 1)
        if not tii:
            tii = self.time_index(daterange[0])
        print("tii", tii)
        cdump = self.concmult * self.cdump_hash[tii]
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

    def get_norm(self, model, r2):
        m_max = np.nanmax(model)
        v_max = np.nanmax(r2.values)
        m_min = np.nanmin(model)
        v_min = np.nanmin(r2.values)
        p_min = np.nanmin([m_min, v_min])
        p_max = np.nanmax([m_max, v_max])
        norm = mpl.colors.Normalize(vmin=p_min, vmax=p_max)
        return norm

    def generate_pairs(self):
        for tii in self.volcat_avg_hash.keys():
            volcat = self.volcat_avg_hash[iii]
            cdump = self.cdump_hash[iii] * self.concmult
            yield volcat, cdump

    def set_bias_correction(self, slope, intercept):
        self.slope = slope
        self.intercept = intercept

    def get_pair(self, tii, coarsen=None, slope=None, intercept=None):
        if isinstance(tii, int):
            iii = tii
        elif isinstance(tii, datetime.datetime):
            iii = self.time_index(tii)
        volcat = self.volcat_avg_hash[iii]
        cdump = self.cdump_hash[iii] * self.concmult

        if coarsen:
            volcat = volcat.coarsen(x=coarsen, boundary="trim").mean()
            volcat = volcat.coarsen(y=coarsen, boundary="trim").mean()
            cdump = cdump.coarsen(x=coarsen, boundary="trim").mean()
            cdump = cdump.coarsen(y=coarsen, boundary="trim").mean()
        if not slope and self.slope:
            slope = self.slope
        if slope:
            # print('multiplying cdump by {}'.format(1-slope))
            cdump = cdump * (1 - slope)
        if not intercept and self.intercept:
            intercept = self.intercept
        if intercept:
            # print('shifting cdump by {}'.format(-1*intercept))
            cdump = cdump - intercept
            # if slope is positive then remove negative values.
            f2 = xr.where(cdump > 0, cdump, 0)
            # if slope is negative then do not add to values that were previously 0.
            f2 = xr.where(cdump < 0.0001, 0, f2)
            cdump = f2
        return volcat, cdump

    def match_volcat(self, forecast):
        time = pd.to_datetime(forecast.time.values)
        tii = self.time_index(time)
        volcat = self.volcat_avg_hash[tii]
        return volcat

    def compare_forecast_dist(self, forecast, thresh=None):
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
        # remove nans
        data2 = data[~np.isnan(data)]
        if thresh:
            data2 = data2[data2 > thresh]
            # apply threshold.
            # vpi = data2 < thresh
            # data2[vpi] = np.nan
            # data2 = data2[~np.isnan(data2)]
        return data2

    def set_sampling_time(self, default="start"):
        if "time description" in self.cdump.attrs.keys():
            if "start" in str.lower(self.cdump.attrs["time description"]):
                return "start"
            if "end" in str.lower(self.cdump.attrs["time description"]):
                return "end"
        else:
            return default

    def add_obs(self, vdata):
        """
        vdata is a pandas data frame with mass,lat,lon,time
        """
        self.vdata = vdata

    def prepare_one_time(
        self, daterange, das=None, htoptions=0, st="start", zvals=None
    ):
        return None, None

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
            a1, a2, b1, b2 = self.clip(temp.fillna(0))
            temp = temp[a1:a2, b1:b2]
            temp.plot.pcolormesh(x="longitude", y="latitude")
            plt.show()
        avra = vset.fillna(0).mean(dim="time")
        avra2 = vset.mean(dim="time")
        a1, a2, b1, b2 = self.clip(avra)
        avra = avra[a1:a2, b1:b2]
        avra2 = avra2[a1:a2, b1:b2]
        print("Average with nans set to 0")
        avra.plot.pcolormesh(x="longitude", y="latitude")
        plt.show()
        print("Average of values above 0")
        avra2.plot.pcolormesh(x="longitude", y="latitude")
        plt.show()
        diff = avra2.fillna(0) - avra
        diff.plot.pcolormesh(x="longitude", y="latitude")

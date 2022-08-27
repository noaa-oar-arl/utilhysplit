import sys
import os
import datetime
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import numpy as np
import numpy.ma as ma
import pandas as pd
import seaborn as sns
from scipy.stats import describe
from utilvolc import volcat
from monetio.models import hysplit
from utilhysplit.evaluation import ensemble_tools
from utilhysplit.evaluation import plume_stat
from utilvolc.ash_inverse import InverseAsh
import utilvolc.ash_inverse as ainv
from utilhysplit.evaluation import statmain
from utilhysplit.plotutils import colormaker
from utilhysplit.evaluation import hysplit_boxplots
from utilhysplit.evaluation import cdf_matching


class matchplots:
    """
    create plots using pandas DataFrame output from cdf_match
    """

    def __init__(self, df, fignum=1):
        self.df = df
        self.fignum = 1
        self.set_axis(fignum)

    def set_axis(self, fignum=1):
        sns.set_style("whitegrid")
        sns.set_context("paper")
        self.fig = plt.figure(fignum)
        self.ax = self.fig.add_subplot(1, 1, 1)

    def plot_residual(self):
        df2 = pd.pivot_table(self.df, values="residual", index="time", columns="ens")
        df2.plot(ax=self.ax, legend=False)
        self.ax.set_ylabel("Residual")
        self.fig.autofmt_xdate()

    def plot_intercept(self):
        df2 = pd.pivot_table(self.df, values="intercept", index="time", columns="ens")
        df2.plot(ax=self.ax, legend=False, alpha=0.8)
        self.ax.set_ylabel("Intercept")
        self.fig.autofmt_xdate()
        mean = pd.pivot_table(self.df, values="intercept", index="time", aggfunc="mean")
        mean.plot(ax=self.ax, alpha=0.5, LineWidth=10, color="k", legend=False)

    def plot_one(self, ens):
        ax = self.set_axis(1)
        df2 = self.df[self.df["ens"] == ens]
        df2["intercept"].plot(ax=ax)
        ax2 = self.set_axis(2)
        df2["slope"].plot(ax=ax2)

    def plot_slope(self):
        df2 = pd.pivot_table(self.df, values="slope", index="time", columns="ens")
        df2.plot(ax=self.ax, legend=False)
        self.ax.set_ylabel("Slope")
        self.fig.autofmt_xdate()
        mean = pd.pivot_table(self.df, values="slope", index="time", aggfunc="mean")
        mean.plot(ax=self.ax, alpha=0.5, LineWidth=10, color="k", legend=False)


def cdf_match(aeval, tiilist, enslist, pfit=1, thresh=0.1):
    """
    aeval :  instance of AshEval class
    tiilist : list of times to use
    enslist : list of ensemble indices to use
    pfit : 1

    RETURNS
    dfcdf : pandas DataFrame with columns of slope, intercept, time, ens, pmthresh, residual.
    """

    elist = []
    slopes = []
    slopes2 = []
    intercepts = []
    pmthresh = []
    polyhash = {}
    residuals = []
    tlist = []
    makeplots = False

    for ens in enslist:
        for tii in tiilist:
            aeval.set_bias_correction(slope=None, intercept=None, dfcdf=pd.DataFrame())
            volcat, forecast = aeval.get_pair(tii)
            forecast = forecast.fillna(0)
            thresh2, poly, fc, res = cdf_matching.cdf_match_volcat(
                forecast,
                volcat,
                thresh=thresh,
                ens=ens,
                scale=1,
                pfit=pfit,
                makeplots=makeplots,
            )
            residuals.append(res["residuals"][0])
            tlist.append(aeval.get_time(tii))
            try:
                ensname = forecast.ens.values[ens]
            except:
                ensname = ens
            elist.append(ensname)
            slopes.append(poly[0])
            slopes2.append(poly[-1])
            intercepts.append(poly[1])
            pmthresh.append(float(thresh2))

    polyhash["slope"] = slopes
    polyhash["slope2"] = slopes2
    polyhash["intercept"] = intercepts
    polyhash["time"] = tlist
    polyhash["ens"] = elist
    polyhash["pmthresh"] = pmthresh
    polyhash["residual"] = residuals
    dfcdf = pd.DataFrame.from_dict(polyhash)
    return dfcdf


class AshEval(InverseAsh):
    def __init__(
        self,
        tdir,
        fname,
        vdir,
        vid,
        configdir=None,
        configfile=None,
        verbose=False,
        ensdim="ens",
    ):
        super().__init__(tdir, fname, vdir, vid, configdir, configfile, verbose, ensdim)
        self.evaldf = pd.DataFrame()

    def calc_massload(self):
        self.massload = hysplit.hysp_massload(self.cdump) * self.concmult

    def calc_max_conc(self):
        levels = self.cdump.z.values
        enslist = self.cdump.ens.values

    # def calc_fss(self,iii,thresh=None):
    #    volcat = self.volcat_avg_hash[iii]
    #    ens_forecast = self.cdump_hash[iii]
    #    dft, dft2 = ensemble_tools.ens_time_fss(ens_forecast, volcat,
    #                                            threshold=threshold,
    #                                            neighborhoods = nb,
    #                                            plot=False,
    #                                            pixel_match=pixel_match)

    def conc_cdf_plot(self, timelist=None, enslist=None, threshold=0):
        conc = self.cdump * self.concmult
        self.cdf_plot(conc, timelist, enslist, threshold)

    def volcat_dist_plot(self, threshold=0):
        date = []
        ens = []
        mean = []
        var = []
        skew = []
        kurt = []
        num = []
        small = []
        big = []

        for jjj, iii in enumerate(self.volcat_avg_hash.keys()):

            volc = self.volcat_avg_hash[iii]
            volc = volc.values.flatten()
            volc = [x for x in volc if x > threshold]
            sts = describe(volc)
            print(jjj, iii)
            date.append(jjj)
            ens.append("obs")
            mean.append(sts.mean)
            var.append(sts.variance)
            skew.append(sts.skewness)
            kurt.append(sts.kurtosis)
            num.append(sts.nobs)
            small.append(sts.minmax[0])
            big.append(sts.minmax[1])
        data = zip(date, ens, mean, var, skew, kurt, num, small, big)
        colnames = [
            "date",
            "ens",
            "mean",
            "variance",
            "skewness",
            "kurtosis",
            "N",
            "min",
            "max",
        ]
        dfout = pd.DataFrame(data)
        dfout.columns = colnames
        self.dfstats = dfout
        return dfout

    def model_cdf_plot(self, ens=0, threshold=0, cii=None):
        step = 5
        clr = ["-m", "-r", "-b", "-c", "-g"]

        for jjj, iii in enumerate(self.cdump_hash.keys()):
            # print(self.cdump.time.values[iii])
            cdump = self.manage_cdump(iii, cii=cii)

            if isinstance(ens, str):
                cdump = cdump.sel(ens=ens)
            elif isinstance(ens, int):
                cdump = cdump.isel(ens=ens)
            # cdump = cdump.isel(ens=ens)
            # print(self.cdump.time.values[iii])
            # try:
            #    volcat = self.cdump_hash[iii].sel(ens=ens) * self.concmult
            # except:
            #    volcat = self.cdump_hash[iii].isel(ens=ens) * self.concmult
            cdump = cdump.values.flatten()
            cdump = [x for x in cdump if x > threshold]
            try:
                sdata, y = statmain.cdf(cdump)
            except:
                print("cannot calculate cdf for {}".format(iii))
                continue
            ax = plt.gca()
            if jjj % 5 == 0:
                lw = 3
            else:
                lw = 1
            ax.step(sdata, y, clr[jjj % len(clr)], LineWidth=lw)
            ax.set_xscale("log")
            ax.set_xlabel("Mass Loading")

    def volcat_boxplot(self, threshold):
        vdata = []
        datelist = []
        for jjj, iii in enumerate(self.volcat_avg_hash.keys()):
            # print(self.cdump.time.values[iii])
            volcat = self.volcat_avg_hash[iii]
            volcat = volcat.values.flatten()
            volcat = [x for x in volcat if x > threshold]
            datelist.append(self.get_time(iii))
            vdata.append(volcat)
        # dj is a pandas dataframe
        dj = hysplit_boxplots.prepare_boxplotdata(datelist, vdata)
        hysplit_boxplots.make_boxplot(dj)

    def volcat_cdf_plot(self, threshold=0):
        step = 5
        clr = ["-m", "-r", "-b", "-c", "-g"]
        for jjj, iii in enumerate(self.volcat_avg_hash.keys()):
            # print(self.cdump.time.values[iii])
            volcat = self.volcat_avg_hash[iii]
            volcat = volcat.values.flatten()
            volcat = [x for x in volcat if x > threshold]
            try:
                sdata, y = statmain.cdf(volcat)
            except:
                print("cannot calculate cdf for {}".format(iii))
                continue
            ax = plt.gca()
            if jjj % 5 == 0:
                lw = 3
                # print('here')
            else:
                lw = 1
            ax.step(sdata, y, clr[jjj % len(clr)], LineWidth=lw)

    def pixel_matching(time, enslist=None, threshold=0):
        mass = self.massload.sel(time=time)
        volcat = self.match_volcat(forecast)

    def mass_pdf_plot(
        self,
        timelist=None,
        enslist=None,
        threshold=0,
        use_pixel_match=False,
        plotdiff=False,
        figname=None,
        xscale="log",
    ):
        mass = self.massload
        ilist = [self.time_index(time) for time in timelist]
        kshash = {}
        for ttt, iii in enumerate(ilist):
            if iii not in self.volcat_avg_hash.keys():
                continue

            # get and process volcat data
            volcat = self.volcat_avg_hash[iii]
            volcat = volcat.values.flatten()
            volcat = [x for x in volcat if x > threshold]
        sns.distplot(mass.values.flatten())

    def mass_boxplot(self, timelist=None, enslist=None, threshold=0, clist=None):
        mass = self.massload
        ensemble_tools.ens_boxplot(mass, threshold=threshold, clist=clist)

    def mass_cdf_plot(
        self,
        timelist=None,
        enslist=None,
        threshold=0,
        use_pixel_match=False,
        plotdiff=False,
        figname=None,
        xscale="log",
        colors=None,
        plotvolcat=True,
    ):
        # mass = self.massload
        ilist = [self.time_index(time) for time in timelist]
        kshash = {}
        cdfhash = {}
        print("running mass_cdf_plot", ilist, timelist)
        for ttt, iii in enumerate(ilist):
            volcat, mass = self.get_pair(iii)
            if plotvolcat:
                if iii not in self.volcat_avg_hash.keys():
                    continue

                # get and process volcat data
                volcat = self.volcat_avg_hash[iii]
                volcat = volcat.values.flatten()
                volcat = [x for x in volcat if x > threshold]
                # CDF of volcat data
                sdata, y = statmain.cdf(volcat)

            # CDF of model data
            if not use_pixel_match:
                cdhash = self.cdf_plot(
                    mass,
                    [timelist[ttt]],
                    enslist,
                    threshold,
                    pixel_match=None,
                    xscale=xscale,
                    colors=colors,
                )
            else:
                cdhash = self.cdf_plot(
                    mass,
                    [timelist[ttt]],
                    enslist,
                    threshold,
                    pixel_match=len(sdata),
                    xscale=xscale,
                    colors=colors,
                )
            # add obs to cdhash
            cdhash[(timelist[ttt], "obs")] = (sdata, y)
            cdfhash.update(cdhash)
            # plot the volcat data
            ax = plt.gca()
            ax.step(sdata, y, "--k", LineWidth=2)
            if not plotdiff:
                continue
            # compute the ks value.
            for key in cdhash.keys():
                cdf = cdhash[key]
                try:
                    kstest = statmain.kstest_sub(sdata, y, cdf[0], cdf[1])
                except:
                    continue
                maxval = np.max(kstest[1])
                minval = np.min(kstest[1])
                if np.abs(minval) > np.abs(maxval):
                    kshash[key] = [minval]
                else:
                    kshash[key] = [maxval]
                if plotdiff:
                    plt.plot(kstest[0], kstest[1])
            if figname:
                print("saving", figname)
                plt.savefig(figname)
            plt.show()

        # now plot ks score as function of time or ensemble member.
        if plotdiff:
            fig = plt.figure(1)
            self.ksdf = self.kshash_to_df(kshash)
            if isinstance(timelist, np.ndarray) or isinstance(timelist, list):
                if len(timelist) > 1:
                    # plot as function of time
                    self.ksdf.plot()
                else:
                    # plot as function of ense member
                    temp = self.ksdf.T
                temp.plot()
        # 1return self.ksdf
        return cdfhash

    def kshash_to_df(self, kshash):
        # index should be date. columns should be ensemble members. values should be kstest output.
        ksdf = pd.DataFrame.from_dict(kshash)
        ksdf = ksdf.T
        ksdf = ksdf.reset_index()
        ksdf.columns = ["date", "ens", "ks"]
        ksdf = ksdf.pivot(values="ks", columns="ens", index="date")
        return ksdf

    def unpack_kshash(self, kshash):
        temp = list(kshash.items())
        test2 = list(zip(*temp))
        values = test2[1]
        temp = list(zip(*test2[0]))
        dates = temp[0]
        ens = temp[1]
        return dates, ens, values

    def cdf_plot(
        self,
        dra,
        timelist=None,
        enslist=None,
        threshold=0,
        pixel_match=None,
        xscale="log",
        colors=None,
        label="time",
    ):
        sns.set_style("whitegrid")
        sourcelist = None
        cdhash = ensemble_tools.ens_cdf(
            dra,
            enslist=enslist,
            sourcelist=sourcelist,
            timelist=timelist,
            threshold=threshold,
            pixel_match=pixel_match,
            plot=True,
            xscale=xscale,
            colors=colors,
            label=label,
        )
        print("RETURNING", type(cdhash))
        return cdhash

    def ks_time_series(self, thresh=0, drange=None):
        # THIS IS OUTDATED.
        dra = self.massload
        kshash = {}
        forecast_area = {}
        volcat_area = {}
        datelist = []
        # each run gets it's own time series.
        for hrun in dra.ens.values:
            kshash[hrun] = {}
            forecast_area[hrun] = []
            volcat_area[hrun] = []
        for sdate in dra.time.values:
            print(type(sdate), type(drange[0]))
            if drange:
                if pd.to_datetime(sdate) < drange[0]:
                    continue
                if pd.to_datetime(sdate) > drange[1]:
                    break
            temp = dra.sel(time=sdate)
            datelist.append(sdate)
            for hrun in dra.ens.values:
                print("working on {} {}".format(sdate, hrun))
                temp2 = temp.sel(ens=hrun)
                ks = self.compare_forecast_dist(temp2, thresh)
                print("got ks")
                kshash[hrun][sdate] = ks
                a1, a2 = self.compare(temp2, thresh)
                print("got area")
                volcat_area[hrun].append(a2)
                forecast_area[hrun].append(a2)
        self.kshash = kshash
        self.forecast_area = forecast_area
        self.datelist = datelist
        self.volcat_area = volcat_area

    def compare_forecast_dist(self, forecast, thresh=None, ax=None):
        from utilhysplit.evaluation import statmain

        volcat = self.match_volcat(forecast)
        print("got volcat")
        evals = forecast.values
        vvals = volcat.values
        # exval,eyval =  statmain.nancdf(evals.flatten(),thresh)
        # vxval,vyval =  statmain.nancdf(vvals.flatten(),thresh)
        # print('got cdf')
        # if not ax:
        #    fig = plt.figure(1)
        #    ax = fig.add_subplot(1,1,1)
        # ax.step(exval, eyval, '-r',label='HYSPLIT')
        # ax.step(vxval, vyval, '-b',label='Volcat')
        ks1 = statmain.kstestnan(evals.flatten(), vvals.flatten(), thresh)
        try:
            print("Kolmogorov-Smirnov Parameter {}".format(np.max(np.abs(ks1))))
        except:
            return 1
        return np.max(np.abs(ks1))

    def compare(self, forecast, thresh=None):
        volcat = self.match_volcat(forecast)
        evals = self.remove_nans(forecast.values.flatten(), thresh)
        vvals = self.remove_nans(volcat.values.flatten(), thresh)
        # number of pixles above threshold in each.
        forecast_parea = len(evals)
        volcat_parea = len(vvals)
        return forecast_parea, volcat_parea

    def compare_forecast(
        self,
        forecast,
        cmap="viridis",
        ptype="pcolormesh",
        vloc=None,
        include="all",
        thresh=0.00001,
        prob=False,
        cmap_prob="cool",
        overlay=False,
        **kwargs
    ):
        """
        possible kwargs:
        xlim
        ylim
        vmin
        vmax
        log
        cmap3
        """
        # forecast should be an xarray in mass loading format with no time dimension.
        if "log" in kwargs.keys():
            logscale = kwargs["log"]
        else:
            logscale = False
        sns.set()
        sns.set_style("whitegrid")

        if include == "all":
            fig = plt.figure(figsize=[15, 4])
            ax1 = fig.add_subplot(1, 3, 1)
            ax2 = fig.add_subplot(1, 3, 2)
            ax3 = fig.add_subplot(1, 3, 3)
        else:
            fig = plt.figure(figsize=[10, 3])
            ax2 = fig.add_subplot(1, 2, 1)
            ax3 = fig.add_subplot(1, 2, 2)

        time = pd.to_datetime(forecast.time.values)
        tii = self.time_index(time)
        volcat = self.volcat_avg_hash[tii]
        evals = forecast.values.copy()
        vpi = evals < thresh
        evals[vpi] = np.nan

        if logscale:
            evals = np.log10(evals)
            vvals = np.log10(volcat.values.copy())
        else:
            vvals = volcat.values.copy()

        if "vmin" in kwargs.keys() and "vmax" in kwargs.keys():
            if logscale:
               norm = mpl.colors.LogNorm(vmin=kwargs["vmin"], vmax=kwargs["vmax"])
            else: 
               norm = mpl.colors.Normalize(vmin=kwargs["vmin"], vmax=kwargs["vmax"])
            vmax = kwargs["vmax"]
        else:
            norm = self.get_norm(vvals, evals)
            if not logscale:
               vmax = np.max(evals) + 10
            else:
               vmax = np.max(evals)

        vpi = vvals < 0.001
        vvals[vpi] = np.nan
        if "clevels" in kwargs.keys():
            clevels = kwargs["clevels"]
        else:
            clevels = [0.02, 0.2, 2, 5, 10, 50]

        dlevels = [0.1, 2.0, 10.0]
        if "cmap3" in kwargs.keys():
            print("HERE HERE HERE", kwargs["cmap3"])
            cmap3a = kwargs["cmap3"][0]
            cmap3b = kwargs["cmap3"][1]

            cm = colormaker.ColorMaker(cmap3a, len(clevels), ctype="rgb")
            cm3 = colormaker.ColorMaker(cmap3b, len(clevels), ctype="rgb")
        else:
            cm = colormaker.ColorMaker(cmap, len(clevels), ctype="rgb")
            cm3 = cm
            cmap3a = cmap

        if prob:
            cm2 = colormaker.ColorMaker(cmap_prob, len(clevels), ctype="rgb")
        else:
            cm2 = colormaker.ColorMaker(cmap, len(clevels), ctype="rgb")

        if ptype == "pcolormesh":
            if include == "all":
                if logscale:
                    cb = ax1.pcolormesh(
                        volcat.longitude,
                        volcat.latitude,
                        vvals,
                        cmap=cmap,
                        shading="nearest",
                     )
                else:
                    cb = ax1.pcolormesh(
                        volcat.longitude,
                        volcat.latitude,
                        vvals,
                        norm=norm,
                        cmap=cmap,
                        shading="nearest",
                     )
            # cb2 = ax2.pcolormesh(
            #    forecast.longitude,
            #    forecast.latitude,
            #    evals,
            # norm=norm,
            #    cmap=cmap,
            #    shading="nearest",
            # )
            if "plevels" in kwargs.keys():
                plevels = kwargs["plevels"]
            else:
                plevels = clevels
            if ptype == "contour":
                cb2 = ax2.contourf(
                    forecast.longitude,
                    forecast.latitude,
                    evals,
                    norm=norm,
                    levels=plevels,
                    colors=cm2()
                )
            else:
                print('here')
                if logscale:
                    cb2 = ax2.pcolormesh(
                        forecast.longitude,
                        forecast.latitude,
                        evals,
                        cmap=cmap,
                        shading="nearest",
                        )
                else:
                    cb2 = ax2.pcolormesh(
                        forecast.longitude,
                        forecast.latitude,
                        evals,
                        norm=norm,
                        cmap=cmap,
                        shading="nearest",
                        ) 
                if not logscale:
                   vals = np.log10(evals)
                else:
                   vals = evals

                cb3 = ax3.contourf(
                    forecast.longitude,
                    forecast.latitude,
                    vals,
                    levels=np.log10(clevels),
                    # colors=cm(),
                    linewidths=5,
                    cmap=cmap3a,
                    alpha=1
                )
            if overlay:
                if "plevels" in kwargs.keys():
                    plevels = kwargs["plevels"]
                if not logscale:
                   vals = np.log10(vvals)
                else:
                   vals = vvals
                cb4 = ax3.contourf(
                    volcat.longitude,
                    volcat.latitude,
                    vals,
                    levels=np.log10(plevels),
                   # colors = cm3(),
                    # linewidths=3,
                    cmap=cmap3b,
                    alpha=1,
                )
                #cb4 = ax3.contour(
                #    volcat.longitude,
                #    volcat.latitude,
                #    vals,
                #    levels=np.log10(plevels),
                #    # colors = cm3(),
                #    linewidths=3,
                #    cmap=cmap3b,
                #    alpha=1,
                #)
        plt.title(time.strftime("%Y %m/%d %H:%M UTC"))
        if include == "all":
            plt.colorbar(cb, ax=ax1)
        if np.max(evals) > vmax:
            plt.colorbar(cb2, ax=ax2, extend="max")
        else:
            plt.colorbar(cb2, ax=ax2)
        cb33 = plt.colorbar(cb3, ax=ax3, shrink=1.00)
        #cb33.ax.tick_params(size=10, width=10)
        ticklabels = [str(x) for x in clevels]
        #cb33.ax.set_yticklabels(['0.01','0.1','1','10'])
        cb33.ax.set_yticklabels(ticklabels)
        if overlay:
            cb44 = plt.colorbar(cb4, ax=ax3, shrink=1.00)
            cb44.ax.set_yticklabels(['','','',''])
        # sets the linewidths of contour lines in the plot.
        # plt.setp(cb4.collections,linewidth=1)
        # cb44=plt.colorbar(cb4, ax=ax3)
        # cb44.ax.tick_params(labelsize=8)
        # sets linewidths in the colorbar.
        # cb44.ax.get_children()[0].set_linewidths(10)
        # if not 'plevels' in kwargs.keys(): cb44.set_ticks([])
        if "xlim" in kwargs.keys():
            if include == "all":
                ax1.set_xlim(kwargs["xlim"])
            ax2.set_xlim(kwargs["xlim"])
            ax3.set_xlim(kwargs["xlim"])
        else:
            if include == "all":
                xlim = ax1.get_xlim()
            ax2.set_xlim(xlim)
            ax3.set_xlim(xlim)
        if "ylim" in kwargs.keys():
            if include == "all":
                ax1.set_ylim(kwargs["ylim"])
            ax2.set_ylim(kwargs["ylim"])
            ax3.set_ylim(kwargs["ylim"])
        else:
            if include == "all":
                ylim = ax1.get_ylim()
            ax2.set_ylim(ylim)
            ax3.set_ylim(ylim)
        if vloc:
            ms=10
            if include == "all":
                ax1.plot(vloc[0], vloc[1], "k^", markersize=ms)
            ax2.plot(vloc[0], vloc[1], "k^",markersize=ms)
            ax3.plot(vloc[0], vloc[1], "k^",markersize=ms)
        return fig

    def compare_plots(self, daterange, levels=None):
        fig = plt.figure(1, figsize=(10, 5))
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)
        tii = self.time_index(daterange[0])
        print("tii", tii)
        cdump = self.cdump_hash[tii]
        volcat = self.volcat_avg_hash[tii]
        csum = cdump.sum(dim="ens")
        volcat.plot.pcolormesh(x="longitude", y="latitude", levels=levels, ax=ax1)
        # cdump.sum(dim='ens').plot.contour(x='longitude',y='latitude',ax=ax2)
        # plt.pcolormesh(csum.longitude, csum.latitude, np.log10(csum),cmap='Reds')
        # plt.pcolormesh(csum.longitude, csum.latitude, csum,cmap='Reds',shading='nearest')
        cb = csum.plot.pcolormesh(x="longitude", y="latitude", cmap="viridis", ax=ax2)
        # cb = plt.pcolormesh(volcat.longitude, volcat.latitude, np.log10(volcat),cmap='Blues',levels=levels)
        # cb = plt.scatter(volcat.longitude, volcat.latitude, c=np.log10(volcat),s=2,cmap='Blues')
        # cb = plt.scatter(volcat.longitude, volcat.latitude, c=volcat.values,s=2,cmap='viridis',levels=levels)
        # cb = plt.contour(volcat.longitude, volcat.latitude, np.log10(volcat),cmap='Blues')
        # plt.colorbar(cb)
        plt.tight_layout()
        return ax1, ax2

    # def get_cdump(self,tii,dfcdf=pd.DataFrame(),cii=None,coarsen=None):
    #    print('here')
    #    if isinstance(tii,int):
    #       iii = tii
    #    elif isinstance(tii,datetime.datetime):
    ##       iii = self.time_index(tii)


#
#        cdump = self.cdump_hash[iii]*self.concmult
#        if not cii: cii = tii
#        if not dfcdf.empty:
#        # create xarray with slope and intercept with coordinate of ens.
#            temp = dfcdf[dfcdf['time']==cii]
#            temp = temp[['slope','intercept','ens']]
#            temp = temp.set_index('ens')
#            xrt = temp.to_xarray()
#            # apply slope and intercept correction
#            cdump = cdump * (1-xrt.slope)
#            cdump = cdump - xrt.intercept
#            f2 = xr.where(cdump>0,cdump,0)
#            f2 = xr.where(cdump<0.0001,0,f2)
#            cdump = f2
#
#        if coarsen:
#           cdump = cdump.coarsen(x=coarsen,boundary='trim').mean()
#           cdump = cdump.coarsen(y=coarsen,boundary='trim').mean()
#
#        return  xrt, cdump


def calc_stats_function(hxr, ashmass, threshold):
    from utilhysplit.evaluation import plume_stat
    from utilhysplit.evaluation import ensemble_tools

    """
        Currently this is just for one time period. 
        """
    if isinstance(threshold, float):
        threshold = [threshold]

    # stack ensemble and source dimensions
    dim = "ens"
    if "source" in hxr.dims and "ens" in hxr.dims:
        hxr = hxr.stack(ens=("ens", "source"))
    elif "source" in hxr.dims:
        dim = "source"

    # want to evaluate ensemble mean.
    hxr_mean = hxr.mean(dim=dim)

    # Creating dummy xarray for merging below
    data = np.zeros(len(hxr.ens.values))
    # statsxr = xr.DataArray(name='dummy', data=data, attrs=attrs,
    statsxr = xr.DataArray(name="dummy", data=data, dims=dim, coords=[hxr.ens.values])
    BSlist = []
    fss_mean_list = []
    for thresh in threshold:
        # Converting to probabilistic (0-1) field for model data.
        atl = ensemble_tools.ATL(hxr, norm=True, thresh=thresh)
        # Converting to VOLCAT binary field for BS calculation
        ashmass_binary = xr.where(ashmass >= thresh, 1.0, 0.0)

        # Calculating Brier Score of each ensemble member
        a = 0
        PClistcent = []
        PClistuncent = []
        PClistcentavg = []
        PClistuncentavg = []

        BS = plume_stat.calc_bs(atl, ashmass_binary)
        BSlist.append(BS.values)
        stats = plume_stat.CalcScores(
            hxr_mean,
            ashmass_binary,
            threshold=thresh,
            verbose=False,
            makeplots=False,
            szra=[3, 5, 7, 9, 11],
        )

        # fssdf = stats.calc_fss(

        # for member in hxr.ens.values:
        #    print(member)
        # Calculating pattern correlation coefficients, centered and uncentered
        #    stats = plume_stat.CalcScores(ashmass,
        #                          hxr.sel(ens=member), threshold=thresh, verbose=False)
        #    PCcent, PCuncent = stats.calc_pcorr()
        #    PClistcent.append(PCcent.values)
        #    PClistuncent.append(PCuncent.values)

        # Creating binary field of hysplit output
        # hxr_binary = xr.where(hxr.isel(ens=member) >= thresh, 1., 0.)
        # Calculating the BS values of each ensemble member
        #    BS = ps.calc_bs(atl, ashmass_binary)
        #    BSlist.append(BS.values)
        # Adding Brier Scores to the netcdf, with dimension source
        # if thresh == 0.1:
        #    thresh = '0p1'
        # else:
        #    thresh = str(thresh)

        # threshstr = str(thresh)+' g/m^2'
        # BSxr = xr.DataArray(BSlist, dims='source').load().rename('BS'+thresh)
        # BSxr.attrs['long name'] = 'Brier Score compared to volcat'
        # BSxr.attrs['threshold'] = threshstr
        # BSavgxr = xr.DataArray(BSlistavg, dims='source').load().rename('BSavg'+thresh)
        # BSavgxr.attrs['long name'] = 'Brier Score compared to 1hr avg volcat'
        # BSavgxr.attrs['threshold'] = threshstr
        # PCxr = xr.DataArray(PClistcent, dims='source').load().rename('PC'+thresh)
        # PCxr.attrs['long name'] = 'Pattern Correlation (centered) compared to volcat'
        # PCxr.attrs['threshold'] = threshstr
        # PCxruc = xr.DataArray(PClistuncent, dims='source').load().rename('PCuc'+thresh)
        # PCxruc.attrs['long name'] = 'Pattern Correlation (uncentered) compared to volcat'
        # PCxruc.attrs['threshold'] = threshstr
        # PCavgxr = xr.DataArray(PClistcent, dims='source').load().rename('PC'+thresh)
        # PCavgxr.attrs['long name'] = 'Pattern Correlation (centered) compared to 1hr avg volcat'
        # PCavgxr.attrs['threshold'] = threshstr
        # PCavgxruc = xr.DataArray(PClistuncent, dims='source').load().rename('PCuc'+thresh)
        # PCavgxruc.attrs['long name'] = 'Pattern Correlation (uncentered) compared to 1hr avg volcat'
        # PCavgxruc.attrs['threshold'] = threshstr

        # statsxr = xr.merge([statsxr, BSxr, BSavgxr, PCxr, PCxruc, PCavgxr,
        #                    PCavgxruc], combine_attrs='drop_conflicts')
    # Dropping dummy variable
    # statsxr = statsxr.drop(labels='dummy')
    return BSlist

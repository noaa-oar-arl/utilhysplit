import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#import seaborn as sns
#import xarray as xr
from utilhysplit.evaluation import  plume_stat
from utilhysplit.evaluation import ensemble_tools

"""
FUNCTIONS
The functions in this file are for ensemble fractional skill score calculations

ens_time_fss  
ens_fss       
plot_ens_fss_ts    
plot_ens_fss
plot_ens_afss
plot_ens_afss_ts


"""

# _______________________________________________________________________

# 2023 Jan 10 amc moved functions in ensemble_tools which depend on plume_stat here.


# calls ens_fss
def ens_time_fss(
    indralist,
    obsralist,
    enslist=None,
    sourcelist=None,
    # timelist=None,
    neighborhoods=[1, 3, 5, 7],
    threshold=0,
    plot=True,
    pixel_match=False,
):
    """
    RETURNS
    pandas dataframe with columns
    Nlen, FBS, FBS_ref, FSS, ens, time

    pandas dataframe with columns
    MSE, MAE, threshold, exclude_zeros, N, ens
    """
    dflist = []
    df2list = []
    for pairs in zip(indralist, obsralist):
        dff, df2 = ens_fss(
            pairs[0],
            pairs[1],
            enslist,
            sourcelist,
            neighborhoods,
            threshold,
            plot,
            return_objects=False,
            pixel_match=pixel_match,
        )
        dflist.append(dff)
        df2list.append(df2)
    return pd.concat(dflist), pd.concat(df2list)


# utilizes plume_stat
def ens_fss(
    indra,
    obsra,
    enslist=None,
    sourcelist=None,
    # timelist=None,
    neighborhoods=[1, 3, 5, 7],
    threshold=0,
    plot=True,
    return_objects=False,
    pixel_match=False,
):
    """
    indra and obsra need to be same time period and grid.

    RETURNS
    dfall: pandas dataframe with columns
    Nlen, FBS, FBS_ref, FSS, ens, time

    dfall2: pandas dataframe with columns
    MSE, MAE, threshold, exclude_zeros, N, ens

    dfall3: pandas dataframe with columns
    contingency table./dft2
    """
    dra, dim = ensemble_tools.preprocess(indra, enslist, sourcelist)
    dflist = []
    df2list = []
    df3list = []
    iii = 0
    # calculate fss for each ensemble member.
    for ens in dra[dim].values:
        if dim == "ens":
            subdra = dra.sel(ens=ens)
        elif dim == "source":
            subdra = dra.sel(source=ens)
        scores = plume_stat.CalcScores(
            obsra, subdra, threshold=threshold, pixel_match=pixel_match
        )
        df1 = scores.calc_fss(makeplots=False, szra=neighborhoods)
        df2 = scores.calc_accuracy_measures(threshold=0)
        df3 = scores.table2csi(scores.get_contingency_table())
        df1["ens"] = ens
        df2["ens"] = ens
        df3["ens"] = ens
        dflist.append(df1)
        df2list.append(df2)
        df3list.append(df3)
        iii += 1

    # calculate fss for ensemble mean
    meanra = dra.mean(dim=dim)
    mean_scores = plume_stat.CalcScores(
        obsra, meanra, threshold=threshold, pixel_match=pixel_match
    )
    if plot:
        mean_scores.binxra2.plot.pcolormesh()
        # print('Mean sum', mean_scores.binxra2.sum())
        plt.show()
    df1 = mean_scores.calc_fss(makeplots=False, szra=neighborhoods)
    df2 = mean_scores.calc_accuracy_measures(threshold=0)
    df3 = mean_scores.table2csi(mean_scores.get_contingency_table())
    df1["ens"] = "mean"
    df2["ens"] = "mean"
    df3["ens"] = "mean"
    dflist.append(df1)
    df2list.append(df2)
    df3list.append(df3)

    # calculate fss for probabilistic output
    prob_scores = plume_stat.CalcScores(
        obsra, dra, threshold=threshold, probabilistic=True, pixel_match=pixel_match
    )
    if plot:
        prob_scores.binxra2.plot.pcolormesh()
        plt.show()
    # print('Prob sum', prob_scores.binxra2.sum())
    df1 = prob_scores.calc_fss(makeplots=False, szra=neighborhoods)
    df1["ens"] = "prob"
    dflist.append(df1)

    # add time to dataframe.
    dfall = pd.concat(dflist)
    dfall2 = pd.concat(df2list)
    dfall3 = pd.concat(df3list)
    if "time" in indra.coords:
        dfall["time"] = pd.to_datetime(indra.coords["time"].values)
        dfall2["time"] = pd.to_datetime(indra.coords["time"].values)
        dfall3["time"] = pd.to_datetime(indra.coords["time"].values)
    if return_objects:
        return mean_scores, prob_scores, dfall, dfall2
    dfall4 = dfall3.merge(dfall2, how="outer", on=["time", "ens"])
    return dfall, dfall4


# calls ensfss
def plot_ens_fss_ts(ensdf, nval=5, clrs=None):
    """
    Plot FSS on y and time on x.
    ensdf : pandas DataFrame output from  ens_time_fss
    nval : neighborhood size to plot.
    """
    if nval:
        if nval not in ensdf["Nlen"].unique():
            pvals = np.array(ensdf["Nlen"].unique())
            idx = np.abs(pvals - nval).argmin()
            nval = pvals[idx]
            print("nval not in possible values. changing to {}".format(nval))
        tempdf = ensdf[ensdf["Nlen"] == nval]

    fig, ax = plt.subplots(1, 1)
    ensfss = tempdf.pivot(columns="ens", values="FSS", index="time")
    uniform = tempdf.pivot(columns="ens", values="uniform", index="time")
    if not clrs:
        ensfss.plot(ax=ax, legend=None, colormap="tab20")
    else:
        ensfss.plot(ax=ax, legend=None, color=clrs)
    colA = uniform.columns[0]
    uniform.plot(ax=ax, y=colA, linestyle="--", legend=None, colormap="winter")
    if "mean" in ensfss.columns:
        ensfss.plot(
            ax=ax, y="mean", linewidth=5, colormap="winter", legend=None, label="mean"
        )
    if "prob" in ensfss.columns:
        ensfss.plot(
            ax=ax,
            y="prob",
            linewidth=3,
            colormap="gist_gray",
            legend=None,
            label="prob",
        )
    ax.set_ylabel("FSS")


def plot_afss_ts(ensdf, clrs=None):
    fig, ax = plt.subplots(1, 1)
    afss = ensdf[["ens", "afss", "time"]]
    afss = afss.drop_duplicates()
    afss = afss.pivot(index="time", columns="ens", values="afss")
    if not clrs:
        afss.plot(ax=ax, legend=None)
    else:
        afss.plot(ax=ax, legend=None, color=clrs)
    if "mean" in afss.columns:
        afss.plot(ax=ax, y="mean", linewidth=5, colormap="winter", legend=None)
    if "prob" in afss.columns:
        afss.plot(ax=ax, y="prob", linewidth=3, colormap="gist_gray", legend=None)
    ax.set_ylabel("AFSS")
    return afss


def plot_afss(ensdf):
    fig, ax = plt.subplots(1, 1)
    afss = ensdf[["ens", "afss"]]
    afss = afss.drop_duplicates()
    afss = afss.set_index("ens")
    afss.plot(ax=ax, linestyle="", Marker="o", legend=None)
    ax.set_ylabel("AFSS")


def plot_ens_fss(
    ensdfin,
    sizemult=1,
    xmult=1,
    #timelist=None,
    #enslist=None,
    plotafss=False,
    clrs=None,
    ax=None,
):
    """
    Plot FSS on y and Neighborhood length on x.

    ensdf : pandas DataFrame output from ens_fss function
    sizemult : set to value other than one to convert to degrees.
    """
    if not ax:
        # if not isinstance(ax, matplotlib.axes._subplots.AxesSubplot):
        fig, ax = plt.subplots(1, 1)
    ensdf = ensdfin.copy()
    ensdf["Nlen"] = ensdf["Nlen"] * xmult
    if sizemult != 1:
        ensdf["length (degrees)"] = ensdf["Nlen"] * sizemult
        ensfss = ensdf.pivot(columns="ens", values="FSS", index="length (degrees)")
    else:
        ensfss = ensdf.pivot(columns="ens", values="FSS", index="Nlen")
    if not clrs:
        ensfss.plot(ax=ax, legend=None, colormap="spring")
    else:
        ensfss.plot(ax=ax, legend=None, color=clrs)
    random = ensdf["random"].unique()
    nmin = float(np.min(ensdf["Nlen"])) * sizemult
    nmax = float(np.max(ensdf["Nlen"])) * sizemult
    if "mean" in ensfss.columns:
        ensfss.plot(ax=ax, y="mean", linewidth=5, colormap="winter", legend=None)
    if "prob" in ensfss.columns:
        ensfss.plot(ax=ax, y="prob", linewidth=3, colormap="gist_gray", legend=None)
    # plot random forecast
    for randomval in random:
        ax.plot([nmin, nmax], [randomval, randomval], "--k")
    uniform = ensdf["uniform"].unique()
    # plot uniform forecast
    for uniformval in uniform:
        ax.plot([nmin, nmax], [uniformval, uniformval], "--r")
    if plotafss:
        afss = ensdf["afss"].unique()
        for afssval in afss:
            ax.plot([nmin, nmax], [afssval, afssval], "--k", alpha=0.5)

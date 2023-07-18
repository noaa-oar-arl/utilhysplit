import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import BoundaryNorm
from utilhysplit.plotutils import vtools
import utilhysplit.evaluation.web_ensemble_plots as wep

from monetio.models import hysplit



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


def plotcdump(
    cdump,
    time_index=0,
    cmap="viridis",
    ptype="pcolormesh",
    vloc=None,
    include="all",
    thresh=0.00001,
    prob=False,
    central_longitude=0,
    logscale=True,
    **kwargs
):
    """ """
    figsize = [10, 10]
    transform = cartopy.crs.PlateCarree(central_longitude=central_longitude)
    vtransform = cartopy.crs.PlateCarree(central_longitude=0)
    sns.set()
    cdump = cdump.isel(time=time_index)
    forecast = hysplit.hysp_massload(cdump)

    ncols = 1
    fig, ax1 = plt.subplots(
        nrows=1,
        ncols=ncols,
        figsize=figsize,
        constrained_layout=True,
        subplot_kw={"projection": transform},
    )

    time = pd.to_datetime(forecast.time.values)
    evals = forecast.values.copy()
    vpi = evals < thresh
    evals[vpi] = np.nan

    if logscale:
        evals = np.log10(evals)

    if "vmin" in kwargs.keys() and "vmax" in kwargs.keys():
        if logscale:
            norm = mpl.colors.LogNorm(vmin=kwargs["vmin"], vmax=kwargs["vmax"])
        else:
            norm = mpl.colors.Normalize(vmin=kwargs["vmin"], vmax=kwargs["vmax"])
        vmax = kwargs["vmax"]
    else:
        norm = get_norm(evals, evals)
        if not logscale:
            vmax = np.max(evals) + 10
        else:
            vmax = np.max(evals)

    cb2 = ax1.pcolormesh(
        forecast.longitude,
        forecast.latitude,
        evals,
        norm=norm,
        cmap=cmap,
        shading="nearest",
    )
    plt.colorbar(cb2, ax=ax1, extend="max")
    # )
    plt.title(time.strftime("%Y %m/%d %H:%M UTC"))
    wep.format_plot(ax1, transform)
    return fig

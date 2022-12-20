import datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from utilhysplit.fixlondf import fixlondf


def frequency_plots_all(df,sdate,dtt,dres=1,vloc=None):
    """
    df : pandas aataframe
    sdate : datetime object
    dtt :
    """
    times = df.time.unique()
    times.sort()
    done = False
    if isinstance(dtt,int):
       dtt = datetime.timedelta(hours=dtt)
    while not done:
          onetime = df[df['time']==sdate]
          onetime = fixlondf(onetime, colname='longitude', neg=False)
          frequency_plot(onetime,dres,vloc)
          plt.title(sdate)
          sdate += dtt
          if sdate > pd.to_datetime(times[-1]): 
             print('done', sdate)
             print(times)
             done = True
          plt.show()

def frequency_plot(df, dres = 1, vloc=None):
    """
    df : pandas dataframe with columns of longitude, latitude.
         weight column is optional. If it does not exist then all points are 
         evenly weighted.
    dres : resolution of bins in degrees.
    """
    sns.set()
    if 'weight' not in df.columns:
        df['weight'] = 1
    dres = dres
    xmin = np.floor(np.min(df.longitude.values))
    xmax = np.ceil(np.max(df.longitude.values))
    ymin = np.floor(np.min(df.latitude.values))
    ymax = np.ceil(np.max(df.latitude.values))

    xlon = np.arange(int(xmin),int(xmax),dres)
    ylat = np.arange(int(ymin),int(ymax),dres)
    onetime = df.copy()
    onetime = fixlondf(onetime,colname='longitude',neg=False)
    hhh, xedges, yedges = np.histogram2d(onetime.longitude.values,onetime.latitude.values,weights=onetime.weight,bins=[xlon,ylat])
    cb = plt.pcolormesh(xedges, yedges, hhh.T, cmap='viridis')
    plt.colorbar(cb)
    if isinstance(vloc,(list,tuple)):
        plt.plot(360+vloc[1],vloc[0], 'y^',markersize=10)
    return hhh, xedges, yedges

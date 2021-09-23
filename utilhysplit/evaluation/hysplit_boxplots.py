import numpy as np
import pandas as pd
import datetime
import matplotlib.pyplot as plt
import seaborn as sns


def make_boxplot(dj, cols=None):
    """
    dj : pandas dataframe
    cols : list of indices of columns to plot box plots for
    """
    sns.set_style('whitegrid')
    fig = plt.figure(1,figsize=(20,5))
    if isinstance(cols,(list,np.ndarray)):
        dj.iloc[:,cols].boxplot()
    else:
        dj.boxplot()
    ax = plt.gca()
    ax.set_yscale('log')
    #ax.set_xlim([d1,d2])
    fig.autofmt_xdate()
    plt.show()  

def prepare_boxplotdata(datelist, vdata):
    """
    datelist : list
    vdata : list of lists
    """

    newdata = []
    mlen=2
    for ra in vdata:
        ra = np.array(ra)
        mlen = np.max([mlen,ra.shape[0]])
    for ra in vdata:
        ra = np.array(ra)
        if ra.shape[0] < mlen:
           newra = np.pad(ra, (0,mlen-ra.shape[0]),'constant',constant_values=(0,-999))
        else:
           newra = ra
        newdata.append(newra)
    newdata = np.array(newdata)
    dj = pd.DataFrame(newdata.T, columns=datelist)
    dj = dj.replace(-999,np.NaN)
    return dj


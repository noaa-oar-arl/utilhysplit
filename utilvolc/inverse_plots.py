import datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import seaborn as sns
import utilvolc.ash_inverse as ai
from utilhysplit.plotutils import colormaker

"""

"""

def getlabels(subdirlist):
    labels = []
    for sdir in subdirlist:
        temp = sdir.split('/')
        temp = temp[-1]
        labels.append(temp)
    #dt = [0,1,2,3,4,5,6,7,8]
    #d1 = datetime.datetime(2020,10,21,23)
    #dlist = [d1+datetime.timedelta(hours=n*1) for n in dt]
    #print(dlist)
    #labels = [x.strftime("%m/%d %H:00 UTC") for x in dlist]
    return labels

def emisplots(subdirlist,fignum=1, tag='A',labels=None,plotmean=False,hunit='km'):
    dflist = []
    #print('tag {}'.format(tag))
    fig2 = plt.figure(fignum,figsize=(10,5))
    clrs = colormaker.ColorMaker('viridis',len(subdirlist),ctype='hex',transparency=None)
    colors=clrs()
    sns.set_style('whitegrid')
    nnn=1
    ccc=2
    ax1 = fig2.add_subplot(nnn,ccc,1)
    ax2 = fig2.add_subplot(nnn,ccc,2)
    if not labels:
        labels = getlabels(subdirlist)

    totalmassb = []
    totalmass2b = []
    #for iii, subdir in enumerate(runtag):
    for iii, subdir in enumerate(subdirlist):
        #print(clrlist[iii])
        print(subdir)
        temp = subdir.split('/')
        cname = os.path.join(subdir, temp[-1] + '.csv')
        if not os.path.isfile(cname): 
           print('could not find {}'.format(cname))
           continue
        df = ai.read_emis_df(cname)
        dflist.append(df)
        ai.plot_outdat_ts_function(df,ax=ax1,clr='#'+colors[iii],marker='.')
        #print(df)
        axr,mass = ai.plot_outdat_profile_function(df,ax=ax2,clr='#'+colors[iii],label=labels[iii],marker='.',unit=hunit)
        totalmassb.append(mass)
      
        temp = pd.concat(dflist)

    temp = temp.reset_index().groupby('0').mean()
    if plotmean:
        axr,mass = ai.plot_outdat_profile_function(temp,ax=ax2,clr='r',label='mean',lw=5,alpha=0.5,unit=hunit)
    ai.plot_outdat_ts_function(temp,ax=ax1,clr='r',lw=5,alpha=0.5)
    handles,labels = axr.get_legend_handles_labels()
    fig2.autofmt_xdate()

    for ax in [ax1,ax2]:
        ax.tick_params(labelsize=15)
    fig2.autofmt_xdate()

    xpos=0.80
    ypos=0.90
    ax1.text(xpos,ypos,'(a)',transform=ax1.transAxes,bbox=dict(facecolor='white'),size=20,color='k')
    ax2.text(xpos,ypos,'(b)',transform=ax2.transAxes,bbox=dict(facecolor='white'),size=20,color='k')
    plt.tight_layout()
    plt.savefig('emissions{}.png'.format(tag))
    plt.show()

    # create the legend
    figlegend = plt.figure()
    sns.set_style('white')
    axg = figlegend.add_subplot(1,1,1)
    axg.legend(handles,labels,loc="center",fontsize=20)
    axg.axis("off")
    plt.tight_layout()
    plt.savefig('emissions{}_legend.png'.format(tag))
    plt.show()

    labels = np.arange(0,len(totalmassb))

    fig = plt.figure(3)
    sns.set_style('whitegrid')
    axm = fig.add_subplot(1,1,1)
    axm.plot(labels, totalmassb,'-k.')
    axm.tick_params(labelsize=15)
    axm.set_ylabel('Total Mass Emitted (Tg)',fontsize=15)
    fig.autofmt_xdate()
    axm.text(xpos,ypos,'(e)',transform=axm.transAxes,bbox=dict(facecolor='white'),size=20,color='k')
    plt.tight_layout()
    plt.savefig('emissions{}_totalmass.png'.format(tag))
    return totalmassb, totalmass2b

 

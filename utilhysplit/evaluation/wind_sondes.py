import datetime
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys

from utilhysplit.evaluation import reliability

# 01/05/2022
# read files for wind sonde comparison from Binyu.

def test1(wmax):
    ddd = 5
    wmin = 320
    wlist = np.arange(wmin,wmax,ddd)
    wlist = [x if x<=360 else x-360 for x in wlist]
    return wlist


def readfile(fname):
    dstr = "%Y%m%d_%H%M00"
    df = pd.read_csv(fname, sep='\s+',parse_dates=['DATE'])
    df['date'] = df.apply(lambda row: datetime.datetime.strptime(row['DATE'],dstr),axis=1)
    return df 

def ccw_diff(start,stop):
    """
    counter clockwise difference
    start : starting angle
    stop : stopping angle
    Returns :
        how many degrees need to move counter-clockwise to get from small to large.
    """
    if start > stop:
       diff = start-stop
    else:
       diff = start + 360 - stop
    return diff

def cw_diff(start, stop):
    """
    clockwise difference
    start : starting angle
    stop : stopping angle
    Returns :
        how many degrees need to move clocwise to get from small to large.
    """
    dres = 0.001
    if stop >= start:
       diff = stop - start
    else:
       diff = stop + 360-start
    if np.abs(diff-360) < dres:
       diff = 0
    return diff    

def diff_wdir(wdir1, wdir2):
    if wdir1 > wdir2:
       big = wdir1
       small = wdir2
    else:
       big = wdir2
       small = wdir1
    diff1 = big - small
    diff2 = small + (360-big)
    if diff1 > diff2: 
       diff = diff2
       # move counter clockwise from smallest to largest.
       direction = 'ccw' #counter clockwise
    else: 
       diff = diff1
       # move clockwise from smallest to largest.
       direction = 'cw'  #clockwise
    return direction, diff

def spread(wlistin,verbose=False):
    """
    wlist: list of wind directions.
    order from least to greatest (0 to 360)
    Take clock-wise difference between each subsequent one.
    Find largest clock-wise gap.
    Set start degree to direction that occurs after this gap.
    Set end degree to direction that occurs before the gap.
    Return:
        start degree
        end degree
        spread 
    """
    wlist = wlistin.copy()
    if isinstance(wlist, pd.Series):
       wlist = wlist.values
    wlist.sort()
    dlist = []
    dmax = 0
    imax = 0
    if verbose: print(wlist)
    for iii, www in enumerate(wlist[1:]):
        diff = cw_diff(wlist[iii],www)
        if verbose: print(www,wlist[iii],diff)
        dlist.append(diff)
        if diff > dmax:
           dmax = diff
           imax = iii
    # check diff between last and first.
    diff = cw_diff(wlist[-1],wlist[0])
    if verbose: print(wlist[-1],wlist[0], diff)
    if verbose: print('current', dmax, imax)
    dlist.append(diff)
    if diff > dmax:
       dmax = diff
       imax = len(wlist)-1
       imin = 0
       if verbose: print('new', dmax, imax)
    else:
       imin = imax + 1
    spread = cw_diff(wlist[imin], wlist[imax])
    return wlist[imin], wlist[imax],  spread

def process_files(flist, addobs=None,wspd=None):
    """
    flist : list of files
    addobs : str. 'ODIR' or 'O_SPEED'
    wspd : float. If not None then removes rows where observed wind speed is below this value.
    RETURNS
    dfall : pandas dataframe.
    """
    dflist = []

    for fname in flist:
        df = readfile(fname)
        ens = fname.split('.')[0]
        if '/' in ens:
            ens  = ens.split('/')[-1]
        df['ens'] = ens
        #df = df[['date','FDIR','ens']].copy()
        dflist.append(df)

    dfall = pd.concat(dflist)
    if wspd:
       dfall = dfall[dfall['O_SPEED']>wspd]

    # creates dataframe with only three columns, date, speed or direction, ensid.
    if isinstance(addobs,str):
        if not addobs in ['ODIR','O_SPEED']:
           print('addobs not recognized {}'.format(addobs))
        else:   
           dfnew = readfile(flist[0])
           obs = dfnew[['date',addobs]].copy()
           obs = obs.drop_duplicates()
           obs['ens'] = 'obs'

           if addobs=='ODIR': new='FDIR'
           elif addobs=='O_SPEED': new = 'F_SPEED'

           obs.columns = ['date',new,'ens']
        dfall = dfall[['date',new,'ens']]
        dfall = pd.concat([dfall,obs])

    return dfall

def get_wspd_df(flist, addobs=True,addspread=False,wspd=None):
    if addobs: aaa = 'O_SPEED'
    else: aaa=None
    dfall = process_files(flist,addobs=aaa,wspd=wspd)
    dfwdir = pd.pivot(dfall, columns = 'ens', values='F_SPEED', index='date')
    return dfwdir

def get_wdir_df(flist, addobs=True,addspread=True,wspd=None):
    """
    flist : list of files
    addobs : boolean
    addspread : boolean
    wspd : float. If not None then removes rows where observed wind speed is below this value.
    RETURNS
    dfwdir : pandas dataframe.
             index is date
             values are wind direction.
             columns are ensemble member, obs
             if addspread then last three columns give starting angle, ending angle and spread.  
    """
    if addobs: aaa = 'ODIR'
    else: aaa = None
    dfall = process_files(flist,addobs=aaa,wspd=wspd)
    dfwdir = pd.pivot(dfall, columns = 'ens', values='FDIR', index='date')
    if addspread:
        dfwdir = add_spread(dfwdir)
    return dfwdir

def add_spread(dfwdirin):
    dfwdir = dfwdirin.copy()
    dftemp = dfwdir.apply(lambda row: spread(row), axis=1)
    dftemp = dftemp.to_frame().reset_index()
    dftemp.columns = ['date','temp']
    dftemp['spread'] = dftemp.apply(lambda row: row['temp'][2],axis=1)
    dftemp['start']  = dftemp.apply(lambda row: row['temp'][0],axis=1)
    dftemp['stop']  = dftemp.apply(lambda row: row['temp'][1],axis=1)
    dftemp = dftemp[['date','spread','start','stop']]
    dfwdir = pd.concat([dfwdirin, dftemp.set_index('date')], axis=1)
    return dfwdir


def recalc_wdir(dfwdir):
    """
    dfwdir : pandas dataframe with wdir for obs, ensemble members.
             Also needs columns start, stop and spread.
    wdir values are recalculated as degrees clock-wise from the starting value.
    """
    df = dfwdir.copy()
    cnames = dfwdir.columns
    cnames = [x for x in cnames if x not in ['spread','start','stop']]

    def calc(xval,start):
        # calculate clockwise distance from start. 
        newval = cw_diff(start,xval) 
        return newval

    for ens in cnames:
        df[ens] = df.apply(lambda row: calc(row[ens],row['start']),axis=1)
    return df 


class Comparison:

    def __init__(self,tdir,tag='A',plist =['P200','P300','P500'],wspd=0):
        """
        tdir : directory where files are located
        tag  : site identifier
        plist : filenames should be ensid.com.p*.txt where p* indicates the pressure                level. plist gives the pressure levels to use.
        """
        self.tag = tag
        self.tdir = tdir
        self.plist = plist
        self.filehash = {}

        self.dfhash = {}  # wind direction with obs added
        self.dfhash2 = {} # wind diection without obs added

        self.spdhash = {} # wind speed with obs added

        self.rankhash = {}  # rank histograms for wind direction
        self.srankhash = {} # rank histograms for wind speed

        self.dfallhash = {} # 

        self.recalc(wspd=wspd)


    def recalc(self,wspd):
        self.dfhash = {}  # wind direction with obs added
        self.dfhash2 = {} # wind diection without obs added

        self.spdhash = {} # wind speed with obs added

        self.rankhash = {}  # rank histograms for wind direction
        self.srankhash = {} # rank histograms for wind speed

        self.dfallhash = {} # 
  
        for pressure in self.plist:
            fff = '/*.com.{}.txt'.format(pressure)
            self.filehash[pressure] = glob.glob(self.tdir + fff) 

            self.dfallhash[pressure] = process_files(self.filehash[pressure])

   
            self.dfhash[pressure] = get_wdir_df(self.filehash[pressure],
                                             addobs=True,
                                             addspread=True,
                                             wspd=wspd) 
            self.dfhash2[pressure] = get_wdir_df(self.filehash[pressure],
                                             addobs=False,
                                             addspread=True, 
                                             wspd=wspd) 

            self.spdhash[pressure] = get_wspd_df(self.filehash[pressure],
                                             addobs=True,
                                             wspd=wspd) 



    def _return_plist(self,pressure=None):
        if not isinstance(pressure,(list,str)):
           pressure = self.plist
        elif isinstance(pressure,str):
           pressure = [pressure]
        return pressure


    def plot_diff(self,fignum=1):
        cmap = 'viridis'
        bins=[25,72]
        yyy = 15 
        pressure = self._return_plist()
        fig = plt.figure(fignum,figsize=[10,15])
        sns.set_style('whitegrid')
        ax1 = fig.add_subplot(3,1,1)
        ax2 = fig.add_subplot(3,1,2)
        ax3 = fig.add_subplot(3,1,3)
        axlist = [ax1,ax2,ax3]
        label_list = ['(a)','(b)','(c)']
        xll = 0
        for iii, prss in enumerate(pressure):
            dfall = self.dfallhash[prss]
            biny = np.arange(-180,185,5)
            xmax = np.max(dfall.O_SPEED)
            if xmax> xll: xll=xmax
            binx = np.arange(0,xmax+5,2)
            bins = [binx,biny]
            cb = axlist[iii].hist2d(dfall.O_SPEED, dfall.DIR_ERR, bins=bins,cmap=cmap) 
            plt.colorbar(cb[3])
            #axlist[iii].colorbar(cb)
            plbl = prss.replace('P','') + ' mb'
            label_list[iii] += '{}'.format(plbl)
        axlist[-1].set_xlabel('Observed wind speed (m/s)')
        axlist[-1].set_ylabel('Wind Direction error', fontsize=10)
        plt.tight_layout()
        for iii, ax in enumerate(axlist):
            ms=10
            lbl = label_list[iii]
            ax.set_ylabel('Wind Direction error', fontsize=10)
            ax.text(0.01,1.01,lbl,va='bottom',ha='left',rotation='horizontal',transform=ax.transAxes,size=ms)
            ax.plot([0,xll],[yyy,yyy],'--w',linewidth=1)
            ax.plot([0,xll],[-yyy,-yyy],'--w',linewidth=1)
            ax.set_xlim(0,xll)
        fname = '{}_directionerr_vs_ospeed.png'.format(self.tag)
        fig.savefig(os.path.join(self.tdir,fname))
        return fig

    def plotspread(self,pressure=None):
        pressure = self._return_plist(pressure)
        for prss in pressure:
            sns.set_style('whitegrid')
            if not prss in self.dfhash.keys():
               print('pressure not found {}'.format(prss)) 
               continue
            plt.plot(self.dfhash2[prss].spread, '--b.', linewidth=0.2, label='with obs')
            plt.plot(self.dfhash[prss].spread, 'r+', markersize=4, label='No obs')
            plt.title(self.tag + ' ' + prss)
            plt.xticks(rotation=45)
            plt.show()




    def plot_hist_spread(self,pressure=None,fname=None,fignum=2):
        pressure = self._return_plist(pressure)
        fig = plt.figure(fignum,figsize=[12,12])
        ax1 = fig.add_subplot(3,3,1)
        ax2 = fig.add_subplot(3,3,2)
        ax3 = fig.add_subplot(3,3,3)
        ax4 = fig.add_subplot(3,3,4)
        ax5 = fig.add_subplot(3,3,5)
        ax6 = fig.add_subplot(3,3,6)
        ax7 = fig.add_subplot(3,3,7)
        ax8 = fig.add_subplot(3,3,8)
        ax9 = fig.add_subplot(3,3,9)
        axlist = [ax1,ax4,ax7,ax2,ax5,ax8,ax3,ax6,ax9]
        label_list = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'] 
        for iii, prss in enumerate(pressure):
            sns.set_style('whitegrid')
            # use the spread calculated without the obs.
            new = self.dfhash2[prss].copy()
            obsdf = self.dfhash[prss].copy()
            new = new[['spread','start','stop']]
            bins = np.arange(0,370,10)
            axlist[iii].hist(new['spread'],bins=bins)
            axlist[iii+3].hist(new['start'],bins=bins,alpha=0.5)
            axlist[iii+3].hist(new['stop'],bins=bins,alpha=0.5)
            axlist[iii+6].hist(obsdf['obs'],bins=bins,color='g',alpha=0.5)

            plbl = prss.replace('P','') + ' mb'
            label_list[iii] += ' {} spread'.format(plbl)
            label_list[iii+3] += ' start (blue) and stop (red)'
            label_list[iii+6] += ' observed'
            #new.boxplot(ax=axlist[iii])
        for iii, ax in enumerate(axlist):
            ms=10
            lbl = label_list[iii]
            ax.text(0.01,1.01,lbl,va='bottom',ha='left',rotation='horizontal',transform=ax.transAxes,size=ms)
        xlbl = 'Wind Direction (degrees)'
        ylbl = 'Counts'
        for ax in [ax7, ax8,ax9]:
            ax.set_xlabel(xlbl)    
        for ax in [ax1,ax4,ax7]:
            ax.set_ylabel(ylbl)       
 
        if isinstance(fname,str):
            if fname=='default':
               fname = '{}_SpreadHist.png'.format(self.tag)
            fig.savefig(os.path.join(self.tdir,fname))
        return fig


    def plot_mindiffhist(self,fignum=10,binwidth=5):
        fig = plt.figure(fignum,figsize=[15,5])
        ax1 = fig.add_subplot(1,3,1)
        ax2 = fig.add_subplot(1,3,2)
        ax3 = fig.add_subplot(1,3,3)
        pressure = self._return_plist(None)
        axlist = [ax1,ax2,ax3]
        label_list = ['(a)','(b)','(c)'] 
        for iii, prss in enumerate(pressure):
            sns.set_style('whitegrid')
            temp = self.dfallhash[prss]
            temp = pd.pivot(temp,columns='ens',values='DIR_ERR',index='date')
            # get absolute value of differences
            temp = np.abs(temp)
            # get minimum in each row 
            mindiff = temp.min(axis=1)
            # default bins are 5 wide 
            step = binwidth
            bins = np.arange(0,180+step,step)
            # set density equal to true
            axlist[iii].hist(mindiff.values,bins=bins,alpha=0.5,density=True)
            plbl = prss.replace('P','') + ' mb'
            label_list[iii] += ' {}'.format(plbl)
        for iii, ax in enumerate(axlist):
            # since bins are 5 wide, multiply the y axis value by 5. 

            yticks = np.array([0,10,20,30,40,50,60,70,80,90,100])
            yticklabels = [str(x) for x in yticks]

            yticks = yticks/100.0/step
            ax.set_yticks(yticks)
            ax.set_ylim(0,yticks[-1])
            ax.set_yticklabels(yticklabels)
            ax.set_xlabel('smallest difference in wind direction')
            ms=10
            lbl = label_list[iii]
            ax.text(0.70,0.90,lbl,va='bottom',ha='left',rotation='horizontal',transform=ax.transAxes,size=ms)
        ax1.set_ylabel('percent in bin')
        fname = '{}_mindiffhist_wspd.png'.format(self.tag)
        fig.savefig(os.path.join(self.tdir,fname))


    def plot_boxB(self,daterange,fname=None,spd=True,fignum=3):
        fig = plt.figure(fignum,figsize=[12,12])
        ax1 = fig.add_subplot(3,1,1)
        ax2 = fig.add_subplot(3,1,2)
        ax3 = fig.add_subplot(3,1,3)
        axlist = [ax1,ax2,ax3]
        label_list = ['(a)','(b)','(c)'] 
        pressure = self._return_plist(None)
        for iii, prss in enumerate(pressure):
            sns.set_style('whitegrid')
            temp = self.spdhash[prss].reset_index()
            temp = temp[temp.date >= daterange[0]]    
            temp = temp[temp.date <= daterange[1]]    
            temp = temp.set_index('date')
            obs = temp['obs']
            temp = temp.drop(['obs'],axis=1)
            temp = temp.T
            sns.boxplot(data=temp,color='grey',saturation=0.5,ax=axlist[iii])
            xxx = [str(pd.to_datetime(x)) for x in obs.index.values]
            axlist[iii].plot(xxx,obs.values,'ro',markersize=5)
            plt.xticks(rotation=45)
            plbl = prss.replace('P','') + ' mb'
            label_list[iii] += ' {}'.format(plbl)
            plt.xticks(np.arange(0,len(xxx)+1,2))
        for iii, ax in enumerate(axlist):
            ms=10
            lbl = label_list[iii]
            ax.text(0.01,1.05,lbl,va='bottom',ha='left',rotation='horizontal',transform=ax.transAxes,size=ms)
            ax.set_ylabel('Wind Speed (m/s)')

        ax1.set_xticklabels([])
        ax2.set_xticklabels([])
        ax1.set_xlabel('')
        ax2.set_xlabel('')
        if isinstance(fname,str):
            if fname=='default':
               sss='WSPD'
               dstr = '%Y%m%d'
               ddd = '{}_{}'.format(daterange[0].strftime(dstr), daterange[1].strftime(dstr))
               fname = '{}_BoxPlotB_{}_{}.png'.format(self.tag,sss,ddd)
            fig.savefig(os.path.join(self.tdir,fname))
        return fig

    def plot_boxA(self,pressure=None,fname=None,spd=True,fignum=4):
        """
        plot boxplots for distribution for each ensemble member and observation.
        """
        pressure = self._return_plist(pressure)
        fig = plt.figure(fignum,figsize=[12,12])
        ax1 = fig.add_subplot(3,1,1)
        ax2 = fig.add_subplot(3,1,2)
        ax3 = fig.add_subplot(3,1,3)
        axlist = [ax1,ax2,ax3]
        label_list = ['(a)','(b)','(c)'] 
        for iii, prss in enumerate(pressure):
            sns.set_style('whitegrid')
            if spd:
                new = self.spdhash[prss].copy()
            else:
                new = self.dfhash[prss].copy()
                new = new.drop(['spread','start','stop'],axis=1)
            new.boxplot(ax=axlist[iii])  
            plt.xticks(rotation=45)
            plbl = prss.replace('P','') + ' mb'
            label_list[iii] += ' {}'.format(plbl)
        ax1.set_xticklabels([])
        ax2.set_xticklabels([])
        for iii, ax in enumerate(axlist):
            ms=10
            lbl = label_list[iii]
            ax.text(0.01,1.05,lbl,va='bottom',ha='left',rotation='horizontal',transform=ax.transAxes,size=ms)
        if isinstance(fname,str):
            if fname=='default':
               if spd: sss = 'WSPD'
               else: sss='WDIR'
               fname = '{}_BoxPlot_{}.png'.format(self.tag,sss)
            fig.savefig(os.path.join(self.tdir,fname))
        return fig

    def get_axlist(self,fig):
        ax1 = fig.add_subplot(2,3,1)
        ax2 = fig.add_subplot(2,3,2)
        ax3 = fig.add_subplot(2,3,3)
        ax4 = fig.add_subplot(2,3,4)
        ax5 = fig.add_subplot(2,3,5)
        ax6 = fig.add_subplot(2,3,6)
        axlist = [ax1,ax2,ax3,ax4,ax5,ax6]
        label_list = ['(a)','(b)','(c)','(d)','(e)','(f)']
        return axlist, label_list

    def calc_rank_diagram(self,pressure=None,fname=None,fignum=5):
        pressure = self._return_plist(pressure)
        # TODO adjust number of axis depending on length of pressure list.
        fig = plt.figure(fignum,figsize=[12,5])
        axlist, label_list = self.get_axlist(fig)
        for iii, prss in enumerate(pressure):
            # recalculates wind direction with reference to starting angle.
            new = recalc_wdir(self.dfhash[prss])
            new = new.drop(['spread','start','stop'],axis=1)
            if prss not in self.rankhash.keys():
                drank = reliability.Talagrand(thresh=0,nbins=32)
                drank.add_data(new)
                srank = reliability.Talagrand(thresh=0,nbins=32)
                srank.add_data(self.spdhash[prss])
                self.rankhash[prss] = drank
                self.srankhash[prss] = srank
            plbl = prss.replace('P','') + ' mb'
            label_list[iii]  += ' direction ' + plbl 
            label_list[iii+3] += ' speed ' + plbl 
            self.rankhash[prss].plotrank(ax=axlist[iii])
            self.srankhash[prss].plotrank(ax=axlist[iii+3])
        plt.sca(axlist[1])
        plt.title(self.tag)
        # add text labels
        for iii, ax in enumerate(axlist):
            ms=10
            lbl = label_list[iii]
            ax.text(0.05,0.8,lbl,va='bottom',ha='left',rotation='horizontal',transform=ax.transAxes,size=ms)

        if isinstance(fname,str):
            if fname=='default':
               fname = '{}_RankHistogram.png'.format(self.tag)
            fig.savefig(os.path.join(self.tdir,fname))
        return fig, axlist

            
if __name__=="__main__":
   
   tdir = sys.argv[1]
   tag = sys.argv[2]
   station = Comparison(tdir,tag=tag)

   fig, axlist = station.calc_rank_diagram(fname='default')  
   fig.clear()

   fig = station.plot_hist_spread(fname='default')
   fig.clear()

   fig = station.plot_boxA(fname='default')
   fig.clear()
 

   df = station.dfallhash['P200']
   dates=df['date'].values 
   d1 = pd.to_datetime(dates[0])
   d2 = d1 + datetime.timedelta(hours=30*24)
   print('MAX date', np.max(dates))
   df2 = station.dfhash['P200']
   print('sIZE', df2.shape)
   station.plot_boxB(daterange=[d1,d2],fname='default')
   fig.clear()

   fig = station.plot_diff()
   fig.clear() 
   station.plot_mindiffhist(binwidth=5)


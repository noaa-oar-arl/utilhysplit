import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from utilhysplit.evaluation import web_ensemble_plots as wep
from utilhysplit.fixlondf import fixlondf
from utilhysplit.plotutils import map_util
from utilvolc import utiltraj
from utilhysplit import geotools
import shapely.geometry as sgeo
from utilvolc.helperinterface import DFInterface

class TrajDF:
    required = ['time','traj_num','met_grid','forecast_hour','traj_age','latitude','longitude',
                'altitude','pressure','run_num','weight','obsalt','area']

    def __init__(self):
        self._edf = pd.DataFrame()
       

    @property
    def edf(self):
        return self._edf

    # to do - add some checks.
    @edf.setter
    def edf(self,edf):
        self._edf = edf

class TrajAnalDF(DFInterface):
    # columns which are required
    required=['run_num','traj_num','dist','Alt Initial','Closest Point']

    def __init__(self):
        self._edf = pd.DataFrame()

    @property
    def required_columns(self):
        return self.required

    @property
    def edf(self):
        return self._edf

    @edf.setter
    def edf(self,edf):
        if self.are_required_columns_present(edf):
            self._edf = edf

    def read(self):
        TODO = -1
        return TODO

    def save(self,overwrite=False):
        TODO = -1
        return TODO

    def add_df(self,df):
        TODO = -1
        return TODO

    def are_required_columns_present(self,df):
        answer=True
        for req in self.required:
            if req not in df.columns:
               logger.info('TRAJAnalDF WARNING: data frame does not contain required columns {}'.format(req))
               answer = False
        return answer 


class TrajAnalysis:
    """
    trajectory analysis

    """

    def __init__(self):
        """
        Analysis of pandas dataframe with many trajectories.
        """

        self.eruption_start = None # datetime
        self.eruption_end   = None # datetime
        self.vloc = None
        self._trajdf = TrajDF()
        self._analdf = TrajAnalDF() 
 
    @property
    def trajdf(self):
        return self._trajdf.edf

    @trajdf.setter
    def trajdf(self, edf):
        self._trajdf.edf = edf 
   
    @property
    def analdf(self):
        return self._analdf.edf

    @analdf.setter
    def analdf(self, edf):
        self._analdf.edf = edf 
   
    @property
    def start_time(self):
        start = self.get_start()
        return start.time.unique()[0]

    def get_start_hull(self,alpha):
        start = self.get_start()
        lon = start.longitude.values
        lat = start.latitude.values
        mpts = geotools.make_multi(lon,lat)
        numpts = len(lon)
        if numpts > 4:
           ch, ep = geotools.concave_hull(mpts,alpha)
        else:
           ch = mpt.convex_hull
           ep = None
        return ch, ep

    def compare(self,other,trajlist):
        linestyle=''
        markersize=8
        shull,ep = other.get_start_hull(alpha=1)
        ctime = other.start_time 
        clow = datetime.datetime(ctime.year, ctime.month,ctime.day,ctime.hour)
        
        dfall2 = self.trajdf
        dfall2 = dfall2[dfall2.time>=clow]
        dfall2 = dfall2[dfall2.time<=clow+datetime.timedelta(hours=0.5)]

        start = self.trajdf[self.trajdf.traj_age==0]
        if isinstance(trajlist,list):
           dfall2 = dfall2[dfall2.run_num.isin(trajlist)]
           start2 = start[start.run_num.isin(trajlist)]
        else:
           start2 = start
        for rnum in dfall2.run_num.unique():
            temp = dfall2[dfall2.run_num==rnum]
            alist = []
            zlist = []
            for trajnum in temp.traj_num.unique():
                #print('working on traj {}'.format(trajnum))
                temp2 = temp[temp.traj_num==trajnum]
                ialt = start[start.run_num == rnum]
                ialt = ialt[ialt.traj_num==trajnum]
                lat = ialt.latitude.values
                lon = ialt.longitude.values
                ialt = ialt.altitude.values[0]
                point = sgeo.Point(temp2.longitude.values[0], temp2.latitude.values[0])
                zval = temp2.altitude.values[0]
                if shull.contains(point): 
                    plt.plot(temp2.longitude, temp2.latitude,linestyle=linestyle,
                         marker='.',color='b',markersize=markersize,label=ialt)
                    alist.append((lon,lat,ialt,point,zval))
                    ax = plt.gca()
                    ax.text(temp2.longitude,temp2.latitude,str(ialt),transform=ax.transAxes)
                else:    
                    plt.plot(temp2.longitude, temp2.latitude,linestyle=linestyle,
                         marker='+',color='r',markersize=markersize,label=ialt)
                    zlist.append((lon,lat,ialt,point,zval))
                    
        return alist,zlist 

    def make_dist_df(self,vloc):
        self.vloc = vloc
        blist = []
        dfall2 = self.trajdf
        dfall2 = dfall2[dfall2.time>=self.eruption_start]
        dfall2 = dfall2[dfall2.time<=self.eruption_end]
        start = self.trajdf[self.trajdf.traj_age==0]
        for rrr in dfall2.run_num.unique():
            tdump = dfall2[dfall2.run_num==rrr] 
            for traj_num in tdump.traj_num.unique():
                temp = tdump[tdump.traj_num==traj_num]
                a, dist = utiltraj.traj_volc_dist2(vloc,temp)
                obsalt = temp.obsalt.values
                ialt = start[start.run_num == rrr]
                ialt = ialt[ialt.traj_num==traj_num].altitude.values[0]
                blist.append((rrr,traj_num,dist,ialt,a))
             
        self.analdf = pd.DataFrame.from_records(blist, 
                                columns=['run_num','traj_num','dist','Alt Initial','Closest Point'])

        return self.analdf
 
    def plotdist(self,run_num,vloc):
        self.vloc = vloc
        blist = []
        sns.set()
        dfall2 = self.trajdf
        dfall2 = dfall2[dfall2.time>=self.eruption_start]
        dfall2 = dfall2[dfall2.time<=self.eruption_end]
        start = self.trajdf[self.trajdf.traj_age==0]
        if isinstance(run_num,(int,float)):
           run_num = [run_num]
        tdump0 = dfall2[dfall2.run_num.isin(run_num)]
        for rrr in tdump0.run_num.unique():
            tdump = tdump0[tdump0.run_num==rrr] 
            alist=[]
            for traj_num in tdump.traj_num.unique():
                temp = tdump[tdump.traj_num==traj_num]
                #print('HERE2', traj_num, temp)
                a, dist = utiltraj.traj_volc_dist2(vloc,temp)
                #print('Here3', dist)
                obsalt = temp.obsalt.values
                ialt = start[start.run_num == rrr]
                ialt = ialt[ialt.traj_num==traj_num].altitude.values[0]

                blist.append((rrr,traj_num,dist,ialt,a))
                alist.append((traj_num,dist,ialt,a))
                #print(a.x, a.y)
                #print(traj_num, dist, ialt)
               
            #print('observed altitude is', set(obsalt)) 
            #print('-----------------')
            zzz = list(zip(*alist))
            plt.plot(zzz[2],zzz[1],'--k.')
            obs = obsalt[0]*1000.0
            plt.plot([obs,obs],[0,10],'-r')
        ax = plt.gca()
        ax.set_xlabel('initial altitude (m)')
        ax.set_ylabel('Closest distance (km)')

        analdf = pd.DataFrame.from_records(blist, 
                                columns=['run_num','traj_num','dist','Alt Initial','Closest Point'])
        self.analdf = analdf
        return self.analdf

    def hists(self,p=1,bins=None,dthresh=5):
        
        adf = pd.pivot_table(self.analdf,index='Alt Initial',columns='run_num',values='dist')
        fig = plt.figure(1)
        ax = fig.add_subplot(1,1,1)
        if p==1:
            ax.hist(self.analdf.dist)
            ax.set_xlabel('closest distance (km)')
            plt.title('All trajectories')
        elif p==2:
            hmin = adf.min(axis=1)
            hmin.plot(marker='.')
            hmean = adf.mean(axis=1)
            hmean.plot(marker='*')
            ax.set_xlabel('Initial Altitude')
            ax.set_ylabel('Closest distance (mean and min)')
        elif p==3:
            run_min = adf.min(axis=0)
            ax.hist(run_min,bins=bins)
            ax.set_xlabel('Closest distance to vent (km)')
            ax.set_ylabel('Number of instances')
            plt.title('One distance per observation point')
        elif p==4:
            nlevels = adf[adf<=dthresh].count()
            #nlevels.plot(marker='.',linestyle='')
            #plt.show()
            plt.hist(nlevels,bins=bins)
            ax=plt.gca()
            ax.set_xlabel('Number of levels that produce trajectories that come within {} km'.format(dthresh))
            ax.set_ylabel('Number of instances')
            plt.tight_layout()


    def get_start(self):
        start = self.trajdf[self.trajdf.traj_age==0]
        return start

    def plotA(self,vloc=None,trajlist=None,skip=25,markersize=1):
          
        sns.set()
        # only look at end points between possible eruption start and end.
        #dfall2 = self.trajdf[pd.to_datetime(self.trajdf.time)>=self.eruption_start]
        dfall2 = self.trajdf
        dfall2 = dfall2[dfall2.time>=self.eruption_start]
        dfall2 = dfall2[dfall2.time<=self.eruption_end]
        start = self.trajdf[self.trajdf.traj_age==0]
        if isinstance(trajlist,list):
           dfall2 = dfall2[dfall2.run_num.isin(trajlist)]
           start2 = start[start.run_num.isin(trajlist)]
        else:
           start2 = start

        for rnum in dfall2.run_num.unique():
            temp = dfall2[dfall2.run_num==rnum]
            if rnum%skip ==0: linestyle='--'
            else: linestyle = '' 
            for trajnum in temp.traj_num.unique():
                #print('working on traj {}'.format(trajnum))
                temp2 = temp[temp.traj_num==trajnum]
            
                ialt = start[start.run_num == rnum]
                ialt = ialt[ialt.traj_num==trajnum].altitude.values[0]
                plt.plot(temp2.longitude, temp2.latitude,linestyle=linestyle,
                         marker='.',markersize=markersize,label=ialt)
                #utiltraj.traj_volc_dist(vloc,temp)
            
            plt.plot(start.longitude, start.latitude,linestyle='',marker='.',alpha=0.5)
            plt.plot(start2.longitude, start2.latitude,linestyle='',marker='*',alpha=0.5)
            plt.plot(vloc[1],vloc[0],'r^',markersize=20)
            ax = plt.gca()
            handles,labels = ax.get_legend_handles_labels()


            
        



def plottraj_map(tdump, clr="-c", ms=0.5, alpha=0.5, central_longitude=180):
    """
    central_longitude should be 180 if trajectories are crossing dateline.
    Otherwise use 0.
    """
    transform = wep.get_transform(central_longitude)
    fig, axra = plt.subplots(
        nrows=1,
        ncols=1,
        figsize=(10, 10),
        constrained_layout=False,
        subplot_kw={"projection": transform},
    )
    plottraj(tdump, axra, clr, alpha, ms)
    wep.format_plot(axra, transform)
    return axra


def plottraj(tdump2, ax, clr="-c", alpha=0.5, ms=0.5):
    sns.set()
    tnum = tdump2.traj_num.unique()
    tlist = tnum
    for tnum in tlist:
        temp = tdump2[tdump2["traj_num"] == tnum]
        temp = temp.sort_values(by="time")
        alt = temp.iloc[0].altitude
        xval = temp["longitude"].values
        xval = wep.shift_xvals(xval, 180)
        ax.plot(xval, temp["latitude"], clr, linewidth=0.1, markersize=ms, alpha=alpha)
    return ax


def frequency_plots_all(df, sdate, dtt, dres=1, vloc=None):
    """
    df : pandas aataframe
    sdate : datetime object
    dtt :
    """
    times = df.time.unique()
    times.sort()
    done = False
    if isinstance(dtt, int):
        dtt = datetime.timedelta(hours=dtt)
    while not done:
        onetime = df[df["time"] == sdate]
        onetime = fixlondf(onetime, colname="longitude", neg=False)
        frequency_plot(onetime, dres, vloc)
        plt.title(sdate)
        sdate += dtt
        if sdate > pd.to_datetime(times[-1]):
            print("done", sdate)
            print(times)
            done = True
        plt.show()


def frequency_plot(df, dres=1, vloc=None):
    """
    df : pandas dataframe with columns of longitude, latitude.
         weight column is optional. If it does not exist then all points are
         evenly weighted.
    dres : resolution of bins in degrees.
    """
    sns.set()
    if "weight" not in df.columns:
        df["weight"] = 1
    xmin = np.floor(np.min(df.longitude.values))
    xmax = np.ceil(np.max(df.longitude.values))
    ymin = np.floor(np.min(df.latitude.values))
    ymax = np.ceil(np.max(df.latitude.values))

    xlon = np.arange(int(xmin), int(xmax), dres)
    ylat = np.arange(int(ymin), int(ymax), dres)
    onetime = df.copy()
    onetime = fixlondf(onetime, colname="longitude", neg=False)
    hhh, xedges, yedges = np.histogram2d(
        onetime.longitude.values,
        onetime.latitude.values,
        weights=onetime.weight,
        bins=[xlon, ylat],
    )
    cb = plt.pcolormesh(xedges, yedges, hhh.T, cmap="viridis")
    plt.colorbar(cb)
    if isinstance(vloc, (list, tuple)):
        plt.plot(360 + vloc[1], vloc[0], "y^", markersize=10)
    return hhh, xedges, yedges

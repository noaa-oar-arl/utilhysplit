import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns        
from utilhysplit.evaluation import ensemble_tools
 
def plot_freq(obs,prob,name='freq', thresh=2.5):
    #plt.figure(figsize=(5,20))
    fs=20
    figsize=(20,3)
    fig = plt.figure(10,figsize=figsize)
    sns.set_style('whitegrid', {"fontsize":fs})
    obscolor = '-k'
    probcolor = '-b'
    ax1 = fig.add_subplot(1,1,1)
    ax2 = ax1.twinx()
    ax2.plot(prob, probcolor  ,linewidth=2, alpha=0.7) 
    ax2.set_yticks([0,0.2,0.4,0.6,0.8,1.0])
    #ax2.set_yticklabels(['0','1','2','3','4','5'])
    ax2.tick_params(axis='y', which='major', labelsize=fs,colors='blue')
    ax1.tick_params(axis='y', which='major', labelsize=fs)
    ax1.tick_params(axis='x', which='major', labelsize=fs)
    ax2.grid(b=None, which='major', axis='y')

    ax1.set_ylabel('Measured \n Concentration \n (ppb)',fontsize=fs)
    ax2.set_ylabel('Frequency > {}ppb'.format(2.5),color='blue',fontsize=fs)
 
    ax1.plot(obs, obscolor, label='obs')  
    plt.savefig(name + '.png')
     

def process_data(df, resample_time, resample_type):
    if resample_type == 'max':
       df = df.resample(resample_time).max()  
    elif resample_type == 'mean':
       df = df.resample(resample_time).mean()  
    elif resample_type == 'sum':
       df = df.resample(resample_time).sum()  
    else:
       print('Warning: resample_type not recognized', resample_type)
    return df


def make_talagrand(dflist,thresh,bnum,
                     resample_time=None,
                     resample_type='max',
                     rname='talagrand',
                     background=0,
                     verbose=False):
    tal = Talagrand(thresh,bnum,background=background,verbose=verbose)
    for df in dflist:
        if resample_time:
           df = process_data(df,resample_time,resample_type)
        tal.add_data(df)
    #tal.plotrank(nbins=bnum)
    return tal

def make_reliability(dflist,thresh,bnum,
                     resample_time=None,
                     resample_type='max',
                     rname='reliability'):
    """
    dflist : list of pandas DataFrames
    """
    # seperate the observation column out.
    # bnum is number of bins in the reliability curve.
    
    # initialize reliability curve obejct.
    # num is number of bins to use.
    # here use one bin per simulation.
    rc = ReliabilityCurve(thresh, num=bnum)
    for df in dflist:
        if resample_time:
           df = process_data(df,resample_time,resample_type)

        obs_col = [x for x in df.columns if 'obs' in x]
        obs_col = obs_col[0]
        #print(obs_col)
        obsdf = df[obs_col]
        df2 = df.drop([obs_col],axis=1)
        # not sure why take the transpose?
        num = len(df2.columns)
        df2 = df2.T
        # num is number of simulations.
         
        # gives percentage of simulations above threshold.
        prob = (df2[df2 > thresh].count())/num 
        plot_freq(obsdf, prob, name=obs_col+'_freq',thresh=thresh)
        plt.savefig('{}_ts.png'.format(rname))
        plt.show()

        # add data to curve
        rc.reliability_add(obsdf, prob)
    # plot curve
    rc.reliability_plot(rname=rname) 
    return rc

class Talagrand:
    # get list of forecasts and order them.
    # find where in the list the observation would go.
    # p372 Wilks
    # if obs smaller than all forecasts then rank is 1
    # if obs is larger than all forecasts then rank is n_ens+1
    # calculate rank for each observation
    # create histogram.

    # simple way when there are not duplicat values of forecasts
    # is to just sort the list and find index of observation.

    # However when there are a lot of duplicate values then must
    # create the bins first and fill them as you go along. 

    # do we only look at observations above threshold?
    # What if multiple forecasts are the same? What rank is obs given?
    # This would occur in case of 0's usually. What if 10 forecasts are 0 and
    # observation is 0 and 5 forecasts are above 0? 
    # ANSWER: create the bins first - one for each forecast.
    # 

    def __init__(self, thresh,nbins,background=0,verbose=False):
        self.verbose = verbose
        self.thresh = thresh
        self.background = background
        self.ranklist = [] 
        self.obsvals = []
        # create bins for histogram.
        self.binra = np.zeros(nbins)
        # when obs are 0, count how many forecasts are 0.
        self.zeronum = []
        # when rank is 27, count how many forecasts 
        # are non-zero.
        # number which are 0 shows how many are completely missed.
        self.nonzero = []
        self.nonzerorow = []
        # value of obs which is higher than all forecasts
        self.obsmax = []
        # consider values with difference less than this the same.
        self.tolerance = 1e-2

    def check1(self,fname=None,fs=10):
        # creates scatter plot of observations whic
        # were higher than all forecasts and the
        # largest forecast
        import matplotlib.colors
        sns.set_style('whitegrid')
        obs = []
        maxf = []
        for vals in self.nonzerorow:
            obs.append(vals[-1])
            maxf.append(np.max(vals[0:-1]))
        yline = np.arange(0,int(np.max(maxf)),1)
        xline = yline
        norm = matplotlib.colors.LogNorm(vmin=0.1,vmax=None,clip=False)
        cb = plt.hist2d(obs,maxf,bins=(50,50),density=False,cmap=plt.cm.BuPu,norm=norm)
        plt.plot(xline,yline,'-k')
        plt.colorbar()
        ax = plt.gca()
        ax.set_xlabel('Measured Value (ppb)',fontsize=fs)
        ax.set_ylabel('Highest forecast(ppb)',fontsize=fs)
        if fname: plt.savefig(fname + '.png')
        #plt.plot(obs,maxf,'k.',markersize=5)
        #plt.show()
        return ax

    def plotrank(self, fname=None,fs=10):
        sns.set_style('whitegrid')
        nbins = self.binra.shape[0] + 1
        xval = np.arange(1,nbins)
        plt.bar(xval,self.binra,alpha=0.5)
        plt.xlabel('Rank',fontsize=fs)
        plt.ylabel('Counts',fontsize=fs)
        plt.tight_layout() 
        if fname:
           plt.savefig('{}.png'.format(fname))
       
    def plotrank_old(self,nbins):
        sns.set()
        plt.hist(self.ranklist,bins=nbins,histtype='bar',
                 color='b',rwidth=0.9,density=True,align='mid')
        #plt.show()
        #plt.hist(self.obsvals,bins=nbins,histtype='bar',
        #         color='g',rwidth=0.9,density=True,align='mid')

    def add_data_temp(self,df):
        # NOT WORKING.
        obs_col = [x for x in df.columns if 'obs' in x]
        obs_col = obs_col[0]
        # this keeps rows which have at least one value above threshold.
        df2 = df[(df>self.thresh).any(1)]
        for iii in np.arange(0,len(df2)):
            # selects row
            row = df2.iloc[iii]
            # sorts values
            row = row.sort_values()
            # creates dataframe with columns of index, name, value
            temp = pd.DataFrame(row).reset_index().reset_index()
            temp.columns = ['iii','name','val']
            # gets index of the observation
            # add one since index starts at 0.
            if self.verbose: print(temp)
             
            rank = float(temp[temp.name==obs_col].iii) 
            obsval = float(temp[temp.name==obs_col].val)

    def sort_wind_direction(self,row):
        # NOT WORKING YET.
        # check for when distance between first and last
        # value is greater than 360-R-1 + R0
        # 0--------------------------------------------360
        #     R0                             R-1
        #
        if (row[-1] - row[0]) > (360-row[-1] + row[0]):
            return -1 
 
    def add_data(self,df, wind_direction=False):
        # Assumes no duplicate values in forecast.
        # ToDO - should be able to apply to row so that the rank becomes
        # another column in the dataframe.
        obs_col = [x for x in df.columns if 'obs' in x]
        obs_col = obs_col[0]
        # this keeps rows which have at least one value above threshold.
        df2 = df[(df>self.thresh).any(1)]
        for iii in np.arange(0,len(df2)):
            # selects row
            row = df2.iloc[iii]
            rowcheck = row.copy()
            # if below background then set to 0.
            # otherwise was running into problems with
            # obs never having lowest rank.
            if row[obs_col] < self.background: 
               row[obs_col] = 0
            for col in row.index.values:
                if row[col] < self.background: 
                   row[col] = 0
            #print('ROW',row)
            # sorts values
            rowcheck = rowcheck.sort_values()
            row = row.sort_values()
            print('ROW',row)
            # creates dataframe with columns of index, name, value
            temp = pd.DataFrame(row).reset_index().reset_index()
            temp.columns = ['iii','name','val']
            #print('TEMP',temp)
            # gets index of the observation
            # add one since index starts at 0.
            
            # do not add +1 because now using it as index of self.binra.
            # which starts at 0. 
            rank = int(float(temp[temp.name==obs_col].iii))
            obsval = float(temp[temp.name==obs_col].val)
            print('Rank', rank, obsval)

            temp['zval'] = np.abs(obsval - temp['val'])
            # should be one zero value where obs matches itself.
            # if multiple forecasts are within the tolerance then
            # point is split among them evenly.
            if len(temp[temp.zval <= self.tolerance]) > 1:
                val2add = 1.0 / len(temp[temp.zval<=self.tolerance]) 
                # indices of where forecasts match obs.
                rankvals = np.where(temp.zval<=self.tolerance)
                # add value to binra.
                #print(rankvals[0], type(rankvals))
                for rrr in rankvals[0]:
                    rrr = int(rrr)
                    #print('Adding to binra', rank, self.binra)
                    self.binra[rrr] += val2add

            # else just add one to position of obs value.
            else:
                #print('Adding to binra', rank)
                self.binra[rank] += 1

            # This is incorrect way. Making rank middle
            # of multiple 0 values. 
            #if obsval < self.tolerance:
            #   rowz = row[row<self.tolerance]
            #   rank = int(len(rowz)/2.0)
            #   self.zeronum.append(len(rowz)-1)
     
            # this is used in check1 method.
            if rank == len(self.binra)-1:
               if self.verbose: print('high val', row)
               rowz = rowcheck[rowcheck>self.tolerance]
               self.nonzero.append(len(rowz)-1)
               self.nonzerorow.append(rowcheck)
               self.obsmax.append(obsval)
 
            #print(rank)
            self.ranklist.append(rank)                
            self.obsvals.append(obsval)          

 
class ReliabilityCurve:

    def __init__(self, thresh, num):
        """
        thresh : float
        num : int or list of floats.

        Attributes
        self.binlist : list of floats from 0 to 1.

        self.yeshash : dictionary
                       key is bin number, value is number of times
                       obs above threshold for this bin.

        self.nohash  : dictionary
                       key is bin number, value is number of times obs
                       below threshold for this bin.
                   
        self.thresh  : float
 

        """
        # use num to define the bins
        if isinstance(num, int):
            self.binlist = list(np.arange(0,1, 1.0/num) + 1.0/(num*10))
            self.binlist.append(1)
        elif isinstance(num, list):
            self.binlist = num
        else: 
            self.binlist = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
     
        self.yeshash = {}
        self.nohash = {}
        for binval in self.binlist:
            self.yeshash[binval] = 0
            self.nohash[binval] = 0
  
        self.thresh = thresh


    def reliability_add(self, obs, prob):
        """
        obs : pandas time series
        prob : pandas time series

        """
        df = pd.concat([obs, prob], axis=1)
        df.dropna(axis=0, inplace=True)
        cols = ['obs','prob']
        df.columns = cols
        self.reliability_add_sub(df)
        #print(type(df))
        #print(df[0:10])


    def reliability_add_xra(self,obs,forecast,fill=True):
        """
        obs      : xarray with observations. Needs to be same size in x,y as forecast.
        forecast : xarray with dimensions of x,y  and 'ens' and/or 'source'
        """
        prob = ensemble_tools.ATL(forecast, thresh= self.thresh, norm=True) 
        modelra = prob.values.flatten()
        obsra = obs.values.flatten()
        if modelra.size != obsra.size:
           print('Cannot add data, inputs not the same size {} {}', modelra.size, obsra.size)
           return -1
        else:
           dfin = pd.DataFrame(zip(obsra,modelra))
           dfin.columns = ['obs','prob']
        if fill:
            self.reliability_add_sub(dfin.fillna(0))
        else:
            self.reliability_add_sub(dfin.dropna())
        return dfin 
          

    def reliability_add_sub(self, dfin):
        """
        dfin : pandas dataframe with 'prob' and 'obs' columns.
        """
        for index, row in dfin.iterrows():
            prob = row['prob'] 
            for pbin in self.binlist:
                if prob <= pbin: 
                   binval = pbin
                   #print('break', binval)
                   break 
            if row['obs'] >= self.thresh:
               self.yeshash[binval] += 1
            else:
               self.nohash[binval] += 1


    def get_binvals(self):
        xxx = [] # the bin values for x axis
        yyy = [] # total number abo
        nnn = [] # total number
        yesval = [] 
        noval = []
        for binval in self.binlist:
            x = binval
            
            n = self.yeshash[binval] + self.nohash[binval]
            if n != 0:
                noval.append(self.nohash[binval])
                yesval.append(self.yeshash[binval])
                y = self.yeshash[binval] / n
                xxx.append(x)
                yyy.append(y)
                nnn.append(n)
            else:
                yesval.append(0)
                noval.append(0) 
                xxx.append(x)
                yyy.append(-0.1) 
                nnn.append(0)
        return xxx, yyy, nnn, yesval, noval

    def reliability_plot(self,rname='reliability'):
        xxx,yyy,nnn,yesval,noval = self.get_binvals()
        sns.set()
        sns.set_style("whitegrid")
        fig = plt.figure(1,figsize=(6,12))

        ax1 = fig.add_subplot(2,1,1)
        ax2 = fig.add_subplot(2,1,2)   
        sub_reliability_plot(self,ax1)
        sub_reliability_number_plot(self,ax2)
       
        plt.tight_layout()
        plt.savefig('{}.png'.format(rname))

        fig2 = plt.figure(2)
        ax3 = fig2.add_subplot(1,1,1)
        #with open('reliability.txt', 'w') as fid:
        #    fid.write(yesval)
        #    print(noval)
        #print('meas > thresh', yesval)
        #print('meas < thresh', noval)
        plt.plot(xxx, yesval, '-g.',label='Obs above')
        plt.plot(xxx, noval, '-r.',label='Obs below')
        plt.show()


def sub_reliability_number_plot(rc,ax, clr='-g.',fs=10,label=''):
    xxx,yyy,nnn,yesval,noval = rc.get_binvals()
    xx2 = []
    nn2 = []
    for xval in zip(xxx,nnn):
        if xval[1] >=0: 
           xx2.append(xval[0])
           nn2.append(xval[1])             
        
    ax.plot(xx2, nn2, clr,label=label)
    ax.set_yscale('log')
    ax.set_xlabel('Fraction model runs',fontsize=fs)
    ax.set_ylabel('Number of points',fontsize=fs)
    ax.tick_params(labelsize=fs)
    #ax.set_yticks([1,10,100,1000])
    #ax.set_yticklabels([1,10,100,1000,10000])
    return ax

def sub_reliability_plot(rc,ax, clr='-g.',fs=10,label=''):
    xxx,yyy,nnn,yesval,noval = rc.get_binvals()
    xx2 = []
    yy2 = []
    for xval in zip(xxx,yyy):
        if xval[1] >=0: 
           xx2.append(xval[0])
           yy2.append(xval[1])             
    ax.plot(xx2,yy2, clr,label=label)
    ax.plot([0,1], [0,1], '-k')
    yes_sum = np.array(yesval).sum()
    no_sum = np.array(noval).sum()
    climate = yes_sum / (yes_sum + no_sum)
    print(climate)
    ax.set_ylabel('Fraction observations',fontsize=fs)
    clr2 = clr.replace('-','--')[0:-1]
    print(clr2)
    ax.plot([0,1],[climate,climate], clr2 )
    ax.tick_params(labelsize=fs)
    return ax

#import matplotlib.pyplot as plt
#import textwrap
import datetime
import os
import sys
import numpy as np
import pandas as pd
from utilhysplit import par2conc



def write_dataA(dfin, c_col='mean', name='model.txt',thresh=1):
    mult=1e12
    model=True
    dataA=False
    data=False
    if dataA:
        reorder = ['Num','Site','Lat','Lon','Yr','Mo','Da','Hr','Meas','Calc']
        rename = ['Num', 'date','Lat','Lon','Site','Meas','mean','max','min','dur']
        droplist = ['date','dur','mean','max','min']
    if model or data:
        reorder = ['Yr','Mo','Da','Hr','dur','Lat','Lon','Calc','Site']
        rename = ['Num', 'date','Lat','Lon','Site','Meas','mean','max','min','dur']
        droplist = ['date','mean','max','min','Meas']
    if  data:
        reorder = ['Yr','Mo','Da','Hr','dur','Lat','Lon','Meas','Site']
        droplist = ['date','mean','max','min','Calc']
    df = dfin.copy()
    df.reset_index(inplace=True)
    df.columns = rename
    def apply_thresh(row):
        mult=1e12
        val = float(row) * mult
        if val < thresh: val=0
        return val
        
    df['Yr'] = df.apply(lambda row: row['date'].year, axis=1)
    df['Mo'] = df.apply(lambda row: row['date'].month, axis=1)
    df['Da'] = df.apply(lambda row: row['date'].day, axis=1)
    df['Hr'] = df.apply(lambda row: row['date'].hour*100, axis=1)
    df['Calc'] = df.apply(lambda row: apply_thresh(row[c_col]),axis=1)
    df = df.drop(droplist,axis=1)
    df = df[reorder]
    df.to_csv(name, index=False, sep=' ') 
    return df

#def par2stn(stndf,pdict,nnn=None,maxht=300,mlist=None,fit='all'):
    #if fit=='all':
#
#    else:
#       mfitlist,outdf = par2stn 
#       write_dataA(outdf)
#    return -1

def par2df(stndf,pardf,nnn=None,ht=10, maxht=300,dd=0.01,dh=0.01,
            buf=[0.05,0.05],mlist=None,method='gmm',method2='fit1'):
    udates = stndf.date.unique()
    udur = stndf.dur.unique()
    outdf = pd.DataFrame()
    iii=0
    mfitlist = []
    for date in udates:
        sdf2 = stndf[stndf['date']==date]        
        for dur in udur:
            sdf = sdf2[sdf2['dur']==dur]        
            if sdf.empty: 
               continue
            tmave = int(dur)/100 * 60
            #jjj, dfnew = combine_pdict(self.pdict, pd.to_datetime(date), tmave)
            d1 = pd.to_datetime(date)
            d2 = d1 + datetime.timedelta(minutes=tmave)
            pardfnew = pardf[pardf.date < d2]
            pardfnew = pardfnew[pardfnew.date >= d1]
            print(pardfnew.date.unique())
            #pdictnew = par2conc.subset_pdict(pdict, pd.to_datetime(date), tmave)
            if mlist: mval = mlist[iii]
            else: mval=None 
            if method2 == 'fit1':
               #fit each time period separately
               submlist = sub(pardfnew,nnn,maxht,mval,method)
            if method2 == 'fit2':
               #fit all particles in averaging time period together.
               submlist = sub2(pardfnew,nnn,maxht,mval,method)
            concdf = get_concdf(submlist, sdf, 
                               ht=ht,dd=dd,dh=dh,buf=buf)
            if not mlist: mfitlist.append(submlist)
            if iii==0: outdf = concdf
            else: 
               try:
                   outdf = pd.concat([outdf, concdf], axis=0) 
               except:
                   print('par2df: outdf', outdf)
                   print('par2df: concdf', concdf)
            iii+=1
    #outdf can be input into write_dataA
    return mfitlist, outdf             

def sub2(pardf,nnn,maxht, mlist=None, method='gmm'):
    # fit all particles in averaging time period at once.
    # mass needs to be divided by averaging time periods.
    massmult  = 1.0 / float(len(pardf.date.unique()))
    jjj=0
    submlist = []
    masslist = []
    pdn = pardf.cop()
    if maxht: pdn = pdn[pdn['ht']< maxht]
    mfit = par2conc.par2fit(pdn,mult=massmult,nnn=nnn, method=method)
    if not mlist: 
       mfit = par2conc.par2fit(pdn,nnn=nnn, method=method)
    else: 
       mfit = mlist[jjj]
    return submlist 

def sub(pardf, nnn, maxht, mlist=None,method='gmm'):
    jjj=0
    submlist = []
    #masslist = []
    #fit each unique date in the period.
    print('in sub')
    for ndate in pardf.date.unique(): 
        pdn = pardf[pardf.date==ndate]
        print('sub', pdn)
        if maxht: pdn = pdn[pdn['ht']< maxht]
        if not mlist:
           if jjj==0:
               mfit = par2conc.par2fit(pdn,nnn=nnn, method=method)
               print('sub first', mfit)
           else:
               mfit = par2conc.par2fit(pdn,nnn=nnn, method='p_bgm',pfit=pfit)
               print('sub mfit with previous', mfit)
        else: 
           mfit = mlist[jjj]
        if not mfit.fit: continue
        pfit = mfit.gfit
        #masslist.append(pdn['pmass'].sum())
        submlist.append(mfit)
        jjj+=1
    return submlist 


def get_concdf( mfitlist, stndf,  
               ht=50,
               dd=0.01,
               dh=0.01,
               buf=[0.05,0.05]):
    """
    ht should be input in meters.
    dd
    dh
    buf gives area around to calculated
    """
    ht = ht/1000.0
    dlist = []
    #dd = 0.01
    #dh = 0.05
    #dh = 0.01
    #buf = [0.05,0.01]
    measname = 'pmch'
    for row in stndf.itertuples(index=True, name='Pandas'): 
        phash={}
        date = getattr(row,'date')
        lat = getattr(row,'lat')
        lon = getattr(row,'lon')
        stn = getattr(row,'stn')
        meas = getattr(row,measname)
        dur = getattr(row,'dur')
        print('get_concdf buf', buf, 'dh', dh)
        concra = par2conc.average_mfitlist(mfitlist,
                          dd=dd,dh=dh,buf=buf,lat=lat,lon=lon,ht=ht)
        concra = par2conc.shift_underground(concra) 

        phash['date']=date
        phash['lat'] = lat                 
        phash['lon'] = lon                
        phash['stn'] = stn
        phash['meas'] = meas
        phash['mean'] = float(concra.mean())
        phash['max'] = float(concra.max())
        phash['min'] = float(concra.min())               
        phash['dur'] = dur
        dlist.append(phash)
    return pd.DataFrame(dlist)       
      


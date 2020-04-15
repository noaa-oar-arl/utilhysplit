#import matplotlib.pyplot as plt
#import textwrap
import datetime
import os
import sys
import numpy as np
import pandas as pd
from utilhysplit import parutils



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

def par2stn(stndf,pdict,nnn=None,maxht=300,mlist=None,fit='all'):
    #if fit=='all':
#
#    else:
#       mfitlist,outdf = par2stn 
#       write_dataA(outdf)
    return -1

def par2df( stndf,pdict,nnn=None,ht=10, maxht=300,dd=0.01,dh=0.01,
            buf=[0.05,0.05],mlist=None,method='gmm'):
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
            print('par2df: adding :', date, dur)
            tmave = int(dur)/100 * 60
            #jjj, dfnew = combine_pdict(self.pdict, pd.to_datetime(date), tmave)
            pdictnew = parutils.subset_pdict(pdict, pd.to_datetime(date), tmave)
            if mlist: mval = mlist[iii]
            else: mval=None 
            print('par2df buf', buf, 'dh', dh)
            masslist, submlist = sub(pdictnew,nnn,maxht,mval,method)
            concdf = get_concdf(submlist, sdf, masslist, 
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

def sub(pdictnew, nnn, maxht, mlist=None,method='gmm'):
    jjj=0
    submlist = []
    masslist = []
    for pdn in pdictnew: 
        if maxht: pdn = pdn[pdn['ht']< maxht]
        if not mlist: 
           mfit = parutils.par2fit(pdn,nnn=nnn, method=method)
        else: 
           mfit = mlist[jjj]
        if not mfit.fit: continue
        masslist.append(pdn['pmass'].sum())
        submlist.append(mfit)
        jjj+=1
    return masslist, submlist 


def par2df_b( stndf,pdict,nnn=None,  maxht=300,mlist=None):
    """
    stndf : Pandas dataframe. Should have columns
            date, dur, pmch
    """
    #sdf = stndf.sort_values(by=['date','dur'],axis=1)
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
                print('EMPTY', date, dur)
                continue
            print('adding', date, dur)
            tmave = int(dur)/100 * 60
            jjj, dfnew = parutils.combine_pdict(pdict, pd.to_datetime(date), tmave)
            if dfnew.empty: 
               print('dfnew empty')
               continue
            if maxht: dfnew = dfnew[dfnew['ht']< maxht]
            if not mlist:
               mfit = parutils.par2fit(dfnew,nnn=nnn)
               mfitlist.append(mfit)
            else:
               mfit = mlist[iii]
            mass = dfnew['pmass'].sum() / jjj 
            #mfitlist.append(mfit)
            concdf = get_concdf([mfit], sdf, mass,ht,dd,dh,buf)
            if iii==0: outdf = concdf
            else: outdf = pd.concat([outdf, concdf], axis=0) 
            iii+=1
    return mfitlist, outdf             

def get_concdf( mfitlist, stndf, mass, 
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
        concra = parutils.average_mfitlist(mfitlist,mass,
                          dd=dd,dh=dh,buf=buf,lat=lat,lon=lon,ht=ht)
        concra = parutils.shift_underground(concra) 

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
      


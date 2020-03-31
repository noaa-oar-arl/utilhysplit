import xarray as xr
#import matplotlib.pyplot as plt
#import textwrap
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import numpy as np
import numpy.ma as ma
import pandas as pd
import datetime
import seaborn as sns
import os
import sys
from Tools import volcMER
import utilhysplit.pardump as pardump
import matplotlib.pyplot as plt
#import reventador as reventador


class ParDat:
    """
    Read in and process a particle.dat file from STILT options. 
    """

    def __init__(self,fname, stime):
        """
        fname : str
        stime : datetime. date of beginning of run.
        
        """
        self.fname = fname
        self.stime = stime #start time of simulation
        self.keepvals = ['dtime','index','site','lat','lon','agl','grdht','foot']
        self.colnames = ['date','sorti','site','lat','lon','agl','grdht','foot']
        self.df = self.read()
        self.df = self.process(self.df)

    def read(self):
        df = pd.read_csv(self.fname, sep='\s+', header=[0])
        # time is minutes since start of run
        df['dtime'] = df.apply(lambda row: self.stime +
                                datetime.timedelta(minutes=row['time']),axis=1) 
        return df
 
    def process(self, df):
        df = df[self.keepvals]
        df.columns = self.colnames
        self.df = df
        return df



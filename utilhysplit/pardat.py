#import matplotlib.pyplot as plt
#import textwrap
import pandas as pd
import datetime
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
        self.keepvals =
                      ['dtime','index','site','lat','lon',
                       'agl','grdht','foot','dmass']
        self.colnames =
                      ['date','sorti','site','lat','lon',
                       'agl','grdht','foot','pmass']
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



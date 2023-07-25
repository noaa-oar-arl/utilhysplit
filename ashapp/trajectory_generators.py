import numpy as np
import pandas as pd

"""
Used as input into RunTrajectory class.
the trajectory generators output a dictionary with height, latitude, longitude.
"""


def generate_traj_from_config(inp):
    outp={}
    lat = inp["latitude"]
    lon = inp["longitude"]
    height = inp["height"]
    if not isinstance(height,(list,np.ndarray)):
       height = [height]
    for hgt in height:
        outp['height'] = hgt
        outp['latitude'] = lat
        outp['longitude'] = lon
        yield outp

def generate_traj_from_df(df):
    for rrr in df.iterrows():
        outp = {}
        row = rrr[1]
        # this is specifically for file from hunga tonga so2 data.
        # heightI uses the interpolated data as well.
        outp["height"] = row.height * 1000
        outp["latitude"] = row.lat
        outp["longitude"] = row.lon
        yield outp 

def generate_traj_from_obsdf(csvname):
    """
    generate back trajectories from a csv file.
  
    RETURNS
    tuple (time, generate_traj_from_df function)
    """
    obsdf = pd.read_csv(csvname, parse_dates=["time"])
    #return obsdf
    outp = {}
    # one trajectory run per time period.
    timelist = obsdf['time'].unique()
    for time in timelist:
        newdf = obsdf[obsdf['time']==time]
        time = pd.to_datetime(time)
        yield time, generate_traj_from_df(newdf)


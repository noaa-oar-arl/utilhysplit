import numpy as np
import pandas as pd
import xarray as xr
import datetime
import os
import math
import glob
import monet
import matplotlib.pyplot as plt
import seaborn as sns
from natsort import natsorted
from monetio.models import hysplit
from monetio.models import hytraj
from ashapp.backtraj import combine_traj
 
import sys


"""
Bavand Sadeghi:

ABSTRACT: functions for calculating the nearest point of HYPLIT Back Trajectories to a volcano vent. The code reads the resolved trajectories, caclulates the distances of every point of the trajectory output to the volcano and finds the characteristics of the nearest point of each trajectory to the volcano. The characteristics of the nearest point will be presented further below for each function. 


"""
def trajvold(volcano, df):
    """
    The function calculates the distances of each trajectory (at each time step) to the volcano.
    
    Inputs
    volcano: array
             lat and lon of the volcano
    df: dataframe
             characteristics of the trajectory

    Outputs
    dist : array
           values of the distances between each point of the trajectory and the volcano     

    """
    dist = []
    for i in range(len(df)):
        span = (df.latitude[i]-volcano[0])**2 + (df.longitude[i]-volcano[1])**2
        span = math.sqrt(span)
        dist = np.append(dist, span)
        dist = np.array(dist)
    return dist

def trajdisthgt(tnames,obs_path,volcano):
    """
    The function calculates the plume thickness through a series of trajectories coming close to the volcano

    Inputs
    tnames: name and path of resolved back trajectory (list)
    obs_path: dataframe
              contains the characteristics of the data measurement  
    volcano:  array
              lat and lon of the volcano

    outputs:
    characteristics of the nearest point of the resolved trajectories to the volcano
    the names are abbreviated:
        dist_len: array
                  the distance length between the closest point of the trajectory and the volcano
        dist_hgt: array
                  height at which the trajectories come closest to the volcano
        dist_lat, dist_lon: array
                  the coordinates of the nearest points to the volcano eruption
        dist_time:
                  the time when the nearest point reaches the volcano
        obs_"variables":
                  the locations of the back trajectory starting points (coordinates and altitude)
    """
    (dist_len, dist_hgt, dist_lat, dist_lon, dist_time, obs_time, obs_lat,
    obs_lon, init_alt, obs_height, dist_weight) = ([] for i in range(11))
    for traj_no in range(500):

        df = combine_traj([tnames[traj_no]], csvfile = obs_path).copy()
        df['init_alt'] = df['altitude'].iloc[0]
        obs_path_csvfile = pd.read_csv(obs_path)
        df['heightI_obs'] = (obs_path_csvfile['heightI_obs'].iloc[0])*1000
        df.loc[df['longitude'] < 0, 'longitude'] = df['longitude'] + 360
        # select the trajectory points within the first 12 hours after the eruption
        df2 = df.loc[df['time'].between('2022-01-15 04:00:00', '2022-01-15 16:00:00')]
        df2 = df2.reset_index(level=None, drop=True, inplace=False)
        dist_out = trajvold(volcano, df2)
        closer_dist_value = np.min(dist_out)
        closer_dist_index = np.argmin(dist_out)

        dist_len = np.append(dist_len, dist_out[closer_dist_index])
        dist_hgt = np.append(dist_hgt, df2.altitude[closer_dist_index])
        dist_lat = np.append(dist_lat, df2.latitude[closer_dist_index])
        dist_lon = np.append(dist_lon, df2.longitude[closer_dist_index])
        dist_time  = np.append(dist_time, df2.time[closer_dist_index])
        obs_time = np.append(obs_time, df.time[0])
        obs_lat = np.append(obs_lat, df.latitude[0])
        obs_lon = np.append(obs_lon, df.longitude[0])
        init_alt = np.append(init_alt, df2.init_alt[closer_dist_index])
        obs_height = np.append(obs_height, df2.heightI_obs[closer_dist_index])
        dist_weight = np.append(dist_weight, df2.weight[closer_dist_index])

    return (dist_len, dist_hgt, dist_lat, dist_lon, dist_time, obs_time, obs_lat, obs_lon, init_alt,
            obs_height, dist_weight)


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
Classes to read output files from the HYSPLIT trajectory. 

"""
# Function of distances between trajectories and volcano
def trajvold(volcano, df):
    dist = []
    for i in range(len(df)):
        span = (df.latitude[i]-volcano[0])**2 + (df.longitude[i]-volcano[1])**2
        span = math.sqrt(span)
        dist = np.append(dist, span)
        dist = np.array(dist)
    return dist

# Function of calculations of the thickness of ash cloud
def trajdisthgt(tnames,obs_path,volcano):
    (dist_len, dist_hgt, dist_lat, dist_lon, dist_time, obs_time, obs_lat,
    obs_lon, init_alt, obs_height, dist_weight) = ([] for i in range(11))
    for traj_no in range(500):

        # read the trajectory file (HYSPLIT output)
        df = combine_traj([tnames[traj_no]], csvfile = obs_path).copy()
        # altitude of trajectory starting point
        df['init_alt'] = df['altitude'].iloc[0]
        # measured altitude of observation point
        obs_path_csvfile = pd.read_csv(obs_path)
        df['heightI_obs'] = (obs_path_csvfile['heightI_obs'].iloc[0])*1000
        # adjust the coordinates of the points with negative longitude
        df.loc[df['longitude'] < 0, 'longitude'] = df['longitude'] + 360
        # filter the trajectory points for the 12 hours period after the volcano eruption
        df2 = df.loc[df['time'].between('2022-01-15 04:00:00', '2022-01-15 16:00:00')]
        df2 = df2.reset_index(level=None, drop=True, inplace=False)
        # find the distances between the trajectories and volcano, nearest point, and its index
        dist_out = trajvold(volcano, df2)
        closer_dist_value = np.min(dist_out)
        closer_dist_index = np.argmin(dist_out)

        # dist_len, the distance length between the closest point of the trajectory and the volcano
        dist_len = np.append(dist_len, dist_out[closer_dist_index])
        # dist_hgt; height at which the trajectories come closest to the volcano
        dist_hgt = np.append(dist_hgt, df2.altitude[closer_dist_index])
        # dist_lat and dist_lon; the coordinates of the nearest points to the volcano
        dist_lat = np.append(dist_lat, df2.latitude[closer_dist_index])
        dist_lon = np.append(dist_lon, df2.longitude[closer_dist_index])
        # dist_time; the time when the nearest point reaches the volcano
        dist_time  = np.append(dist_time, df2.time[closer_dist_index])
        # the locations of the back trajectory starting points (coordinates and altitude)
        obs_time = np.append(obs_time, df.time[0])
        obs_lat = np.append(obs_lat, df.latitude[0])
        obs_lon = np.append(obs_lon, df.longitude[0])
        init_alt = np.append(init_alt, df2.init_alt[closer_dist_index])
        obs_height = np.append(obs_height, df2.heightI_obs[closer_dist_index])
        dist_weight = np.append(dist_weight, df2.weight[closer_dist_index])

    return (dist_len, dist_hgt, dist_lat, dist_lon, dist_time, obs_time, obs_lat, obs_lon, init_alt,
            obs_height, dist_weight)


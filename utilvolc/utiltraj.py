import pandas as pd
import os
import datetime
import math
import sys
import glob
import datetime
import monet
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import numpy.ma as ma
import xarray as xr
import shapely.geometry as sgeo
from shapely.ops import nearest_points


#from natsort import natsorted
from math import pi
from monetio.models import hytraj
from monetio.models import hysplit
#from utilvolc.utiltraj import combine_traj
from utilhysplit import emitimes
from utilvolc import volcat
from utilvolc.volcat import VolcatName
from utilvolc import get_area
from utilhysplit import geotools

from utilvolc.make_data_insertion import make_1D_sub, EmitName

# 2023 Dec 04 (amc)  added make_btraj_filename function
# 2023 Dec 04 (amc)  modified combine_traj to handle tdump files with multiple trajectories.
# 2023 Dec 04 (amc)  modified combine_traj to add observed time, altitude, and area to file
# 2023 Dec 04 (amc)  modified trajectory_input_csv function for different naming convention.
# 2023 Dec 04 (amc)  added modified functions for computing distance from trajectory to volcano.
# 2023 Dec 04 (amc)  removed dependence on natsorted




def make_btraj_filename(volcat_fname):
    suffix = ''
    enn = EmitName(None)
    return enn.make_filename(volcat_fname, prefix='btraj', suffix=suffix) + '.csv'


def trajectory_input_csv(dataset, data_dir, layer_height=None):
    """
    Functions reads the xarray datasets and generates csv files which can be used by ash_main.py.
    Input:
        xarray dataset representing volcat data.
    Output:
        csv file used by ash_main.py to create a set of back trajectory runs from the observation 
        points. 
    """
    obs_data_orig = pd.DataFrame(make_1D_sub(dataset), columns=['lat','lon','mass','height','area'])
    obs_data_orig['time'] = dataset.time.values
    if layer_height:
       obs_data_orig["heightI"] = layer_height
       fname = f'btraj{"%02d" %layer_height}km.csv'
    else:
       fname = make_btraj_filename(dataset.dataset_name)
    out = os.path.join(data_dir,fname)
    print('writing', out)
    obs_data_orig.to_csv(out, index = False)
    return obs_data_orig


def combine_traj(fnames, csvfile=None):
    """
    fnames  : list of str. trajectory file names. full path.
    csvfile : csv file output by sample_and_write which contains weighting information.
    combined trajectories in different files into one dataframe.
    """

    # 01 November 2023. changed to handle tdump files which may have multiple trajectories.
    #                   do not over-write the traj_num column. instead add a run_num column.

    #                   Also add the observed altitude and area to the new dataframe.
    trajlist = []
    if csvfile:
        weightcsv = pd.read_csv(csvfile)
    for iii, fnn in enumerate(fnames):
        try:
            df1 = hytraj.open_dataset(fnn)
        except Exception as eee:
            print('Failed to open {}'.format(fnn))
            print(eee)
           
        # get trajectory number from the file name
        temp = fnn.split(".")
        trajnum = int(temp[-1])
        # add new column to dataframe with run number
        # the run may have multiple trajectories in it.
        #df1["traj_num"] = trajnum
        df1["run_num"] = trajnum
        #print('TRAJNUM', trajnum)
        # add weight information from csvfile to the dataframe
        if csvfile:
            temp = weightcsv.loc[trajnum]
        #    weight = temp.massI
            weight = temp.mass
            obsalt = temp.height
            obsarea = temp.area
        else:
            weight = 1
            obsalt = None
            obsarea = None
        df1["weight"] = weight
        df1["obsalt"] = obsalt
        df1["area"] = obsarea
 
        trajlist.append(df1.copy())
    # concatenate the trajectories into one dataframe.
    trajdf = pd.concat(trajlist)
    return trajdf


def read_traj_output(data_dir, fname, num_layer):
    """
    Function reads the trajectories and measured observations. These paths of the trajectories will
    be inputs for calculations of plume heights
    """
    tnames_path = []
    obs_path = []
    trajectory_data = {}
    for ii in range(1,num_layer+1):
        tdump_files = glob.glob(data_dir + f'{fname}_{"%02d" %ii}km/' + 'tdump.ashtest_btraj*')
        trajectory_data[f'tdump{ii}'] = tdump_files
        tnames_path.append(trajectory_data[f'tdump{ii}'])
        obs_path.append(data_dir + f'{fname}_{"%02d" %ii}km/' + f'btraj{"%02d" %ii}km.csv')
    return (tnames_path, obs_path)


def traj_volc_dist2(vloc, df):
    # shapely only does computations in 2d. 
    # although you can define point and line with a z coordinate.
    vpoint = sgeo.Point((vloc[1],vloc[0]))
    x=df.longitude
    y=df.latitude
    xy=list(zip(x,y))
    # creates line segments from trajectory points.
    tline = sgeo.LineString(xy)
    # finds closest point on the linestring to volcano.
    # do this because it may be between trajectory points.
    a,b = nearest_points(tline,vpoint)
    # returns distance in km
    distance = geotools.distance(a,b)
    return a, distance


def traj_volc_dist_km(vloc, df):
    from utilhysplit.geotools import calculate_distance
    """
    The function calculates the distances between each trajectory (at each time step) and the volcano vent.
    Similar to traj_volc_dist but output is in km.
    Inputs:
        volcano: array; lat and lon of the volcano
        df: dataframe, trajectory characteristics
    Outputs:
        dist : array; values of the distances between each point of the trajectory and the volcano.
    """
    distlist = []
    for iii, row in df.iterrows():
        dist = calculate_distance(row.latitude,row.longitude,vloc[0],vloc[1])
        distlist.append(dist)
        dist = np.array(distlist)
    return dist

def traj_volc_dist(vloc, df):
    """
    The function calculates the distances between each trajectory (at each time step) and the volcano vent.
    Inputs:
        volcano: array; lat and lon of the volcano
        df: dataframe, trajectory characteristics
    Outputs:
        dist : array; values of the distances between each point of the trajectory and the volcano.
    """
    dist = []
    deg2km = 111.111
    for i in range(len(df)):
    # distance calculation to meters
        span = ((df.latitude[i] - vloc[0])*deg2km)**2 + ((df.longitude[i] - vloc[1])*deg2km)**2
        span = math.sqrt(span)
        dist = np.append(dist, span)
        dist = np.array(dist)
    return dist


def traj_layer_cal(tnames_path,obs_path,vloc):
    """
    The function calculates the characteristics of trajectories coming close to the volcano over entire observations.
    It performs calculations for the number of measurements.

    Inputs:
        tnames: name and path of resolved back trajectory (list)
        obs_path: dataframe; contains the characteristics of the data measurement
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
    for traj_no in range(716):

        # read the trajectory file (HYSPLIT output)
        df = combine_traj([tnames_path[traj_no]], csvfile = obs_path).copy()
        # altitude of trajectory starting point
        df['init_alt'] = df['altitude'].iloc[0]
        # measured altitude of observation point
        obs_path_csvfile = pd.read_csv(obs_path)
        df['height'] = (obs_path_csvfile['height'].iloc[traj_no])*1000
        df.loc[df['longitude'] < 0, 'longitude'] = df['longitude'] + 360
        # filter the trajectory points for the 12 hours period after the volcano eruption. check for each volcano.
        #df2 = df.loc[df['time'].between('2022-01-15 04:00:00', '2022-01-15 16:00:00')]
        df2 = df
        df2 = df2.reset_index(level=None, drop=True, inplace=False)
        # find the distances between the trajectories and volcano, nearest point, and its index
        dist_out = traj_volc_dist(vloc, df2)
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
        obs_height = np.append(obs_height, df2.height[closer_dist_index])
        dist_weight = np.append(dist_weight, df2.weight[closer_dist_index])

    return (dist_len, dist_hgt, dist_lat, dist_lon, dist_time, obs_time, obs_lat, obs_lon, init_alt,
            obs_height, dist_weight)



def traj_layer_dataframe(tnames_path, obs_path, vloc):
    """
    This function outputs the characteristics of the trajectories' nearest point to the volcano.
    """
    var1 = {}

    for i in range(0, 30):
        print(obs_path[i])

        dist_func = traj_layer_cal(tnames_path[i], obs_path[i],vloc)

        dist_len    = dist_func[0]
        dist_hgt    = dist_func[1]
        dist_lat    = dist_func[2]
        dist_lon    = dist_func[3]
        dist_time   = dist_func[4]
        obs_time    = dist_func[5]
        obs_lat     = dist_func[6]
        obs_lon     = dist_func[7]
        init_alt    = dist_func[8]
        obs_height  = dist_func[9]
        dist_weight = dist_func[10]

        obs_data = {'dist_len': dist_len, 'dist_hgt': dist_hgt, 'dist_lat': dist_lat, 'dist_lon': dist_lon,
                    'dist_time': dist_time, 'obs_time': obs_time, 'obs_lat':obs_lat, 'obs_lon': obs_lon,
                    'init_alt': init_alt, 'obs_height': obs_height, 'dist_weight': dist_weight}
        var1[i] = obs_data

    df = {}
    for ii in range(716):
        df_variable = pd.DataFrame()
        for j in range(0,30):
            df_variable = pd.concat([df_variable, pd.DataFrame(var1[j])[ii:ii+1]])
            df[ii] = df_variable

    return (df)



def plume_thick_cal(df):
    """
    This function selects the layers at which the shortest distance between trajectories and volcano points occurs.
    It sets a cut-off value and identifies the base and top of the cloud, allowing the code to calculate the cloud thickness.
    Outputs: a dataframe containing the cloud characteristics (columns) for entire observations (rows)
    """
    df_dist_min = pd.DataFrame()
    df_dist_min_criteria_pass = pd.DataFrame()
    df_closest_row = pd.DataFrame()
    cutoff = 0.08
    critera_pass_traj_num = 0

    for i in range(716):

        # finds the layers where the trajectories come within a certain distance to the volcano
        selected_rows = df[i].loc[df[i]['dist_len'] < cutoff]

        # If there was no single layer where the trajectory comes within this threshold distance of the vent, the code will choose
        # the initial altitude of the nearest trajectory layer and set it as the top of the cloud.
        if selected_rows.empty:
            df_closest_row = pd.DataFrame(df[i].iloc[df[i]['dist_len'].argmin()]).transpose()
            df_closest_row['cloud_top'] = df_closest_row['init_alt']
            df_closest_row['cloud_bottom'] = df_closest_row['cloud_top'] - 1000
            df_closest_row['thickness'] = 1000
            df_dist_min = pd.concat([df_dist_min, df_closest_row])

        # For the trajectories that come within the threshold, the code determine both a bottom and top level for the plume
        # and display the number of them multiplied by 1 km as the thickness.
        else:
            selected_values = selected_rows['init_alt'].values
            range_of_values = (np.min(selected_values), np.max(selected_values))
            df_closest_row = selected_rows.nsmallest(1, 'dist_len')
            df_closest_row['cloud_top'] = np.max(selected_values)
            df_closest_row['cloud_bottom'] = np.min(selected_values)
            df_closest_row['thickness'] = len(selected_rows)*1000 # np.max(selected_values) - np.min(selected_values)
            df_closest_row['thickness'] = df_closest_row['thickness'].clip(lower=1000)
            df_dist_min = pd.concat([df_dist_min, df_closest_row])
            df_dist_min_criteria_pass = pd.concat([df_dist_min_criteria_pass, df_closest_row])
            critera_pass_traj_num = critera_pass_traj_num + 1
    print(critera_pass_traj_num)
    return(df_dist_min, df_dist_min_criteria_pass, critera_pass_traj_num)


def conc_emitimes_data(df):
    """"
    This function prepares the EMITIMES input file for Hysplit Concentration Calculations.
    """
    df_dist_fwd_data = pd.DataFrame()
    df_dist_fwd = {}
    dg_dist_fwd = {}
    cutoff = 0.08
    critera_pass_traj_num = 0
    for i in range(716):
        selected_rows = df[i].loc[df[i]['dist_len'] < cutoff]
        if selected_rows.empty:
            df_dist_fwd[i] = pd.DataFrame(df[i].iloc[df[i]['dist_len'].argmin()]).transpose()
            print('length of the variable here is', len(df_dist_fwd[i]))
            # convert the DU to g/m^2
            #df_dist_fwd[i]['tload'] = (df_dist_fwd[i]['dist_weight']) * (2.6867E20) * (64.066) * (1/(6.022E23))
            df_dist_fwd[i].loc[:, 'tload'] = (df_dist_fwd[i]['dist_weight']) * (2.6867E20) * (64.066) * (1/(6.022E23))
            # convert the g/m^2 to g (7.0 km ~ 7000 m)
            # df_dist_fwd[i]['tload'] = df_dist_fwd[i]['tload'] * 49000000
            # will add a new column showing the number of source emission layers at the point
            df_dist_fwd[i] = df_dist_fwd[i].assign(count = len(df_dist_fwd[i]))
            # will add a new column for the emission rate (1/hr): df_dist_fwd[i]['tload'] * (1/[emission time (hr)])
            df_dist_fwd[i]['rate'] = df_dist_fwd[i]['tload'] / (1/12)
            df_dist_fwd[i]['YYYY'] = df_dist_fwd[i]['obs_time'].dt.year
            df_dist_fwd[i]['MM']   = df_dist_fwd[i]['obs_time'].dt.month
            df_dist_fwd[i]['DD']   = df_dist_fwd[i]['obs_time'].dt.day
            df_dist_fwd[i]['HH']   = df_dist_fwd[i]['obs_time'].dt.hour
            # df_dist_fwd[i]['dist_weight'] = df_dist_fwd[i]['dist_weight']/len(df_dist_fwd[i])
            df_dist_fwd_data = pd.concat([df_dist_fwd_data, df_dist_fwd[i]])

        else:
            dg_dist_fwd[i] = selected_rows
            df_dist_fwd[i] = selected_rows
            # convert the DU to g/m^2
            # df_dist_fwd[i]['tload'] = (df_dist_fwd[i]['dist_weight']) * (2.6867E20) * (64.066) * (1/(6.022E23))
            df_dist_fwd[i].loc[:, 'tload'] = (df_dist_fwd[i]['dist_weight']) * (2.6867E20) * (64.066) * (1/(6.022E23))
            # convert the g/m^2 to g (7.0 km ~ 7000 m)
            # df_dist_fwd[i]['tload'] = df_dist_fwd[i]['tload'] * 49000000
            df_dist_fwd[i].loc[:, 'tload'] = df_dist_fwd[i]['tload'] / len(df_dist_fwd[i])
            df_dist_fwd[i] = df_dist_fwd[i].assign(count = len(df_dist_fwd[i]))
            # will add a new column for the emission rate (1/hr): df_dist_fwd[i]['tload'] * (1/[emission time (hr)])
            df_dist_fwd[i]['rate'] = df_dist_fwd[i]['tload'] / (1/12)
            df_dist_fwd[i]['YYYY'] = df_dist_fwd[i]['obs_time'].dt.year
            df_dist_fwd[i]['MM']   = df_dist_fwd[i]['obs_time'].dt.month
            df_dist_fwd[i]['DD']   = df_dist_fwd[i]['obs_time'].dt.day
            df_dist_fwd[i]['HH']   = df_dist_fwd[i]['obs_time'].dt.hour
            # print(df_dist_fwd[i])
            df_dist_fwd_data = pd.concat([df_dist_fwd_data, df_dist_fwd[i]])
    return(df_dist_fwd_data)




def read_traj_output(data_dir, fname, num_layer):
    """
    Function reads the trajectories and measured observations. These paths of the trajectories will
    be inputs for calculations of plume heights
    """
    tnames_path = []
    obs_path = []
    trajectory_data = {}
    for ii in range(1,num_layer+1):
        tdump_files = glob.glob(data_dir + f'{fname}_{"%02d" %ii}km/' + 'tdump.ashtest_btraj*')
        trajectory_data[f'tdump{ii}'] = tdump_files
        tnames_path.append(trajectory_data[f'tdump{ii}'])
        obs_path.append(data_dir + f'{fname}_{"%02d" %ii}km/' + f'btraj{"%02d" %ii}km.csv')
    return (tnames_path, obs_path)




def traj_layer_cal(tnames_path,obs_path,vloc):
    """
    The function calculates the characteristics of trajectories coming close to the volcano over entire observations.
    It performs calculations for the number of measurements.

    Inputs:
        tnames: name and path of resolved back trajectory (list)
        obs_path: dataframe; contains the characteristics of the data measurement
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
    for traj_no in range(716):

        # read the trajectory file (HYSPLIT output)
        df = combine_traj([tnames_path[traj_no]], csvfile = obs_path).copy()
        # altitude of trajectory starting point
        df['init_alt'] = df['altitude'].iloc[0]
        # measured altitude of observation point
        obs_path_csvfile = pd.read_csv(obs_path)
        df['height'] = (obs_path_csvfile['height'].iloc[traj_no])*1000
        df.loc[df['longitude'] < 0, 'longitude'] = df['longitude'] + 360
        # filter the trajectory points for the 12 hours period after the volcano eruption. check for each volcano.
        #df2 = df.loc[df['time'].between('2022-01-15 04:00:00', '2022-01-15 16:00:00')]
        df2 = df
        df2 = df2.reset_index(level=None, drop=True, inplace=False)
        # find the distances between the trajectories and volcano, nearest point, and its index
        dist_out = traj_volc_dist(vloc, df2)
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
        obs_height = np.append(obs_height, df2.height[closer_dist_index])
        dist_weight = np.append(dist_weight, df2.weight[closer_dist_index])

    return (dist_len, dist_hgt, dist_lat, dist_lon, dist_time, obs_time, obs_lat, obs_lon, init_alt,
            obs_height, dist_weight)



def traj_layer_dataframe(tnames_path, obs_path, vloc):
    """
    This function outputs the characteristics of the trajectories' nearest point to the volcano.
    """
    var1 = {}

    for i in range(0, 30):
        print(obs_path[i])

        dist_func = traj_layer_cal(tnames_path[i], obs_path[i],vloc)

        dist_len    = dist_func[0]
        dist_hgt    = dist_func[1]
        dist_lat    = dist_func[2]
        dist_lon    = dist_func[3]
        dist_time   = dist_func[4]
        obs_time    = dist_func[5]
        obs_lat     = dist_func[6]
        obs_lon     = dist_func[7]
        init_alt    = dist_func[8]
        obs_height  = dist_func[9]
        dist_weight = dist_func[10]

        obs_data = {'dist_len': dist_len, 'dist_hgt': dist_hgt, 'dist_lat': dist_lat, 'dist_lon': dist_lon,
                    'dist_time': dist_time, 'obs_time': obs_time, 'obs_lat':obs_lat, 'obs_lon': obs_lon,
                    'init_alt': init_alt, 'obs_height': obs_height, 'dist_weight': dist_weight}
        var1[i] = obs_data

    df = {}
    for ii in range(716):
        df_variable = pd.DataFrame()
        for j in range(0,30):
            df_variable = pd.concat([df_variable, pd.DataFrame(var1[j])[ii:ii+1]])
            df[ii] = df_variable

    return (df)



def plume_thick_cal(df):
    """
    This function selects the layers at which the shortest distance between trajectories and volcano points occurs.
    It sets a cut-off value and identifies the base and top of the cloud, allowing the code to calculate the cloud thickness.
    Outputs: a dataframe containing the cloud characteristics (columns) for entire observations (rows)
    """
    df_dist_min = pd.DataFrame()
    df_dist_min_criteria_pass = pd.DataFrame()
    df_closest_row = pd.DataFrame()
    cutoff = 0.08
    critera_pass_traj_num = 0

    for i in range(716):

        # finds the layers where the trajectories come within a certain distance to the volcano
        selected_rows = df[i].loc[df[i]['dist_len'] < cutoff]

        # If there was no single layer where the trajectory comes within this threshold distance of the vent, the code will choose
        # the initial altitude of the nearest trajectory layer and set it as the top of the cloud.
        if selected_rows.empty:
            df_closest_row = pd.DataFrame(df[i].iloc[df[i]['dist_len'].argmin()]).transpose()
            df_closest_row['cloud_top'] = df_closest_row['init_alt']
            df_closest_row['cloud_bottom'] = df_closest_row['cloud_top'] - 1000
            df_closest_row['thickness'] = 1000
            df_dist_min = pd.concat([df_dist_min, df_closest_row])

        # For the trajectories that come within the threshold, the code determine both a bottom and top level for the plume
        # and display the number of them multiplied by 1 km as the thickness.
        else:
            selected_values = selected_rows['init_alt'].values
            range_of_values = (np.min(selected_values), np.max(selected_values))
            df_closest_row = selected_rows.nsmallest(1, 'dist_len')
            df_closest_row['cloud_top'] = np.max(selected_values)
            df_closest_row['cloud_bottom'] = np.min(selected_values)
            df_closest_row['thickness'] = len(selected_rows)*1000 # np.max(selected_values) - np.min(selected_values)
            df_closest_row['thickness'] = df_closest_row['thickness'].clip(lower=1000)
            df_dist_min = pd.concat([df_dist_min, df_closest_row])
            df_dist_min_criteria_pass = pd.concat([df_dist_min_criteria_pass, df_closest_row])
            critera_pass_traj_num = critera_pass_traj_num + 1
    print(critera_pass_traj_num)
    return(df_dist_min, df_dist_min_criteria_pass, critera_pass_traj_num)


def conc_emitimes_data(df):
    """"
    This function prepares the EMITIMES input file for Hysplit Concentration Calculations.
    """
    df_dist_fwd_data = pd.DataFrame()
    df_dist_fwd = {}
    dg_dist_fwd = {}
    cutoff = 0.08
    critera_pass_traj_num = 0
    for i in range(716):
        selected_rows = df[i].loc[df[i]['dist_len'] < cutoff]
        if selected_rows.empty:
            df_dist_fwd[i] = pd.DataFrame(df[i].iloc[df[i]['dist_len'].argmin()]).transpose()
            print('length of the variable here is', len(df_dist_fwd[i]))
            # convert the DU to g/m^2
            #df_dist_fwd[i]['tload'] = (df_dist_fwd[i]['dist_weight']) * (2.6867E20) * (64.066) * (1/(6.022E23))
            df_dist_fwd[i].loc[:, 'tload'] = (df_dist_fwd[i]['dist_weight']) * (2.6867E20) * (64.066) * (1/(6.022E23))
            # convert the g/m^2 to g (7.0 km ~ 7000 m)
            # df_dist_fwd[i]['tload'] = df_dist_fwd[i]['tload'] * 49000000
            # will add a new column showing the number of source emission layers at the point
            df_dist_fwd[i] = df_dist_fwd[i].assign(count = len(df_dist_fwd[i]))
            # will add a new column for the emission rate (1/hr): df_dist_fwd[i]['tload'] * (1/[emission time (hr)])
            df_dist_fwd[i]['rate'] = df_dist_fwd[i]['tload'] / (1/12)
            df_dist_fwd[i]['YYYY'] = df_dist_fwd[i]['obs_time'].dt.year
            df_dist_fwd[i]['MM']   = df_dist_fwd[i]['obs_time'].dt.month
            df_dist_fwd[i]['DD']   = df_dist_fwd[i]['obs_time'].dt.day
            df_dist_fwd[i]['HH']   = df_dist_fwd[i]['obs_time'].dt.hour
            # df_dist_fwd[i]['dist_weight'] = df_dist_fwd[i]['dist_weight']/len(df_dist_fwd[i])
            df_dist_fwd_data = pd.concat([df_dist_fwd_data, df_dist_fwd[i]])

        else:
            dg_dist_fwd[i] = selected_rows
            df_dist_fwd[i] = selected_rows
            # convert the DU to g/m^2
            # df_dist_fwd[i]['tload'] = (df_dist_fwd[i]['dist_weight']) * (2.6867E20) * (64.066) * (1/(6.022E23))
            df_dist_fwd[i].loc[:, 'tload'] = (df_dist_fwd[i]['dist_weight']) * (2.6867E20) * (64.066) * (1/(6.022E23))
            # convert the g/m^2 to g (7.0 km ~ 7000 m)
            # df_dist_fwd[i]['tload'] = df_dist_fwd[i]['tload'] * 49000000
            df_dist_fwd[i].loc[:, 'tload'] = df_dist_fwd[i]['tload'] / len(df_dist_fwd[i])
            df_dist_fwd[i] = df_dist_fwd[i].assign(count = len(df_dist_fwd[i]))
            # will add a new column for the emission rate (1/hr): df_dist_fwd[i]['tload'] * (1/[emission time (hr)])
            df_dist_fwd[i]['rate'] = df_dist_fwd[i]['tload'] / (1/12)
            df_dist_fwd[i]['YYYY'] = df_dist_fwd[i]['obs_time'].dt.year
            df_dist_fwd[i]['MM']   = df_dist_fwd[i]['obs_time'].dt.month
            df_dist_fwd[i]['DD']   = df_dist_fwd[i]['obs_time'].dt.day
            df_dist_fwd[i]['HH']   = df_dist_fwd[i]['obs_time'].dt.hour
            # print(df_dist_fwd[i])
            df_dist_fwd_data = pd.concat([df_dist_fwd_data, df_dist_fwd[i]])
    return(df_dist_fwd_data)





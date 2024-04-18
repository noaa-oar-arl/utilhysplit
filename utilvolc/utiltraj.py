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
from utilhysplit.fixlondf import fixlondf
from utilvolc import volcat
from utilvolc import volcat_so2
from utilvolc.volcat import VolcatName
from utilvolc import get_area
from utilhysplit import geotools
from utilvolc.make_data_insertion import make_1D_sub, EmitName

from sklearn.cluster import KMeans

"""
This function reads the VOLCAT data and performs K-Means clustering that partition a dataset into "K"
non-overlapping clusters.
"""

def volcat_kmeans_clustering(vname, num_samples, time_obs_threshold):
    """
    vname: the path of the VOLCAT retrieval data.
    num_samples: number of samples that we seek to transform the total observations into. # of clusters.
    time_obs_threshold: minimum count of measurements for each time period.
    """
    # read data and preprocess
    vvolc = volcat_so2.volcatSO2L3(fname=vname)
    vframe = fixlondf(vvolc.points2frame(), colname='lon',neg=False)
    vframe_df = vframe[['lat', 'lon', 'massI', 'heightI', 'time']]
    vframe_df.rename(columns={"massI": "mass", "heightI": "height"}, inplace = True)
    vframe_df = vframe_df.dropna(subset=['height','mass'])
    vframe_df = vframe_df.sort_values(by='time', ascending=True)

    # count the number of observations at each distinct time
    time_counts = vframe_df['time'].value_counts()

    # exclude times with observation counts below the threshold
    exclude_times = [str(time) for time in time_counts[time_counts < time_obs_threshold].index]
    df_filtered = vframe_df[~vframe_df['time'].isin(exclude_times)]

    sample_data = pd.DataFrame(columns=['lat', 'lon', 'mass', 'height', 'time'])

    # perform clustering
    for distinct_time in df_filtered['time'].unique():
        df_filtered_time = df_filtered[df_filtered['time'] == distinct_time]
        num_samples_per_time = math.floor((len(df_filtered_time)/len(df_filtered)) * num_samples)
        # skip if the number of samples is below 1
        if num_samples_per_time < 1:
            continue

        # extract features for clustering and perform K-means clustering
        features = df_filtered_time[['lat','lon','height']].values
        kmeans = KMeans(n_clusters=num_samples_per_time, random_state=42)
        cluster_labels = kmeans.fit_predict(features)
        cluster_centers = kmeans.cluster_centers_
        cluster_counts = np.bincount(cluster_labels)

        # extract cluster information and calculate sample data
        for cluster_id in range(num_samples_per_time):
            cluster_points = df_filtered_time[cluster_labels == cluster_id]
            total_mass = cluster_points['mass'].sum()

            weighted_lat = cluster_centers[cluster_id, 0]        
            weighted_lon = cluster_centers[cluster_id, 1]
            weighted_count = cluster_counts[cluster_id]
            weighted_height = cluster_centers[cluster_id, 2]

            sample_data = sample_data.append({
                'lat': weighted_lat,
                'lon': weighted_lon,
                'mass': total_mass,
                'count': weighted_count,
                'height': weighted_height,
                'time': distinct_time,
            }, ignore_index=True)


    return sample_data

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

#____________________________________________________________
# updated function trajectory_volc_dist by 12/05/2023 !Bavand
def trajectory_volc_dist(obs_data_point, tdumpfiles, vlocation, trange_start_time, trange_end_time):
    """
    Functions reads the xarray datasets of observations and tdump files and finds nearest trajectories,
    the output is a dataframe.
    Input:
        xarray dataset representing volcat data (obs_data_point).
        vlocation has longitude and latitude of the volcano.
        tdump files are the results of back trajectories calculated for a layer.
        trange_start_time and trange_end_time are start time and end time of each trajectory for calculation
        of nearest point to vent.
    Outputs:
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

    traj_dataframe = hytraj.combine_dataset(tdumpfiles, renumber=True)
    traj_dataframe['adjusted_longitude'] = np.where(traj_dataframe['longitude'] < 0, traj_dataframe['longitude'] + 360, traj_dataframe['longitude'])
    # add new column of initial altitude to mark the starting height of each trajectories
    traj_dataframe['init_alt']    = traj_dataframe['altitude'].iloc[0]
    # sort the dataframe with the traj_num -from beginning to the end- and time from the trajectory start point to the last point 
    traj_dataframe = traj_dataframe.sort_values(by=['traj_num','time'], ascending=[True, False]).copy()
    traj_dataframe.reset_index(drop=True, inplace=True)
    # select those trajectories that fall within a certain time range after the eruption event
    #trange_start_time = '2022-01-15 04:00:00'
    #trange_end_time   = '2022-01-15 16:00:00'
    traj_dataframe = traj_dataframe[(traj_dataframe['time'] >= trange_start_time) & (traj_dataframe['time'] <= trange_end_time)].copy()
    deg2km = 111.111
    # find the distances between the trajectories and volcano, nearest point and its index
    traj_dataframe['dist_len'] = np.sqrt(((traj_dataframe['latitude']-vlocation[0])*(np.cos(np.radians(0.5*(traj_dataframe['latitude']+vlocation[0]))))*deg2km)**2 + ((traj_dataframe['adjusted_longitude']-vlocation[1])*deg2km)**2) 

    min_dist_rows = traj_dataframe.groupby('traj_num')['dist_len'].idxmin()
    result_df = traj_dataframe.loc[min_dist_rows]
    result_df.reset_index(drop=True, inplace=True)

    result_df['obs_time'] = pd.to_datetime(obs_data_point['time'], format='%Y-%m-%d %H:%M:%S', errors='coerce')
    result_df['obs_lat'] = obs_data_point['lat']
    result_df['obs_lon'] = obs_data_point['lon']
    result_df['obs_height'] = obs_data_point['height']*1000
    result_df['dist_weight'] = obs_data_point['mass']
    result_df['obs_count'] = obs_data_point['count']
    result_df['col_conc'] = result_df['dist_weight']/result_df['obs_count']
    result_df.rename(columns={'latitude':'dist_lat', 'adjusted_longitude':'dist_lon', 'altitude':'dist_hgt', 'time':'dist_time'}, inplace = True)
    
    return result_df

def traj_layer_dataframe2(datadir, fname, suffix, num_lyr, vlocation, trange_start_time, trange_end_time):
    """
    This function outputs the characteristics of the trajectories' nearest point to the volcano.
    """
    var1 = {}

    for i in range(1, num_lyr):

        print(i)
        tdump_files1 = glob.glob(datadir +  f'{fname}_{"%02d" %i}km/' + f'tdump.ashtest_btraj{"%02d" %i}{suffix}.*')
        obs_data_point1 = pd.read_csv(datadir +  f'{fname}_{"%02d" %i}km/' + f'btraj{"%02d" %i}{suffix}.csv')
#        print(tdump_files1)
#        print('***********************************')
#        print(obs_data_point1)
#        print('***********************************')
        print('****************end first step********************')
        var1[i] = trajectory_volc_dist(obs_data_point1, tdump_files1, vlocation, trange_start_time, trange_end_time)

#    return (var1)

#        obs_data = {'dist_len': dist_len, 'dist_hgt': dist_hgt, 'dist_lat': dist_lat, 'dist_lon': dist_lon,
#                    'dist_time': dist_time, 'obs_time': obs_time, 'obs_lat':obs_lat, 'obs_lon': obs_lon,
#                    'init_alt': init_alt, 'obs_height': obs_height, 'dist_weight': dist_weight}
#        var1[i] = obs_data
    print('**************************************************')
#    print('****************end second step********************')
    df = {}

    for ii in range(len(obs_data_point1)):
        print(ii)
        df_variable = pd.DataFrame()
        for j in range(1, num_lyr):

            df_variable = pd.concat([df_variable, pd.DataFrame(var1[j])[ii:ii+1]])
            df[ii] = df_variable[['obs_time', 'obs_lat', 'obs_lon', 'obs_height', 'init_alt', 'dist_time', 'dist_hgt', 
                                'dist_lat', 'dist_lon', 'dist_len', 'dist_weight','obs_count', 'col_conc']]

    return (df)
#____________________________________________________________

def combine_traj_old(fnames, csvfile=None):
    """
    fnames  : list of str. trajectory file names. full path.
    csvfile : csv file output by sample_and_write which contains weighting information.
    combined trajectories in different files into one dataframe.
    """
    trajlist = []
    if csvfile:
        weightcsv = pd.read_csv(csvfile)
    for iii, fnn in enumerate(fnames):
        try:
            df1 = hytraj.open_dataset(fnn)
        except:
            print('Failed {}'.format(fnn))
            continue
        # get trajectory number from the file name
        temp = fnn.split(".")
        trajnum = int(temp[-1])
        # add new column to dataframe with trajectory number
        df1["traj_num"] = trajnum
        #print('TRAJNUM', trajnum)
        # add weight information from csvfile to the dataframe
        if csvfile:
            temp = weightcsv.loc[trajnum]
        #    weight = temp.massI
            weight = temp.mass
        else:
            weight = 1
        df1["weight"] = weight
 
        trajlist.append(df1.copy())
    # concatenate the trajectories into one dataframe.
    trajdf = pd.concat(trajlist)
    return trajdf


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
        except:
            print('Failed {}'.format(fnn))
            continue
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
        #    obsarea = temp.area
        else:
            weight = 1
            obsalt = None
            obsarea = None
        df1["weight"] = weight
        df1["obsalt"] = obsalt
        #df1["area"] = obsarea

        trajlist.append(df1.copy())
        print(fnn)
        print('bavand')
    # concatenate the trajectories into one dataframe.
    trajdf = pd.concat(trajlist)
    return trajdf



def combine_traj2(fnames, csvfile=None):
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
#03/28/2024
#        weightcsv.rename(columns={weightcsv.columns[6]: 'bavand'}, inplace=True)

#03/28/2024
    for iii, fnn in enumerate(fnames):
        try:
            df1 = hytraj.open_dataset(fnn)
        except:
            print('Failed {}'.format(fnn))
            continue
        # get trajectory number from the file name
        temp = fnn.split(".")
        trajnum = int(temp[-1])
        # add new column to dataframe with run number
        # the run may have multiple trajectories in it.
        #df1["traj_num"] = trajnum

#03/28/2024 BS7
        df1["run_num"] = trajnum
        #df1["traj_num"] = trajnum
#03/28/2024 BS7


        #print('TRAJNUM', trajnum)
        # add weight information from csvfile to the dataframe
        if csvfile:
            temp = weightcsv.loc[trajnum]
        #    weight = temp.massI
            weight = temp.mass
            obsalt = temp.height
            #obsarea = temp.area

        else:
            weight = 1
            obsalt = None
            obsarea = None
        df1["weight"] = weight
        df1["obsalt"] = obsalt
        #df1["area"] = obsarea

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
    for ii in range(1,num_layer):
        tdump_files = natsorted(glob.glob(data_dir + f'{fname}_{"%02d" %ii}km/' + 'tdump.ashtest_btraj*'))
        trajectory_data[f'tdump{ii}'] = tdump_files
        tnames_path.append(trajectory_data[f'tdump{ii}'])
        obs_path.append(data_dir + f'{fname}_{"%02d" %ii}km/' + f'btraj{"%02d" %ii}kmDES.csv')
    return (tnames_path, obs_path)


def traj_layer_dataframe(tnames_path, obs_path, vloc, tdump_num):
    """
    This function outputs the characteristics of the trajectories' nearest point to the volcano.
    """
    var1 = {}

    for i in range(0, len(obs_path)):

        print(obs_path[i])
        dist_func = traj_layer_cal(tnames_path[i], obs_path[i], vloc, tdump_num)

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

    for ii in range(tdump_num):
        df_variable = pd.DataFrame()
        for j in range(0, len(obs_path)):

            df_variable = pd.concat([df_variable, pd.DataFrame(var1[j])[ii:ii+1]])
            df[ii] = df_variable

    return (df)

# the function was modified by Bavand 12/05/2023
def plume_thick_cal(df):
    """
    This function selects the layers at which the shortest distance between trajectories and volcano points occurs.
    It sets a cut-off value and identifies the base and top of the cloud, allowing the code to calculate the cloud thickness.
    Outputs: a dataframe containing the cloud characteristics (columns) for entire observations (rows)
    """
    df_dist_min = pd.DataFrame()
    df_dist_min_criteria_pass = pd.DataFrame()
    df_closest_row = pd.DataFrame()

    critera_pass_traj_num = 0
    
#    cutoff = 500
    for i in range(len(df)):
#        cutoff = 1.5*df[i]['dist_len'].min()
        cutoff = 800
        # finds the layers where the trajectories come within a certain distance to the volcano
        selected_rows = df[i].loc[df[i]['dist_len'] <= cutoff]


        # If there was no single layer where the trajectory comes within this threshold distance of the vent, the code will choose
        # the initial altitude of the nearest trajectory layer and set it as the top of the cloud.
#        if len(selected_rows) == 1:
        if len(selected_rows) == 0 or len(selected_rows) == 1:
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
            # df_dist_min is the dataframe of plume characteristics over entire observation points
            df_dist_min = pd.concat([df_dist_min, df_closest_row])
            # df_dist_min is the dataframe of plume characteristics over observation points that come within the cut-off value.
            df_dist_min_criteria_pass = pd.concat([df_dist_min_criteria_pass, df_closest_row])
            critera_pass_traj_num = critera_pass_traj_num + 1

    return(df_dist_min, df_dist_min_criteria_pass, critera_pass_traj_num)

# This function was modified by a dynamic method of selecting the cutoff, 02/21/2023
def conc_emitimes_data(df):
    """"
    This function prepares the characteristics of plume emission for EMITIMES needed for forward dispersion run.
    """
    df_dist_fwd_data = pd.DataFrame()
    df_dist_fwd = {}
    deg2km = 111.111

    critera_pass_traj_num = 0
    for i in range(1,len(df)):

        cutoff = 1.5*df[i]['dist_len'].min()
        # finds the layers where the trajectories come within a certain distance to the volcano
        selected_rows = df[i].loc[df[i]['dist_len'] < cutoff]

        # if there was no single layer where the trajectory comes within this threshold distance of the vent, the code will
        # choose the initial altitude of the nearest trajectory layer and set it as the top of the cloud
        if len(selected_rows) == 1:

            df_dist_fwd[i] = pd.DataFrame(df[i].iloc[df[i]['dist_len'].argmin()]).transpose()
            #print('length of the variable here is', len(df_dist_fwd[i]))
            # (I) check if it is needed to convert the DU to g/m^2
            # df_dist_fwd[i].loc[:, 'tload'] = (df_dist_fwd[i]['dist_weight']) * (2.6867E20) * (64.066) * (1/(6.022E23))
            # (II) check if it is needed to convert the g/m^2 to g (7.0 km ~ 7000 m)
            # df_dist_fwd[i]['tload'] = df_dist_fwd[i]['tload'] * 49000000
            # (III) check what convertion the unit of mass requires. Units will remain unchanged
            df_dist_fwd[i] = df_dist_fwd[i].copy()
        
            # to add a new column showing the number of source emission layers at the point
            df_dist_fwd[i] = df_dist_fwd[i].assign(count = len(df_dist_fwd[i]))
            obs_lat_radians = np.radians(df_dist_fwd[i]['obs_lat'].iloc[0])
            area = (np.abs(np.cos(obs_lat_radians)) * 0.045 * deg2km * 1000) * (0.045 * deg2km * 1000)
            print('area')
            print(area)

            df_dist_fwd[i].loc[:, 'tload'] = df_dist_fwd[i]['col_conc'] * df_dist_fwd[i]['obs_count'] * area
#            df_dist_fwd[i].loc[:, 'tload'] = df_dist_fwd[i]['col_conc'] * df_dist_fwd[i]['obs_count'] * 25000000
#            df['altitude'].iloc[0]
            df_dist_fwd[i].loc[:, 'tload'] = df_dist_fwd[i]['tload'] / len(df_dist_fwd[i])
            # will add a new column for the emission rate (1/hr): df_dist_fwd[i]['tload'] * (1/[emission time (hr)])
            # rate unit (1/hr)
            df_dist_fwd[i]['rate'] = df_dist_fwd[i]['tload'] / (1/12)

            # will set the date
            df_dist_fwd[i]['YYYY'] = df_dist_fwd[i]['obs_time'].dt.year
            df_dist_fwd[i]['MM']   = df_dist_fwd[i]['obs_time'].dt.month
            df_dist_fwd[i]['DD']   = df_dist_fwd[i]['obs_time'].dt.day
            df_dist_fwd[i]['HH']   = df_dist_fwd[i]['obs_time'].dt.hour
            #df_dist_fwd[i]['MIN']  = df_dist_fwd[i]['obs_time'].dt.minute
            df_dist_fwd[i]['MIN']  = (df_dist_fwd[i]['obs_time'].dt.minute // 5) *5

            # will add up the emission to data file 
            df_dist_fwd_data = pd.concat([df_dist_fwd_data, df_dist_fwd[i]])

        else:

            df_dist_fwd[i] = selected_rows
            # (I) check if it is needed to convert the DU to g/m^2
            # df_dist_fwd[i].loc[:, 'tload'] = (df_dist_fwd[i]['dist_weight']) * (2.6867E20) * (64.066) * (1/(6.022E23))
            # (II) check if it is needed to convert the g/m^2 to g (7.0 km ~ 7000 m)
            # df_dist_fwd[i]['tload'] = df_dist_fwd[i]['tload'] * 49000000
            # (III) check what convertion the unit of mass requires. Units will remain unchanged
            df_dist_fwd[i] = df_dist_fwd[i].copy()
            # to add a new column showing the number of source emission layers at the point
            df_dist_fwd[i] = df_dist_fwd[i].assign(count = len(df_dist_fwd[i]))
            obs_lat_radians = np.radians(df_dist_fwd[i]['obs_lat'].iloc[0])
            area = (np.abs(np.cos(obs_lat_radians)) * 0.045 * deg2km * 1000) * (0.045 * deg2km * 1000)
            print('area')
            print(area)

#            df_dist_fwd[i].loc[:, 'tload'] = (df_dist_fwd[i]['col_conc'] * df_dist_fwd[i]['obs_count'] * 25000000)
            df_dist_fwd[i].loc[:, 'tload'] = df_dist_fwd[i]['col_conc'] * df_dist_fwd[i]['obs_count'] * area
            # distribute the total mass over the layers
            df_dist_fwd[i].loc[:, 'tload'] = df_dist_fwd[i]['tload'] / len(df_dist_fwd[i])

            # will add a new column for the emission rate (1/hr): df_dist_fwd[i]['tload'] * (1/[emission time (hr)])
            # rate unit (1/hr)
            df_dist_fwd[i]['rate'] = df_dist_fwd[i]['tload'] / (1/12)

            # will set the date
            df_dist_fwd[i]['YYYY'] = df_dist_fwd[i]['obs_time'].dt.year
            df_dist_fwd[i]['MM']   = df_dist_fwd[i]['obs_time'].dt.month
            df_dist_fwd[i]['DD']   = df_dist_fwd[i]['obs_time'].dt.day
            df_dist_fwd[i]['HH']   = df_dist_fwd[i]['obs_time'].dt.hour
            #df_dist_fwd[i]['MIN']  = df_dist_fwd[i]['obs_time'].dt.minute
            df_dist_fwd[i]['MIN']  = (df_dist_fwd[i]['obs_time'].dt.minute // 5) *5
        
            # will add up the emission to data file 
            df_dist_fwd_data = pd.concat([df_dist_fwd_data, df_dist_fwd[i]])

    return(df_dist_fwd_data, df_dist_fwd)

# This function was modified by a dynamic method of selecting the cutoff, 02/21/2023
def conc_emitimes_data_raw(df):
    """"
    This function prepares the characteristics of plume emission for EMITIMES needed for forward dispersion run.
    """
    df_dist_fwd_data = pd.DataFrame()
    df_dist_fwd = {}
    deg2km = 111.111

    critera_pass_traj_num = 0
    for i in range(1,len(df)):

        cutoff = 1.5*df[i]['dist_len'].min()
        # finds the layers where the trajectories come within a certain distance to the volcano
        selected_rows = df[i].loc[df[i]['dist_len'] < cutoff]

        # if there was no single layer where the trajectory comes within this threshold distance of the vent, the code will
        # choose the initial altitude of the nearest trajectory layer and set it as the top of the cloud
        #if len(selected_rows) == 1:

        df_dist_fwd[i] = pd.DataFrame(df[i].iloc[df[i]['dist_len'].argmin()]).transpose()
            #print('length of the variable here is', len(df_dist_fwd[i]))
            # (I) check if it is needed to convert the DU to g/m^2
            # df_dist_fwd[i].loc[:, 'tload'] = (df_dist_fwd[i]['dist_weight']) * (2.6867E20) * (64.066) * (1/(6.022E23))
            # (II) check if it is needed to convert the g/m^2 to g (7.0 km ~ 7000 m)
            # df_dist_fwd[i]['tload'] = df_dist_fwd[i]['tload'] * 49000000
            # (III) check what convertion the unit of mass requires. Units will remain unchanged
        df_dist_fwd[i] = df_dist_fwd[i].copy()
            # to add a new column showing the number of source emission layers at the point
        df_dist_fwd[i] = df_dist_fwd[i].assign(count = len(df_dist_fwd[i]))

        obs_lat_radians = np.radians(df_dist_fwd[i]['obs_lat'].iloc[0])
        area = (np.abs(np.cos(obs_lat_radians)) * 0.045 * deg2km * 1000) * (0.045 * deg2km * 1000)
        print('area')
        print(area)
     
        df_dist_fwd[i].loc[:, 'tload'] = df_dist_fwd[i]['col_conc'] * df_dist_fwd[i]['obs_count'] * area
#            df['altitude'].iloc[0]
        df_dist_fwd[i].loc[:, 'tload'] = df_dist_fwd[i]['tload'] / 1
            # will add a new column for the emission rate (1/hr): df_dist_fwd[i]['tload'] * (1/[emission time (hr)])
            # rate unit (1/hr)
        df_dist_fwd[i]['rate'] = df_dist_fwd[i]['tload'] / (1/12)

            # will set the date
        df_dist_fwd[i]['YYYY'] = df_dist_fwd[i]['obs_time'].dt.year
        df_dist_fwd[i]['MM']   = df_dist_fwd[i]['obs_time'].dt.month
        df_dist_fwd[i]['DD']   = df_dist_fwd[i]['obs_time'].dt.day
        df_dist_fwd[i]['HH']   = df_dist_fwd[i]['obs_time'].dt.hour
        #df_dist_fwd[i]['MIN']  = df_dist_fwd[i]['obs_time'].dt.minute
        df_dist_fwd[i]['MIN']  = (df_dist_fwd[i]['obs_time'].dt.minute // 5) *5

            # will add up the emission to data file
        df_dist_fwd_data = pd.concat([df_dist_fwd_data, df_dist_fwd[i]])

#else:

#            df_dist_fwd[i] = selected_rows

            # (I) check if it is needed to convert the DU to g/m^2
            # df_dist_fwd[i].loc[:, 'tload'] = (df_dist_fwd[i]['dist_weight']) * (2.6867E20) * (64.066) * (1/(6.022E23))
            # (II) check if it is needed to convert the g/m^2 to g (7.0 km ~ 7000 m)
            # df_dist_fwd[i]['tload'] = df_dist_fwd[i]['tload'] * 49000000
            # (III) check what convertion the unit of mass requires. Units will remain unchanged
#            df_dist_fwd[i] = df_dist_fwd[i].copy()

            # to add a new column showing the number of source emission layers at the point
#            df_dist_fwd[i] = df_dist_fwd[i].assign(count = len(df_dist_fwd[i]))


#            df_dist_fwd[i].loc[:, 'tload'] = (df_dist_fwd[i]['col_conc'] * df_dist_fwd[i]['obs_count'] * 25000000 * 1000)
            # distribute the total mass over the layers
#            df_dist_fwd[i].loc[:, 'tload'] = df_dist_fwd[i]['tload'] / 1


            # will add a new column for the emission rate (1/hr): df_dist_fwd[i]['tload'] * (1/[emission time (hr)])
            # rate unit (1/hr)
#           df_dist_fwd[i]['rate'] = df_dist_fwd[i]['tload'] / (1/12)

            # will set the date
#            df_dist_fwd[i]['YYYY'] = df_dist_fwd[i]['obs_time'].dt.year
#            df_dist_fwd[i]['MM']   = df_dist_fwd[i]['obs_time'].dt.month
#            df_dist_fwd[i]['DD']   = df_dist_fwd[i]['obs_time'].dt.day
#            df_dist_fwd[i]['HH']   = df_dist_fwd[i]['obs_time'].dt.hour
#            df_dist_fwd[i]['MIN']  = df_dist_fwd[i]['obs_time'].dt.minute

            # will add up the emission to data file
#            df_dist_fwd_data = pd.concat([df_dist_fwd_data, df_dist_fwd[i]])

    return(df_dist_fwd_data, df_dist_fwd)



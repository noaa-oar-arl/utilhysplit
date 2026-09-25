from ashapp.level_setter import set_qva_levels
import numpy as np
import pandas as pd

"""
Used as input into RunTrajectory class.
the trajectory generators output a dictionary with height, latitude, longitude.
the time generators output a time and a trajectory generator.

"""

# 2023 Dec 04 (amc) added generate_traj_from_config 
# 2023 Dec 04 (amc) added generate_height_traj_from_series
# 2023 Dec 04 (amc) added timegenerate_height_traj_from_obsdf 
# 2023 Dec 04 (amc) changed generate_traj_from_obsdf to timegenerate_traj_from_obsdf


#def generate_traj_from_config(inp):
#    outp={}
#    height = inp["height"]
#    if not isinstance(height,(list,np.ndarray)):
#       height = [height]
#    for hgt in height:
#        outp['height'] = hgt
#        outp['latitude'] = lat
#        outp['longitude'] = lon
#        yield outp


def generate_qva_traj_from_config(inp):
    bottom = int(inp['bottom']*3.28084/100)
    top = int(inp['top']*3.28084/100)
    levels = set_qva_levels(bottom,top,dz=50)
    inp['height'] = levels[0]
    for traj in generate_traj_from_config(inp):
        yield traj


def generate_traj_from_config(inp):
    for key in ['height','top','bottom']:
        if key in inp.keys():
           hkey = key
           break
    outp={}
    lat = inp["latitude"]
    lon = inp["longitude"]
    height = inp[hkey]
    if not isinstance(height,(list,np.ndarray)):
       height = [height]
    for hgt in height:
        outp['height'] = hgt
        outp['latitude'] = lat
        outp['longitude'] = lon
        yield outp

def generate_height_traj_from_series(series: pd.core.series.Series,
                              minz : float,
                              maxz : float,
                              dz   : float):
    # for one location generates trajectories from different heights.
    for zzz in np.arange(minz,maxz,dz):
        outp = {}
        # this is specifically for file from hunga tonga so2 data.
        # heightI uses the interpolated data as well.
        outp["height"] = zzz * 1000
        outp["latitude"] = series.lat
        outp["longitude"] = series.lon
        yield outp 


def generate_traj_from_df(df):
    for rrr in df.iterrows():
        outp = {}
        row = rrr[1]
        # this is specifically for file from hunga tonga so2 data.
        # heightI uses the interpolated data as well.
        if 'heightI' in df.columns:
            outp["height"] = row.heightI * 1000
        else:
            outp["height"] = row.height * 1000
        outp["latitude"] = row.lat
        outp["longitude"] = row.lon
        yield outp 

def timegenerate_traj_from_obsdf(csvname):
    """
    generate back trajectories from a csv file. This will create tdump files grouped by time.
    RETURNS
    tuple (time, generate_traj_from_df function)
    """
    obsdf = pd.read_csv(csvname, parse_dates=["time"])
    #return obsdf
    # one trajectory run per time period.
    timelist = obsdf['time'].unique()
    for time in timelist:
        newdf = obsdf[obsdf['time']==time]
        time = pd.to_datetime(time)
        yield time, generate_traj_from_df(newdf)

def timegenerate_height_traj_from_obsdf(csvname,minht=1,maxht=None,dh=1):
    """
    generate back trajectories from a csv file. 
    this will create tdump files for each observation point with multiple heights.

    RETURNS
    tuple (time, generate_height_traj_from_series function)
    """
    if isinstance(csvname,str):
        obsdf = pd.read_csv(csvname, parse_dates=["time"])
    if isinstance(csvname,pd.DataFrame):
        obsdf = csvname
    # change this to change what the maximum height that is used.
    if maxht is None:
        maxht = np.max(obsdf.height) + 2
    # trajectory run for each location.
    # all heights run in the same run for that location.
    timelist = obsdf['time'].unique()
    for time in timelist:
        newdf = obsdf[obsdf['time']==time]
        for ii,row in newdf.iterrows():
            yield time, generate_height_traj_from_series(row,minht,maxht,dh)



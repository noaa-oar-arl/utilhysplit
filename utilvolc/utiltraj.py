import pandas as pd
from monetio.models import hytraj



def combine_traj(fnames, csvfile=None):
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




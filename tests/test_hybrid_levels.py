
import pandas as pd

def era5_levels():
    fname = '/hysplit-users/alicec/era5/model_levels/SIG_ERA52ARL.MESSAGE'
    a = pd.read_csv(fname)
    return a


def gfs_levels():
    # file from Jianping Huang. For GFS 128 layers.
    fname = '/hysplit-users/alicec/projects/testing/global_hyblev.l128.txt'

    iii = 0
    sp = 1013
    plist = []
    with open(fname) as fid:
         for line in fid.readlines():
             if iii==0: 
                iii+=1
                continue
             temp = line.split()
             pressure = float(temp[0])/100.0 + float(temp[1])*sp
             plist.append(pressure)
             iii+=1
    plt.plot(plist,'k.')
    print(plist[0:10])
    print(plist[-20:-1])

# plume_stat.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
#from math import *
#import sys
#import string
#import datetime
#from itertools import permutations
#import xarray as xr
#from math import pi, cos
"""
routines to help calculates critical success index by a point by point comparison.
data must be gridded the same.
Functions
---------
calc_csi (critical success index, aka figure of merit in space)
calc_fss (fraction skill score)
calc_bs (brier score)
calc_bss (brier skill score)
"""


def calc_fss(ra1, ra2, threshold1=0, threshold2=0, szra=[1, 3, 5, 7], makeplots=False, verbose=0):
    from scipy.signal import convolve2d
    """Calculates the fraction skill score (fss)
    See Robers and Lean (2008) Monthly Weather Review
    and Schwartz et al (2010) Weather and Forecasting
    for more information.

    ra1: observations/satellite
    ra2: the forecast/model

    Can plot fractions if desired (double check calculations)
    threshold1: value for data threshold
    threshold2: value for model threshold
    szra: a list of the number of pixels (neightborhood length) to use in fractions calculation default is to use 1, 3, 5, 7 pixels size squares
    makeplots: boolean
    verbose:

    Return
       df : pandas dataframe
    """
    # Creating FSS dictionary
    fss_dict = {}
    bigN = ra1.size

    # create binary fields
    mask1 = mask_threshold(ra1, threshold=threshold1)
    mask2 = mask_threshold(ra2, threshold=threshold2)

    # loop for the convolutions
    for sz in szra:
        if sz == 1:
            filter_array = np.zeros((3, 3))
            filter_array[1, 1] = 1.
            conv_array = (1, 1)
        else:
            filter_array = np.ones((sz, sz))
            filter_array = filter_array * (1/np.sum(filter_array))
            conv_array = np.shape(filter_array)

        if verbose:
            print('Convolution array size: ', np.shape(filter_array))
        if verbose:
            print('Convolution array: ', filter_array)

        start = time.time()
        frac_arr1 = convolve2d(mask1, filter_array, mode='same')
        frac_arr2 = convolve2d(mask2, filter_array, mode='same')
        end = time.time()
        if verbose:
            print('Calculation time: ', end - start)

        # Calculate the Fractions Brier Score (FBS)
        fbs = np.power(frac_arr1 - frac_arr2, 2).sum() / float(bigN)
        print('FBS ', fbs)
        # Calculate the worst possible FBS (assuming no overlap of nonzero fractions)
        fbs_ref = (np.power(frac_arr1, 2).sum() + np.power(frac_arr2, 2).sum()) / float(bigN)
        print('FBS reference', fbs_ref)
        # Calculate the Fractional Skill Score (FSS)
        fss = 1 - (fbs / fbs_ref)
        print('FSS ', fss)
        fss_tmp = dict({'Nlen': conv_array[0], 'FBS': fbs, 'FBS_ref': fbs_ref, 'FSS': fss})
        if (makeplots == True):
            fig = plt.figure(1)
            ax1 = fig.add_subplot(2, 1, 1)
            ax2 = fig.add_subplot(2, 1, 2)
            ax1.imshow(frac_arr1)
            ax2.imshow(frac_arr2)
            plt.show()
        print(' ')
        fss_dict[sz] = fss_tmp
    df = pd.DataFrame.from_dict(fss_dict, orient='index')
    return df


def calc_csi(ra1, ra2, area=None, nodata_value='', threshold=0, verbose=0):
    """ CSI equation: hits / (hits + misses + false alarms)
    ra1 and ra2 need to be on the same grid.
    See monet.remap_nearest to remap arrays.

    ra1: observations/satellite xarray
    ra2: the forecast/model xarray
    area: optional array of grid areas (default = None)
    nodata_value: flag for expanded array grid cells created if ash near boundary
    threshold: value for data threshold
    verbose: boolean
    """
    area = np.array(area)
    # Creating a csihash disctionary
    csihash = {}
    # Converting ra1 and ra2 to arrays of 0's and 1's (1 with values, 0 no values)
    mask1 = mask_threshold(ra1, threshold=0.)
    mask2 = mask_threshold(ra2, threshold=threshold)

    # Calculating hits (matchra), misses (onlyra1), false alarms (onlyra2) for CSI calculation
    matchra = mask2 * mask1  # returns array with ones where the two arrays both have valid values.
    onlyra1 = mask1 - matchra  # returns array with ones where ra1 has valid values and ra2 does not.
    onlyra2 = mask2 - matchra  # returns array with ones where ra2 has valid values and ra1 does not.
    # Assigning a, b, and c arrays to csihash dictionary
    csihash['matched'] = matchra
    csihash['ra1'] = onlyra1
    csihash['ra2'] = onlyra2

    allra = matchra + onlyra1 + onlyra2

    vpi = np.where(allra == 1)  # Why is this necessary?
    # Find pattern correlation (from Zidikheri and Potts) doesn't make sense to me.
    totalpts = matchra.sum() + onlyra1.sum() + onlyra2.sum()
    totalpts = mask2.shape[0] * mask2.shape[1]
    ra1ave = (onlyra1.sum() + matchra.sum()) / float(totalpts)
    ra2ave = (onlyra2.sum() + matchra.sum()) / float(totalpts)

    ra1corr = (mask1 - ra1ave)
    ra2corr = (mask2 - ra2ave)
    norm = ((ra1corr * ra1corr).sum())**0.5 * ((ra2corr*ra2corr).sum())**0.5
    pcorr = (ra1corr * ra2corr).sum() / norm
    print('PCORR', pcorr.values)
    print('ra1ave (obs)',  ra1ave.values)
    print('ra2ave (calc)', ra2ave.values, end=' ')
    print('NORM', norm.values, 'ra1corr*ra2corr',  (ra1corr * ra2corr).sum().values)
    print(mask1.sum().values, mask2.sum().values, np.max(ra1corr), np.min(ra1corr))

    ra1ave = 0
    ra2ave = 0
    ra1corr = (mask1 - ra1ave)
    ra2corr = (mask2 - ra2ave)
    norm = ((ra1corr * ra1corr).sum())**0.5 * ((ra2corr*ra2corr).sum())**0.5
    pcorr = (ra1corr * ra2corr).sum() / norm
    print('PCORR (uncentered)', pcorr.values)

    # Find where data arrays have no information.
    if nodata_value != '':
        vpND = np.where(logical_or(ra1 == nodata_value, ra2 == nodata_value))
        maskND = ra2[:] * 0
        maskND[vpND] = 1
        vp_ra1ND = np.where(logical_and(maskND == 1, onlyra1 == 1))
        if vp_ra1ND[0] != []:
            onlyra1[vp_ra1ND] = 0
        vp_ra2ND = np.where(logical_and(maskND == 1, onlyra2 == 1))
        if vp_ra2ND[0] != []:
            onlyra2[vp_ra2ND] = 0

    if area.shape == matchra.shape:
        csihash['CSI'] = (matchra*area).sum() / ((matchra*area).sum() +
                                                 (onlyra1*area).sum() + (onlyra2*area).sum())
        csihash['POD'] = (matchra*area).sum() / ((matchra*area).sum() + (onlyra1*area).sum())
        csihash['FAR'] = (onlyra2*area).sum() / ((matchra*area).sum() + (onlyra2*area).sum())
        if verbose == 1:
            print('used area')
            print((matchra*area).sum().values, (onlyra1*area).sum().values, (onlyra2*area).sum().values)
            print('CSI POD FAR', csihash['CSI'].values, csihash['POD'].values, csihash['FAR'].values)
    else:
        csihash['CSI'] = matchra.sum() / (matchra.sum() + onlyra1.sum() + onlyra2.sum())
        # hit rate or probability of detection (p 310 Wilks)
        csihash['POD'] = matchra.sum() / (matchra.sum() + onlyra1.sum())
        # false alarm ratio (p 310 Wilks)
        csihash['FAR'] = onlyra2.sum() / (matchra.sum() + onlyra2.sum())
        csihash['PCORR'] = pcorr
        if verbose == 1:
            # print('HERE' , matchra.size.values, matchra.shape.values, onlyra2.shape.values)
            print(matchra.sum().values, onlyra1.sum().values, onlyra2.sum().values)
    print('CSI POD FAR', csihash['CSI'].values, csihash['POD'].values, csihash['FAR'].values)

    return csihash


def calc_bs(xra1, xra2):
    """
    Calculating the Brier Score 
    BS = 1/N (sum of (probability - reference)^2)
    Inputs:
         xra1: probability forecast array (values from 0 to 1) (xarray.DataArray)
         xra2: binary reference or observation array (xarray.DataArray)
    Outputs:
         BS: Brier score (float)
    """
    # tmp = probability - actual
    tmp = xra1 - xra2
    tmp2 = tmp * tmp
    # N = size of tmp2 array
    N = np.size(tmp2)
    # Sum of tmp2 divided by N
    BS = np.sum(tmp2) / N
    return BS


def calc_bss(BS, BSref):
    """
    Calculating the Brier Skill Score
    BSS = 1 - BS/BSref
    Inputs:
         BS: Brier Score of the probabilistic forecast compared to observations (float)
         BSref: Brier Score of probabilistic forecast compared to reference forecast (float)
    Outputs:
         BSS: Brier Skill Score of probabilistic forecast (float)
    """
    BSS = 1 - (BS / BSref)
    return BSS

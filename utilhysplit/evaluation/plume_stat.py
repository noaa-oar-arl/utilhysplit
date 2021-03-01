# plume_stat.py

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import time

#from math import *
#import sys
#import string
#import datetime
#from itertools import permutations
#from math import pi, cos
"""
routines to help calculates critical success index by a point by point comparison.
data must be gridded the same.
---------
Class:
    CalcScores
         calc_csi (critical success index, aka figure of merit in space)
         calc_fss (fraction skill score)
         calc_pcorr (pattern correlation coefficient)
Functions:
    calc_bs (brier score)
    calc_bss (brier skill score)
"""


class CalcScores:

    def __init__(self, xra1, xra2, threshold=0., szra=[1, 3, 5, 7], area=None, verbose=True):
        """ Class of tools for calculating various Scores and Skill Scores, 
        relying on binary arrays and the 2x2 contingency table
        xra1 and xra2 need to be on the same grid.
        See monetio.remap_nearest to remap arrays.
        ----------------------------
        Inputs:
        xra1: observation/satellite array (xarray DataArray)
        xra2: hysplit forecast array (xarray DataArray)
        threshold: data threshold for creating binary fields, default = 0. (float)
        szra: sizes for fractions skill score calculation, default = [1, 3, 5, 7] (list)
        area: optional array of grid areas, must be the same size as xra1 and xra2
        verbose: boolean
        ----------------------------
        Functions:
        calc_csi: calculates the critical success index (gilbert score), probability of detection, false alarm rate, etc.
        calc_fss: calculates the fractions skill score
        calc_pcorr: calculates pattern correlation coefficient (Pearson correlation coefficient) 
        """

        self.xra1 = xra1
        self.xra2 = xra2
        self.threshold = threshold
        self.szra = szra
        self.area = area
        self.verbose = verbose

        self.binxra1 = xr.where(self.xra1 >= self.threshold, 1., 0.)
        self.binxra2 = xr.where(self.xra2 >= self.threshold, 1., 0.)
        self.match = self.binxra1 * self.binxra2
        self.arr1 = self.binxra1 - self.match
        self.arr2 = self.binxra2 - self.match
        self.totalpts = self.binxra1.shape[0] * self.binxra1.shape[1]

    def calc_csi(self):
        """ CSI equation: hits / (hits + misses + false alarms) - aka Gilbert Score

        Inputs:
        match: hits
        arr1: misses
        arr2: false alarms
        area: optional array of grid areas, must be same size as xra1 and xra2 (default = None)
        nodata: flag for expanded array grid cells created if ash near boundary (string)
        verbose: boolean

        Output:
        csihash: Dictionary of values pertaining to Critical Success Index calculation
        contains: hits, misses, false_alarms, CSI, POD, FAR
        """
        area = np.array(self.area)
        # Creating a csihash disctionary
        csihash = {}

        # Assigning a, b, and c arrays to csihash dictionary
        csihash['hits'] = self.match
        csihash['misses'] = self.arr1
        csihash['false_alarms'] = self.arr2

        if area.shape == self.match.shape:
            csihash['CSI'] = (self.match*self.area).sum() / ((self.match*self.area).sum() +
                                                             (self.arr1*self.area).sum() + (self.arr2*self.area).sum())
            csihash['POD'] = (self.match*self.area).sum() / \
                ((self.match*self.area).sum() + (self.arr1*area).sum())
            csihash['FAR'] = (self.arr2*self.area).sum() / \
                ((self.match*self.area).sum() + (self.arr2*area).sum())

            if self.verbose == True:
                print('used area')
                print((self.match*self.area).sum().values, (self.arr1 *
                                                            self.area).sum().values, (self.arr2*self.area).sum().values)
                print('CSI: ', csihash['CSI'].values, 'POD: ',
                      csihash['POD'].values, 'FAR: ', csihash['FAR'].values)
        else:
            csihash['CSI'] = self.match.sum() / (self.match.sum() + self.arr1.sum() + self.arr2.sum())
            # hit rate or probability of detection (p 310 Wilks)
            csihash['POD'] = self.match.sum() / (self.match.sum() + self.arr1.sum())
            # false alarm ratio (p 310 Wilks)
            csihash['FAR'] = self.arr2.sum() / (self.match.sum() + self.arr2.sum())

            if self.verbose == True:
                print('Match: ', self.match.sum().values, 'Misses: ',
                      self.arr1.sum().values, 'FalseAlarms: ', self.arr2.sum().values)
            print('CSI: ', csihash['CSI'].values, 'POD: ',
                  csihash['POD'].values, 'FAR: ', csihash['FAR'].values)

        return csihash

    def calc_fss(self, makeplots=False):
        from scipy.signal import convolve2d
        """Calculates the fraction skill score (fss)
        See Robers and Lean (2008) Monthly Weather Review
        and Schwartz et al (2010) Weather and Forecasting
        for more information.
         
        Can plot fractions if desired (double check calculations)
        szra: a list of the number of pixels (neightborhood length) to use 
        in fractions calculation default is to use 1, 3, 5, 7 pixels size squares
        makeplots: boolean
         
        Return
        df : pandas dataframe """

        # Creating FSS dictionary
        fss_dict = {}
        bigN = self.binxra1.size

        # loop for the convolutions
        for sz in self.szra:
            if sz == 1:
                filter_array = np.zeros((3, 3))
                filter_array[1, 1] = 1.
                conv_array = (1, 1)
            else:
                filter_array = np.ones((sz, sz))
                filter_array = filter_array * (1/np.sum(filter_array))
                conv_array = np.shape(filter_array)

            if self.verbose == True:
                print('Convolution array size: ', np.shape(filter_array))
                print('Convolution array: ', filter_array)

            start = time.time()
            frac_arr1 = convolve2d(self.binxra1, filter_array, mode='same')
            frac_arr2 = convolve2d(self.binxra2, filter_array, mode='same')
            end = time.time()

            # Calculate the Fractions Brier Score (FBS)
            fbs = np.power(frac_arr1 - frac_arr2, 2).sum() / float(bigN)

            # Calculate the worst possible FBS (assuming no overlap of nonzero fractions)
            fbs_ref = (np.power(frac_arr1, 2).sum() + np.power(frac_arr2, 2).sum()) / float(bigN)

            # Calculate the Fractional Skill Score (FSS)
            fss = 1 - (fbs / fbs_ref)
            if self.verbose == True:
                print('Calculation time: ', end - start)
                print('FBS ', fbs)
                print('FBS reference', fbs_ref)
                print('FSS ', fss)
                print(' ')
            fss_tmp = dict({'Nlen': conv_array[0], 'FBS': fbs, 'FBS_ref': fbs_ref, 'FSS': fss})
            if makeplots == True:
                fig = plt.figure(1)
                ax1 = fig.add_subplot(2, 1, 1)
                ax2 = fig.add_subplot(2, 1, 2)
                ax1.imshow(frac_arr1)
                ax2.imshow(frac_arr2)
                plt.show()

            fss_dict[sz] = fss_tmp

        df = pd.DataFrame.from_dict(fss_dict, orient='index')
        return df

    def calc_pcorr(self):
        """ Calculates Pattern Correlation between two arrays
        binxra1 and binxra2 need to be on the same grid. 
        See monetio.remap_nearest or monetio.remap_xesmf to remap arrays.

        Pattern Correlation (uncentered) = (binxra1 * binxra2) / sqrt((binxra1^2)*(binxra2^2))
        Pattern Correlation (centered) = ((binxra1 - arr1)(binxra2 - arr2)) / sqrt(((binxra1-arr1)^2) * ((binxra2 - arr2)^2))

        Outputs:
        pcorr, pcorruncent: pattern correlation (centered), pattern correlation (uncentered) """

        # Space averaged values
        arr1avg = (self.arr1.sum() + self.match.sum()) / float(self.totalpts)
        arr2avg = (self.arr2.sum() + self.match.sum()) / float(self.totalpts)

        # Calculating centered pattern correlation - subtracts space averaged values
        arr1corr = self.binxra1 - arr1avg
        arr2corr = self.binxra2 - arr2avg
        norm = ((arr1corr * arr1corr).sum()) ** 0.5 * ((arr2corr * arr2corr).sum())**0.5
        pcorr = (arr1corr * arr2corr).sum()/norm

        # Calculating uncentered pattern correlation
        norm = ((self.binxra1 * self.binxra1).sum()) ** 0.5 * ((self.binxra2 * self.binxra2).sum())**0.5
        pcorruncent = (self.match).sum()/norm

        if self.verbose == True:
            print('PCORR (centered)', pcorr.values)
            print('PCORR (uncentered)', pcorruncent.values)

        return pcorr, pcorruncent

#FUNCTIONS#


def calc_bs(xra1, xra2):
    """
    Calculating the Brier Score 
    BS = 1/N (sum of (probability - reference)^2)

    Inputs:
    xra1: binary reference or observation array (xarray.DataArray)
    xra2: probability (or binary) forecast array (values from 0 to 1) (xarray.DataArray)

    Outputs:
    BS: Brier score (float)
    """
    # tmp = probability - actual
    tmp = xra2 - xra1
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


def calc_csi(xra1, xra2, area=None, nodata_value='', threshold=0, verbose=0):
    """ CSI equation: hits / (hits + misses + false alarms) - aka Gilbert Score
    ra1 and ra2 need to be on the same grid.
    See monetio.remap_nearest to remap arrays.

    xra1: observations/satellite xarray
    xra2: the forecast/model xarray
    area: optional array of grid areas (default = None)
    nodata_value: flag for expanded array grid cells created if ash near boundary
    threshold: value for data threshold
    verbose: boolean

    Output:
    Dictionary of values pertaining to Critical Success Index calculation
    """
    area = np.array(area)
    # Creating a csihash disctionary
    csihash = {}
# Go into class variables
    # Convert xra1 and xra2 to binary arrays
    binxra1 = xr.where(xra1 >= threshold, 1., 0.)
    binxra2 = xr.where(xra2 >= threshold, 1., 0.)

    # Determining hits or matches (obs yes, forecast yes)
    match = binxra1 * binxra2
    # Determining misses for forecast array (obs yes, forecast no)
    arr1 = binxra1 - match
    # Determining false alarm for forecast array (forecast yes, obs no)
    arr2 = binxra2 - match
    # xra1 and xra2 are the same shape, determining total points in each array
    totalpts = binxra1.shape[0] * binxra1.shape[1]
# End class variables
    # Assigning a, b, and c arrays to csihash dictionary
    csihash['hits'] = match
    csihash['misses'] = arr1
    csihash['false_alarms'] = arr2

    # Find where data arrays have no information.
    if nodata_value != '':
        vpND = np.where(logical_or(xra1 == nodata_value, xra2 == nodata_value))
        maskND = xra2[:] * 0
        maskND[vpND] = 1
        vp_ra1ND = np.where(logical_and(maskND == 1, arr1 == 1))
        if vp_ra1ND[0] != []:
            arr11[vp_ra1ND] = 0
        vp_ra2ND = np.where(logical_and(maskND == 1, arr2 == 1))
        if vp_ra2ND[0] != []:
            arr2[vp_ra2ND] = 0

    if area.shape == match.shape:
        csihash['CSI'] = (match*area).sum() / ((match*area).sum() +
                                               (arr1*area).sum() + (arr2*area).sum())
        csihash['POD'] = (match*area).sum() / ((match*area).sum() + (arr1*area).sum())
        csihash['FAR'] = (arr2*area).sum() / ((match*area).sum() + (arr2*area).sum())

        if verbose == 1:
            print('used area')
            print((match*area).sum().values, (arr1*area).sum().values, (arr2*area).sum().values)
            print('CSI: ', csihash['CSI'].values, 'POD: ',
                  csihash['POD'].values, 'FAR: ', csihash['FAR'].values)
    else:
        csihash['CSI'] = match.sum() / (match.sum() + arr1.sum() + arr2.sum())
        # hit rate or probability of detection (p 310 Wilks)
        csihash['POD'] = match.sum() / (match.sum() + arr1.sum())
        # false alarm ratio (p 310 Wilks)
        csihash['FAR'] = arr2.sum() / (match.sum() + arr2.sum())

        if verbose == 1:
            print('Match: ', match.sum().values, 'Arr1: ', arr1.sum().values, 'Arr2: ', arr2.sum().values)
        print('CSI: ', csihash['CSI'].values, 'POD: ', csihash['POD'].values, 'FAR: ', csihash['FAR'].values)

    return csihash


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


def calc_pcorr(xra1, xra2, threshold=0., verbose=False):
    """
    Calculates Pattern Correlation between two arrays
    xra1 and xra2 need to be on the same grid. 
    See monetio.remap_nearest or monetio.remap_xesmf to remap arrays.

    Pattern Correlation (uncentered) = (binxra1 * binxra2) / sqrt((binxra1^2)*(binxra2^2))
    Pattern Correlation (centered) = ((binxra1 - arr1)(binxra2 - arr2)) / sqrt(((binxra1-arr1)^2) * ((binxra2 - arr2)^2))

    Inputs: 
    xra1:  observations/satellite (xarray DataArray)
    xra2:  the forecast/model (xarray DataArray)
    threshold: value for data threshold, default = 0. (float)
    verbose: boolean (True = print statements)

    Outputs:
    pcorr, pcorruncent: pattern correlation (centered), pattern correlation (uncentered)
    """

    # Convert xra1 and xra2 to binary arrays
    binxra1 = xr.where(xra1 >= threshold, 1., 0.)
    binxra2 = xr.where(xra2 >= threshold, 1., 0.)

    # Determining hits or matches (obs yes, forecast yes)
    match = binxra1 * binxra2
    # Determining misses for forecast array (obs yes, forecast no)
    arr1 = binxra1 - match
    # Determining false alarm for forecast array (forecast yes, obs no)
    arr2 = binxra2 - match
    # xra1 and xra2 are the same shape, determining total points in each array
    totalpts = binxra1.shape[0] * binxra1.shape[1]

    # Everything above this coud be in class as variables - not specific to pcorr calculation
    # Space averaged values
    arr1avg = (arr1.sum() + match.sum()) / float(totalpts)
    arr2avg = (arr2.sum() + match.sum()) / float(totalpts)

    # Calculating centered pattern correlation - subtracts space averaged values
    arr1corr = binxra1 - arr1avg
    arr2corr = binxra2 - arr2avg
    norm = ((arr1corr * arr1corr).sum()) ** 0.5 * ((arr2corr * arr2corr).sum())**0.5
    pcorr = (arr1corr * arr2corr).sum()/norm

    # Calculating uncentered pattern correlation
    norm = ((binxra1 * binxra1).sum()) ** 0.5 * ((binxra2 * binxra2).sum())**0.5
    pcorruncent = (match).sum()/norm

    if verbose == True:
        print('PCORR (centered)', pcorr.values)
        print('PCORR (uncentered)', pcorruncent.values)

    return pcorr, pcorruncent

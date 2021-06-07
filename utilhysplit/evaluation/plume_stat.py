# plume_stat.py

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from utilhysplit.evaluation import ensemble_tools as tools
import time

# from math import *
# import sys
# import string
# import datetime
# from itertools import permutations
# from math import pi, cos
"""
Routines to calculate various statistics like Critical Success Index, Fractions Skill Score,
Pattern Correlation, Brier Scores, and ensemble weighting schemes.
---------
Class:
    CalcScores
         calc_csi (critical success index, aka figure of merit in space)
         calc_fss (fraction skill score)
         calc_pcorr (pattern correlation coefficient)
Functions:
    calc_bs (brier score)
    calc_bss (brier skill score)
    calc_weightsBS (BS weights for ensemble members)
    calc_weightsPC (PC weights for ensemble members)
"""


class CalcScores:

    def __init__(self, xra1, xra2, threshold=0., szra=[1, 3, 5, 7], area=None, verbose=False):
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
        self.allpts = (self.xra1.shape[0] * self.xra1.shape[1])

        self.binxra1 = xr.where(self.xra1 >= self.threshold, 1., 0.)
        self.binxra2 = xr.where(self.xra2 >= self.threshold, 1., 0.)
        self.match = self.binxra1 * self.binxra2
        self.arr1 = self.binxra1 - self.match
        self.arr2 = self.binxra2 - self.match
        self.totalpts = (self.binxra1.shape[0] * self.binxra1.shape[1])

    def calc_csi(self):
        """ CSI equation: hits / (hits + misses + false alarms) - aka Threat Score
        GSS equation: (hits - chance) / ((hits + misses + false alarms) - chance)

        Inputs:
        match: hits
        arr1: misses
        arr2: false alarms
        area: optional array of grid areas, must be same size as xra1 and xra2 (default = None)
        nodata: flag for expanded array grid cells created if ash near boundary (string)
        verbose: boolean

        Output:
        csihash: Dictionary of values pertaining to Critical Success Index calculation
        contains: hits, misses, false_alarms, CSI, POD, FAR, events, total, freq
        """
        area = np.array(self.area)
        # Creating a csihash disctionary
        csihash = {}

        # Assigning a, b, and c arrays to csihash dictionary
        # Made these single values, rather than arrays - AMR 6/4/2021
        csihash['hits'] = self.match.sum().values
        csihash['misses'] = self.arr1.sum().values
        csihash['false_alarms'] = self.arr2.sum().values

        if area.shape == self.match.shape:
            csihash['CSI'] = ((self.match*self.area).sum() / ((self.match*self.area).sum() +
                                                              (self.arr1*self.area).sum() + (self.arr2*self.area).sum())).values
            csihash['POD'] = ((self.match*self.area).sum() /
                              ((self.match*self.area).sum() + (self.arr1*area).sum())).values
            csihash['FAR'] = (self.arr2*self.area).sum() / \
                ((self.match*self.area).sum() + (self.arr2*area).sum())
            # Added 6/3/21 - AMR
            csihash['Events'] = (self.match*self.area).sum() + (self.arr1*self.area).sum()
            csihash['Total'] = self.allpts
            csihash['Freq'] = csihash['Events'].values / csihash['Total'].values
            csihash['Posit'] = (self.match*self.area).sum() + (self.arr2*self.area).sum()
            csihash['Chance'] = csihash['Posit'].values*csihash['Freq'].values
            csihash['GSS'] = ((self.match*self.area).sum()-csihash['Chance'].values) / (csihash['Posit'].values +
                                                                                        csihash['Events'].values - (self.match*self.area).sum()) - csihash['Chance'].values
##

            if self.verbose == True:
                print('used area')
                print((self.match*self.area).sum().values, (self.arr1 *
                                                            self.area).sum().values, (self.arr2*self.area).sum().values, (self.totalpts*self.area).sum().values)
                print('CSI: ', csihash['CSI'].values, 'POD: ',
                      csihash['POD'].values, 'FAR: ', csihash['FAR'].values, 'Frequency: ', csihash['Freq'].values, 'Gilbert Skill Score: ', csihash['GSS'].values)
        else:
            csihash['CSI'] = (self.match.sum() / (self.match.sum() +
                                                  self.arr1.sum() + self.arr2.sum())).values
            # hit rate or probability of detection (p 310 Wilks)
            csihash['POD'] = (self.match.sum() / (self.match.sum() + self.arr1.sum())).values
            # false alarm ratio (p 310 Wilks)
            csihash['FAR'] = (self.arr2.sum() / (self.match.sum() + self.arr2.sum())).values
            # Added 6/3/21 - AMR
            csihash['Events'] = (self.match.sum() + self.arr1.sum()).values
            csihash['Total'] = self.allpts
            csihash['Freq'] = csihash['Events'] / csihash['Total']
            csihash['Posit'] = (self.match.sum() + self.arr2.sum()).values
            csihash['Chance'] = csihash['Posit']*csihash['Freq']
            csihash['GSS'] = ((self.match.sum()-csihash['Chance']) / (self.arr1.sum() +
                                                                      self.arr2.sum() + self.match.sum() - csihash['Chance'])).values

            if self.verbose == True:
                print('Match: ', self.match.sum().values, 'Misses: ',
                      self.arr1.sum().values, 'FalseAlarms: ', self.arr2.sum().values)
                print('CSI: ', csihash['CSI'].values, 'POD: ',
                      csihash['POD'].values, 'FAR: ', csihash['FAR'].values, 'Freq: ', csihash['Freq'].values, 'Gilbert Skill Score: ', csihash['GSS'].values)

        return csihash

    def calc_fss(self, szra=None, makeplots=False):
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
        # 2021 Jun 3 amc added random, uniform and afss to dataframe.

        # Creating FSS dictionary
        fss_dict = {}
        bigN = self.binxra1.size

        if not szra:
            szra = self.szra
        else:
            self.szra = szra

        # calculate frequency of observations and forecast.
        fobs = float(self.binxra1.sum()) / bigN
        fmod = float(self.binxra2.sum() / bigN)
        # random forecast has fss equal to frequency of observations
        random_fss = fobs
        # uniform forecast has fss equal to 0.5 + random
        uniform_fss = 0.5 + random_fss
        # measure of frequency bias
        # fss will asymptote at this value as neighborhood size
        # approaches domain size.
        afss = 2*fobs*fmod/(fobs**2+fmod**2)

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
                fig = plt.figure(1, figsize=(10, 15))
                ax1 = fig.add_subplot(1, 2, 1)
                ax2 = fig.add_subplot(1, 2, 2)
                cb = ax1.imshow(frac_arr1)
                cb2 = ax2.imshow(frac_arr2)
                # plt.colorbar(cb)
                # plt.colorbar(cb2)
                plt.show()

            fss_dict[sz] = fss_tmp

        df = pd.DataFrame.from_dict(fss_dict, orient='index')
        df['random'] = random_fss
        df['uniform'] = uniform_fss
        df['afss'] = afss
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


def calc_weightsBS(xra, scores):
    """
    Calculating the weighting scheme based on brier score values.
    Note, xra and scores should have the same source dimension.
    The scores should be for the correct applied threshold level.
    Inputs:
        xra: binary xarray of ensemble members (source, x, y)
        scores: scores for each ensemble member (xarray) from statistics netcdf
    Output:
        xraprob: ensemble relative frequency (xarray DataArray)
    """
    a = 0
    W = []
    xra2 = xra
    while a < len(xra.source):
        weight = 1-scores[a].values
        W.append(weight)
        xra2[a, :, :] = xra[a, :, :] * weight
        a += 1

    xraprob = xra2.sum(dim='source') / sum(W)
    return xraprob


def calc_weightsPC(xra, scores):
    """Kat62020

    Calculating the weighting scheme based on pattern correlation values.
    Note, xra and scores should have the same source dimension.
    The scores should be for the correct applied threshold level.
    Inputs:
        xra: binary xarray of ensemble members (source, x, y)
        scores: scores for each ensemble member (xarray DataArray) from statistics netcdf
    Output:
        xraprob: ensemble relative frequency (xarray DataArray)
    """
    a = 0
    W = []
    xra2 = xra
    while a < len(xra.source):
        weight = scores[a].values
        W.append(weight)
        xra2[a, :, :] = xra[a, :, :] * weight
        a += 1

    xraprob = xra2.sum(dim='source') / sum(W)
    return xraprob

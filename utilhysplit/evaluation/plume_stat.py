# plume_stat.py

import numpy as np
import pandas as pd
import xarray as xr
import scipy
import matplotlib.pyplot as plt
import time
from utilhysplit.evaluation import ensemble_tools
from utilhysplit.evaluation import statmain

"""
Routines to calculate various statistics like Critical Success Index, Gilbert Skill Score, Fractions Skill Score,
Pattern Correlation, Brier Scores, and ensemble weighting schemes.
---------
Class:
    CalcScores
         calc_accuracy_measures: calculates the scalar accuracy measures for continuous predictands
         calc_basics:
         calc_roc: calculates the relative operating characteristics
         get_contingency_table: puts contingency table together for further calculations
         table2csi: 
         calc_csi (critical success index, aka figure of merit in space)
         calc_fss (fraction skill score)
         calc_pcorr (pattern correlation coefficient)
   
Functions:
    calc_bs (brier score)
    calc_bss (brier skill score)
    calc_weights (calculates list of normalized weights for use in ensemble_tools.ATL)
    calc_weightsBS (BS weights for ensemble members)
    calc_weightsPC (PC weights for ensemble members)
    plot_roc (plots output from CalcScores.calc_roc)
"""
# 2021 Jun 9 amc added plot_roc function to plot output from the calc_roc method in CalcScores.
# AMC - note that we are using ATL from ensemble_tools and ensemble_tools also imports
# from plume_stat.py. Some of the functions need to be rearranged to avoid this situation.
# To go to multicategory forecasts need to add a threshold dimension onto binxra1 and binxra2

# 2021 Jul 2 amr added calc_weights to calculate a normalized list of weights for use with ensemble_tools.ATL function
# 2021 Jul 28 amc added calculation of bias and fractional bias to calc_accuracy_measures method.

def test_fssA():
    """
    create test inputs for CalcScores.
    specifically for testing calc_fss.
    """
    # just a 3x3 array
    sz = (3, 3)
    a = np.zeros(sz)
    b = np.zeros(sz)
    # only 2 observed values
    yval = [2, 0]
    xval = [2, 0]
    for val in zip(xval, yval):
        try:
            a[int(val[0]), int(val[1])] = 1
        except:
            pass
    # could use 2 modeled values.
    # yval2 = [1,0]
    # xval2 = [1,1]
    # only 1 modeled values
    yval2 = [1]
    xval2 = [1]
    for val in zip(xval2, yval2):
        try:
            b[int(val[0]), int(val[1])] = 1
        except:
            pass
    axra = xr.DataArray(a, dims={'x': xval, 'y': yval})
    bxra = xr.DataArray(b, dims={'x': xval2, 'y': yval2})
    cs = CalcScores(axra, bxra, threshold=0.1, pixel_match=False)
    df = cs.calc_fss(szra=[1, 3, 5, 9])
    return df, axra, bxra


def test_fssB():
    """
    create test inputs for CalcScores.
    by drawing from a gaussian distribution.
    """
    sz = (51, 51)
    mean = sz[0]/2.0
    f0 = 0.5  # frequency of obs
    fm = 0.1  # frequency of model
    a = np.zeros(sz)
    b = np.zeros(sz)
    n0 = int(f0 * sz[0]*sz[1])
    nm = int(fm * sz[0]*sz[1])
    yval = np.random.normal(mean, mean/4.0, n0)
    xval = np.random.normal(mean, mean/4.0, n0)
    for val in zip(xval, yval):
        try:
            a[int(val[0]), int(val[1])] = 1
        except:
            pass
    yval2 = np.random.normal(mean, mean/4.0, nm)
    xval2 = np.random.normal(mean, mean/4.0, nm)
    for val in zip(xval2, yval2):
        try:
            b[int(val[0]), int(val[1])] = 1
        except:
            pass
    axra = xr.DataArray(a, dims={'x': xval, 'y': yval})
    bxra = xr.DataArray(b, dims={'x': xval2, 'y': yval2})
    print(axra.sum(), bxra.sum())
    return axra, bxra


class CalcScores:

    def __init__(self, xra1, xra2, threshold=0., szra=[1, 3, 5, 7], area=None, verbose=False,
                 probabilistic=False, pixel_match=False, multi=False, clip=False):
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
        probabilistic : boolean. if True checks for 'ens' or 'source' dimension and creates probabilstic field instead of binary.
        pixel_match : boolean. if True calculate threshold for xra2 by matching number of pixels that are above
                      input threshold in xra1.
        multi: boolean. If True, calculates contingency table for all ensemble members as well as total ensemble.
        clip: boolean. If True, clips around forecast array, reducing number of correct no-forecasts
        ----------------------------
        Functions:
        calc_accuracy_measures: calculates the scalar accuracy measures for continuous predictands
        calc_basics:
        calc_roc: calculates the relative operating characteristics
        get_contingency_table: puts contingency table together for further calculations
        table2csi: 
        calc_csi: calculates the critical success index (gilbert score), probability of detection, false alarm rate, etc.
        calc_fss: calculates the fractions skill score
        calc_pcorr: calculates pattern correlation coefficient (Pearson correlation coefficient)
        """
        # 2021 Jun 3 amc if ens dimension present convert to probabilistic (0-1) field.
        # 2021 Jun 9 amc add pixel matching option for threshold.
        # 2021 Jun 9 amc add new function calc_basics
        # 2021 Jun 9 amc add new function calc_roc
        # 2021 Jun 9 amc add new function get_contingency_table. similar to calc_csi.
        # 2021 Jun 9 amc add self.arr3 which shows correctly forecast 0 values.
        # 2021 Jun 9 amc add 'd' to the csihash
        # 2021 Jun 17 amr added if statement to calculate binxra1/binxra2 if threshold == 0.
        # 2021 Jun 21 amr made get_contingency_table produce values for each ensemble member and total ensemble.
        # 2021 Jun 21 amr Added capability of calc_csi to use get_contingency_table pandas dataframe rather than dictionary.

        self.xra1 = xra1
        self.xra2 = xra2
        self.threshold = threshold
        self.szra = szra
        self.area = area
        self.verbose = verbose

        self.pm_threshold = None  # threshold from pixel matching, if any.

        # TO DO: Depends on input array
        # might need to specify which dimensions to use. Not a problem if using regridded volcat netcdf
        self.allpts = (self.xra1.shape[0] * self.xra1.shape[1])

        # process the  observational data
        if 'ens' in xra1.dims and probabilistic:
            self.binxra1 = ensemble_tools.ATL(xra1, thresh=threshold, norm=True)
        elif 'source' in xra1.dims and probabilistic:
            self.binxra1 = ensemble_tools.ATL(xra1, thresh=threshold, norm=True)
        else:
            if threshold == 0.:
                self.binxra1 = xr.where(self.xra1 > self.threshold, 1., 0.)
            else:
                self.binxra1 = xr.where(self.xra1 >= self.threshold, 1., 0.)

        # process the model data if pixel matching is desired
        if pixel_match:
            # if input is ensemble. Here assume probablistic output is wanted.
            if 'ens' in xra2.dims or 'source' in xra2.dims:
                self.pm_threshold, matchra = ensemble_tools.get_pixel_match(
                    xra2, xra1, threshold, return_binary=True)
                self.binxra2 = ensemble_tools.ATL(matchra, thresh=0.1, norm=True)
            # if input is deterministic
            else:
                self.pm_threshold = statmain.get_pixel_matching_threshold(xra1, xra2, threshold)
                self.binxra2 = xr.where(self.xra2 >= self.pm_threshold, 1., 0.)

        # process the model data with same threshold as observed.
        else:
            # if input is ensemble
            if 'ens' in xra2.dims and probabilistic:
                self.binxra2 = ensemble_tools.ATL(xra2, thresh=threshold, norm=True)
            elif 'source' in xra2.dims and probabilistic:
                self.binxra2 = ensemble_tools.ATL(xra2, thresh=threshold, norm=True)
            # if input is deterministic
            else:
                if threshold == 0.:
                    self.binxra2 = xr.where(self.xra2 > self.threshold, 1., 0.)
                else:
                    self.binxra2 = xr.where(self.xra2 >= self.threshold, 1., 0.)

        self.calc_basics()

    def calc_accuracy_measures(self, threshold=0, exclude_zeros=True):
        """
        scalar accuracy measures for continuous predictands.

        Computes statistics over the domain where either the model or the 
        observations is above the threshold. Does not use points where both are below threshold. 
        """
        # set below threshold values to zero.
        if threshold == 0:
            xra1 = xr.where(self.xra1 > threshold, self.xra1, 0)
            xra2 = xr.where(self.xra2 > threshold, self.xra2, 0)
        else:
            xra1 = xr.where(self.xra1 >= threshold, self.xra1, 0)
            xra2 = xr.where(self.xra2 >= threshold, self.xra2, 0)

        if exclude_zeros:
            # don't use 0,0 pairs.
            xra1 = xr.where((xra1 == 0) & (xra2 == 0), np.nan, xra1)
            xra2 = xr.where((xra1 == 0) & (xra2 == 0), np.nan, xra2)

        # array with all not nan values as 1.
        num = xr.where(~np.isnan(xra1), 1, 0)
        # number of non nan values.
        num = num.sum().values

        # modelmean
        mean1 = xra1.sum() / num
        mean2 = xra2.sum() / num

        # Bias
        bias = (xra1-xra2)
        bias = bias.sum()/num
        # See Stohl 1998 Stohl, A.; Hittenberger, M.; Wotawa, G. 
        # Validation of the Lagrangian particle dispersion model 
        # FLEXPART against large-scale tracer experiment data. 
        # Atmos. Environ. 1998, 32, 4245–4264.
        fracbias = 2*bias / (mean1+mean2)

        # Mean Absolute Error
        # perfect forecast 0.
        # typical magnitude for forecast error in a given verification data set
        maera = np.abs(xra1 - xra2)
        mae = maera.sum()/num

        # Mean Squared Error
        # more sensitive to outliers than MAE.
        msera = (xra1-xra2)**2
        mse = msera.sum() / num
        thash = {'MSE': [float(mse.values)], 'MAE': [float(mae.values)]}
        thash['threshold'] = threshold
        thash['exclude zeros'] = exclude_zeros
        thash['N'] = num
        thash['bias'] = [float(bias.values)]
        thash['fractional_bias'] = [float(fracbias.values)]
        tframe = pd.DataFrame.from_dict(thash)
        return tframe

    def calc_basics(self, probthresh=None, clip=False):
        """
        probthresh : int or float
        The probthresh can be used to convert probabilistic forecasts back to
        deterministic forecasts.
        This can be used for creating things like ROC diagrams.
        """

        if clip:
            # remove all x or y rows that are all 0's.
            temp = xr.concat([self.binxra1, self.binxra2], dim='temp')
            temp = xr.where(temp == 0, np.nan, temp)
            temp = temp.dropna(dim='x', how='all')
            temp = temp.dropna(dim='y', how='all')
            binxra2 = temp.isel(temp=1).fillna(0)
            binxra1 = temp.isel(temp=0).fillna(0)
        else:
            binxra1 = self.binxra1

        if isinstance(probthresh, (int, float)):
            # convert probabilistic forecast to deterministic using threshold.
            binxra2 = xr.where(self.binxra2 >= probthresh, 1.0, 0)
        else:
            binxra2 = self.binxra2
        self.match = binxra1 * binxra2
        self.arr1 = binxra1 - self.match
        self.arr2 = binxra2 - self.match
        self.arr3 = xr.where(self.match > 0, 0, 1)
        self.totalpts = binxra1.shape[0] * binxra1.shape[1]

    def calc_roc(self, clip=True, multi=False, problist = np.arange(0.05,1,0.10)):
        """
        For probabilistic forecasts.
        calculate the ROC (relative operating characteristic)

        Convert probabilistic forecast to binary using a probability threshold.
        Then compute False Alarm Rate and Hit Rate for various probability thresholds.
        ROC curve is the Hit Rate (y-axis) vs. False Alarm Rate (x-axis).
        See Wilks Chapter 8 (7 in 2nd edition).

        One issue with plume forecasts is the number of correclty forecast 0 values
        can be quite large and can be increased simply by increasing domain. This leads
        to a very small value of F (values on x axis).

        Should this be ameliorated by clipping the domain tightly around the plume area?
        clip=True will do this.

        problist : list of floats. Give probability thresholds for the points in the ROC curve.

        Outputs:
        xlist: list of x values for plotting
        ylist: list of y values for plotting
        """
        # Starts at (1,1) ends at (0.,0.) point
        xlist = [1.]
        ylist = [1.]
        # calculate False Alarm Rate (x axis) and
        # Hit Rate (y axis) for each probability threshold.
        for prob in problist:
            #self.calc_basics(prob, clip=clip)
            tframe = self.get_contingency_table(probthresh=prob, clip=clip, multi=False)
            tframe['F'] = tframe.apply(lambda row: row['b'] / (row['b'] + row['d']), axis=1)
            tframe['POD'] = tframe.apply(lambda row: row['a'] / (row['a'] + row['c']), axis=1)
            #csihash = self.calc_csi()
            # xlist.append(csihash['F'])
            # ylist.append(csihash['POD'])
            xlist.append(tframe['F'].values)
            ylist.append(tframe['POD'].values)
        xlist.append(0.)
        ylist.append(0.)
        xlist.reverse()
        ylist.reverse()
        area = scipy.integrate.trapz(y=ylist, x=xlist)
        return xlist, ylist, area[0]

    def get_contingency_table(self, probthresh=None, clip=False, multi=False, verbose=False):

        thash = []
        self.calc_basics(probthresh, clip)

        aval = self.match.sum().values
        cval = self.arr1.sum().values
        bval = self.arr2.sum().values
        dval = self.arr3.sum().values
        if verbose:
            print('a(hits) forecast yes, obs yes : {}'.format(aval))
            print('b(false alarm) forecast yes, obs no  : {}'.format(bval))
            print('c(misses) forecast no, obs yes : {}'.format(cval))
            print('d(correct no) forecast no, obs no : {}'.format(dval))
        thash = [{'a': aval, 'b': bval,
                  'c': cval, 'd': dval}]
        if multi:
            dimens = ['source', 'ens']
            for x in dimens:
                if x in self.binxra2.dims:
                    tval = str(x)
                    num = len(self.binxra2[tval])
                    # AMC - may want to change to using isel/sel for slicing.
                    # AMR - Making sure dimensions are ordered the same for all arrays
                    match = self.match.transpose(tval, 'y', 'x')
                    arr1 = self.arr1.transpose(tval, 'y', 'x')
                    arr2 = self.arr2.transpose(tval, 'y', 'x')
                    arr3 = self.arr3.transpose(tval, 'y', 'x')
                    for i in range(num):
                        element = str(self.binxra2[tval][i].values)
                        aval = float(match[i, :, :].sum().values)
                        cval = float(arr1[i, :, :].sum().values)
                        bval = float(arr2[i, :, :].sum().values)
                        dval = float(arr3[i, :, :].sum().values)
                        tmphash = {tval: element, 'a': aval, 'b':
                                   bval, 'c': cval, 'd': dval}
                        thash.append(tmphash)

        tframe = pd.DataFrame.from_dict(thash)
        if multi:
            tframe[tval] = tframe[tval].fillna('All')
        if isinstance(probthresh, (int, float)):
            tframe['probthresh'] = probthresh
        tframe['threshold'] = self.threshold
        if isinstance(self.pm_threshold, (int, float)):
            tframe['pm_threshold'] = self.pm_threshold
        return tframe

    def table2csi(self, tframe):
        tframe['CSI'] = tframe.apply(lambda row: row['a'] / (row['a'] + row['b'] + row['c']), axis=1)
        # false alarm ratio (p 310 Wilks) b/(a+b)
        # proportion of positive forecasts which were wrong.
        tframe['FAR'] = tframe.apply(lambda row: row['b'] / (row['a'] + row['b']), axis=1)
        tframe['F'] = tframe.apply(lambda row: row['b'] / (row['b'] + row['d']), axis=1)
        tframe['POD'] = tframe.apply(lambda row: row['a'] / (row['a'] + row['c']), axis=1)
        tframe['N'] = tframe.apply(lambda row: row['a'] + row['b'] + row['c'] + row['d'], axis=1)
        tframe['aref'] = tframe.apply(lambda row: (row['a'] + row['b'])*(row['a']+row['c'])/row['N'], axis=1)
        tframe['GSS'] = tframe.apply(lambda row: (row['a'] - row['aref']) /
                                     (row['a'] - row['aref'] + row['b'] + row['c']), axis=1)
        return tframe

    def calc_csi(self, verbose=False):
        """ CSI equation: hits / (hits + misses + false alarms) - aka Gilbert Score
        Inputs:
        match: hits
        arr1: misses
        arr2: false alarms
        area: optional array of grid areas, must be same size as xra1 and xra2 (default = None)
        nodata: flag for expanded array grid cells created if ash near boundary (string)
        multi: Boolean - calculations for each ensemble member and total ensemble
        dframe: Boolean - output as pandas dataframe if desired
        verbose: boolean
        Output:
        csihash: Dictionary of values pertaining to Critical Success Index calculation - can output dataframe if desired
        contains: hits, misses, false_alarms, CSI, POD, FAR
        """
        area = np.array(self.area)
        # Getting contingency table from get_contingency_table()
        #tframe = self.get_contingency_table(multi=multi)
        csihash = {}
        # Assigning a, b, and c arrays to csihash dictionary
        # Made these single values, rather than arrays - AMR 6/4/2021
        csihash['hits'] = self.match.sum().values
        csihash['misses'] = self.arr1.sum().values
        csihash['false_alarms'] = self.arr2.sum().values
        csihash['d'] = self.arr3.sum().values      # d correctly forecast no ash.

        if area.shape == self.match.shape:
            csihash['CSI'] = ((self.match*self.area).sum() / ((self.match*self.area).sum() +
                                                              (self.arr1*self.area).sum() + (self.arr2*self.area).sum())).values
            csihash['POD'] = ((self.match*self.area).sum() /
                              ((self.match*self.area).sum() + (self.arr1*area).sum())).values
            csihash['FAR'] = ((self.arr2*self.area).sum() /
                              ((self.match*self.area).sum() + (self.arr2*area).sum())).values
            # Added 6/3/21 - AMR
            csihash['Events'] = (self.match*self.area).sum() + (self.arr1*self.area).sum()
            csihash['Total'] = self.allpts
            csihash['Freq'] = csihash['Events'].values / csihash['Total'].values
            csihash['Posit'] = (self.match*self.area).sum() + (self.arr2*self.area).sum()
            csihash['Chance'] = csihash['Posit'].values*csihash['Freq'].values
            csihash['GSS'] = ((self.match*self.area).sum()-csihash['Chance'].values) / ((csihash['Posit'].values +
                                                                                         csihash['Events'].values - (self.match*self.area).sum()) - csihash['Chance'].values)

            if self.verbose == True:
                print('used area')
                print((self.match*self.area).sum().values, (self.arr1 *
                                                            self.area).sum().values, (self.arr2*self.area).sum().values)
                print('CSI: ', csihash['CSI'].values, 'POD: ',
                      csihash['POD'].values, 'FAR: ', csihash['FAR'].values, 'Frequency: ', csihash['Freq'].values, 'Gilbert Skill Score: ', csihash['GSS'].values)
        else:
            csihash['CSI'] = self.match.sum() / (self.match.sum() + self.arr1.sum() + self.arr2.sum())
            # hit rate or probability of detection (p 310 Wilks) a/(a+c)
            csihash['POD'] = self.match.sum() / (self.match.sum() + self.arr1.sum())
            # false alarm ratio (p 310 Wilks) b/(a+b)
            # proportion of positive forecasts which were wrong.
            csihash['FAR'] = self.arr2.sum() / (self.match.sum() + self.arr2.sum())
            # false alarm rate (p 311 Wilks) b/(d+b)
            # ratio of false alarms to total number of non-occurences.
            csihash['F'] = self.arr2.sum() / (self.arr2.sum() + self.arr3.sum())
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
            if verbose:
                print('CSI: {:.3f}'.format(csihash['CSI'].values),
                      'POD: {:.3f}'.format(csihash['POD'].values),
                      'FAR: {:.3f}'.format(csihash['FAR'].values),
                      'F  : {:.3f}'.format(csihash['F'].values),
                      'GSS  : {:.3f}'.format(csihash['GSS']))
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

        Note: according to Roberts 2008, the FSS should appproach AFSS as n-> 2N-1 where N is
        the domain size. It does. However the FSS at n=N can be larger than this value
        because of the inclusion of the zero padding at the edges. See the test_fssA function for an example.

        Return
        df : pandas dataframe """
        # 2021 Jun 3 amc added random, uniform and afss to dataframe.

        # Creating FSS dictionary
        fss_dict = {}
        bigN = self.binxra1.size

        if isinstance(szra, (int, float)):
            self.szra = [szra]
        elif isinstance(szra, (list, np.ndarray)):
            self.szra = szra

        # calculate frequency of observations and forecast.
        fobs = float(self.binxra1.sum()) / bigN
        fmod = float(self.binxra2.sum()) / bigN
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
            frac_arr1 = convolve2d(self.binxra1, filter_array, mode='same', fillvalue=0, boundary='fill')
            frac_arr2 = convolve2d(self.binxra2, filter_array, mode='same', fillvalue=0, boundary='fill')
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
                print('AFSS ', afss)
                print('N ', bigN, frac_arr1.size, frac_arr2.size)
                print('fobs ', fobs)
                print('fmod ', fmod)
                print('size ', sz)
            fss_tmp = dict({'Nlen': conv_array[0], 'FBS': fbs, 'FBS_ref': fbs_ref, 'FSS': fss})
            if makeplots == True:
                fig = plt.figure(1, figsize=(10, 15))
                ax1 = fig.add_subplot(1, 2, 1)
                ax2 = fig.add_subplot(1, 2, 2)
                cb = ax1.imshow(frac_arr1)
                cb2 = ax2.imshow(frac_arr2)
                plt.colorbar(cb, ax=ax1)
                plt.colorbar(cb2)
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


def calc_weights(scores, types='BS'):
    """
    Calculating the weighting scheme based on brier score or pattern correlation values.
    The scores xarray should be 1 dimension, either 'ens' or 'source' and
    should be for the desired applied threshold level.
    Inputs:
        scores: scores for each ensemble member (xarray) from statistics netcdf
        types: 9string) 'BS' for Brier Score or 'PC' for Pattern Correlation
    Output:
        wgtscore: normalized list of weights to apply to hysplit ensemble
    """
    a = 0
    scorelist = []
    while a < len(scores):
        if types == 'BS':
            weight = 1-scores.values[a]
        if types == 'PC':
            weight = scores.values[a]
        scorelist.append(weight)
        a += 1
    wgtscore = scorelist / sum(scorelist)
    return wgtscore


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
        weight = 1-scores.values[a]
        W.append(weight)
        xra2[a, :, :] = xra[a, :, :] * weight
        a += 1

    xraprob = xra2.sum(dim='source') / sum(W)
    return xraprob


def calc_weightsPC(xra, scores):
    """
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
        weight = scores.values[a]
        W.append(weight)
        xra2[a, :, :] = xra[a, :, :] * weight
        a += 1

    xraprob = xra2.sum(dim='source') / sum(W)
    return xraprob


def plot_roc(xlist, ylist, ax=None,clr='--ko',label=''):
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
    ax.plot(xlist, ylist, clr, label=label)
    ax.plot([0, 1], [0, 1], '-b')
    ax.set_xlabel('False Alarm Rate')
    ax.set_ylabel('Hit Rate')

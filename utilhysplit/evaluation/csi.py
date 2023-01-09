# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import datetime
import string
import sys
import time
from itertools import permutations
from math import *
from math import cos, pi

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#from xarray import DataArray as da
import xarray as xr


"""
routines to help calculates critical success index by a point by point comparison.
data must be gridded the same.
Functions
---------
match_arrays
get_area
find_threshold
mask_threshold
expand_array
expand_both
trim_array
trim_both
calc_csi (critical success index, aka figure of merit in space)
calc_fss  (fraction skill score)
"""


def match_arrays(ra1, lat1, lon1, ra2, lat2, lon2, missing_value=0, verbose=0):
    """inputs two arrays - adds padding zeros to make arrays fully overlap
    calls for regridding - MONET tool
    output arrays can be input into calc_csi, calc_fss
    ra1: xarray (data)
    lat1: xarray (latitude)
    lon1: xarray (longitude)
    ra2: xarray (data)
    lat2: xarray (latitude)
    lon2: xarray (longitude)
    missing value: 0 (default)
    verbose: 0 (default)
    return ra1new, lat1new, lon1new, ra2new, lat2new, lon2new
    """
    if verbose == 1:
        print('RA1 ', ra1.shape)
        print('RA2 ', ra2.shape)

       # check to make sure data array not already equal.
        if not (np.array_equal(lon1, lon2) and np.array_equal(lat1, lat2)):

            ###MATCH LATITUDES########
            lat_space = round(lat1[1][0] - lat1[0][0], 3)
            lat_space2 = round(lat2[1][0] - lat2[0][0], 3)
            if verbose == 1:
                print('Lat spacing', lat_space, lat_space2)
            # print lat1
            if lat_space == 0:
                print('WARNING. latitude spacing is zero!')
            if lat_space != lat_space2:
                print('WARNING: spacing between two arrays is not the same')
            minlat = np.amin(np.append(lat1, lat2))
            maxlat = np.amax(np.append(lat1, lat2))
            if verbose == 1:
                print('Min latitude ', minlat, ' ', np.amin(lat1),  ' ', np.amin(lat2))
                print('Max latitude ',  maxlat, ' ', np.amax(lat1),  ' ', np.amax(lat2))
            latlist = np.arange(minlat, maxlat+lat_space, lat_space)

            ##Add padding to ra1 ########
            la1 = int(round((np.amin(lat1)-np.amin(latlist))/lat_space))
            la2 = int(round((np.amax(latlist)-np.amax(lat1))/lat_space))
            tempra = np.zeros(len(lat1[0, :])).reshape(1, len(lat1[0, :]))
            tempra = tempra + missing_value
            for x in range(la2):
                ra1 = np.concatenate((ra1, tempra), axis=0)
            for x in range(la1):
                ra1 = np.concatenate((tempra, ra1), axis=0)
            # print 'RA1 \n' , ra1

            ##Add padding to ra2 ########
            la1 = int(round((np.amin(lat2)-np.amin(latlist))/lat_space))
            la2 = int(round((np.amax(latlist)-np.amax(lat2))/lat_space))
            # print la1, la2
            tempra = np.zeros(len(lat2[0, :])).reshape(1, len(lat2[0, :]))
            tempra = tempra + missing_value
            for x in range(la2):
                if verbose:
                    print('la2', x)
                ra2 = np.concatenate((ra2, tempra), axis=0)
            for x in range(la1):
                if verbose:
                    print('la1', x)
                ra2 = np.concatenate((tempra, ra2), axis=0)
            # print 'RA2 \n' , ra2

            ###MATCH LONGITUDES########
            lon_space = round(lon1[0][1] - lon1[0][0], 3)
            lon_space2 = round(lon2[0][1] - lon2[0][0], 3)
            if verbose == 1:
                print('LON space', lon_space, lon_space2)
            if lon_space == 0:
                print('WARNING. longitude spacing is zero!')
            if lon_space != lon_space2:
                print('WARNING: longitude spacing for two arrays is not the same.')
            minlon = np.amin(np.append(lon1, lon2))
            maxlon = np.amax(np.append(lon1, lon2))
            lonlist = np.arange(minlon, maxlon+lon_space, lon_space)
            if verbose == 1:
                print('Min longitude ', minlon, 'Max longitude ', maxlon)

            ##add padding to ra1 ################
            la1 = int(round((np.amin(lon1)-np.amin(lonlist))/lon_space))
            la2 = int(round((np.amax(lonlist)-np.amax(lon1))/lon_space))
            # print la1, la2
            tempra = np.zeros(len(ra1[:, 0])).reshape(1, len(ra1[:, 0]))
            tempra = tempra + missing_value
            # print 'ra shape'   , ra1.shape
            for x in range(la2):
                ra1 = np.concatenate((ra1, tempra.T), axis=1)
            for x in range(la1):
                ra1 = np.concatenate((tempra.T, ra1), axis=1)
            # print 'RA1 \n' , ra1

            ##add padding to ra2 ################
            la1 = int(round((np.amin(lon2)-np.amin(lonlist))/lon_space))
            la2 = int(round((np.amax(lonlist)-np.amax(lon2))/lon_space))
            # print la1, la2
            tempra = np.zeros(len(ra2[:, 0])).reshape(1, len(ra2[:, 0]))
            tempra = tempra + missing_value
            # print 'ra shape'   , ra2.shape
            for x in range(la2):
                ra2 = np.concatenate((ra2, tempra.T), axis=1)
            for x in range(la1):
                ra2 = np.concatenate((tempra.T, ra2), axis=1)
            # print 'RA2 \n' , ra2
            latlongrid = np.meshgrid(lonlist, latlist)
        else:
            latlongrid = np.meshgrid(lon1, lat1)
            # print latlongrid[0]
            # print latlongrid[1]

    else:
        print("Error: shapes of input arrays are incorrect.")
    if verbose == 1:
        print('Shapes of output arrays ra1, ra2, lat, lon', ra1.shape,
              ra2.shape, latlongrid[0].shape, latlongrid[1].shape)
    return latlongrid, ra1, ra2


def get_area_domain(dset, directory='./', write=0, radius=6378.137):
    import xarray as xr
    import numpy as np
    from math import pi, cos
    """Input is a VOLCAT dataset. Use after volcat.open_dataset(fname). 
       Provide directory for where area netcdf should be located. Default is cwd.
       Default to NOT write area array to netcdf file.
       Output is xarray of same size with corresponding area (km2) of each grid cell.
       Converts degrees to meters using a radius of 6378.137 km."""

    # d2r =  pi/180.0              #convert degrees to radians
    # d2km = radius * d2r     #convert degree latitude to kilometers.

    ash_mass = dset.ash_mass  # Pulls out ash mass array
    ash_mass = ash_mass[0, :, :]  # Removes time dimension
    lat = ash_mass.latitude
    lon = ash_mass.longitude
    fname = os.path.join(directory, 'area_whole_domain.nc')
    return compute_area(lat, lon, radius, write, fname)


def compute_area(lat, lon, radius=6378.137, write=0, fname=None):
    # need a faster way to do this.
    # possibly by shifting the arrays and subtracting?
    # it looks like numpy.roll will do this.
    d2r = pi/180.0  # convert degrees to radians
    d2km = radius * d2r  # convert degree latitude to kilometers.

    latrad = lat * d2r  # Creating latitude array in radians
    coslat = np.cos(latrad) * d2km * d2km  # Grouping constant multiplication outside of loop
    # Creates an array copy of ash_mass filled with the fill value
    area = xr.full_like(lat, 0)
    shape = np.shape(area)
    # Begins looping through each element of array
    i = 0
    while i < (shape[0] - 1):
        j = 0
        while j < (shape[1] - 1):
            area[i, j] = abs(lat[i, j] - lat[i+1, j]) * abs(abs(lon[i, j]) - abs(lon[i, j+1])) * coslat[i, j]
            j += 1
        i += 1
    if write == 1:
        # Reformatting array attributes before writing to netcdf
        area.name = 'area'
        area.attrs['long_name'] = 'area of each lat/lon grid box'
        area.attrs['units'] = 'km^2'
        area.to_netcdf(fname)  # Writes area array to netcdf
    return area


def get_area(dset, dset2):
    import xarray as xr
    import numpy as np
    from monet.util import volcat
    """MUST BE RUN AFTER GENERATING AREA FILE FOR FULL DOMAIN!
       Input is a VOLCAT dataset. Use after volcat.open_dataset(fname). 
       Input is also area array. Use xr.open_dataset(fname).
       Output is xarray sized to ash mass data present in VOLCAT file.
       Output is in km^2, will need conversion to m^2 if necessary."""
    area = dset2.area  # Pulling out area array
    ash_mass = dset.ash_mass  # Pulls out ash mass array
    ash_mass = ash_mass[0, :, :]  # Removes time dimension
    area2 = area.where(ash_mass != ash_mass._FillValue, drop=True)
    return area2


def calc_fss(ra1, ra2,
             threshold1=0, threshold2=0,
             szra=[1, 3, 5, 7],
             makeplots=False,
             verbose=0):
    from scipy.signal import convolve2d
    """Calculates the fraction skill score (fss) 
       See Robers and Lean (2008) Monthly Weather Review
       and Schwartz et al (2010) Weather and Forecasting
       for more information.
       
       ra1 : observations/satellite
       ra2 : the forecast/model
        
       Can plot fractions if desired (double check calculations)
       threshold1 = value for data  threshold
       threshold2 = value for model threshold
       szra is  a list of the number of pixels (neightborhood length) to use 
            in fractions calculation
            default is to use 1, 3, 5, 7 pixels size squares

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


def find_threshold(ra1, ra2, nodata_value=None):
    """
       Implement pixel matching.
       Base threshold on matching number of pixels in observations.
       ra1 is the satellite data.
       ra2 is model data"""
    mask1 = mask_threshold(ra1, threshold=0)
    #numpixels = mask1.size - mask1.sum()
    numpixels = mask1.sum()
    list2 = np.copy(ra2)
    if nodata_value:
        vpND = np.where(ra1 == nodata_value)
        # print 'vpND' , len(vpND[0]) , vpND
        vpNDB = np.where(list2[vpND] > 0)
        # print len(vpNDB[0])
        list2[vpND] = 0

    list2 = list2.reshape(list2.shape[0] * list2.shape[1])
    list2.sort()
    list2 = list2[::-1]
    np.set_printoptions()
    print('match treshold', numpixels, list2[numpixels], list2[numpixels + 1000], list2[numpixels-1000])
    print(mask1.size, mask1.sum(), numpixels)
    mask2 = mask_threshold(list2, threshold=0)
    vpi = np.where(list2 > 0)
    print(list2.size, mask2.sum(), np.min(list2[vpi]))
    return list2[numpixels]


def mask_threshold(ra1, threshold=0.):
    """input array. 
        Returns array with 1's where data above threshold
        and 0. if below or equal to threshold (nan = 0.)  
    """
    mask = ra1.fillna(0.)
    mask1 = mask.where(mask <= threshold, 1.)
    return mask1


def expand_array(array, pad=5):
    """Input xarray data array and desired radius of padding.
    Default padding is 10 indices on all 4 sides of array.
    Returns data array with latitudes and longitudes
    of padded array."""
    # Calculating total number of padded indexes
    pads = pad * 2
    # Extracting 1D arrays of latitude and longitude values
    rows = np.array(array.latitude[:, 0])
    cols = np.array(array.longitude[0, :])
    ydim = len(rows)
    xdim = len(cols)
    # Calculating lat and lon deltas, start and end lat/lon values
    # To generate new 1D arrays of expanded dimensions
    dlat = rows[1] - rows[0]
    slat = rows[0] - (dlat * pads)
    elat = rows[-1] + (dlat * (pads + 1))
    newrows = np.arange(slat, elat, dlat)
    dlon = cols[1] - cols[0]
    slon = cols[0] - (dlon * pads)
    elon = cols[-1] + (dlon * (pads + 1))
    newcols = np.arange(slon, elon, dlon)
    # Creating new lat/lon arrays
    newlat, newlon = np.meshgrid(newrows, newcols, indexing='ij')
    # Creating temporary array of zeros (expanded dimensions)
    data = np.zeros((ydim + (pads * 2), xdim + (pads * 2)))
    # Inserting data into new expanded array
    data[pads:(pads + ydim), pads:(pads + xdim)] = array
    # Converting new array to xarray (assigning coordinates)
    expanded = da(data)
    expanded = expanded.rename({'dim_0': 'y'})
    expanded = expanded.rename({'dim_1': 'x'})
    newlat = da(newlat)
    newlat = newlat.rename({'dim_0': 'y'})
    newlat = newlat.rename({'dim_1': 'x'})
    newlon = da(newlon)
    newlon = newlon.rename({'dim_0': 'y'})
    newlon = newlon.rename({'dim_1': 'x'})
    expanded = expanded.assign_coords(latitude=newlat)
    expanded = expanded.assign_coords(longitude=newlon)
    return expanded


def expand_both(array1, array2, pad=3):
    """Input xarray data arrays (1 and 2), padding radius.
    Determines largest common grid for the two arrays.
    Adds designated padding around data.
    Necessary step before regridding can occur.
    Returns two arrays with latitudes and longitudes."""
    pads = pad * 2
    # Determining 1D lat/lon values for each array
    row1 = np.array(array1.latitude[:, 0])
    col1 = np.array(array1.longitude[0, :])
    row2 = np.array(array2.latitude[:, 0])
    col2 = np.array(array2.longitude[0, :])
    # Creating arrays of min/max lat/lon values for both arrays
    yedge = [min(row1), min(row2), max(row1), max(row2)]
    xedge = [min(col1), min(col2), max(col1), max(col2)]
    dlat1 = row1[1] - row1[0]
    dlon1 = col1[1] - col1[0]
    dlat2 = row2[1] - row2[0]
    dlon2 = col2[1] - col2[0]
    # Determine expanded longitudes for both arrays
    slon1 = min(col1) - abs(dlon1 * (abs(ceil((min(xedge) - min(col1)) / dlon1)) + pads))
    elon1 = max(col1) + abs(dlon1 * (abs(ceil((max(xedge) - max(col1)) / dlon1)) + pads))
    slon2 = min(col2) - abs(dlon2 * (abs(ceil((min(xedge) - min(col2)) / dlon2)) + pads))
    elon2 = max(col2) + abs(dlon2 * (abs(ceil((max(xedge) - max(col2)) / dlon2)) + pads))
    # Determine expanded latitudes for both arrays
    slat1 = min(row1) - abs(dlat1 * (abs(ceil((min(yedge) - min(row1)) / dlat1)) + pads))
    elat1 = max(row1) + abs(dlat1 * (abs(ceil((max(yedge) - max(row1)) / dlat1)) + pads))
    slat2 = min(row2) - abs(dlat2 * (abs(ceil((min(yedge) - min(row2)) / dlat2)) + pads))
    elat2 = max(row2) + abs(dlat2 * (abs(ceil((max(yedge) - max(row2)) / dlat2)) + pads))
    # Determine orientation of original array  - maintain start and end for new array
    if abs(col1[-1] - slon1) < abs(col1[0] - slon1):
        tmp = slon1
        slon1 = elon1
        elon1 = tmp
    if abs(col2[-1] - slon2) < abs(col2[0] - slon2):
        tmp = slon2
        slon2 = elon2
        elon2 = tmp
    if abs(row1[-1] - slat1) < abs(row1[0] - slat1):
        tmp = slat1
        slat1 = elat1
        elat1 = tmp
    if abs(row2[-1] - slat2) < abs(row2[0] - slat2):
        tmp = slat2
        slat2 = elat2
        elat2 = tmp
    newcols1 = np.arange(slon1, elon1, dlon1)
    newcols2 = np.arange(slon2, elon2, dlon2)
    newrows1 = np.arange(slat1, elat1, dlat1)
    newrows2 = np.arange(slat2, elat2, dlat2)
    # Creating new lat/lon arrays
    newlat1, newlon1 = np.meshgrid(newrows1, newcols1, indexing='ij')
    newlat2, newlon2 = np.meshgrid(newrows2, newcols2, indexing='ij')
    # Creating temporary array of zeros (expanded dimensions)
    data1 = np.zeros((len(newrows1), len(newcols1)))
    data2 = np.zeros((len(newrows2), len(newcols2)))
    # Inserting data into new expanded array
    stlat1 = (abs(newrows1 - row1[0])).argmin()
    stlon1 = (abs(newcols1 - col1[0])).argmin()
    stlat2 = (abs(newrows2 - row2[0])).argmin()
    stlon2 = (abs(newcols2 - col2[0])).argmin()
    data1[stlat1:(stlat1 + len(row1)), stlon1:(stlon1 + len(col1))] = array1
    data2[stlat2:(stlat2 + len(row2)), stlon2:(stlon2 + len(col2))] = array2
    # Converting new array to xarray (assigning coordinates)
    expanded1 = da(data1)
    expanded2 = da(data2)
    expanded1 = expanded1.rename({'dim_0': 'y'})
    expanded1 = expanded1.rename({'dim_1': 'x'})
    expanded2 = expanded2.rename({'dim_0': 'y'})
    expanded2 = expanded2.rename({'dim_1': 'x'})
    newlat1 = da(newlat1)
    newlon1 = da(newlon1)
    newlat2 = da(newlat2)
    newlon2 = da(newlon2)
    newlat1 = newlat1.rename({'dim_0': 'y'})
    newlat1 = newlat1.rename({'dim_1': 'x'})
    newlon1 = newlon1.rename({'dim_0': 'y'})
    newlon1 = newlon1.rename({'dim_1': 'x'})
    newlat2 = newlat2.rename({'dim_0': 'y'})
    newlat2 = newlat2.rename({'dim_1': 'x'})
    newlon2 = newlon2.rename({'dim_0': 'y'})
    newlon2 = newlon2.rename({'dim_1': 'x'})
    expanded1 = expanded1.assign_coords(latitude=newlat1)
    expanded1 = expanded1.assign_coords(longitude=newlon1)
    expanded2 = expanded2.assign_coords(latitude=newlat2)
    expanded2 = expanded2.assign_coords(longitude=newlon2)
    return expanded1, expanded2


def trim_array(array, pad=5):
    """Trims down array to array size that encompasses 
    data and minor padding for input array.
    Inputs: array to trim, radius of padding around array."""
    pad = pad * 2
    rows = da.any(array, axis=1)
    cols = da.any(array, axis=0)
    ymin, ymax = np.where(rows)[0][[0, -1]]
    xmin, xmax = np.where(cols)[0][[0, -1]]
    # Trims array
    trim = array[(ymin - pad):(ymax + 1 + pad), (xmin - pad):(xmax + 1 + pad)]
    return trim


def trim_both(mask1, mask2, pad=5):
    """Trims two arrays to array size 
    that encompasses data and minor padding for both input arrays. 
    Requires pixel radius (pad) value to determine padding to keep in array. 
    Inputs: array1, array2, radius of padding around array."""
    pad = pad * 2
    print(pad)
    rows1 = np.any(mask1, axis=1)
    cols1 = np.any(mask1, axis=0)
    rows2 = np.any(mask2, axis=1)
    cols2 = np.any(mask2, axis=0)
    ymin1, ymax1 = np.where(rows1)[0][[0, -1]]
    xmin1, xmax1 = np.where(cols1)[0][[0, -1]]
    ymin2, ymax2 = np.where(rows2)[0][[0, -1]]
    xmin2, xmax2 = np.where(cols2)[0][[0, -1]]
    trim1 = mask1[(ymin1 - pad):(ymax1 + 1 + pad), (xmin1 - pad):(xmax1 + 1 + pad)]
    trim2 = mask2[(ymin2 - pad):(ymax2 + 1 + pad), (xmin2 - pad):(xmax2 + 1 + pad)]
    # Determines minimum yarr and xarr between the two arrays
    yarr = [ymin1, ymin2, ymax1, ymax2]
    xarr = [xmin1, xmin2, xmax1, xmax2]
    # Trims each array to size
    trimboth1 = mask1[(min(yarr) - pad):(max(yarr) + 1 + pad), (min(xarr) - pad):(max(xarr) + 1 + pad)]
    trimboth2 = mask2[(min(yarr) - pad):(max(yarr) + 1 + pad), (min(xarr) - pad):(max(xarr) + 1 + pad)]
    return trimboth1, trimboth2


def calc_csi(ra1, ra2, area=None, nodata_value='', threshold=0, verbose=0):
    """ra1 and ra2 need to be on the same grid. 
        See monetio.remap_nearest or monetio.remap_xesmf to remap arrays.
        ra1 = observations/satellite
        ra2 = the forecast/model
        area = optional array of grid areas
        nodata_value = flag for expanded array grid cells created if ash near boundary
        threshold = value for data threshold
        CSI equation: hits / (hits + misses + false alarms)"""
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
            #print('HERE' , matchra.size.values, matchra.shape.values, onlyra2.shape.values)
            print(matchra.sum().values, onlyra1.sum().values, onlyra2.sum().values)
    print('CSI POD FAR', csihash['CSI'].values, csihash['POD'].values, csihash['FAR'].values)

    return csihash

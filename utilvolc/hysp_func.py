# hysp_func.py
# Functions for manipulating HYSPLIT data
# For use with MONET
"""Functions for manipulating HYSPLIT data.
-------------
Functions:
-------------
hysp_heights: determines ash top height from HYSPLIT
hysp_massload: determines total mass loading from HYSPLIT
calc_MER: determines Mass Eruption Rate from HYSPLIT
calc_aml: determines ash mass loading for each altitude layer  from HYSPLIT
hysp_thresh: calculates mask array for ash mass loading threshold from HYSPLIT
"""
from monetio.models import hysplit
from utilhysplit.utilvolc import volcMER
import xarray as xr
import numpy as np


def get_latlongrid(dset, xindx, yindx):
    llcrnr_lat = dset.attrs['Concentration Grid']['llcrnr latitude']
    llcrnr_lon = dset.attrs['Concentration Grid']['llcrnr longitude']
    nlat = dset.attrs['Concentration Grid']['Number Lat Points']
    nlon = dset.attrs['Concentration Grid']['Number Lon Points']
    dlat = dset.attrs['Concentration Grid']['Latitude Spacing']
    dlon = dset.attrs['Concentration Grid']['Longitude Spacing']

    lat = np.arange(llcrnr_lat, llcrnr_lat + nlat * dlat, dlat)
    lon = np.arange(llcrnr_lon, llcrnr_lon + nlon * dlon, dlon)
    print(nlat, nlon, dlat, dlon)
    print('lon shape', lon.shape)
    print('lat shape', lat.shape)
    print(lat)
    print(lon)
    lonlist = [lon[x - 1] for x in xindx]
    latlist = [lat[x - 1] for x in yindx]
    mgrid = np.meshgrid(lonlist, latlist)
    return mgrid


def getlatlon(dset):
    """
    dset : xarray returned by hysplit.open_dataset function
    RETURNS
    lat : 1D array of latitudes
    lon : 1D array of longitudes
    """
    llcrnr_lat = dset.attrs['Concentration Grid']['llcrnr latitude']
    llcrnr_lon = dset.attrs['Concentration Grid']['llcrnr longitude']
    nlat = dset.attrs['Concentration Grid']['Number Lat Points']
    nlon = dset.attrs['Concentration Grid']['Number Lon Points']
    dlat = dset.attrs['Concentration Grid']['Latitude Spacing']
    dlon = dset.attrs['Concentration Grid']['Longitude Spacing']
    lat = np.arange(llcrnr_lat, llcrnr_lat + nlat * dlat, dlat)
    lon = np.arange(llcrnr_lon, llcrnr_lon + nlon * dlon, dlon)
    return lat, lon


def hysp_massload(dset, threshold):
    """ Calculate mass loading from HYSPLIT xarray
    Inputs: xarray, ash mass loading threshold (threshold = xx)
    Outputs: total ash mass loading (summed over all layers), ash mass loading
    Units in (unit mass / m^2)"""
    aml_alts = calc_aml(dset)
    total_aml = aml_alts.sum(dim='z')
    # Calculate conversion factors
    unitmass, mass63 = calc_MER(dset)
    # Calculating the ash mass loading
    total_aml2 = total_aml * mass63
    # Calculating total ash mass loading, accounting for the threshold
    # Multiply binary threshold mask to data
    total_aml_thresh = hysp_thresh(dset, threshold)
    total_aml = total_aml2 * total_aml_thresh
    return total_aml


def hysp_heights(dset, threshold):
    """ Calculate ash top-height from HYSPLIT xarray
    Inputs: xarray, ash mass loading threshold (threshold = xx)
    Outputs: ash top heights, altitude levels """
    # Applying ash mass threshold when calculating ash mass loading
    aml_alts = calc_aml(dset)
    # Create array of 0 and 1 (1 where data exists)
    heights = aml_alts.where(aml_alts == 0., 1.)
    # Multiply each level by the altitude
    alts = dset.coords['z']
    height = _alt_multiply(heights, alts)
    height = height/1000.  # convert to km
    # Determine top height: take max of heights array along z axis
    top_hgt = height.max(dim='z')
    # Apply ash mass loading threshold mask array
    total_aml_thresh = hysp_thresh(dset, threshold)
    top_height = top_hgt * total_aml_thresh
    return top_height


def calc_MER(dset):
    """ Calculate Mass Eruption Rate (kg/s) using Mastin's 2009 equation
    Then converting to units of g/hr then g/m^2 """
    # Extract starting locations and plume height - plume height needed for Mass calc.
    stlocs = dset.attrs['Starting Locations']
    plume = (stlocs[-1][-1] - stlocs[0][-1])/1000  # Put height in km
    # Calculate mass erution rate (Mastin 2009)
    MER = volcMER.mastinMER(plume)
    # Calculating gram equivalent of 1 unit mass
    unitmass, mass63 = volcMER.MER2unit(MER)  # unit: g/hr
    return unitmass, mass63


def calc_aml(dset):
    """ Calculates the ash mass loading at each altitude for the dataset
    Input: xarray
    Output: total ash mass loading """
    # Totals values for all particles
    total_par = _add_species(dset)
    alts = total_par.coords['z']
    # Multiplies the total particles by the altitude layer
    # to create a mass loading for each altitude layer
    aml_alts = _delta_multiply(total_par, alts)
    return aml_alts


def hysp_thresh(dset, threshold):
    """ Calculates a threshold mask array based on the
    ash mass loading from HYSPLIT xarray
    Inputs: xarray, ash mass loading threshold (threshold = xx)
    Outputs: ash mass loading threshold mask array """
    # Calculate ash mass loading for xarray
    aml_alts = calc_aml(dset)
    total_aml = aml_alts.sum(dim='z')
    # Calculate conversion factors
    unitmass, mass63 = calc_MER(dset)
    # Calculating the ash mass loading
    total_aml2 = total_aml * mass63
    total_aml_thresh = total_aml2.where(total_aml2 > threshold, 0.)
    total_aml_thresh = total_aml_thresh.where(total_aml_thresh <= threshold, 1.)
    return total_aml_thresh


def _add_species(dset):
    # Calculate sum of particles
    species = dset.attrs["Species ID"]
    s = 0
    tmp = []
    # Looping through all species in dataset
    while s < len(species):
        tmp.append(dset[species[s]].fillna(0))
        s += 1  # End of loop through species
    total_par = tmp[0]
    p = 1
    # Adding all species together
    while p < len(tmp):
        total_par = total_par + tmp[p]
        p += 1  # End of loop adding all species
    return total_par


def _delta_multiply(pars, alts):
    # Calculate the delta altitude for each layer
    x = 1
    delta = []
    delta.append(alts[0])
    while x < (len(alts)):
        delta.append(alts[x] - alts[x-1])
        x += 1
    # Multiply each level by the delta altitude
    y = 0
    while y < len(delta):
        pars[:, y, :, :] = pars[:, y, :, :] * delta[y]
        y += 1  # End of loop calculating heights
    return pars


def _alt_multiply(pars, alts):
    # For calculating the top height
    # Multiply "1s" in the input array by the altitude
    y = 0
    while y < len(alts):
        pars[:, y, :, :] = pars[:, y, :, :] * alts[y]
        y += 1  # End of loop calculating heights
    return pars

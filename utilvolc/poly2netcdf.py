import xarray as xr
import shapely.geometry as sgeo
from affine import Affine
import numpy as np
from rasterio import features
import matplotlib.pyplot as plt


def transform_from_latlon(lat, lon):
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    trans = Affine.translation(lon[0], lat[0])
    scale = Affine.scale(lon[1]-lon[0], lat[1]-lat[0])
    return trans * scale


class VaacPoly():

    def __init__(self, pts, time, value, tag, height):
        """
        Helper class for poly2xra
        pts : set of tuples defining a polygon
              (longitude,latitude)
        time : datetime object
               time that polygon is valid for.
        value : float or int
        tag  : indicates confidence 
        height : in km or FL. 
                'FL100' or 15
        """
        self.polygon = sgeo.Polygon(pts)
        self.time = time
        self.tag = tag
        self.height = self.checkht(height)

    def plot(self):
        x,y = self.polygon.exterior.xy
        plt.plot(x,y)
 

    def checkht(self, height):
        if isinstance(height, str) and 'fl' in height.lower():
            flight_level = int(height.lower().replace('fl', ''))
            height = flight_level * 0.03048  # convert to km
        elif isinstance(height, str):
            height = float(height)
        return height


def match_dates(dlist1, dlist2):
    # get dates in both lists
    matched = [x for x in dlist1 if x in dlist2]
    # dates in dlist1 that are not matched
    not1 = [x for x in dlist1 if x not in m1]
    # dates in dlist2 that are not matched
    not2 = [x for x in dlist2 if x not in m1]
    return matched, not1, not2


def process_xra(cxra, globalhash={}, atthash={'jobid': '999'}):
    # take an xarray with source and ens dimensions and make a dataset
    # with variables.
    newxra = cxra.stack(newvar=('ens', 'source'))
    newxra = newxra.to_dataset('newvar')
    newvars = newxra.data_vars.keys()
    newnames = {}
    attr1 = {}
    # create attributes for the variable
    iii = 0
    for var in newvars:
        #str1 = str(var[1])
        #print(var, type(var[1]))
        if not isinstance(var[1], str):
            str1 = 'DataInsertion{}'.format(iii)
            str2 = 'DataInsertion'
            atthash['InsertionTime'] = str(var[1])
        else:
            str1 = str(var[1])
        # if isinstance(var[1], np.datetime64):
        #   str1 = 'DataInsertion{}'.format(iii)
        #   atthash['InsertionTime'] = str(var[1])
        #   print('here ', str1)
        # elif isinstance(var[1],str):
        #   str1 = var[1]
        #   print('Bhere ', str1)
        # else:
        #   str1 = ''
        #   print('Chere ', str1)
        newnames[var] = str.join('_', [var[0], str1])
        # print(newnames)
        newxra[var].attrs.update({'Institution': 'NOAA/ARL'})
        newxra[var].attrs.update({'Met_id': var[0]})
        newxra[var].attrs.update({'Source_id': str1})
        newxra[var].attrs.update({'Model': 'HYSPLIT'})
        newxra[var].attrs.update({'Units': 'g/m^3'})
        newxra[var].attrs.update({'FillValue': 'nan'})
        newxra[var].attrs.update({'Run_Type': str2})
        newxra[var].attrs.update({'Creator_Name': 'Alice Crawford'})
        newxra[var].attrs.update({'Creator_Email': 'alice.crawford@noaa.gov'})
        newxra[var].attrs.update(atthash)
        iii += 1

    # rename the stacked variables
    newxra = newxra.rename(newnames)
    newxra.attrs.update(globalhash)
    return newxra


class Poly2xra:

    def __init__(self):
        self.sp1 = 'x'
        self.sp2 = 'y'

        # create these in create_transform
        self.outshape = None
        self.transform = None

    def create3(self,vaac_poly_list, cxra,globalhash={},atthash={}):
        """
        vaac_poly_list : list of VaacPoly objects
        cxra : xarray dataset
        """

        ns2 = self.create1(vaac_poly_list, cxra)
        dset = cxra
        newra = xr.merge([dset,ns2])
        newra.attrs.update(globalhash)
        return newra 

    def create2(self,vaac_poly_list, cxra,globalhash={},atthash={}):
        """
        vaac_poly_list : list of VaacPoly objects
        cxra : xarray data-array
        """
        ns2 = self.create1(vaac_poly_list, cxra)
        dset = process_xra(cxra, globalhash, atthash)
        newra = xr.merge([dset, ns2])
        newra.attrs.update(globalhash)
        return newra

    def create_transform(self, cxra):
        """
        cxra : DataArray with desired lat, lon shape.
        """
        sp1 = 'x'
        sp2 = 'y'
        #spatial_coords = {sp1:cxra.coords[sp1],sp2:cxra.coords[sp2]}

        latitude = cxra.latitude
        longitude = cxra.longitude
        #poly1= sgeo.Polygon(pts)
        outshape = latitude.shape

        latitude = cxra.latitude[:, 0]
        longitude = cxra.longitude[0]
        transform = transform_from_latlon(latitude, longitude)
        self.outshape = outshape
        self.transform = transform
        return outshape, transform, latitude, longitude

    def create1(self,
                vaac_poly_list,
                cxra,
                fill=0,
                tag='poly100',
                coordhash={}):
        """
        vaac_poly_list : list of VaacPoly objects
        times should match the times in the cxra
        cxra : xarray object with model data

        fill : value to fill in area around polygon

        RETURNS:
        ns2 : DataArray with VAAC polygon.
        """
        sp1 = self.sp1
        sp2 = self.sp2
        spatial_coords = {sp1: cxra.coords[sp1], sp2: cxra.coords[sp2]}
        self.create_transform(cxra)
        iii = 0

        # create the DataArray
        for vpoly in vaac_poly_list:

            poly1 = vpoly.polygon
            raster = features.rasterize([poly1], out_shape=self.outshape,
                                     fill=fill, transform=self.transform,
                                     dtype=float)
            newxra = xr.DataArray(raster,coords=spatial_coords, dims=(sp2,sp1))
            print('adding', vpoly.time,iii)
            coordhash['time'] = vpoly.time
            ns2 = newxra.assign_coords(coordhash)
            print('max value', np.max(ns2))
            for coordval in coordhash.keys():
                ns2 = ns2.expand_dims(coordval)

            # expand in the z dimension
            ns2 = self.add_height(ns2, cxra, vpoly.height)

            if iii == 0:
                polyra = ns2
            else:
                polyra = xr.concat([polyra, ns2], 'time')
            iii += 1

        polyra.name = tag
        return polyra.transpose('time', 'z', 'y', 'x')

    def add_height(self, ns2, cxra, polyheight):
        iii = 0
        concatlist = []
        for zval in cxra.z.values:
            ntemp = ns2.assign_coords({'z': zval})
            # if above top heigh tof polygon set values to 0
            # if zval > polyheight:

            #   ntemp = ntemp * 0
            concatlist.append(ntemp)
        # concatenate along z dimension.
        return xr.concat(concatlist, 'z')


"""
    def create1(self,vaac_poly_list, cxra, fill=0, tag='poly100', sourceval='human', ensval='E1'):
    #vaac_poly_list : list of VaacPoly objects
    #times should match the times in the cxra

    #cxra : xarray object with model data

    # pts: list of tuples with polygon points.
    sp1 = 'x'
    sp2 = 'y'
    spatial_coords = {sp1:cxra.coords[sp1],sp2:cxra}

    latitude = cxra.latitude  
    longitude = cxra.longitude
    poly1= sgeo.Polygon(pts) 
    outshape = latitude.shape  
    
    latitude = cxra.latitude[:,0]
    longitude = cxra.longitude[0]
    transform = transform_from_latlon(latitude,longitude)
    iii=0 

    # take the cxra and make a dataset
    for vpoly in  vaac_poly_list:
        
        poly1 =  vpoly.polygon 
        raster = features.rasterize([poly1], out_shape=outshape,
                                 fill=fill, transform=transform,
                                 dtype=float)
        newxra = xr.DataArray(raster,coords=spatial_coords, dims=(sp2,sp1))

        coordhash = {'source':soureval,
                     'ens' : ensval,
                     'time': vpoly.time} 
        ns2 =
        newxra.assign_coords(coordhash)
        ns2 = ns2.expand_dims('ens')
        ns2 = ns2.expand_dims('source')
        ns2 = ns2.expand_dims('time')
        # concatenate along time axis.
        if iii==0:
           polyra = ns2
        else:
           polyra = xr.concat([polyra,ns2],'time')
        iii+=1
    ns2.name = tag

    # concatenate along source dimension
    combined_xra = xr.concat([cxra,ns2],'source')

    # could possibly split into a DataSet 
    dset = combined_xra.to_dataset(dim='source')
 
 
"""

import xarray as xr
import shapely.geometry as sgeo
from affine import Affine

def transform_from_latlon(lat,lon):
    lat = np.asarray(lat)
    lon = np.asaray(lon)
    trans = Affine.translation(lon[0],lat[0])
    scale = Affine.scale(lon[1]-lon[0],lat[1]-lat[0])
    return trans * scale


class VaacPoly():

    def __init__(self, pts, time, value, tag, height):
        """
        Helper class for poly2xra
        pts : set of tuples defining a polygon
              (latitude, longitude)
        time : datetime object
               time that polygon is valid for.
        value : float or int
        tag  : indicates confidence 
        height : in km or FL. 
                'FL100' or 15
        """
        self.polygon = sgeo.Polygon(pts)
        self.time = time
        self.tag =  tag
        self.height = self.checkht(height)

    def checkht(self, height):
        if isinstance(height, str) and 'fl' in height.lower():
           flight_level = int(height.lower().replace('fl',''))
           height = flight_level * 0.03048  #convert to km
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

def process_xra(cxra):
    # take an xarray with source and ens dimensions and make a dataset
    # with variables.
    newxra = cxra.stack(newvar = ('ens','source'))
    newxra = newxra.to_dataset('newvar')
    newvars = newxra.data_vars.keys()
    newnames = {}
    attr1 = {}
    # create attributes for the variable
    for var in newvars:
        newnames[var] = str.join('_',var)
        newxra[var].attrs.update('metid',var[0])
        newxra[var].attrs.update('sourceid',var[1])
        newxra[var].attrs.update('model','HYSPLIT')
    # rename the stacked variables
    newxra = newxra.rename(newnames) 
    return newxra 


class Poly2xra:

    def __init__(self):
        self.sp1 = 'x'
        self.sp2 = 'y'
  
        # create these in create_transform 
        self.outshape = None
        self.transform = None
 
    def create_transform(self,cxra):
        """
        cxra : DataArray with desired lat, lon shape.
        """
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
        self.create_transform(cxra)
        iii=0 

        # create the DataArray
        for vpoly in  vaac_poly_list:
            
            poly1 =  vpoly.polygon 
            raster = features.rasterize([poly1], out_shape=outshape,
                                     fill=fill, transform=transform,
                                     dtype=float)
            newxra = xr.DataArray(raster,coords=spatial_coords, dims=(sp2,sp1))

            coordhash['time'] = vpoly.time
            ns2 = newxra.assign_coords(coordhash)

            for coordval in coordhash.keys():
                ns2 = ns2.expand_dims(coordval)
          
            # expand in the z dimension 
            ns2 = self.add_height(ns2,cxra, vpoly.height)        

            if iii==0:
               polyra = ns2
            else:
               polyra = xr.concat([polyra,ns2],'time')
            iii+=1

        ns2.name = tag
        return ns2

    def add_height(self, ns2, cxra, polyheight):
        iii=0
        concatlist = []
        for zval in cxra.z.values:
            ntemp = ns2.assign_coords({'z':zval})
            # if above top heigh tof polygon set values to 0
            if zval > polyheight:
               ntemp = ntemp * 0 
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

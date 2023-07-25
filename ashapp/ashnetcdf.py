import datetime
import numpy as np
import os
import logging
import xarray as xr
from monetio.models import hysplit
from ashapp import utils
from utilvolc.runhelper import Helper
logger = logging.getLogger(__name__)
utils.setup_logger()

class AshAttributes:

    def __init__(self, attr={},dstr="%Y-%m-%d %H:%M:%S"):
        self.dstr = dstr
        self.attr = attr
        #self._attr = check_attributes(attr,dstr=self.dstr)
        self._keylist = ['mult',
                        'meteorological data',
                        'jobname',
                        'durationOfSimulation',
                        'latitude',
                        'longitude',
                        'eflag']
       
     

    def __call__(self):
        return self._attr

    def __str__(self):
        return str(self._attr)

    @property
    def dstr(self):
        return self._dstr

    @dstr.setter
    def dstr(self,indstr):
        testdate = datetime.datetime(2020,2,1,0)
        default = indstr
        self._dstr = indstr
        # make sure it is valid
        try:
           tstr = testdate.strftime(indstr)
           test = datetime.datetime.strptime(tstr,self._dstr) 
        except ValueError:
           self._dstr = default  
           raise ValueError('invalid input for dstr') 

    @property
    def attr(self):
        return self._attr

    @attr.setter
    def attr(self,inp):
        if not isinstance(inp,dict):
           logger.warning('input should be a dictionary. converting to dictionary.')
           temp = {}
           temp['INPUT1'] = inp
           inp = temp
        inp2 = check_attributes(inp,self.dstr)
        self._attr=inp2

    def update(self,inp):
        if not isinstance(inp,dict):
           logger.warning('input should be a dictionary. converting to dictionary.')
           temp = {}
           keys = self.attr.keys()
           iii = 1
           newkey = 'INPUT1'
           while newkey in keys:
                 newkey = 'INPUT{}'.format(iii)
                 iii+=1
           temp[newkey] = inp
           inp = temp
        inp2 = check_attributes(inp,self.dstr)
        self._attr.update(inp2) 

    def reset(self):
        self._attr = {}

    @property
    def keylist(self):
        return self._keylist

    @keylist.setter
    def keylist(self,klist):
        if isinstance(klist,list):
           self._keylist = klist
        elif isinstance(klist,(str,np.ndarray)):
           self._keylist = list(klist)

    def check(self):
        nflist = []
        for key in self._keylist:
            if key not in self._attr.keys():
               logger.warning('key not found in attributes {}'.format(key))
               nflist.append(key)
        return nflist 

    #def inp2attr(self):
    #    """
    #    return dictionary where all values are strings.
    #    This is so it can be written as attributes to netcdf file.
    #    """
    #    atthash = {}
    #    for key in self.inp.keys():
    #        try:
    #            val = str(self.inp[key])
    #            logger.info("attribute can be written as string {}".format(key))
    #            logger.info("attribute can be written as string {}".format(val))
    #        except:
    #            logger.warning("attribute cannot be written as string {}".format(key))
    #            val = "skip"
    #        if "/" in val:
    #            val = "skip"
    #        if val != "skip":
    #            atthash[key] = val
    #    return atthash

def check_attributes(atthashin,dstr="%Y-%m-%d %H:%M:%S"):
    """
    Make sure all attributes have forms that can be written to netcdf well.
    convert all np.ndarray to lists.
    convert all datetimes to strings.
    convert booleans to strings.
    """
    atthash = atthashin.copy()
    for key in atthash.keys():
        val = atthash[key]
        if isinstance(val, np.ndarray):
            newval = list(val)
            atthash[key] = newval
        elif isinstance(val,datetime.datetime):
            newval = val.strftime(dstr)
            atthash[key] = newval  
        elif isinstance(val,bool):
            newval = str(val)
            atthash[key] = newval  
        elif isinstance(val,type(None)):
            newval = 'None'
            atthash[key] = newval  
        elif isinstance(val,dict):
            newval = check_attributes(val,dstr)
            atthash[key] = newval  
    return atthash

    

class HYSPLITAshNetcdf:

    def __init__(self, fname):
        self.fname = fname

        self._cxra = xr.DataArray()
        self._attr = AshAttributes()
        self._datasetname = 'HYSPLIT'

        # keeps track of status.
        # does file exists?
        self._fexists = self.read()
        # has the multiplication factor
        # been changd through changemult method.
        self._changed = False       
 

    @property
    def cxra(self):
        return self._cxra

    @property
    def attr(self):
        return self._attr

    def empty(self):
        if self._cxra.ndim == 0: return True
        else: return False

    def assign_attrs(self,inp):
        self.attr.update(inp)
        logger.debug('adding attributes {}'.format(self.attr))
        self._cxra = self._cxra.assign_attrs(self.attr.attr)

    def make_cdump_xra(self, blist, century=None, species=None, inp={}):
        self._cxra = hysplit.combine_dataset(
            blist, century=century, sample_time_stamp="start", species=species
        )
        self.assign_attrs(inp)

    def close(self):
        try:
            self._cxra.close()
        except Exception as eee:
            logger.warning(eee)
            print('could not close')
        self._cxra = xr.DataArray()

    def read(self):
        if os.path.isfile(self.fname):
            logger.info("netcdf file exists. Opening {}".format(self.fname))
            cxra = xr.open_dataset(self.fname)
            logger.info('data variable {}'.format(list(cxra.data_vars.keys())[0]))
            cxra = cxra[list(cxra.data_vars.keys())[0]]
            self._cxra = cxra
            self._fexists=True
        else:
            self._fexists=False
        return  self._fexists

    def changemult(self,mult):
        if 'mult' not in self._cxra.attrs.keys():
           logger.warning('mult not in attributes')
           self._changed=False
           return self._changed
        pmult = self._cxra.attrs["mult"]
        if np.isclose(pmult, mult, atol=0, rtol=0.1):
           self._changed = False
        else:
           self._cxra = mult * self._cxra / pmult
           self.assign_attrs({'mult':mult})
           self._changed = True
        return self._changed

    def remove(self):
        if os.path.isfile(self.fname):
            Helper.remove(self.fname)
        self._fexists=False

    def write_with_compression(self,overwrite=False):
        if self._fexists:
           # remove the old file if it has been changed or
           # overwrite is set.
           if self._changed or overwrite: 
              self.remove()
        if not self._fexists:
            cxra2 = self._cxra.to_dataset(name=self._datasetname)
            ehash = {"zlib": True, "complevel": 9}
            vlist = [x for x in cxra2.data_vars]
            vhash = {}
            for vvv in vlist:
                vhash[vvv] = ehash
            cxra2.to_netcdf(self.fname, encoding=vhash)


# create a source term ensemble from TCM.
import datetime
import numpy as np
import xarray as xr
from utilvolc import volcMER
from utilhysplit.evaluation import ensemble_tools
from utilhysplit.evaluation import vaa_atl_montage

# take a tcm matrix and create a source term ensemble.



def sval2info(sval,year):
    temp = sval.split('_')
    ht = int(temp[-1])
    dt = temp[0]
    month = int(dt[0:2])
    day = int(dt[2:4])
    hour = int(dt[4:])
    date = datetime.datetime(year,month,day,hour)
    return date, ht


class SourceEns:

    def __init__(self, tcmra):
        self.sources = tcmra.source.values
        self.tcmra = tcmra       
        self.forecast_list = []
        self.evector_list = []
        self._attrs = {}
        self.attrs = tcmra.attrs

    @property
    def attrs(self):
        return self._attrs

    @attrs.setter
    def attrs(self, atthash):
        if isinstance(atthash,dict):
           self._attrs.update(atthash) 
 
    @property
    def vent_height(self):
        return self._ventht

    @vent_height.setter
    def vent_height(self, vht, unit='m'):
        # unit should be meters.
        if unit.lower() == 'ft':
           vht = vht * 0.3048 
        self._ventht = vht

    @property
    def start_range(self):
        return self.earliest_start, self.latest_start

    @start_range.setter
    def start_range(self,erange):
        self.earliest_start = erange[0]
        self.latest_start = erange[1]
        
    @property
    def end_range(self):
        return self.earliest_end, self.latest_end

    @end_range.setter
    def end_range(self,erange):
        self.earliest_end = erange[0]
        self.latest_end = erange[1]

    @property
    def top_range(self):
        return self.lowest_top, self.highest_top

    @top_range.setter
    def top_range(self,trange):
        self.lowest_top = trange[0]
        self.highest_top = trange[1]

    def return_samples(self, nnn):
        hours = self.latest_start - self.earliest_start
        hours = int(hours.seconds/3600.0)
        if hours > 0:
           starts = np.random.randint(0,hours,nnn,dtype=int)
        startlist = []
        for sss in starts:
            startlist.append(self.earliest_start + datetime.timedelta(hours=int(sss)))

        hours = self.latest_end - self.earliest_end
        hours = int(hours.seconds/3600.0)
        endlist = []
        if hours > 0:
           ends = np.random.randint(0,hours,nnn,dtype=int)
        for eee in ends:
            endlist.append(self.earliest_end + datetime.timedelta(hours=int(eee)))
       
        tops = np.random.randint(self.lowest_top, self.highest_top,nnn,dtype=int)

        return list(zip(startlist, endlist, tops))

    def get_evector(self,start,end,top,year=2023):
        evector = []
        totmass = self.get_totmass(top)
        bottom = top - 5000
        nlev = 0
        if bottom < self.vent_height: bottom = self.vent_height
        for source in self.tcmra.source.values:
            date, ht = sval2info(source,year) 
            if date >= start and date <= end:
               if ht >= bottom and ht <= top:
                   evector.append(1.0) 
                   nlev += 1 
               else:
                   evector.append(0.0)
            else:
               evector.append(0.0)
        umass = totmass / nlev
        evector = np.array(evector) * umass
        emission_vector = xr.DataArray(evector,coords={'source':self.tcmra.source.values})
        return emission_vector


    def get_totmass(self,topheight):
        totmass = volcMER.HT2unit((topheight-self.vent_height)/1000.0,M63=0.1,verbose=False)
        # output is en g/h
        totmass = totmass * 1e3  # convert to amount of mg in 1 h 
        return totmass


    def generate_ensemble(self,nnn=10):
        # dot product in xarray is slow.
        # this is known issue.
        # see https://stackoverflow.com/questions/47180126/xarray-too-slow-for-performance-critical-code
        samples = self.return_samples(nnn)
        for sss in samples:
            print('working on', sss)
            emission_vector = self.get_evector(sss[0],sss[1],sss[2])
            self.evector_list.append(emission_vector)
           # def udot(tcm,ev):
           #     func = lambda x,y: np.dot(x,y)
           #     return xr.apply_ufunc(func,tcm,ev)
            #print('dot product')
            #self.forecast_list.append(udot(self.tcmra.isel(ens=0), emission_vector))
                

            #print('dot product')
            self.forecast_list.append(self.tcmra.dot(emission_vector))

    def plot(self,ilist=None,cmap='viridis',tlist=None):
        temp = self.get_ensemble(ilist)
        if tlist is not None:
           temp = temp.isel(time=tlist) 
        probvm = vaa_atl_montage.VAAMontageATL(temp,cmap=cmap)
        probvm.vaathresh=0.01
        #fig = probvm.plotpage()
        return probvm

    def write_ensemble(self,name):
        temp = get_ensemble(ilist=None)
        temp.to_netcdf(name) 

    def get_ensemble(self,ilist=None):
        # not working for unknown reason.
        # preprocess is messing up the coordinates.
        # seems to work ok on the jupyter notebook but not here.
        if ilist is None:
           fcl = self.forecast_list
        else:
           fcl = self.forecast_list[ilist]
        new = xr.concat(self.forecast_list,dim='source')
        slist = new.source.values
        new = new.assign_coords(source=('source',slist))
        #print(new.source.values)
        #print(new)
        temp,dim= ensemble_tools.preprocess(new)
        #print('---------------------')
        print(type(temp))
        temp = temp.rename({"ens":"metens"})
        temp = temp.rename({"source":"ens"})
        temp = temp.assign_attrs(self.attrs)
        return temp
    


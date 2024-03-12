# create a source term ensemble from TCM.
import datetime
import numpy as np
import xarray as xr
from utilvolc import volcMER

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
        sources = tcmra.source.values
        self.tcmra = tcmra       
 
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
        samples = self.return_samples(nnn)
        self.forecast_list = []
        self.evector_list = []
        for sss in samples:
            print('working on', sss)
            emission_vector = self.get_evector(sss[0],sss[1],sss[2])
            self.evector_list.append(emission_vector)
            print('dot product')
            self.forecast_list.append(self.tcmra.dot(emission_vector))


    


# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import numpy as np
import datetime
import pandas as pd
from pylab import matrix


"""
PGRMMR: Alice Crawford ORG: ARL/CICS
PYTHON 3
ABSTRACT: classes and functions for creating HYSPLIT control and setup files.

   CLASSES
   ModelBin class for parsing binary HYSPLIT CDUMP file


   CHANGES for PYTHON 3
     For python 3 the numpy char4 are read in as a numpy.bytes_ class and need to be converted to a python
     string by using decode('UTF-8').


"""


class ModelBin(object):
  """represents a binary cdump (concentration) output file from HYSPLIT
     methods:
     get_concentration - returns concentrations or mass loadings as either list or array.
     get_latlon - returns latitude longitude positions as either list or array
     define_struct - static method storing structure of cdump binary file in numpy dtypes.
     __init__
     _col_name - creates concentration column name describing pollutant and level
     _mcol_name- creates mass loading column name describing pollutant and level
     thicknesses - calculates self.depth - list of thicknesses corresponding to each top height in self.level
     _readfile - opens and reads contents of cdump file into pandas dataframe

     not working
     add_conc (method to add concentration grids together)
     xtract_point (method to extract concentration at lat-lon point
     ##It isn't clear if we want to keep the restriction that the sample stop must be less than drange[1]. 
     ##TO DO - be able to just get a summary of what is in the file - sample output times and so forth
     ##        without reading in the whole file?
     ##TO DO - method to write a cdump file. 
     ##TO DO - write add_conc method to add concentration grids together.
  """

  def __init__(self, filename, dir='', drange=[], missing=(), century=0, verbose=False, readwrite='r'):
     """
       drange should be a list of two datetime objects. 
        The read method will store data from the cdump file for which the sample start is greater thand drange[0] and less than drange[1] 
        for which the sample stop is less than drange[1]. 
        
     """
     ##It isn't clear if we want to keep the restriction that the sample stop must be less than drange[1]. 
     self.drange = drange
     self.filename = dir + filename
     self.missing = missing
     self.century = century
     self.verbose = verbose
     self.zeroconcdates = []        #list of tuples (date1, date2)  of averaging periods with zero concentrations
     self.nonzeroconcdates = []     #list of tuples  of averaging periods with nonzero concentrtations]
     #The following are appended to in the  __readfile__ method
     self.conc_names=[]  #names of the concentration columns
     self.species=[]     #list of pollutant identifier codes from cdump file.
     #The following are defined in the __readfile__ method
     #self.concframe   #pandas Data frame with lat, lon and concentration for each level and pollutant
     #self.sdate       #first sampling time start dates
     #self.edate       #last sampling time end dates
     #self.levels      #list of levels  (top height of level as stored in cdump file)
     if readwrite == 'r': 
         self.dataflag = self._readfile(dir+filename, drange, verbose, century)
         if self.dataflag:
             self.thicknesses() 
             #self.get_concentration()
             #self.get_latlon()
     
  def add_conc(self, concframe):
      ##TO DO - method to add concentration grids together.
      #self.concframe = concframe
      return -1 

  def _col_name(self, poll, lev):
      """concentration column name for panda dataframe describing pollutant and level"""
      if poll=='':
         poll = 'all'
      if lev=='':
         lev = 'all'
      col_name = 'conc' + '_' + str(poll) + '_' + str(lev)
      return col_name

  def _mcol_name(self, poll, lev):
      """mass loading column name for panda dataframe describing pollutant and level"""
      if poll=='':
         poll = 'all'
      if lev=='':
         lev = 'all'
      col_name = 'mass' + '_' + poll + '_' + str(lev)
      return col_name

  def thicknesses(self):
      """given list of levels as in HYSPLIT CONTROL FILE, return thickness of each level"""
      self.depth = {}
      levp = 0
      for lev in self.levels:
          self.depth[str(lev)] = (lev - levp)
          levp = lev
 

 

  def get_concentration(self,pdate=None, species=None, levels=None, mass_loading =1, units='', cthresh=0, grid=0, multx=1, verbose=False):
      """Returns concentration  for  pollutants specified in list of species and levels specified in list of levels 
         and for the date specified.
 
         This method creates a copy of self.concframe and operates on it. 
         This method does not alter self.concframe or any other class attributes.
  
         Inputs:
         If no date is specified then it will return concentrations for first time period.
         if species=None, then returns average concentration for all species. Should be input as list of species identifier codes.
         TO DO - do we want to be able to specify pollutants by number 1,2,3,4 rather than code?

         if levels=None then returns average concentration over all levels. Should be input as list of level top heights (as stored in cdump file)
         TO DO - do we want to be able to specify levels by index number 1,2,3,4 etc?

         if mass_loading = 1 returns column mass loading else returns concentration.
         if grid==0 then returns list of points at which mass loading  or concentration is non-zero. 
         if grid==1 then returns 2d array of mass loading or concentration over the entire grid.


         TO DO - Do we want to be able to return concentrations for multiple time periods or simply call this function multiple times?
         TO DO - do we want to keep track of total depth and top height for each point?


         TO DO - implement units which would allow unit transformation and unit meta data to be added.
         TO DO - implement cthresh which would implement a threshold on the concentration.
         """

      self.concframe.fillna(0, inplace=True)
      ##pandas.DataFrame.xs - returns a cross section (row(s) or column(s) from the dataframe)
      if pdate is None:                                    #if date is specified then
         concframe = self.concframe.xs(self.pdates[0][0]).copy()   #look at first date only.
         if verbose:
            print('Returning concentration for sampling time' , self.pdates[0][0] , ' to ' ,  self.pdates[0][1])
      else:
         try:
           concframe = self.concframe.xs(pdate).copy()    #takes slice of panda dataframe only for those dates
         except:
           lat2d, lon2d = self.get_latlon(grid=1)
           return np.zeros_like(lat2d) + -1               #if dates don't exist return array of zeros.
      #print concframe.info()
      if species is None:
         species = self.species
      if levels is None:
         levels = self.levels 
      lnames =[]

      ##Add species loop
      ##This loop adds concentrations from different species to get conc_all_lev columns
      for lev in levels:   #This loop adds concentrations from different species to get conc_all_lev columns
          cnames = []
          for sp in species:
              cnames.append(self._col_name(sp, lev))
          cnames2 = list(set(cnames) & set(self.conc_names))  #remove levels/ pollutants with no concentrations
          #creates new column which sums concentration for species listed in cnames2.
          concframe[self._col_name('', lev)] = concframe[cnames2].sum(axis=1)  
          lnames.append(self._col_name('', lev))
          if verbose: print('conc all for each species' , cnames2)
      ##END Add species loop

      
      ##Add concentration * depth loop
      ##This loop adds concentration * depth for each level      
      m_cols=[]
      thickness = 0
      for lev in levels:  
           #The deposition level has 0 depth so it will not contribute to the column mass loading.
           if lev != 0:
               cname = self._col_name('', lev)
               mname = self._mcol_name('',lev)
               m_cols.append(mname)
               concframe[mname] = concframe[cname] * self.depth[str(lev)] * multx #computes a mass loading for each level
               thickness += self.depth[str(lev)]
      ##END Add concentration * depth loop

      if verbose: print('columns to sum' , m_cols)
      concframe['mass_loading'] = concframe[m_cols].sum(axis=1)  #sums the mass loading of all levels
      if verbose:
         print('HERE')
         print(self.species)
         print(cnames)
         print(concframe.describe()) 

      if grid == 0:
         if mass_loading != 1:
            return list(concframe['mass_loading'])              #return mass loading list
         else:
            return list(concframe['mass_loading'] / thickness)  #return concentration list
      else:
         lat2d, lon2d = self.get_latlon(grid=1)
         conc2d = np.zeros_like(lat2d)
         indx = list(concframe['indx'])
         jndx = list(concframe['jndx'])
         ml = list(concframe['mass_loading'])
         k = 0
         for i in indx:
             j =jndx[k]
             conc2d[j-1][i-1]  = ml[k]   #need to subtract one because fortran arrays start at index 1 while python arrays start at index 0.
             k+=1
         if mass_loading != 1:
            conc2d = conc2d / thickness
         return  conc2d                  #returns 2d array of mass loading or concentration.



  def xtract_point(self, latpt, lonpt, drange=[]):
      #look at latlon grid. Find concentration at only a location
      # Find i,j, indx corresponding to latlon location
      ##not functional yet.
      lat, lon = get_latlon(grid=0)
      lat = np.array(lat)
      lon = np.array(lon)
      latdiff = np.abs(lat-latpt)
      londiff = np.abs(lon-lonpt)
      iii = np.where(latdiff == np.min(latdiff))
      jjj = np.where(londiff == np.min(londiff))
      return -1 


      
  def get_latlon(self, grid=0):
      """ returns latitude and longitude. if grid=1 then returns 2d arrays of latittude and longitude. Otherwise returns a list.""" 
      lat = np.arange(self.llcrnr_lat, self.llcrnr_lat+ self.nlat * self.dlat, self.dlat)
      lon = np.arange(self.llcrnr_lon, self.llcrnr_lon+ self.nlon * self.dlon, self.dlon)
      #if 'lat' not in self.concframe.columns:
      #    flat = lambda x : lat[x-1]
      #    flon = lambda x : lon[x-1]
      #    self.concframe['lat'] = self.concframe['jndx'].apply(flat)
      #    self.concframe['lon'] = self.concframe['indx'].apply(flon)
      if grid ==0:
           return list(self.concframe['lat']) , list(self.concframe['lon'])
      else:
           temp1 = np.zeros((lat.shape[0], lon.shape[0]))
           temp2 = np.zeros((lat.shape[0], lon.shape[0]))
           
           temp1[:,] = lon
           lon2d = temp1

           latm = matrix(lat)
           temp2[:,] = latm.T
           lat2d = temp2
           return lat2d, lon2d


  @staticmethod
  def define_struct():
     """Each record in the fortran binary begins and ends with 4 bytes which specify the length of the record.
     These bytes are called pad below. They are not used here, but are thrown out. 
     The following block defines a numpy dtype object for each record in the binary file. """
     real4 ='>f'
     int4  ='>i'
     int2  ='>i2'
     char4 ='>a4'

     rec1 = np.dtype([('pad1'      , int4),
                       ('model_id'  , char4),  #meteorological model id
                       ('met_year'  , int4),   #meteorological model starting time
                       ('met_month' , int4),
                       ('met_day'   , int4),
                       ('met_hr'    , int4),
                       ('met_fhr'   , int4),  #forecast hour
                       ('start_loc' , int4),  #number of starting locations
                       ('conc_pack' , int4),  #concentration packing flag (0=no, 1=yes)
                       ('pad2'      , int4),
                      ])

     #start_loc in rec1 tell how many rec there are.
     rec2 = np.dtype([('pad1'    ,int4),
                       ('r_year'  ,int4),     #release starting time
                       ('r_month' ,int4),
                       ('r_day'   ,int4),
                       ('r_hr'    ,int4),
                       ('s_lat'   ,real4),    #Release location
                       ('s_lon'   ,real4),
                       ('s_ht'    ,real4),
                       ('r_min'   ,int4),     #release startime time (minutes)
                       ('pad2'    ,int4),
                      ])

     rec3 = np.dtype([('pad1'    ,int4),
                      ('nlat'    ,int4),
                      ('nlon'    ,int4),
                      ('dlat'    ,real4),
                      ('dlon'    ,real4),
                      ('llcrnr_lat'    ,real4),
                      ('llcrnr_lon'    ,real4),
                       ('pad2'    ,int4),
                      ])

     rec4a = np.dtype([('pad1'    ,int4),
                      ('nlev'    ,int4),     #number of vertical levels in concentration grid
                      ])

     rec4b = np.dtype([('levht'   ,int4),    #height of each level (meters above ground)
                      ])

 
     rec5a = np.dtype([('pad1'    ,int4),
                       ('pad2'    ,int4),
                       ('pollnum'    ,int4),  #number of different pollutants
                      ])

     rec5b = np.dtype([('pname'   ,char4),    #identification string for each pollutant
                      ])
   
     rec5c = np.dtype([('pad2'   ,int4),
                      ])

     rec6 = np.dtype([('pad1'    ,int4),
                      ('oyear'   ,int4),      #sample start time.
                      ('omonth'  ,int4),
                      ('oday'    ,int4),
                      ('ohr'     ,int4),
                      ('omin'    ,int4),
                      ('oforecast'  ,int4),
                      ('pad3'    ,int4),
                      ])

     #rec7 has same form as rec6.            #sample stop time.

    #record 8 is pollutant type identification string, output level.

     rec8a = np.dtype([('pad1'  ,int4),  
                      ('poll'   ,char4),  
                      ('lev'    ,int4),
                      ('ne'     ,int4), #number of elements
                      ])
    
     rec8b = np.dtype([('indx'   ,int2),
                       ('jndx'   ,int2),
                       ('conc'   ,real4),
                      ])
     
     rec8c = np.dtype([('pad2'   ,int4),
                      ])
 
     return rec1, rec2, rec3, rec4a , rec4b, rec5a, rec5b, rec5c, rec6, rec8a, rec8b, rec8c

  def _readfile(self, filename, drange, verbose, century): 
     """Data from the file is stored in a pandas dataframe.
        returns False if all concentrations are zero else returns True.
        INPUTS
        filename - name of cdump file to open
        drange - [date1, date2] - range of dates to load data for. if [] then loads all data.
                 date1 and date2  should be datetime ojbects.
        verbose - turns on print statements
        century - if 0 or none will try to guess the century by looking at the last two digits of the year.
        For python 3 the numpy char4 are read in as a numpy.bytes_ class and need to be converted to a python
        string by using decode('UTF-8').

     """
        ##8/16/2016 moved species=[]  to before while loop. Added print statements when verbose.

     self.pdates=[]  #list of tuples giving the (sample start date, sample end date)
     fp = open(filename, 'rb') 

     ##each record in the fortran binary begins and ends with 4 bytes which specify the length of the record.
     ##These bytes are called pad1 and pad2 below. They are not used here, but are thrown out. 
     ##The following block defines a numpy dtype object for each record in the binary file.
     rec1 , rec2, rec3, rec4a, rec4b, rec5a, rec5b, rec5c, rec6, rec8a, rec8b, rec8c = self.define_struct()
     rec7 = rec6

     #start_loc in rec1 tell how many rec there are.
     tempzeroconcdates =[]   
     #Reads header data. This consists of records 1-5.
     hdata1=np.fromfile(fp,dtype=rec1, count=1)
     nstartloc = hdata1['start_loc'][0]
     ##in python 3 np.fromfile reads the record into a list even if it is just one number.
     ##so if the length this record is greater than one something is wrong.
     if len(hdata1['start_loc']) != 1:
        print('WARNING in ModelBin _readfile - number of starting locations incorrect')
        print(hdata1['start_loc'])
     hdata2=np.fromfile(fp,dtype=rec2, count=nstartloc)
     hdata3=np.fromfile(fp,dtype=rec3, count=1)
     self.nlat = hdata3['nlat'][0]
     self.nlon = hdata3['nlon'][0]
     self.dlat = hdata3['dlat'][0]
     self.dlon = hdata3['dlon'][0]
     self.llcrnr_lon = hdata3['llcrnr_lon'][0]
     self.llcrnr_lat = hdata3['llcrnr_lat'][0]
     self.metmodel = hdata1['model_id'][0].decode('UTF-8')
     self.sourcedate=[]
     self.slat=[]
     self.slon=[]
     self.sht=[]
     
     ##Loop through starting locations
     for n in  range(0,nstartloc):
        #create list of starting latitudes, longitudes and heights.
        self.slat.append(hdata2['s_lat'][n])
        self.slon.append(hdata2['s_lon'][n])
        self.sht.append(hdata2['s_ht'][n])

        #try to guess century if century not given
        if century == 0:
           if hdata2['r_year'][0] < 50:
              century = 2000
           else:
              century= 1900
        ##add sourcedate which is datetime.datetime object
        self.sourcedate.append(datetime.datetime(century+hdata2['r_year'][n] , \
                            hdata2['r_month'][n] , hdata2['r_day'][n] , \
                            hdata2['r_hr'][n], hdata2['r_min'][n]))
     
     #read record 4 which gives information about vertical levels.
     hdata4a=np.fromfile(fp,dtype=rec4a, count=1)                  #gets nmber of levels
     hdata4b=np.fromfile(fp,dtype=rec4b, count=hdata4a['nlev'][0]) #reads levels, count is number of levels.
     self.levels = hdata4b['levht']

     #read record 5 which gives information about pollutants / species.
     hdata5a=np.fromfile(fp,dtype=rec5a, count=1)
     hdata5b=np.fromfile(fp,dtype=rec5b, count=hdata5a['pollnum'][0])
     hdata5c=np.fromfile(fp,dtype=rec5c, count=1)

     if verbose:
         print('REC1 ***************************************')
         print('pad, MetId, Met starting time (year, month, day, hour, forecast-hour), number starting loc., packing flag')
         print(hdata1)
         print('REC 2 **************************************')
         print('Release start time (year, month, day, hour), start location (lat, lon, height), release start time (minutes)')
         print(hdata2) 
         print('REC 3 **************************************')
         print(hdata3) 
         print('REC 4 **************************************')
         print(hdata4a)
         print(hdata4b)
         print('nlev', hdata4a['nlev'])
         print('height of levels' , hdata4b['levht'])
         print('REC 5 **************************************')
         print(hdata5a)
         print(hdata5b)
         print(hdata5c)
         print('number of pollutants', hdata5a['pollnum'])

     #Loop to reads records 6-8. Number of loops is equal to number of output times.
     #Only save data for output times within drange. if drange=[] then save all.
     ##Loop to go through each sampling time
     ii=0
     iii=0
     imax = 1e3 #Safety valve - will not allow more than 1000 loops to be executed.
     testf=True
     while testf:
        hdata6=np.fromfile(fp,dtype=rec6, count=1)
        hdata7=np.fromfile(fp,dtype=rec6, count=1)
        if len(hdata6) == 0:       #if no data read then break out of the while loop.
           break
        if verbose:
            print('REC 6 & 7 ***************')
            print(hdata6)
            print(hdata7)
        #pdate1 is the sample start
        #pdate2 is the sample stop
        pdate1 = datetime.datetime(century+hdata6['oyear'], hdata6['omonth'], hdata6['oday'], hdata6['ohr'])
        pdate2 = datetime.datetime(century+hdata7['oyear'], hdata7['omonth'], hdata7['oday'], hdata7['ohr'])
        savedata=True
        pdatea = pd.Timestamp(pdate1)
        pdateb = pd.Timestamp(pdate2)
        tdelta = pdateb -pdateb

        #if pdate1 is within drange then save the data.
        #AND if pdate2 is within drange then save the data.
        #if drange[0] > pdate1 then stop looping to look for more data
        ##this block sets savedata to true if data within specified time range or time range not specified
        if drange ==[]:
           savedata=True
           ii=0
        elif pdate1 >= drange[0] and pdate1 <= drange[1] and pdate2 <= drange[1]:
           savedata=True
           ii=0
        elif pdate1 > drange[1] or pdate2 > drange[1]:
           testf=False
           savedata=False
        else:
           savedata=False
        ##END block

          
        if verbose: 
           print(savedata , 'DATES :' , pdate1 , pdate2) 
        if savedata:
           self.pdates.append((pdate1,pdate2))  #add sample start and sample stop time to pdates list.

        datelist = [] 
        inc_iii=False
 
        ##LOOP to go through each level
        for lev in range(hdata4a['nlev'][0]):
            ##LOOP to go through each pollutant
            for pollutant in range(hdata5a['pollnum'][0]):
                ##record 8a has the number of elements (ne). If number of elements greater than 0 than there are concentrations.
                hdata8a=np.fromfile(fp,dtype=rec8a, count=1)
                
                if hdata8a['ne'] >= 1:  #if number of elements is nonzero then
                   hdata8b = np.fromfile(fp,dtype=rec8b, count=hdata8a['ne'][0]) #get rec8 - indx and jndx
                   self.nonzeroconcdates.append(pdate1)                          #add sample start time to list of start times with non zero conc
                else:
                   tempzeroconcdates.append(pdate1)                              #or add sample start time to list of start times with zero conc.
                hdata8c=np.fromfile(fp,dtype=rec8c, count=1)

                ##if savedata is set and nonzero concentrations then save the data in a pandas dataframe
                if savedata and hdata8a['ne'] >=1: 
                   #self.nonzeroconcdates.append(pdate1)
                   inc_iii=True  #set to True to indicate that there is data to be saved.
                   #create column name for data
                   col_name = self._col_name(hdata8a['poll'][0].decode('UTF-8'), hdata8a['lev'][0]) 
                   if col_name not in self.conc_names:
                       self.conc_names.append(col_name)
                       if verbose: print('appending column name' , col_name)
                   if hdata8a['poll'][0].decode('UTF-8') not in self.species:
                       if verbose: print('appending species' , hdata8a['poll'][0])
                       self.species.append(hdata8a['poll'][0].decode('UTF-8'))
                          
                   ndata = hdata8b.byteswap().newbyteorder()  #otherwise get endian error.
                   if ii==0:  #if this is the first time through. create dataframe for first level and pollutant.
                      concframe = pd.DataFrame.from_records(ndata)
                      concframe.rename(columns={'conc':col_name}, inplace=True)
                   else:     #create dataframe for level and pollutant and then merge with main dataframe.
                      concframe2 = pd.DataFrame.from_records(ndata)
                      concframe2.rename(columns={'conc':col_name}, inplace=True)
                      concframe =  pd.merge(concframe, concframe2, how='outer', on=['jndx','indx'])
                   ii+=1
                if verbose:
                    print('REC 8 **************************************')
                    print(hdata8a)
            ##END LOOP to go through each pollutant
        ##END LOOP to go through each level
        if ii > imax:  #safety check - will stop sampling time while loop if goes over imax iterations.
           testf=False
        if inc_iii:
            if iii==0:  
                #if first time through sampling time loop then create final dataframe.
                concframe3 = concframe
                concframe3['sdate'] = pdatea #create column sdate which is sampling start time
                concframe3['edate'] = pdateb #create column edate which is sampling stop time
            else:
                ##add the next sampling time onto the final dataframe.
                #print '******concat*************' , pdatea , '  time' , str(iii)
                concframe['sdate'] = pdatea
                concframe['edate'] = pdateb
                concframe3 = pd.concat([concframe3, concframe])
            if inc_iii: 
                iii+=1

     ##END OF Loop to go through each sampling time


     if iii==0:
        print('Warning: ModelBin class _readfile method: no data in the date range found')
        return False
        #If concframe3 does not exist then file had no concentrations in date range specified
     concframe3['idx'] = list(range(0,concframe3.shape[0]))
     concframe3.set_index(['sdate', 'idx'],inplace=True)
     self.concframe = concframe3
     self.sdate = self.pdates[0][0] 
     self.edate = self.pdates[-1][1]

     ##add latitude longitude columns
     lat = np.arange(self.llcrnr_lat, self.llcrnr_lat+ self.nlat * self.dlat, self.dlat)
     lon = np.arange(self.llcrnr_lon, self.llcrnr_lon+ self.nlon * self.dlon, self.dlon)
     flat = lambda x : lat[x-1]
     flon = lambda x : lon[x-1]
     self.concframe['lat'] = self.concframe['jndx'].apply(flat)
     self.concframe['lon'] = self.concframe['indx'].apply(flon)

     for dt in tempzeroconcdates:
         if dt not in self.nonzeroconcdates:
            self.zeroconcdates.append(dt)

     return True 



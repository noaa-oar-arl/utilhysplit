# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
#import sys 
import numpy as np
import datetime


"""class EmitTimes 
   class EmitCycle
   class EmitLine
"""

##Cannot handle multiple species.


class EmitTimes(object):

   def __init__(self, filename='EMITIMES.txt'):
       self.filename=filename
       self.cycle_list = []  #list of EmitCycle objects.
       self.ncycles = 0
       self.chash = {}
       self.sdatelist =[] 

   def header_str(self):
      returnval =  'YYYY MM DD HH    DURATION(hhhh) #RECORDS \n'
      returnval += 'YYYY MM DD HH MM DURATION(hhmm) LAT LON HGT(m) RATE(/h) AREA(m2) HEAT(w)  \n'
      return returnval 


   def findmaxrec(self):
       maxrec = 0
       for ec in self.cycle_list:
           if ec.nrecs > maxrec:
              maxrec = ec.nrecs
       return maxrec


   def write_new(self, filename):
       maxrec = self.findmaxrec()
       with open(filename, 'a') as fid:
           fid.write(self.header_str())
       for ecycle in self.cycle_list:
           for iii in range(0, maxrec - ecycle.nrecs):
              ecycle.add_dummy_record()
           ecycle.write_new(filename)

   def read_file(self):
       with open(self.filename, 'r') as fid:
            fid.readline()
            fid.readline()
            done=False
            while not done:
               ec = EmitCycle()
               check = ec.read_cycle(fid)
               if not check: 
                  done=True
               else:
                  self.cycle_list.append(ec)
                  self.ncycles += 1

   def add_cycle(self, sdate, duration):
        self.ncycles +=1 
        ec = EmitCycle(sdate, duration)
        self.cycle_list.append(ec)
        d1 = sdate
        dt = datetime.timedelta(hours=int(duration))        
        d2 = sdate + dt
        self.chash[self.ncycles-1] = (d1, d2)
        self.sdatelist.append(sdate)
 
   def filter_records(self,llcrnr, urcrnr):
       for ec in self.cycle_list:
           ec.filter_records(llcrnr, urcrnr)

 
   def add_record(self, date, duration, lat, lon,
                   height, rate, area, heat):
        """
        adds a record to a cycle based on the date of the record.
        """
        ###This block determines which cycle the record goes into
        ###based on the date.
        cycle_number = -1
        for ccc in self.chash.keys():
            if date >= self.chash[ccc][0] and date < self.chash[ccc][1]:
               cycle_number = ccc
        if cycle_number == -1:
            rvalue = False
        else:
            self.cycle_list[cycle_number].add_record(date, duration, lat, lon,
                                                 height, rate, area, heat)
            rvalue = True 
        return rvalue 

class EmitCycle(object):
   """Helper class for EmitTimes
   This represents a cycle in an EmitTimes file.
   Each cycle begins with a line which has the start date, duration
   and number of records. Then the records follow.
   """
   #def __init__(self, filename='EMITIMES.txt'):
   def __init__(self, sdate=None, duration=None):
       self.sdate = sdate
       self.duration = duration  #duration of the cycle.
       self.recordra=[]
       self.nrecs = 0
       ##all cycles in a file must have same number of records.
       ##so some cycles may need to have dummy records
       ##with zero emissions.
       self.dummy_recordra = []
       self.drecs = 0

   def parse_header(self, header):
       """
       read header in the file.
       """
       temp = header.split()
       year = int(temp[0])
       month = int(temp[1])
       day = int(temp[2])
       hour = int(temp[3])
       #minute = int(temp[4])
       dhour = int(temp[4])
       nrecs = int(temp[5])
       self.sdate = datetime.datetime(year, month, day, hour)
       self.duration = datetime.timedelta(hours=dhour) 
       self.nrecs = nrecs
       return nrecs

   def write_new(self, filename):
       """
       write new emittimes file.
       """
       maxrec = self.nrecs + self.drecs 
       datestr = self.sdate.strftime('%Y %m %d %H ')
       #print('FILENAME EMIT', filename)
       with open(filename, 'a') as fid:
            #fid.write(self.header_str())
            fid.write(datestr + ' ' + self.duration + ' ' + str(maxrec) + '\n')
            for record in self.recordra:
                fid.write(str(record)) 
            for record in self.dummy_recordra:
                fid.write(str(record)) 

   def parse_record(self, record):
       temp = record.split()
       year = int(temp[0])
       month = int(temp[1])
       day = int(temp[2])
       hour = int(temp[3])
       dhour = int(temp[5][0:2])
       dmin = int(temp[5][-2:])
       duration = temp[5]
       sdate = datetime.datetime(year, month, day, hour)
       lat = float(temp[6])
       lon = float(temp[7])
       ht = float(temp[8])
       rate = float(temp[9])
       try:
          area = float(temp[10])
       except:
          area=0
       try:
          heat = float(temp[11])
       except:
          heat=0
       return EmitLine(sdate, duration, lat, lon, ht, rate, area, heat) 

   def add_dummy_record(self):
       """uses last record in the recordra to get date and position"""
       rc = self.recordra[-1]
       eline = EmitLine(rc.date, "0100", rc.lat, rc.lon, 0, 0, 0, 0)
       self.dummy_recordra.append(eline)
       self.drecs +=1

   def add_record(self, sdate, duration, lat, lon, ht, rate, area, heat):
       """Inputs
       sdate
       duration
       lat
       lon
       height
       rate
       area
       heat
       """
       eline = EmitLine(sdate, duration, lat, lon, ht, rate, area, heat)
       self.recordra.append(eline)
       self.nrecs +=1

   def read_cycle(self, fid):
       check=True
       recordra=[]
       header = fid.readline()
       if not header: 
         check=False
       else:
           nrecs =  self.parse_header(header)
           for line in range(0,nrecs):
               temp = fid.readline()
               recordra.append(self.parse_record(temp))
           self.recordra.extend(recordra)
           self.nrecs += nrecs         
       return check


   def filter_records(self,llcrnr, urcrnr):
       """ removes records which are outside the box
           described by llcrnr = (lat, lon) lower left corner
                        urcrnr = (lat, lon) upper right corner
       """
       iii=0
       rrr=[]
       for record in self.recordra:
           if record.lat < llcrnr[1] or record.lat > urcrnr[1]:
              rrr.append(iii)
           elif record.lon< llcrnr[0] or record.lon > urcrnr[0]:
              rrr.append(iii)
           iii+=1
       for iii in sorted(rrr, reverse=True):
           self.recordra.pop(iii)
           self.nrecs -= 1


class EmitLine(object):
   """

   """
   def __init__(self, date, duration, lat, lon, height, rate, area=0, heat=0):
       self.date = date
       self.duration = duration
       self.lat = lat
       self.lon = lon
       self.height=height
       self.rate = rate
       self.area = area
       self.heat = heat

   def __str__(self):
       returnstr = self.date.strftime("%Y %m %d %H %M ")
       returnstr += self.duration + ' '
       returnstr += str(self.lat) + ' '
       returnstr += str(self.lon) + ' '
       returnstr += str(self.height) + ' '
       returnstr += '{:1.2e}'.format(self.rate) + ' '
       returnstr += '{:1.2e}'.format(self.area) + ' '
       returnstr += '{:1.2e}'.format(self.heat) + ' \n'
       return returnstr


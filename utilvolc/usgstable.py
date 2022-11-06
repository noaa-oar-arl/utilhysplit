#!/opt/Tools/anaconda3/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
#from math import *
import numpy as np
import string
import datetime

##CLASSEs
## USGStable - object which represents USGS table of volcanoes. Attribute is name and location of csv file.
##             method to obtain data for one volcano in the table.
################################################################################################
class USGStable(object):
      """represents usgs preliminary spreadshet of eruption source parameters for volcanoes of the world"""
      def __init__(self, dir="/pub/Scratch/alicec/Plumeria/" , fname="usgs_table.csv"):
          self.fname = dir + fname

      def volclist(self, firstletter, accents=1):
          """Input - name of volcano - must match name in usgs file. Not case senstive.
             returns information from the USGS file for the given volcano"""
          #fid = codecs.open(self.fname, encoding='utf-8', mode='r')
          fid = open(self.fname, 'r')
          firstletter = firstletter.lower().strip()
          i = 0
          vhash = {}
          while i== 0:
                try:
                    wline = fid.readline()
                except:
                    pass
                if wline == '':
                   i=1
                   vhash=-1
                   print("EOF reached.")
                else:
                   if '"'in wline:
                       temp = wline.split('"')
                       name = temp[1]
                       name = name.split(',')
                       name = name[0]
                       wline = temp[0] + name + temp[2]
                       #print wline 
                       #print name
                   wline = wline.split(',')
                   vname = self.remove_accents(wline[4].lower().strip(), accents)
                   vname = vname.replace(' ', '')
                   #vname = vname.replace('-', '')
                   #print vname
                   if firstletter == vname[0]:
                      print(vname)
          fid.close()

      def volcano(self, volcname, accents=0):
          """Input - name of volcano - must match name in usgs file. Not case senstive.
             returns information from the USGS file for the given volcano"""
          #fid = codecs.open(self.fname, encoding='utf-8', mode='r')
          fid = open(self.fname, 'r')
          volcname = volcname.lower().strip()
          i = 0
          vhash = {}
          while i== 0:
                wline = fid.readline()
                if wline == '':
                   i=1
                   vhash=-1
                   print("EOF reached." , volcname , "not found in USGS table." , self.fname)
                else:
                   if '"'in wline:
                       temp = wline.split('"')
                       name = temp[1]
                       name = name.split(',')
                       name = name[0]
                       wline = temp[0] + name + temp[2]
                   wline = wline.split(',')
                   vname = self.remove_accents(wline[4].lower().strip(), accents)
                   vname = vname.replace(' ', '')
                   if vname == volcname:
                      i = 1
                      vhash= self.parseline(wline, accents=accents)
          fid.close()
          return vhash

      def parseline(self, wline, accents):
              vhash={}
              vhash['number'] = wline[0].strip() 
              vhash['rn'] = wline[1].strip() 
              vhash['sn'] = wline[2].strip() 
              vhash['vn'] = wline[3].strip() 
              vhash['name'] = self.remove_accents(wline[4].strip(), accents)
              vhash['namea'] = wline[4].strip()
              vhash['location'] = self.remove_accents(wline[5].strip(), accents)
              vhash['locationa'] = wline[5].strip(), accents
              vhash['status'] = wline[6].strip() 
              vhash['latitude'] = float(wline[7])
              if wline[8].strip().upper() == "S":
                   vhash['latitude'] = -1 * vhash['latitude']
          
              vhash['vf'] = wline[9]

              vhash['longitude'] = float(wline[10])
              if wline[11].strip().upper() == "W":
                   vhash['longitude'] = -1 * vhash['longitude']

              vhash['elev'] = float(wline[12]) / 1000.0
              vhash['type'] = wline[13].strip()
              vhash['timeframe'] = wline[14].strip()
              vhash['eruption type'] = wline[15].strip()
              return vhash 


      def remove_accents(self, word, accents ):
          """The USGS table writes all volcano and placenames with appropriate accents.
             if accents==1 then this will remove accents. if accents !=1 then it will do nothing."""
          if accents==1:
             word = word.replace("\xe9", "e")
             word = word.replace("\xe8", "e")
             word = word.replace("\xeA", "e")
             word = word.replace("\xeB", "e")
             
             word = word.replace("\xF9", "u")
             word = word.replace("\xFA", "u")
             word = word.replace("\xFB", "u")
             word = word.replace("\xFC", "u")
             
             word = word.replace("\xEC", "i")
             word = word.replace("\xED", "i")
             word = word.replace("\xEE", "i")
             word = word.replace("\xEF", "i")
         
             word = word.replace("\xFD", "y")
             word = word.replace("\xFF", "y")


             word = word.replace("\xF2", "o")
             word = word.replace("\xF3", "o")
             word = word.replace("\xF4", "o")
             word = word.replace("\xF5", "o")
             word = word.replace("\xF6", "o")
             word = word.replace("\xF8", "o")
           
             word = word.replace("\xE1", "a")
             word = word.replace("\xE2", "a")
             word = word.replace("\xE3", "a")
             word = word.replace("\xE4", "a")
             word = word.replace("\xE5", "a")
             
             word = word.replace("\xE6", "ae")
             
             word = word.replace("\xE7", "c")
             
             word = word.replace("\xF1", "n")
          return word

      #def list(self, location='', 


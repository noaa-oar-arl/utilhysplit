#!/n-home/alicec/anaconda/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
from math import *
#import sys 
#from scipy.io import netcdf
#from pylab import *
import numpy as np
import pandas as pd
#import numpy.ma as ma
#import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import string
import datetime
#from matplotlib.path import Path
#import shapely.geometry as sgeo
from scipy import interpolate
#from read_arason import *
#import subprocess
#import sys, getopt
#from operator import itemgetter
#import codecs
import os.path
##mytools is a module which contains conversions and other useful routines.
#from mytools import *  
#import random
#from pymet.meteorology import *
#import mastin_eq



## meteo_profile - object which contains information from a profile.txt file created by the fortran program profile.
#                  method to write meteorological input for plumeria
#                  method to compute wind direction
#                  method getvar which returns an array with the value requested as well as the units.



def sphu2relh(sphu, temp, pres ):
    """converts specific humidity, sphu, to relative humidity.
       esat is saturated vapor pressure, pres is pressure, temp is temperature (Kelvin) """
    esat = sat_vap_density(temp)
    #Pc=220.640  #hPa
    #Tc=647.096  #kelvin
    #C1 = -7.859
    #C2 = 1.844
    #C3 = -11.7866
    #C4 = 2.6807
    #C5 = -15.9618
    #C6 = 1.80123
    #theta = 1 - temp/Tc
    #M = Tc/T * (C1*theta + C2 * theta**1.5 + C3 * theta **3  + C4*theta**3.5 + C5*theta**4+C6*theta**7.5)
    #Pws = Pc * np.exp(M)
    return sphu*(pres) / (0.622*esat)

def sat_vap_density(t, type=1):
    """calculates saturation vapour pressure used to convert sphu to relh"""
    esat = np.exp(21.4-5351.0/t)
    return esat

def relh2sphu(relh, temp, pres):
    esat = sat_vap_density(temp)
    evap = relh * esat
    sphu = 0.622 * evap / pres * 1000
    return sphu

                    
def callprofile(hdir, mdir, meteo, lat, lon, dt, mfname='', nstop=24, offset=0):
    import subprocess
    """ Calls the fortran program profile.
        mdir directory of meteorological file
        meteo name of meteorological file
        lat, lon latitude longitude 
        dt time step
        mfname - name of file to move profile.txt to """
    callstr = os.path.join(hdir + "profile")
    p1 = "-d" + mdir
    p2 = "-f" + meteo
    p3 = "-y" + str(lat).strip()
    p4 = "-x" + str(lon).strip()
    p5 = "-t" + str(dt)
    p6 = "-n" + str(nstop)
    p7 = "-o" + str(offset)
    print(p1,p2,p3,p4,p5,p6,p7)
    subprocess.call([callstr, p1, p2, p3, p4, p5, p6, p7])
    if mfname !='':
        subprocess.call(["mv", "profile.txt", mfname])
  


class MeteoProfile(object):
    """Information read from a profile.txt file created by the fortran program profile.
       There may be information about multiple times in the profile.txt file.
       fname - name of profile.txt file. Default is profile.txt
       dir - directory in which profile.txt can be found. Default is current directory.
       type is the type of meteorological file the profile was created from (gdas, wrf, narr...) 
       datesvalid gives range of dates to retrieve data for.  
    """

    def __init__(self, fname='profile.txt', pdir='', tpe='', datesvalid=[]):
        """reads info from a profile file."""
        #fname="/pub/Scratch/alicec/KASATOCHI/meteo/wrf.profile.txt"
        #fname="/pub/Scratch/alicec/KASATOCHI/meteo/narr.profile.txt"
        #fname="/pub/Scratch/alicec/KASATOCHI/meteo/gdas1.profile.txt"
        self.pdir = pdir
        self.fname=pdir + fname
        self.valra = {}
        self.valra2d = {}
        self.date_ra=[]


        self.uwind_ra=[]
        self.vwind_ra=[]
        self.wwind_ra=[]
        self.temp_ra =[]
        self.height_ra =[]
        self.tpot_ra =[]
        self.humidity_ra =[]
        self.magwind_ra=[]
        self.pressure_ra=[]

        self.time_ra=[]
        self.descrip_list=[]
       
        height_rat=[] 

        self.var3d = []
        self.var2d = [] 

        #count=0
        #time_index=0
        #firstday=[]
        #txtfile= open(self.fname, "r")
        #dvalid=0
        self.readnew(os.path.join(pdir,fname), datesvalid)
     
    def readnew(self, fname, datesvalid=[]):
        txtfile = open(fname, "r")
        cnt=0
        dvalid = True
        iii = 0
        init=True
        for line in txtfile:
            #print 'COUNT', cnt, line
            if '________' in line:
               cnt=0
            elif cnt==10 and dvalid:
               self.valra[time].append(line.strip())
            elif cnt == 3 and dvalid:
               self.valra2d[time].append(line.strip())
               cnt = 0
            elif cnt == 2:
               self.dunits2 = line
               cnt = 3
            elif cnt == 1:
               self.var2d = self.columnnames(line)
               cnt = 2
            elif cnt == 5:
               self.dunits3 = line
               cnt = 10
            elif cnt == 4:
               self.var3d = self.columnnames(line)
               cnt = 5
            elif "2D Fields" in line:
                cnt=1
            elif "3D Fields" in line:
                cnt=4
            elif "Profile Time" in line:
                time = self.gettime(line)     
                dvalid = self.testdate(time, datesvalid)
            #    print 'valid time', dvalid
                if dvalid:
                    self.valra[time] = []
                    self.valra2d[time] = []
                    self.date_ra.append(time)
            iii+=1
            if iii> 999990:
               break
        txtfile.close()
        #print self.date_ra
          
    def columnnames(self, line):
        return(line.split())

    def get_3dvar_df(self):
        """
        retunrs pandas dataframe with 3d variables in it.
        """
        varra = []
        dates = self.date_ra
        for tm in dates:
           for line in self.valra[tm]:
               temp2 = line.split()
               temp2 = list(map(float,temp2))
               # TO DO. 
               # skip lines which contain missing values.
               # need to implement fixed column width reader
               # cannot assume spaces are delimiters.
               if len(temp2) != len(self.var3d)+1: 
                  #print('warning: line contains missing values ', line)
                  continue
               temp = [tm] 
               temp.extend(temp2)
               #print(temp)
               try:
                   varra.append(temp)
               except:
                   pass
        cols = ['time','PRES1']
        cols.extend(self.var3d)
        df = pd.DataFrame(varra,columns=cols)
        return df

    def get_var(self, var, dates=[], hts=[]):
        varra = []
        if dates ==[]:
           dates = self.date_ra
        if var in  self.var3d:
           iii = self.var3d.index(var) + 1
           for tm in dates:
               for line in self.valra[tm]:
                   temp = line.split()
                   print(tm)
                   try:
                       varra.append(float(temp[iii]))
                   except:
                       pass
        elif var in self.var2d:
           iii = self.var2d.index(var) + 1
           for tm in dates:
               temp = self.valra2d[tm][0].split()
               varra.append(float(temp[iii]))
        else:
           return -1
        return varra     

    def testdate(self, time, datesvalid, verbose=False):
        if datesvalid == []:
           return True
        else:
           try:
             datesvalid[0]
           except:
             datesvalid = [datesvalid, datesvalid]
           try:
             datesvalid[1]
           except:
             datesvalid.append(datesvalid)
           if time >= datesvalid[0] and time <= datesvalid[1]:
              return True 
           else:
              return False
             

    def gettime(self,line, century=''):
                
            time=line.split()
            year = int(time[2])
            if century == '':
               #if no century provided then guess century.
               if year > datetime.datetime.now().year - 2000:
                  year += 1900
               else:
                  year += 2000
            month= time[3]
            day  = time[4]
            hour = time[5]
            hr = time[5]
            minute=0    #for input into datetime object. Minutes are always 00.
            testdate = datetime.datetime(int(year),int(month),int(day),int(hour),int(minute))
            return testdate


    #def __str__(self):
    #    """prints out in format suitable for entry into plume rise model"""
    #    ht = self.height_ra
    #    descrip="descrip"
    #    return descrip        

    def wind_direction(self, vwind, uwind):
            #vwind is magnitude of wind going from south to north
            #uwind is magnitude of wind going from West to East
            #print(len(vwind), len(uwind))
            uwind = np.array(uwind)
            vwind = np.array(vwind)
            wind_dir = np.arctan(vwind/uwind)*180/pi
            upos = np.where(uwind >= 0)
            wind_dir[upos]= 270 - wind_dir[upos] 
            uneg = np.where(uwind < 0)
            wind_dir[uneg] = 90 - wind_dir[uneg] 
            #print(wind_dir.shape , wind_dir[2])
            #print('U' , uwind.shape, uwind[2])
            #print('V' , vwind.shape, vwind[2])
            return wind_dir

class Radiosonde(MeteoProfile):
   """radiosonde file
   """

   def __init__(self, fname , datesvalid=[]):
       fid = open(fname, "r")
       datelist = []
       date = datesvalid[0] - datetime.timedelta(36000)

       self.magwind_ra =[]
       self.humidity_ra = []
       self.temp_ra = []
       self.height_ra = []
       self.pressure_ra=[]
       self.winddir_ra=[]
 
       self.date_ra = []
       self.time_ra = []

       ##not used. but must be compatible with profile files
       self.uwind_ra=[]
       self.vwind_ra=[]
       self.tpot_ra=[]
       magwind =[]
       humidity = []
       temp = []
       height = []
       pressure=[]
       winddir=[]

       for line in fid:
           wl = line.split()
           lintype = wl[0].strip()
           if lintype == '254':
              #This is an indentification line
              datestr = wl[4].strip() + '-' + wl[3] + '-' + wl[2] + ' ' + wl[1] + 'z'
              date = datetime.datetime.strptime(datestr , '%Y-%b-%d %Hz')
              if magwind != []:            #append last sets of data
                  self.magwind_ra.append(magwind)
                  self.humidity_ra.append(humidity)
                  self.temp_ra.append(temp)
                  self.height_ra.append(height)
                  self.pressure_ra.append(pressure)
                  self.winddir_ra.append(winddir)
                  magwind =[]
                  humidity = []
                  temp = []
                  height = []
                  pressure=[]
                  winddir=[]

              if date <= datesvalid[1] and date >=datesvalid[0]:
                      #print date
                  
                      self.date_ra.append(date)

           if date <= datesvalid[1] and date >=datesvalid[0]:
               if lintype in ['4','5','6','7','8','9']:       #dataline
                  #print(wl)
                  pressure.append(float(wl[1]))
                  height.append(float(wl[2]))
                  temp.append(float(wl[3]))
                  humidity.append(float(wl[4]))
                  winddir.append(float(wl[5]))
                  mw = float(wl[6])
                  if mw != 99999:
                      mw = mw / 10.0
                  magwind.append(mw)
                  #print 'height' , np.array(height).shape
               elif lintype == '1':       #station identification line
                  wban = wl[1].strip()
                  wmo  = wl[2].strip()
                  lat =  wl[3].strip()
                  lon =  wl[4].strip()
                  elev = wl[5].strip()
                  rtime = wl[6].strip()
               elif lintype == '2':       #sounding checks line
                  hyrdo = wl[1].strip()
                  maxwd = wl[2].strip()   #pressure of level having max wind.
                  tropl = wl[3].strip()   #pressure of level containing the tropopause. 
                  line = wl[4].strip()    #number of levels in the sounding.
                  tindex = wl[5].strip()  #7 if good estimation of tropopause. 11 if tropopause is 'suspected'
                  source = wl[6].strip()  #source.
               elif lintype == '3':       #station identifier
                  staid   = wl[1].strip() #station identifier
                  sonde   = wl[2].strip() #type of radiosonde code from TTBB. 
                  wsunits = wl[3].strip() #wind speed units. ms is tenths of meters per second.
                       
       fid.close()






#######################################################################################################



def crop_contour_plot(x,y,z, xr=[], yr=[]):
    """crops 2d arrays so that x data lies in range given by xr and y
       data lies in range given by yr. xr and yr need to be list/array with [xmin, xmax]"""
    print('YR', yr)
    print('XR', xr)
    if yr != []:
        vp = np.where(np.logical_and(y> yr[0] , y < yr[1]))
        min_yi = np.min(vp[1])
        max_yi = np.max(vp[1]) 
        print(min_yi , max_yi)
        x=x[0:,min_yi:max_yi]
        y=y[0:,min_yi:max_yi]
        z=z[0:,min_yi:max_yi]
    if xr != []:
        vp = np.where(np.logical_and(x> xr[0] , x < xr[1]))
        min_xi = np.min(vp[0])
        max_xi = np.max(vp[0]) 
        print(min_xi , max_xi)
        x=x[min_xi:max_xi,]
        y=y[min_xi:max_xi,]
        z=z[min_xi:max_xi,]
    return x , y , z 

def write_inputs(p, fname="input.txt", humidity="normal", datevalid=[]):
    namehash={"wrf":"/pub/Scratch/alicec/KASATOCHI/meteo/wrf.profile.txt",
            "narr":"/pub/Scratch/alicec/KASATOCHI/meteo/narr.profile.txt",
            "eyja":"/pub/Scratch/alicec/WRF/wrf.eyjap3.profile.txt",
            "gdas":"/pub/Scratch/alicec/KASATOCHI/meteo/gdas1.profile.txt",
            "ecmwf":"/pub/Scratch/alicec/KASATOCHI/meteo/ecmwf.profile.txt"}
    data=meteo_profile(namehash[p], datesvalid=[])
    data.write_plumeria_input(fname=fname , humidity=humidity, datevalid=datevalid)

def meteo_stats(fname, v, clevs=[], datesvalid=[], htrange=[]):
    """htrange should be a list of tuples. given a profile.txt file , a variable, valid dates and height ranges,
       returns statistics over the given height and date ranges for the variable."""
    meanlist = []
    maxlist= []
    minlist = []
    medianlist=[] 

    fig = plt.figure(fignum)
    ax = fig.add_subplot(splt)
    data=MeteoProfile(fname, datesvalid=datesvalid)
    #time=data.time_ra
    var, units = data.getvar(v)
    h2plot=np.array(data.height_ra) / 1000.0
    i=0
    for ht in h2plot:
        htr = htrange[i]
        xy = list(zip(ht, var[i])) 
        xy = list(filter(lambda x: x[0] < htr[1] and x[0] > htr[0]))
        xy = np.array(list(zip(*xy)))
        meanlist.append(np.mean(xy[1]))
        maxlist.append(np.amax(xy[1]))
        minlist.append(np.amin(xy[1]))
        medianlist.append(np.median(xy[1]))
        i+=1 
        
        
    shash['mean'] = meanlist
    shash['max'] = maxlist
    shash['min'] = minlist
    shash['median'] = medianlist 
    shahs['date'] =  data.date_ra
    return shash


def meteo_contours(fname, v, clevs=[], fignum=1, splt='111',  yr=[], datesvalid=[], rad = 0):
    """fname is the name of the meteo profile txt file to use
       v is the name of the variable to plot
    """
    fig = plt.figure(fignum)
    ax = fig.add_subplot(splt)
    if rad ==1:
       data=radiosonde(fname, datesvalid=datesvalid)
       var, units = data.getvar(v)
       i=0
       for ht in data.height_ra:
              data.height_ra[i] , var[i] = removefromset(data.height_ra[i], var[i], 99999)
              data.height_ra[i] , var[i] = removefromset(data.height_ra[i], var[i], 9999)
              data.height_ra[i] = np.array(data.height_ra[i]) / 1000
              i+=1
       h2plot = data.height_ra 
    else:
       data=meteo_profile(fname, datesvalid=datesvalid)
       var, units = data.getvar(v)
       h2plot=np.array(data.height_ra) / 1000.0
    #time=data.time_ra
    t2p=data.date_ra
    var, units = data.getvar(v)
    #print np.array(data.height_ra).shape
    #h2plot=np.array(data.height_ra) / 1000.0
    #h2plot=data.height_ra
    h2p=np.array(h2plot[0])
    print('HEIGHTS' , h2p)
    if yr != []:
       h2p = h2p[np.where(h2p<= yr[1])]
       h2p = h2p[np.where(h2p>= yr[0])]
    if h2p == []:
       print("WARNING: y range is invalid")
    var=np.array(var) 
    vshape=var.shape
    var2plot = []
    ##This creates output at standard heights.
    for i in range(0,vshape[0]):
        fail = 0
        var[i] = np.array(var[i])
        try:
           tspl = spline_v(h2plot[i], var[i], k=1)
        except:
           print('Interpolation failed for' , t2p[i])
           print(var[i])
           print('*******************************\n')
           print(h2plot[i])
           plt.figure(i)
           plt.plot(var[i], h2plot[i], '-b.')
        varint = interpolate.splev(h2p, tspl[0], der=0)
        var2plot.append(varint) 
    #tspl = spline_v(h2plot, var, k=1)
    #for t in tspl:
    #     varint = interpolate.splev(h2p, t , der=0)
    #     var2plot.append(varint)

    var2plot = np.array(var2plot)

    plt.figure(fignum)
    if clevs==[]:
       cs = plt.contourf(t2p,h2p,var2plot.T)
    else:
       cs = plt.contourf(t2p,h2p,var2plot.T, levels=clevs)
    cbar = plt.colorbar()
    cbar.set_label(units)
    set_date_axis(ax, t2p)
    #ax.set_xlabel('Hours  since ' +  data.startime)
    #ax.set_ylabel('Height (m)')
    #plt.title(data.meteofile + ' at grid point closest to ' + 'Lat ' + data.latitude + ' Lon ' + data.longitude )


def plot_contours(p,v, clevs=[], fignum=1, splt='111' , xr=[], yr=[],datesvalid=[]):
    """plots contour of  variable v as a function of height and time.
       """
    fig=plt.figure(fignum, figsize=(8,5))
    ax=fig.add_subplot(splt) 
    namehash={"wrf":"/pub/Scratch/alicec/KASATOCHI/meteo/wrf.profile.txt",
            "narr":"/pub/Scratch/alicec/KASATOCHI/meteo/narr.profile.txt",
            "ecmwf":"/pub/Scratch/alicec/KASATOCHI/meteo/ecmwf.profile.txt",
            "eyja":"/pub/Scratch/alicec/WRF/wrf.eyjap3.profile.txt",
            "gdas":"/pub/Scratch/alicec/KASATOCHI/meteo/gdas1.profile.txt",
            "popo":"/pub/Scratch/alicec/POPO/gdasprofile_dec14.txt"}
    typehash={"wrf":"wrf",
            "narr":"narr",
            "eyja":"wrf",
            "gdas":"gdas",
            "popo":"gdas",
            "ecmwf":"gdas"}
    data=MeteoProfile(namehash[p], datesvalid=[], type=typehash[p])
    time=data.time_ra
    var2plot, units = data.getvar(v)
    h2plot=np.array(data.height_ra)
    var2plot=np.array(var2plot) 
    vshape=var2plot.shape
    print(vshape , h2plot.shape , np.array(time).shape)
    ntile=vshape[1]
    t2p=np.tile(time,(ntile,1))
    t2p=t2p.T
    h2p=np.array(h2plot) 
    print(t2p.shape , h2p.shape, var2plot.shape)
    if xr!=[] or yr!=[]:
       t2p, h2p, var2plot = crop_contour_plot(t2p, h2p, var2plot, xr=xr, yr=yr)
    if clevs==[]:
       cs = plt.contourf(t2p,h2p,var2plot)
    else:
       cs = plt.contourf(t2p,h2p,var2plot, levels=clevs)
    cbar = plt.colorbar()
    cbar.set_label(units)
    ax.set_xlabel('Hours  since ' +  data.startime)
    ax.set_ylabel('Height (m)')
    plt.title(data.meteofile + ' at grid point closest to ' + 'Lat ' + data.latitude + ' Lon ' + data.longitude )
 
def plot_profiles(plist,tlist,vlist, clrs=['r','b','g','k','m'], syms=['.','o','*','3','+','x','d','1','.'], fignum=1, splt = '111'):
    """plots variables in vlist at times in tlist for profile files in plist"""
    namehash={"wrf":"/pub/Scratch/alicec/KASATOCHI/meteo/wrf.profile.txt",
            "narr":"/pub/Scratch/alicec/KASATOCHI/meteo/narr.profile.txt",
            "eyja":"/pub/Scratch/alicec/WRF/wrf.eyjap3.profile.txt",
            "ecmwf":"/pub/Scratch/alicec/KASATOCHI/meteo/ecmwf.profile.txt",
            "gdas":"/pub/Scratch/alicec/KASATOCHI/meteo/gdas1.profile.txt",
            "popo":"/pub/Scratch/alicec/POPO/gdasprofile_dec14.txt",
            "popo2":"/pub/Scratch/alicec/POPO/gp5.txt"}
    phash={}   #key=descrip of meteo file (wrf,narr or gdas), value is meteo_profile object.
    vhash={}   #key=descrip of meteo file (wrf,narr or gdas), value is array of variable of interest.
    print(vlist)
    if vlist==["all"]:
       vlist=["uwind","vwind","magwind","temp","tpot","humidity","wind_direction"]

    numplots=len(vlist)+1
    fig=[0]*numplots
    ax=[0]*numplots
    for p in plist:
        try:
          fname=namehash[p]
        except:
          fname = p
        phash[p]=MeteoProfile(fname)
    fignum_i=0
    for v in vlist:
        clr_i=0
        fig[fignum_i]=plt.figure(fignum)
        ax[fignum_i]=fig[fignum_i].add_subplot(splt)
        for p in plist:
            pclr=clrs[clr_i] 
            vhash[p], units= phash[p].getvar(v)
            #print vhash[p]
            sym_i=0
            print(phash[p].date_ra)
            for t in tlist:
                #tsym='-' + pclr + syms[sym_i]
                tsym=syms[sym_i]
                try:
                    tindex = phash[p].date_ra.index(t)
                except:
                    tindex=-1
                    print('TINDEX NOT FOUND' , phash[p].date_ra)
                    print(t,  tlist)
                if tindex !=-1:  
                    ltxt=p.replace('.txt','').upper() + ' '+ t.strftime('%m/%d %H:%MZ')
                    var2plot = vhash[p][tindex]  
                    ht2plot=phash[p].height_ra[tindex]
                    linetype='-' + pclr+tsym
                    ax[fignum_i].plot(var2plot, ht2plot, linetype, label=ltxt)      
                    #plt.plot(var2plot, ht2plot, pclr )      
                sym_i+=1   
            clr_i+=1
        ax[fignum_i].set_xlabel(units, fontsize=18)
        ax[fignum_i].set_ylabel('Height (m)', fontsize=18)
        print("SETTING AXIS LABELS")
        #plt.ylim(ymin=0)
        plt.grid(b=True, which='both', color='0.65',linestyle='-')
        handles , labels = ax[fignum_i].get_legend_handles_labels()
        plt.legend(handles, labels, loc='best', prop={'size':12})       
        plt.ylim((0,20000)) 
        fignum+=1
        fignum_i+=1


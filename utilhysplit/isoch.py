#!/opt/Tools/anaconda3/bin/python
import datetime
import os.path
import warnings
from optparse import OptionParser

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap

from pyhysplit.hysplit import HycsControl, ModelBin

mpl.use('Agg')

#import time
#import sys
#import re

#pylint:  disable=C0103
## TO DO - used warnings label to catch and suppress a warning when contourf was
##called. Was getting a basemap/__init__.py:???? VisibleDeprecationWarning:
##using a non-integer number instead of an integer will result in an error in the future.
##This might be a bug in basemap. One change in python3 is that an integer divided by an integer returns a float.

"""
Python3.5 (actually matplotlib >= v1.2)
Also needs pandas module and pyhysplit.hysplit
 -h for help on input options.
 Reads a cdump file.
 Reads a control file if available
       determines rate and duration from input arguments. If no input argument then
       will get rate and duration from first pollutant defined in the CONTROL file.
TO DO?
Should height of release be written on the plots?
handling multiple pollutants.

Do we need to change from basemap to cartopy as basemap no longer supported?

Written by Alain Wain to plot arrival times.
Modified by Alice Crawford.
            added option parser.
            Read cdump file directly rather than using con2asc to create ascii files from the cdump.
            Get information from CONTROL file.
"""

parser = OptionParser()

parser.add_option("-i", type="string", dest="fname", default='cdump',
                  help="Path and name of input cdump binary file")
parser.add_option("-o", type="string", dest="oname", default='toa',
                  help="name for output .png files")
parser.add_option("-t", type="string", dest="mytv", default='1e-19',
                  help="threshold concentration")
parser.add_option("-b", type="string", dest="bnds", default='',
                  help="latitutde and longitude bounds for plot. \
                  Format is latmin:latmax:lonmin:lonmax \
                  if none given, then the program will make an educated guess")
parser.add_option("--title", type="string", dest="title", default='NOAA HYSPLIT Simulation',
                  help="Main title to put on plots. Default is NOAA HYSPLIT Simulation")
parser.add_option("--loc", type="string", dest="location", default='Test Location',
                  help="Name of release location")
parser.add_option("-c", type="string", dest="controlfile", default='CONTROL',
                  help="name of control file to use for looking up release duration")
parser.add_option("--r", type="string", dest="myrate", default='',
                  help="Rate of release in mass / hour. \
                        Will try to read from CONTROL file if not given.")
parser.add_option("-d", type="string", dest="myrel", default='',
                  help="Duration of release in hours. \
                        Will try to read from CONTROL file if not given")
parser.add_option("--verbose", action="store_true", dest="verbose", default=False,
                  help="Verbose")

(options, args) = parser.parse_args()

fname = options.fname
scale = 1.77778
verbose = options.verbose
if os.path.isfile(fname):
    #Get source data from cdump file file. ModelBin is a class representing cdump file.
    hrun = ModelBin(fname, verbose=verbose)

    slat = hrun.slat
    slon = hrun.slon
    latlon = list(set(zip(slat, slon)))
    slon = latlon[0][0]  #if multiple, use the first longitude for putting on plot
    slat = latlon[0][1]  #if multiple, use the first latitude for putting on chart
    if len(latlon) != 1:
        print('Warning: more than one source location')
    sht = str(hrun.sht)  #currently the height of release is not used.
    mystart = hrun.sourcedate[0]                       #get date from cdump file.
    mystart = mystart.strftime('%d %b %Y %H%M UTC')    #convert datetime to a string

    ##if bounds are input then use those.
    obounds = False
    if options.bnds != '':
        berr = 0
        bounds = options.bnds.strip().split(':')
        try:
            blatmin = float(bounds[0])
        except:
            berr = 1
        try:
            blatmax = float(bounds[1])
        except:
            berr = 1
        try:
            blonmin = float(bounds[2])
        except:
            berr = 1
        try:
            blonmax = float(bounds[3])
        except:
            berr = 1
        if berr == 1:
            print('Warning: Input lat lon bounds does not parse. Format should be \
                   latmin:latmax:lonmin:lonmax. ', options.bnds)
            obounds = False
        else:
            obounds = True

    ##To get the name of location, rate of release and duration of release -
    ##first look at if they are input as as an option.
    myloc = options.location
    myrel = options.myrel
    myrate = options.myrate
    #This is the time that this program is run at.
    exec_time = datetime.datetime.now()
    mysrc = str(latlon[0][0]) + ' ' + str(latlon[0][1])
    isstime = exec_time.strftime('%d %b %Y  %H%M UTC')

    #get the 4 letter species / pollutant code(s) from cdump file.
    myspec = ''
    for msp in hrun.species:
        myspec = str(msp) +  ' '

    #get the concentration threshold from a system input argument.
    mytv = options.mytv
    thold = float(mytv)

    ##Meteorological model code from cdump file
    mymodel = hrun.metmodel


    ##determine output interval by subtracting first output time from second output time
    dts = list(set(hrun.concframe.index.get_level_values(0)))
    dts.sort()
    tint = dts[1] - dts[0]
    tint = tint.seconds / 3600.0
    if verbose:
        print('---------------------------------')
        print('output interval ', tint, ' hrs')
        print('sample output times: ')
        for idts in dts:
            print(idts)
        print('---------------------------------')

    ##Get lat lon grids from cdump file.
    lats, lons = hrun.get_latlon(grid=1)

    ##What information do we need from the control file?
    controlfile = options.controlfile
    if os.path.isfile(controlfile):
        controlfile = HycsControl(fname=controlfile)
        controlfile.read(verbose=verbose)
        ##If myrate and myrel were not input from the command line
        ##then get them from the control file.
        ##TO DO - if there are multiple species this will use the rate
        ## and duration only for the first species.
        if myrate == '':
            myrate = str(controlfile.species[0].rate)
        if myrel == '':
            myrel = str(controlfile.species[0].duration)
    else:
        print('Warning: control file ', controlfile, ' does not exist.')

    mylevs = ''
    if len(hrun.levels) > 1 and hrun.levels[0] != 0:
        mylevs = str(hrun.levels[0]) + ' to ' + str(hrun.levels[-1]) + ' m'
    else:
        mylevs = 'surface to ' + str(hrun.levels[-1]) + ' m'

    ## Use a list of arrays. 
    conc = []
    for dt in hrun.pdates:
        conc.append(hrun.get_concentration(pdate=dt[0], grid=1, mass_loading=0))

    ##TO DO: This part could be generalized to accomodate
    ##more or less days or different time periods.
    numdays = 3

    day = []
    #create an array the same size as the latitude array with all 0's for each day. 
    for ii in range(0, numdays):
        day.append(np.zeros_like(lats))
    dt = 0
    latmin = 90
    latmax = -90
    lonmin = 180
    lonmax = -180
    #loop through each concentration array.
    for concarray in conc:
        dt += tint  #concentration arrays are separated in time by tint. dt keeps track of time of concentration array being analyzed in loop.
        if dt <= 24:  #Add to the day1 array
            ##np.where returns index of points where the day[0] array is 0 (no pollutant has arrived previously.
            ##and where concarray is now greater than threshold (pollutant has arrived at this time, for the first time.
            vpi = np.where(np.logical_and(day[0] == 0, concarray >= thold))
            ##set the points with index vpi equal to the arrival time dt.
            day[0][vpi] = dt
            ##This part finds min and max lat's and lons for determining map boundaries.
            temp = np.min(lats[vpi])
            if temp < latmin:
                latmin = temp
            maxtemp = np.max(lats[vpi])
            if maxtemp > latmax:
                latmax = maxtemp

            temp = np.min(lons[vpi])
            if temp < lonmin:
                lonmin = temp
            maxtemp = np.max(lons[vpi])
            if maxtemp > lonmax:
                lonmax = maxtemp

        if dt > 24 and dt <= 48:  #Now find arrival times that occur in day2.
            vpi = np.where(np.logical_and(day[1] == 0, concarray >= thold))
            day[1][vpi] = dt
            #These two lines set places which had arrival times on day 0 to 0.
            vpi2 = np.where(day[0] != 0)
            day[1][vpi2] = 0
            temp = np.min(lats[vpi])
            if temp < latmin:
                latmin = temp
            maxtemp = np.max(lats[vpi])
            if maxtemp > latmax:
                latmax = maxtemp
            temp = np.min(lons[vpi])
            if temp < lonmin:
                lonmin = temp
            maxtemp = np.max(lons[vpi])
            if maxtemp > lonmax:
                lonmax = maxtemp

        if dt > 48 and dt <= 72: #find arrival times that occur in day3.
            vpi = np.where(np.logical_and(day[2] == 0, concarray >= thold))
            day[2][vpi] = dt
            #These two lines set places which had arrival times on day 0 to 0.
            vpi2 = np.where(day[0] != 0)
            day[2][vpi2] = 0
            #These two lines set places which had arrival times on day 0 to 0.
            vpi2 = np.where(day[1] != 0)
            day[2][vpi2] = 0
            temp = np.min(lats[vpi])
            if temp < latmin:
                latmin = temp
            maxtemp = np.max(lats[vpi])
            if maxtemp > latmax:
                latmax = maxtemp
            temp = np.min(lons[vpi])
            if temp < lonmin:
                lonmin = temp
            maxtemp = np.max(lons[vpi])
            if maxtemp > lonmax:
                lonmax = maxtemp


    ## Do General Map things ############################################
    if obounds:
        lllat = blatmin
        lllon = blonmin
        urlat = blatmax
        urlon = blonmax

        if blatmin > latmin or blatmax < latmax or blonmin > lonmin or blonmax < latmax:
            print('Warning: some portions of the data may be cutoff. \
                  suggested boundaries (latmin:latmax:lonmin:lonmax)', \
                 latmin, ':', latmax, ':', lonmin, ':', lonmax)
    else: #if bounds not specified then find them.
#       print('latmin, latmax, lonmin, lonmax ', \
#            latmin, ':', latmax, ':', lonmin, ':', lonmax)
        latpad = 5
        lllat = latmin - latpad
        lllon = lonmin - latpad
        urlat = latmax + latpad
        urlon = lonmax + latpad
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_axes([0.2,0.2,0.6,0.5])
    hatch = ['.', '\\', '-', '|']

    ####Following block is to preserve aspect ratio of the plot. If this is not done
    ### then text may overlap plot.
    aratio = (urlon - lllon) / (urlat-lllat)
#   print('latmin, latmax, lonmin, lonmax, aratio ', \
#        latmin, ':', latmax, ':', lonmin, ':', lonmax, ':', aratio)
    loncen = 0.5*(lllon + urlon)
    latcen = 0.5*(lllat + urlat)
    if aratio < scale:
        latdiff = urlat - lllat
#       londiff = latdiff * 16.0 / 9.0
#       londiff = (latdiff * 16.0 / 9.0) - (urlon - lllon)
#       if verbose:
#           print('adding ', londiff, 'to longitude range')
#       urlon = lllon + londiff
#       urlon += 0.5*londiff
#       lllon -= 0.5*londiff
        delx = latdiff*scale
        lllon = loncen-delx/2.0
        urlon = loncen+delx/2.0
    if aratio > scale:
        #latitude difference is too small
        londiff = urlon - lllon
#       latdiff = (londiff *  9 / 16.00) - (urlat - lllat)
#       if verbose:
#           print('adding ', latdiff, ' to latitude range')
#       urlat += 0.5 * latdiff
#       lllat -= 0.5 * latdiff
        dely = londiff/scale
        lllat = latcen-dely/2.0
        urlat = latcen+dely/2.0
    #####

    ##check to make sure min and max lat and lon are not out of bounds.
    if urlat > 90:
        urlat -= 90
    if lllat < -90:
        lllat = -90

    clevs = []
    clevsa = []

    ##To do - should contours be changeable?
    clevs.append([0.0001, 6.0001, 12.0001, 18.0001, 24.0001])
    clevsa.append([0.0001, 24.0001])
    clevs.append([25, 30.0001, 36.001, 42.001, 48.001])
    clevsa.append([0.0001, 48.001])
    clevs.append([49, 54, 60, 66, 72])
    #clevs = [6,12,18,24,30,36,42,48,54,60,66,72,78]

    ############################################## AMC - Commented out shape file
    ## SHAPEFILE INFO
    #SHAPEDIR="./GIS/"
    #SHAPEDIR2="/g/sc/ophome/nmoc_share/bin/hysplit/shapes/"
    #cmapz="cntry08"
    tmap = Basemap(llcrnrlon=lllon, llcrnrlat=lllat, urcrnrlon=urlon, urcrnrlat=urlat,
                   resolution='l', projection='cyl')

    x, y = tmap(lons, lats)
    ############# Annotation #########
    main_title = options.title

    for ii in range(0, numdays):
        #inputs to be passed/read
        if ii == 0:
            label = "0-24 hours"
        if ii == 1:
            label = "25-48 hours"
        if ii == 2:
            label = "49-72 hours"
        ###########################################

        ncols = ['#ccebc5', '#a8ddb5', '#7bccc4', '#43a2ca', '#0868ac']
        #ncols2 = ['#ccebc5', '#a8ddb5', '#7bccc4', '#43a2ca', '#0868ac']
        #ncols3 = ['#ccebc5', '#a8ddb5', '#7bccc4', '#43a2ca', '#0868ac']
        plt.figtext(0.45, 0.95, main_title, fontsize=16, ha='center', weight='bold')
        plt.figtext(0.45, 0.92, "Time Of Arrival Chart: "+label, color="k", fontsize=14,
                    ha='center', weight='bold')
        plt.figtext(0.45, 0.88, "Simulation begins "+mystart, color="b", fontsize=14,
                    ha='center', weight='bold')
        plt.figtext(0.45, 0.85, "Level displayed: " + mylevs, color="b", fontsize=14,
                    ha='center', weight='bold')

        plt.figtext(0.10, 0.15, "Source Location: "+myloc, color="b", fontsize=14, weight='bold')
        plt.figtext(0.10, 0.12, "Source Coordinates : "+mysrc, color="b", fontsize=14, weight='bold')
        plt.figtext(0.10, 0.09, "Release Begins: "+mystart, color="k", fontsize=14, weight='bold')
        plt.figtext(0.10, 0.06, "Release Rate: "+myrate+"/hr", color="k", fontsize=14, weight='bold')
        plt.figtext(0.10, 0.03, "NWP Data Source : "+str(mymodel), color="b", fontsize=14, weight='bold')

        plt.figtext(0.52, 0.09, "Release Duration: "+myrel+" hrs", color="k", fontsize=14,
                    weight='bold')
        plt.figtext(0.30, 0.06, "Radionuclide: "+myspec, ha='left', color="r", fontsize=14,
                    weight='bold')
        plt.figtext(0.52, 0.06, "Threshold for Arrival Time: "+mytv, color="k", fontsize=14,
                    weight='bold')

        plt.figtext(0.60, 0.03, "Chart Issued : "+isstime, color="g", fontsize=11, weight='bold')

        tmap.drawcoastlines(linewidth=1)
        tmap.drawcountries(linewidth=1, color='#0099FF')
        tmap.drawstates(linewidth=1)
        tmap.drawparallels(np.arange(-80, 81, 10), labels=[1, 1, 0, 0])
        #shpinfo1=tmap.readshapefile(SHAPEDIR2+cmapz,'coast',drawbounds=True,color='black',linewidth=1.0)
        tmap.drawmeridians(np.arange(-180, 180, 10), labels=[0, 0, 0, 1])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
#           cs = tmap.contourf(x, y, day[ii], clevs[ii], cmap=mpl.colors.ListedColormap(ncols),
#                          hatches=hatch, extend='neither')
            if ii == 0:
                cs11 = tmap.contour(x, y, day[ii],clevs[ii], linewidths=0.5, colors='k' )
                cs1 = tmap.contourf(x, y, day[ii], clevs[ii], cmap=mpl.colors.ListedColormap(ncols),
                           hatches=hatch, extend='neither')
            if ii == 1:
                cs21 = tmap.contour(x, y, day[ii-1],clevsa[ii-1], linewidths=0.5, colors='k' )
                cs2 = tmap.contourf(x, y, day[ii-1], clevsa[ii-1], colors='#bdbdbd',
                                linestyles='-', extend='neither')
                cs11 = tmap.contour(x, y, day[ii],clevs[ii], linewidths=0.5, colors='k' )
                cs1 = tmap.contourf(x, y, day[ii], clevs[ii], cmap=mpl.colors.ListedColormap(ncols),
                           hatches=hatch, extend='neither')
            if ii == 2:
                cs31 = tmap.contour(x, y, day[ii-2],clevsa[ii-2], linewidths=0.5, colors='k' )
                cs3 = tmap.contourf(x, y, day[ii-2], clevsa[ii-2], colors='#bdbdbd',
                                linestyles='-', extend='neither')
                cs21 = tmap.contour(x, y, day[ii-1],clevsa[ii-1], linewidths=0.5, colors='k' )
                cs2 = tmap.contourf(x, y, day[ii-1], clevsa[ii-1], colors='#bdbdbd',
                                linestyles='-', extend='neither')
                cs11 = tmap.contour(x, y, day[ii],clevs[ii], linewidths=0.5, colors='k' )
                cs1 = tmap.contourf(x, y, day[ii], clevs[ii], cmap=mpl.colors.ListedColormap(ncols),
                           hatches=hatch, extend='neither')
        c = plt.colorbar(orientation='vertical', format='%.d')
#       if ii > 0:  #This part shades the map gray for areas covered by previous days.
#           with warnings.catch_warnings():
#               warnings.simplefilter("ignore")
#               cs1 = tmap.contourf(x, y, day[ii-1], clevsa[ii-1], colors='#bdbdbd',
#                               linestyles='-', extend='neither')

        ###################################
        #plot the source point
        x1, y1 = tmap(float(slat), float(slon))
        tmap.plot(x1, y1, marker='p', color='r', markersize=9)
        pstr = 'toa_day' + str(ii+1) + '.png'
#       pstr2 = 'toa_day' + str(ii+1) + '.pdf'
#       pstr3 = 'toa_day' + str(ii+1) + '.jpg'
        if verbose:
            print("put marker at source point ", slat, " ", slon)
            print("plotting ", pstr)
        plt.savefig(pstr)
#       plt.savefig(pstr2)
#       plt.savefig(pstr3)
        plt.clf()

else:
    print('Warning: cdump file  ', fname, ' does not exist')

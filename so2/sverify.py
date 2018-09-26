import os
import subprocess
import pandas as pd
import numpy as np
import pickle as pickle
from optparse import OptionParser
import datetime
import sys
import seaborn as sns
import matplotlib.pyplot as plt
#from monet.obs import cems_mod
#from monet.obs import aqs_mod
#from monet.obs import airnow
#from monet.obs import ish_mod
#import monet.obs.obs_util as obs_util
#from arlhysplit import runh
#from arlhysplit.runh import date2dir
#from arlhysplit.runh import source_generator
#from arlhysplit.runh import create_plume
#from arlhysplit.tcm import TCM
#from arlhysplit.models import emittimes
#from arlhysplit.models.datem import writedatem_sh
#from arlhysplit.models.datem import frame2datem
#from arlhysplit.models.datem import mk_datem_pkl
#from monet.obs.epa_util import convert_epa_unit

"""
INPUTS: Dates to run
        Area to consider emissions from
        list of states to consider emissions from
        top level directories for
            1. code executables
            2. where output should be located.
            3. (optional) where any csv files are with emissions data.

STEPS
A. Preliminary.
1. Find Emission sources.
2. Find measurement stations in area 
3. Produce map of sources and measurements.
4. Create plots of emissions vs. time
5. Create plots of measurements vs. time

B.1 trajectory runs

B.2 Dispersion runs 
Create HYSPLIT input files and run HYSPLIT
1. Use emissions source to create HYSPLIT sources
TODO - need to be able to specify stack height for each source.
       Stack height for sources not available from current data sources.
       These must be looked up individually.
TODO - resolve ambiguity in local / daylight savings time in the CEMS.
2. Use measurement stations to create concenctration grids
   OR just use a large concentration grid that covers alls stations.
3. Run HYSPLIT.
TODO - should these be unit runs so a TCM can be creatd later or do runs in
some kind of chunks. 

C. Evaluate HYSPLIT output and measurements.
1. create dfile (datem file) for the measurements in the HYSPLIT output directories.
TODO check that averaging time of measurements and HYSPLIT matches 
2. Run c2datem and collect output from cfiles in each directory.
3. Create a modeled measurements vs. time using the output from c2datem
4. Compare the modeled and measured concentrations.
TODO what kinds of statistics?
5. create plume plots (TO DO)

EPA method code 60  = EQSA-0486-060

"""


def create_map(fignum):
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    fig = plt.figure(fignum)
    proj = ccrs.PlateCarree()
    ax = plt.axes(projection=proj)
    gl = ax.gridlines(draw_labels=True, linewidth=2, color='gray')
    gl.ylabels_right=False
    gl.xlabels_top = False
    states = cfeature.NaturalEarthFeature(category='cultural',
             name='admin_1_states_provinces_lines', scale='50m',
             facecolor='none')
    ax.add_feature(states, edgecolor='gray')
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.LAKES)
    ax.add_feature(cfeature.RIVERS)
    ax.add_feature(cfeature.COASTLINE)
    return(ax)



parser = OptionParser()

parser.add_option('-a', type="string", dest="area", default="ND",\
                  help='two letter state code (ND)')
parser.add_option('-d', type="string", dest="drange", \
                  default="2106:1:1:2016:2:1", \
                  help='daterange YYYY:M:D:YYYY:M:D')
parser.add_option('--cems', action="store_true", dest="cems", default=False)
parser.add_option('--obs', action="store_true", dest="obs", default=False)
##-----##
parser.add_option('--run', action="store_true", dest="runh", default=False)
parser.add_option('--map', action="store_true", dest="emap", default=False)
parser.add_option('--plote', action="store_true", dest="eplot", default=False, \
                  help='plot emissions')
parser.add_option('--datem', action="store_true", dest="datem", default=False)
parser.add_option('--rundatem', action="store_true", dest="rundatem", default=False)
parser.add_option('--pickle', action="store_true", dest="pickle", default=False)
parser.add_option('--tcm', action="store_true", dest="tcm", default=False)
parser.add_option('--test', action="store_true", dest="runtest", default=False)
#parser.add_option('--obs', action="store_true", dest="findobs", default=False)
parser.add_option('-x', action="store_false", dest="opkl", default=True)
(options, args) = parser.parse_args()

opkl = options.opkl

temp = options.drange.split(':')
try:
    d1 = datetime.datetime(int(temp[0]), int(temp[1]), int(temp[2]), 0)
except:
    print('daterange is not correct ' + options.drange)
try:
    d2 = datetime.datetime(int(temp[3]), int(temp[4]), int(temp[5]), 0)
except:
    print('daterange is not correct ' + temp)

if options.area.lower().strip() == 'nd':
    area = [-105.0, -97.0, 44.5, 49.5]
    state=['nd']
else:
    area = None
    state=[options.area.strip()]

#sv = SO2Verify([d1,d2], area, state)

##emissions are on order of 1,000-2,000 lbs (about 1,000 kg)
##10,000 particles - each particle would be 0.1 kg or 100g.
##0.05 x 0.05 degree area is about 30.25 km^2. 30.25e6 m^2.
##50 meter in the vertical gives 1.5e9 m^3.
##So 1 particle in a 0.05 x 0.05 degree area is 0.067 ug/m3.
##Need a 100 particles to get to 6.7 ug/m3.
##This seems reasonable.

rfignum=1
if options.cems:
    from semissions import SEmissions
    ef = SEmissions([d1,d2], area, state)
    ef.find()
    ef.plot()
    ef.create_emittimes(ef.d1)
    rfignum = ef.fignum 
    if not options.obs:
        mapfig = plt.figure(fignum)
        axmap = create_map(rfignum)
        ef.map(axmap)
        plt.show()

if options.obs:
    from sobs import SObs
    obs = SObs([d1,d2], area, state)
    obs.fignum=rfignum
    obs.find(pload=opkl)
    obs.obs2datem() 
    obs.plot()
    fignum = obs.fignum 
    axmap = create_map(fignum)
    obs.map(axmap)
    if options.cems:
       ef.map(axmap)
    plt.show()


##------------------------------------------------------##
if options.runtest:
   sv.testsources()

#if options.plume:
#   sv.plotplume()

#if options.findobs:
#    sv.find_obs(pload=opkl)

if options.runh:
    sv.find_emissions()
    sv.runHYSPLIT()

if options.emap:
    ax = create_map(fignum=1)
    sv.find_emissions()
    sns.set()
    sv.map_emissions(ax)
    sv.find_obs(pload=opkl)
    sv.map_obs(ax)
    sv.obs2datem() 
    plt.show()

#if options.oplot:
#    from sobs import SObs
#    obs = SObs([d1,d2], area, state)
#    obs.find(pload=opkl)
#    obs.plot()

if options.datem:
    sv.find_obs(pload=opkl)
    sv.obs2datem() 
    #sv.rundatem() 

if options.rundatem:
    sv.rundatem() 

if options.pickle:
   sv.mkpkl()

if options.tcm:
   sv.mktcm()
   sv.find_obs(pload=opkl)
   sv.plot_both()

# write_Raikoke_emitimes_newarea.py
# Writing necessary variables from VOLCAT file to HYSPLIT emitimes file
# Using new area calculation for mass emission rate
# Written by: Allison M. Ring
import xarray as xr
import numpy as np
import numpy.ma as ma
from utilvolc import volcat
#from math import pi, cos, sin, radians
from monet.util import csi_updated as csi
from datetime import datetime
from os import path

#hh = ['00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23']
#mm = ['00','10','20','30','40','50']
dd = ['22', '23', '24']
mm = ['40']
hh = ['02']
#hh = ['19','20','21','22','23']
#dd = ['21']
directory = '/pub/ECMWF/JPSS/VOLCAT/Raikoke/Ash/'
direct = '/hysplit-users/allisonr/Raikoke/EMIT_TIMES/Data_Insertion/VOLCAT/'

# Hard coded values
# number of particles, particle distrubution, heat, simulation run time, duration of emission
#par = 4
par = 1
#par_dist = [0.008, 0.068, 0.25, 0.67]
par_dist = [1.0]
Heat = '0.00e+00'
#run_time = 36
duration = '0010'

# Looping over days, hours, and minutes for each VOLCAT file
k = 0
while k < len(dd):
    i = 0
    while i < len(hh):
        j = 0
        while j < len(mm):
            # Creating datetime object
            dt = datetime(2019, 6, int(dd[k]), int(hh[i]), int(mm[j]))
            print('k: '+str(k)+' i: '+str(i)+' j: '+str(j))

            # Reading in VOLCAT file
            fname = 'SCOPE_NWC_ASH-L2-ASH_PRODUCTS-HIMAWARI8_NOAA-RAIKOKE-201906' + \
                dd[k]+'-'+hh[i]+mm[j]+'00-fv2.nc'
            # Check if file exists
            if path.exists(directory+fname):
                dset = volcat.open_dataset(directory+fname)
            # Extracting ash top height and mass loading
                height = volcat.get_height(dset)
                height2 = height.values * 1000.  # Put into meters
                massload = volcat.get_mass(dset).values

            # Creating lat and lon arrays
                latitude = height.latitude.values
                longitude = height.longitude.values

            # Calculating pixel area
                fname2 = 'area_whole_domain.nc'
                dir2 = '/hysplit-users/allisonr/Raikoke/'
                dset2 = xr.open_dataset(dir2+fname2)
                # Need to drop area array indices that are not within mass array
                area = dset2.area
                area2 = csi.get_area(dset, dset2).values

            # Calculating mass - rate is (mass/hr) in HYSPLIT
                # area needs to be in m^2 not km^2!!
                ashmass = massload * area2 * 1000. * 1000.

            # Checking to see if height and ashmass arrays are the same size
            # If not, determines which array is smaller and trims the other to same size
                hgtsz = np.size(height2)
                ashsz = np.size(ashmass)
                # If height array is LARGER than ash array
                if hgtsz != ashsz:
                    diff1 = hgtsz - ashsz
                    m2 = dset.ash_mass
                    m2 = m2[0, :, :]
                    if diff1 > 0.:
                        height2 = height.where(m2 != m2._FillValue, drop=True).values * 1000
                    # If height array is SMALLER than ash array
                    if diff1 < 0.:
                        h2 = dset.ash_cth
                        h2 = h2[0, :, :]
                        area2 = area.where(h2 != h2._FillValue, drop=True).values
                        massload = m2.where(h2 != h2._FillValue, drop=True).values
                        ashmass = massload * area2 * 1000. * 1000.
            # Creating compressed 1-d arrays (removing nans) to write to a file
                hgt_nan = ma.compressed(height2)  # array with nan values
                hgt = hgt_nan[~np.isnan(hgt_nan)]  # removing nan values
                mass = ma.compressed(ashmass)[~np.isnan(hgt_nan)]  # removing nan values
                lat = ma.compressed(latitude)[~np.isnan(hgt_nan)]  # removing nan values
                lon = ma.compressed(longitude)[~np.isnan(hgt_nan)]  # removing nan values
                area3 = ma.compressed(area2)[~np.isnan(hgt_nan)]  # removing nan values
                area3 = area3 * 1000000  # Moving from km^2 to m^2

            # Creating data file in emitimes format
                f = open(direct+'VOLCAT_newarea.Raikoke_201906'+dd[k]+'.'+hh[i]+mm[j]+'00.par'+str(par), 'w')
                f.write('YYYY MM DD HH        DURATION(hhhh) #RECORDS \n')
                f.write('YYYY MM DD HH MM Duration(hhmm) LAT LON HGT(m) RATE(g/hr) AREA (m2) HEAT(w) \n')
                f.write('{:%Y %m %d %H  } {} {}\n'.format(dt, duration, np.shape(hgt)[0]*par))
                h = 0
                while h < len(hgt):
                    l = 0
                    while l < len(par_dist):
                        f.write('{:%Y %m %d %H %M} {} {:9.6f} {:10.6f} {:8.2f} {:.2E} {:.2E} {} \n'.format(
                            dt, duration, lat[h], lon[h], hgt[h], mass[h]*6*par_dist[l], area3[h], Heat))
                        l += 1
                    h += 1

                f.close()
            j += 1  # mm loop
        i += 1  # hh loop
    k += 1  # dd loop

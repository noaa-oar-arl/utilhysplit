# make_L1L2_netcdf.py
# Makes netcdf of ash mass loading, ash cloud height, ash effective radius, lat, lon
# From L1 and L2 VOLCAT files
from utilvolc import volcat
from glob import glob

vdir = '/pub/ECMWF/JPSS/VOLCAT/Vincent/Ash/'
fname1 = glob(vdir+'*L1*2021104.02*')
fname2 = glob(vdir+'*L2*2021104.02*')

a = 0
while a < len(fname1):
    volcat.create_netcdf(fname1[a], fname2[a])
    a += 1

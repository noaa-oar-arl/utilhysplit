# write_cylsource_emitimes.py
# Writes a HYSPLIT EMITIMES.txt file for a cylindrical volcanic source
from utilhysplit import cylsource
from utilhysplit import emitimes
from datetime import datetime

# Write cylinder EMITIMES.txt file
# Consistent Inputs
wdir = '/hysplit-users/allisonr/Raikoke/EMIT_TIMES/Cylinder/'
lat = 48.292
lon = 153.385
# Radius of cylinder around volcano vent (m)
radius = 20000
#radius = 10000
# dr should be evenly divisible into radius
#dr = 10000
dr = 5000
# Duration of ash emission
duration = '1200'
pollnum = 1
pollpercents = [1]
altitude = [2000, 16000]
# Need to change if doing Umbrella (2: emissions in upper half of column, 3: emissions in upper third of column)
umbrella = 1

# Creating date loops
hrs = [18]
mins = [10]
day = 21

i = 0
while i < len(hrs):  # Loop over hours
    j = 0
    while j < len(mins):  # Loop over mins
        # Varying Inputs
        dt = datetime(2019, 6, day, hrs[i], mins[j])
        fname = 'EMIT_Raikoke_cyl_'+dt.strftime("%Y%m%d%H%M")+'_12hrsUmbrella.txt'
        #fname = 'EMIT_Raikoke_cyl_'+dt.strftime("%Y%m%d%H%M")+'_3hrsUmbrella.txt'
        #fname = 'EMIT_Raikoke_cyl_'+dt.strftime("%Y%m%d%H%M")+'_10hrs.txt'
        #fname = 'EMIT_Raikoke_cyl_'+dt.strftime("%Y%m%d%H%M")+'.txt'
        filename = wdir+fname

        # Creating emitimes files
        efile = cylsource.EmitTimes(filename=filename)
        latlist, lonlist = efile.calc_cyl(lat, lon, radius, dr)
        nrecords = efile.calc_nrecs(latlist, pollnum=pollnum, umbrella=umbrella)
        efile.write_data(dt, latlist, lonlist, nrecords, pollnum, duration=duration,
                         pollpercents=pollpercents, height=altitude)

        j += 1
    i += 1

# Read number of records
#efile = emitimes.EmiTimes(filename)
#efile.read_file(num_species = 1)
#recs = efile.findmaxrec()
# print(recs)

import datetime

from utilhysplit import hcontrol, metfiles

# Read a default CONTROL file.
control = hcontrol.HycsControl(fname='CONTROL.0',
                               working_directory = './',
                               )
control.read(verbose=True)

# change the runtime
runtime = 24
control.add_duration(runtime)

# change the start date
d1 = datetime.datetime.now()
control.add_sdate(d1)

# change the meteorological files using the MetFiles class.
metfmt = '/pub/forecast/%Y%m%d/hysplit.t00z.href.m01'
control.remove_metfile(rall=True)
met_files = metfiles.MetFiles(metfmt)
mfiles = met_files.get_files(control.date, runtime) 
print(mfiles)
for mf in mfiles:
    control.add_metfile(mf[0],mf[1])

# change the CONTROL file name
suffix = 'test'
control.fname = 'CONTROL.{}'.format(suffix)

# write the CONTROL file
control.write(metgrid=1,query=False,overwrite=True)
print('writing control {} {}'.format(control.wdir, control.fname)) 

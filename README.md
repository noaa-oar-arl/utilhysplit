python code for reading and creating HYSPLIT input files as well as processing HYSPLIT output files.


pardump.py contains class for reading and writing HYSPLIT binary particle (PARDUMP) file.

hcontrol.py contains classes and functions for writing and reading HYSPLIT input files, CONTROL, SETUP.CFG and ASC.CFG

isoch.py creates time of arrival charts. Needs a concentration (cdump) file as input. Optional input is CONTROL file.

UPDATES
3/25/2019 A lot of development was moved to other repositories including monet. Now, some of it has been moved back. The models directory was
changed to the utilhysplit directory and most of the functional code is now there.

python code for reading and creating HYSPLIT input files as well as processing HYSPLIT output files.


* ashapp directory contains code for creating HYSPLIT runs. Dependencies on utilvolc and utilhysplit.
* utilvolc directory contains code specific to volcanic ash applications. Dependencies on utilhysplit.
* utilhysplit directory contains general utilities for working with HYSPLIT. 
* utilml directory contains code which has dependency on the sklearn module
* utiltesting contains tests

hcontrol.py contains classes and functions for writing and reading HYSPLIT input files, CONTROL, SETUP.CFG and ASC.CFG

isoch.py creates time of arrival charts. Needs a concentration (cdump) file as input. Optional input is CONTROL file.

UPDATES
3/25/2019 A lot of development was moved to other repositories including monet. Now, some of it has been moved back. The models directory was changed to the utilhysplit directory and most of the functional code is now there.

5/8/2020 pardump.py was moved to monetio

6/11/2020 Added utilvolc folder: contains utilities for volcano hysplit applications and volcat data manipulation

11/5/2022 start with new conda environment to look at dependencies.
install xarray, cartopy, shapely, netCDF4, lxml, seaborn
dask dependency in monetio.
removed everything in the the __init__.py in monetio.obs directory so don't need dask.

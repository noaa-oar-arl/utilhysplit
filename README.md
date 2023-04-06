python code for reading and creating HYSPLIT input files as well as processing HYSPLIT output files.

hcontrol.py contains classes and functions for writing and reading HYSPLIT input files, CONTROL, SETUP.CFG and ASC.CFG

isoch.py creates time of arrival charts. Needs a concentration (cdump) file as input. Optional input is CONTROL file.

UPDATES

3/25/2019 A lot of development was moved to other repositories including monet. Now, some of it has been moved back. The models directory was changed to the utilhysplit directory and most of the functional code is now there.

5/8/2020 pardump.py was moved to monetio

6/11/2020 Added utilvolc folder: contains utilities for volcano hysplit applications and volcat data manipulation

11/5/2022 start with new conda environment to look at dependencies.
install xarray, cartopy, shapely, netCDF4, lxml, seaborn
xrarray needs h5netcdf to support encoding zlib when writing to a netcdf file.:w

4/6/2023 Update the code

dependency on monet. correcting parallax uses the nearest_ij function in the monet accessor for the xarray.
dask dependency in monetio.
removed everything in the the __init__.py in monetio.obs directory so don't need dask.

11/7/2022 dependency.
utilhysplit should contain only dependency on external approved modules.
utilvolc can contain dependency on utilhysplit
ashapp can contain dependency on utilvolc and utilhysplit


11/10/2022
code which has dependency on sklearn moved to utilml subdirectory.


DEPENDENCIES
ashapp directory

ash_run.py
   ashensemble.py
   ashbase.py
   ashinverse.py
   ashdatainsertion.py
   ashtrajectory.py 
   metfiles.py (gefs_suffix_list)
   runhelper 
   utils (setup_logger)


ashbase.py 
   utilvolc.ashapp.metfiles
   utilhysplit.hcontrol
   monetio.models hysplit
   utilvolc.ashapp.runhelper
   utilvolc.volcMER (HT2Unit)

ashtrajectory.py
   utilvolc.ashapp.runhelper (Helper)
   utilvolc.ashapp.ashbase (AshRun)

ashinverse.py
   utilvolc.ashapp.runhelper (Helper)
   utilvolc.ashapp.runhandler (ProcessList)
   utilvolc.ashapp.ashbase (AshRun)
   utilvolc.ashapp.ensemble_tools (massload_plot)
   monetio.models hysplit

ashensemble.py
   utilvolc.ashapp.ashbase (AshRun)
   utilvolc.ashapp.metfiles (gefs_suffix_list)
   utilvolc.ashapp.runhelper (Helper)
   utilvolc.ashapp.runhandler (ProcessList)
   ensemble_tools.py 
   utilvolc.ashapp.cdump2xml (HysplitKml)

ashdatainsertion.py
   utilvolc.ashapp.ashbase (AshRun)
   utilvolc.ashapp.runhelper (AshDINameComposer)
   utilvolc.ashapp.runhandler (ProcessList)
   monetio.models hysplit
   
cdump2netcdf.py
   monetio.models (hysplit)

cdumpxml.py
   None
   datetime, loggin, os, zipfile, matplotlib, pandas, lxml

ensemble_tools.py
   NOT tracked?
   cartopy
   None

metfiles.py (also in utilhysplit)
   None

emitimes.py (also in utilhysplit)
   None
   numpy, datetime

runhandler.py
   None

runhelper.py
   utilhysplit.hcontrol (NameList)
   datetime, glob, loggin, os, pathlib, shutil, subprocess, sys

so2base.py
   # NOT finished
   # NOT tracked

volc_inverse.py
   not tracked. old version of ashinverse.py?

utils.py
   None
   logging, sys 



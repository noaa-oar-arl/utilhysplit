python code for reading and creating HYSPLIT input files as well as processing HYSPLIT output files.


* utilhysplit directory contains general utilities for working with HYSPLIT.
   * hcontrol.py contains classes and functions for reading and writing HYSPLIT CONTROL, SETUP.CFG and ASC.CFG files.
       * example.py shows example of how to use classes in hcontrol.py
   * emitimes.py contains classes and functions for reading and writing emit-times files.
       * emitimes_example.py shows examples of use of classes in emitimes.py 
   * metfiles.py contains classes and functions to automate adding meteorological data to CONTROL files.
       * metexamples.py example use of the getmetfiles function.
       * tests/test_metdata.py has tests for the MetFileFinder class which can also be used as examples.
   * message_parse.py will read in the MESSAGE file and can create some basic plots of things like time step, emrise, 
   * hysplit_gridutil.py contains some functions which can be useful for combining different xarray datasets of cdump data.
   * vmixing.py classes and functions for reading in and plotting vmixing output. 
   * xtrct_stn.py classes and functions for reading in and plotting output from xrct_stn.py
   * hysplit_profile.py classes and functions for reading and plotting ouput of profile program
   * runhandler.py class for handling multiple runs. Can be used for instance if need to do many HYSPLIT runs but only want to have 10 running at a time. runhandler can help start new runs as old ones finish.
   * geotools.py contains some functions for 
       * computing concave hull around a set of points
       * finding distance and bearing between two points
   * arlmet.py class for reading ARL formatted files. Does not currently work.
* utilhysplit/particlesize
   * contains some functions for computing fall velocity by various methods.
* utilhysplit/plotutils
   * colormaker.py class to create list of colors. can produce colors suitable for kml files. 
* utilhysplit/evaluation 
   * plume_stat.py calculate contingency table and derived statistics.
   * reliability.py classes for calculating reliability diagrams and rank histograms.
   * statmain.py can compute most of the statistics in statmain function. contains a MatchedData class.
   * hysplit_boxplots.py
   * ensemble_tools.py and ensemble_stats.py for working with xarrays which have dimensions of latitude,longitude,height,time,ensemble_number
   * datem.py for handling data in datem format. 
* utilml directory contains code which has dependency on the sklearn module
* utilvolc directory contains code specific to volcanic ash applications. Dependencies on utilhysplit.
* ashapp directory contains code for creating automated HYSPLIT runs. Dependencies on utilvolc and utilhysplit.
* utiltesting contains tests

UPDATES

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


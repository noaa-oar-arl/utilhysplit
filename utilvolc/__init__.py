"""Volcanic Emissions Forecasts and Data

VOLCAT ash
    volcat.py         opening and manipulating volcat netcdf files with ash data
    volcat_2dplots.py quick plotting of volcat data
    volcat_plots.py   plots summarizing many volcat files
    get_area.py       function for calculating area of grid cells
    process_volcat.py executable for writing the parallax corrected files
    volcalert.py      reads volcat alert files

Workflow
    qva_logic.py      workflow for creating HYSPLIT forecasts with data fusion
    dtwin.py          digitial twin. in progress

Inversion Algorithm
    ash_inverse.py    create TCMS and run inversion algorithm
    inverse_plots.py  plots of emissions from inversion algorithm

Data Insertion
    make_data_insertion.py creates emitimes files from VOLCAT netcdf files.

Evaluation and Verification
    ash_eval.py       evaluation and verification of ash forecasts

SO2
    volcat_so2.py     opening and manipulating volcat netcdf files with SO2 data.
    so2_inverse.py    create TCMS and run inversion using SO2 volcat data
    trajectory_analysis.py  under development
    plot_CrIS_so2.py

Volcano data 
    volcMER.py        functions using Mastin equation.
    usgstable.py      parses USGS table of volcano information

Volcanic ash advisories
    iwxxmVAA.py      for getting and reading VAA in IWXXM format from Washington VAAC webpage.
    poly2netcdf.py   can rasterize polygon into xarray grid. 

Other satellite data
    caliophdf.py     for reading HDF files with CALIOP data
    tropomi.py       for reading tropomi data

Other data
    sondes.py        for analyzing sonde data for Turrialba from Garry Morris.


Possible Legacy files
make_netcdf.py
vaa.py           reads text vaa files. very fragile. hasn't been updated in awhile.


LEGACY files
volcat_legacy.py  functions no longer maintained.
make_cylinder_source.py no longer needed now that HYSPLIT will initialize in circular area.

"""

from . import volcalert, volcat, volcMER, volcat_so2, read_csv

__all__ = ["volcalert", "volcat", "volcMER", "volcat_so2", "read_csv"]

__name__ = "utilvolc"

#!/opt/Tools/anaconda3/envs/hysplit/bin/python
# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# ash_run.py - run HYSPLIT model on web and create plots
#
# 01 JUN 2020 (AMC) - adapted from locusts-run.py
# 09 DEC 2022 (AMC) - changed latspan to 180 and lonspan to 360 for global grid.
# 09 DEC 2022 (AMC) - changed numpar to 10000. this should not be hardwired.
# 28 FEB 2023 (AMC) - change to write_cxra so if netcdf file exists then set it as self.cxra
# 28 FEB 2023 (AMC) - check attributes added as stand alone function but still kept as static method.
#
# -----------------------------------------------------------------------------
# To run in offline mode use python ash_run.py -999
#
#
# -----------------------------------------------------------------------------


# from abc mport ABC, abstractmethod
import datetime
import glob
import logging
import os
import shutil
import subprocess
import sys
import zipfile

import numpy as np
import requests
import xarray as xr

from monetio.models import hysplit
from utilhysplit import hcontrol
import utilhysplit.metfiles as metfile
from utilvolc.runhelper import  Helper, JobFileNameComposer
from utilhysplit.plotutils.concplotutil import ConcplotColors

# from runhandler import ProcessList
from utilvolc.volcMER import HT2unit

# from cdump2xml import HysplitKml

# from ashbase import AshRun
# from ashensemble import EnsembleAshRun

# from locusts import LocustsFileNameComposer, Helper, LocustsSetUpParser, FlightPlannerFactory, \
#                    SingleHeightSource, LocustSwarm
logger = logging.getLogger(__name__)

def check_attributes(atthash):
    # when writing to netcdf file, attributes which are numpy arrays do not write properly.
    # need to change them to lists.
    for key in atthash.keys():
        val = atthash[key]
        if isinstance(val, np.ndarray):
            newval = list(val)
            atthash[key] = newval
    return atthash


def print_usage():
    print(
        """\
USAGE: ash_run.py JOBID

The following environment variables must be set prior to calling this script:
    RUN_API_KEY         - secret key to access Locusts web APIs.
    RUN_URL             - URL to the Locusts web application."""
    )


# Base class is AshRun
# EnsembleAshRun inherits from AshRun
# TrajectoryAshRun inherits from AshRun


# def meters2FL(meters):
#    flight_level = meters * 3.28084 / 100.0
#    return int(np.ceil(flight_level / 10.0) * 10)


def FL2meters(flight_level):
    meters = flight_level * 100 / 3.28084
    return int(meters)
    # return int(np.ceil(flight_level / 10.0) * 10)


def make_chemrate(wdir):
    fname = "CHEMRATE.TXT"
    fstr = "1 2 0.01 1.0"
    fpath = os.path.join(wdir, fname)
    with open(fpath, "w") as fid:
        fid.write(fstr)


# class SO2Run(AshRun):
#
#    def __init__(self, JOBID):
#        super().__init__(JOBID)
#        self.default_control = 'CONTROL.so2'
#        self.default_setup =  'SETUP.so2'

#    def add_inputs(self, inp):
#        super().add_inputs(inp)
#        make_chemrate(self.inp['WORK_DIR'])


# def ncfile_encoding(dset):
#    encoding = {}
#    ekeys = {"_FillValue", "dtype","scale_factor",'add_offset','grid_mapping')
#    for dvar in dset.data_vars:
#        encoding[data_var] = {key: v:q


class AshRun:
    def __init__(self, JOBID):
        self.JOBID = str(JOBID)  # string
        self.inp = {}  # dictionary from JobSetUP

        self.default_control = "CONTROL.default"
        self.default_setup = "SETUP.default"

        self.apistr = None
        self.urlstr = None
        self.headerstr = None

        self.metfilefinder = None
        self.filelocator = None
        self.maptexthash = {}
        self.awips = True

        # read cdump file into this xarray.
        self.cxra = xr.DataArray()

        self.so2 = False

    def get_cdump_xra(self):
        if self.cxra.ndim == 0:
            blist = []
            cdumpname = self.filelocator.get_cdump_filename(stage=0)
            source_tag = "Line to {:1.0f} km".format(self.inp["top"] / 1000.0)
            met_tag = self.inp["meteorologicalData"]
            blist = [(cdumpname, source_tag, met_tag)]
            century = 100 * (int(self.inp["start_date"].year / 100))
            # If SO2 then only put SO2 in the xarray
            if self.so2:
                species = ["pS02"]
            else:
                species = None
            cdumpxra = hysplit.combine_dataset(
                blist, century=century, sample_time_stamp="start", species=species
            )
            self.cxra = cdumpxra
        return self.cxra

    def inp2attr(self):
        """
        return dictionary where all values are strings.
        This is so it can be written as attributes to netcdf file.
        """
        atthash = {}
        for key in self.inp.keys():
            try:
                val = str(self.inp[key])
                logger.info("attribute can be written as string {}".format(key))
                logger.info("attribute can be written as string {}".format(val))
            except:
                logger.warning("attribute cannot be written as string {}".format(key))
                val = "skip"
            if "/" in val:
                val = "skip"
            if val != "skip":
                atthash[key] = val
        return atthash

    def write_cxra(self):
        # convert to mg /m3
        fname = "xrfile.{}.nc".format(self.JOBID)
        mult = self.get_conc_multiplier()

        # check if file exists
        if os.path.isfile(fname):
            logger.info("netcdf file exists. Opening {}".format(fname))
            cxra = xr.open_dataset(fname)
            # get first data variable
            logger.info('data variable {}'.format(list(cxra.data_vars.keys())[0]))
            cxra = cxra[list(cxra.data_vars.keys())[0]]
            pmult = cxra.attrs["mult"]

            # if it exists and multiplication factor hasn't changed then doing nothing more.
            if np.isclose(pmult, mult, atol=0, rtol=0.1):
               self.cxra = cxra
               return []
          
            # if multiplication has changed then remove the old file and write new one.
            logger.info('changing mult factor {} {}'.format(pmult, mult))
            cxra = mult * cxra / pmult
            self.cxra = cxra
            # remove netcdf file so new one can be written.
            Helper.remove(os.path.join(self.inp["WORK_DIR"], fname))
        # if file does not already exist then create it.
        else:
            logger.info("netcdf file does not exist. Creating {}".format(fname))
            self.cxra = mult * self.get_cdump_xra()

        self.cxra = self.cxra.assign_attrs({"mult": mult})
        self.cxra = self.cxra.assign_attrs(self.inp2attr())
        logger.info("writing nc file {}".format(fname))
        # cxra.attrs.update(self.inp2attr())

        # need to change to dataset to add compression
        self.write_with_compression(self.cxra, fname)

        if not self.cxra.size > 1:
            logger.info(
                "make_awips_netcdf: cxra empty. cannot create awips\
                         files"
            )
            return []

    # when writing to netcdf file, attributes which are numpy arrays do not write properly.
    # need to change them to lists.
    @staticmethod
    def check_attributes(atthash):
        for key in atthash.keys():
            val = atthash[key]
            if isinstance(val, np.ndarray):
                newval = list(val)
                atthash[key] = newval
        return atthash

    def write_with_compression(self, cxra, fname):
        atthash = self.check_attributes(cxra.attrs)
        cxra = cxra.assign_attrs(atthash)
        cxra2 = cxra.to_dataset(name='POL')
        ehash = {"zlib": True, "complevel": 9}
        vlist = [x for x in cxra2.data_vars]
        vhash = {}
        for vvv in vlist:
            vhash[vvv] = ehash
        cxra2.to_netcdf(fname, encoding=vhash)

    def make_awips_netcdf(self):
        import cdump2netcdf

        ghash = {}
        ghash["source_latitude"] = self.inp["latitude"]
        ghash["source_longitude"] = self.inp["longitude"]
        ghash["source_name"] = self.inp["VolcanoName"]
        ghash["emission_start"] = self.inp["start_date"]
        ghash["emission_duration_hours"] = self.inp["emissionHours"]
        # in mg
        mult = self.get_conc_multiplier()
        mer = mult / 1e6 / 3600  # kg released in a second.
        ghash["MER"] = mer
        ghash["MER_unit"] = "kg/s"
        logger.debug("MULT value for awips {:3e}".format(mult))
        awipsname = self.filelocator.get_awips_filename(stage=0)
        c2n = cdump2netcdf.Cdump2Awips(
            self.cxra, awipsname, munit="mg", jobid=self.JOBID, globalhash=ghash
        )
        awips_files = c2n.create_all_files()
        # returns list of awips files that were created.
        return awips_files

    # def postprocessing(self):
    #    return -1

    # def preprocessing(self):
    #    # no preprocessing for this workflow
    #    return -1

    def add_api_info(self, apistr, urlstr, headerstr):
        self.apistr = apistr
        self.urlstr = urlstr
        self.headerstr = headerstr

    def add_inputs(self, inp):
        logger.info("adding inputs")
        self.inp = inp
        self.filelocator = JobFileNameComposer(
            self.inp["WORK_DIR"], self.JOBID, self.inp["jobname"]
        )
        self.metfilefinder = metfile.MetFileFinder(self.inp["meteorologicalData"])
        self.metfilefinder.set_forecast_directory(self.inp["forecastDirectory"])
        self.metfilefinder.set_archives_directory(self.inp["archivesDirectory"])
        self.maptexthash = self.get_maptext_info()
        if "control" in self.inp.keys():
            self.default_control = self.inp["control"]
        if "setup" in self.inp.keys():
            self.default_setup = self.inp["setup"]
        logger.warning("Using control {}".format(self.default_control))
        logger.warning("Using setup {}".format(self.default_setup))

    def update_run_status(self, jobId, status):
        if self.apistr:
            API_KEY = os.environ[self.apistr]
            RUN_URL = os.environ[self.urlstr]
            statusUrl = "{}/status/{}/{}".format(RUN_URL, jobId, status)
            r = requests.put(statusUrl, headers={self.headerstr: API_KEY})
        else:
            logger.info("Running in offline test mode")
        logger.info("Posted status {} for job {}".format(status, jobId))

    def handle_crash(self, stage=0):
        self.update_run_status(self.JOBID, "CRASHED")
        logger.info("The model has crashed for job {}.".format(self.JOBID))
        logger.info(datetime.datetime.now())

    def setup_basic_control(self, stage=0, rtype="dispersion"):
        """
        inp : dictionary with inputs
        used for both dispersion and trajectory runs.
        """
        # logger.debug("Setting up control stage {}".format(stage))
        duration = self.inp["durationOfSimulation"]
        stime = self.inp["start_date"]
        # c.jobid_str = self.JOBID
        # work_dir = self.WORK_DIR

        # read the default CONTROL file
        control = hcontrol.HycsControl(
            fname=self.default_control,
            working_directory=self.inp["DATA_DIR"],
            rtype=rtype,
        )
        control.read()
        # rename the file
        newname = self.filelocator.get_control_filename(stage)
        control.rename(newname, self.inp["WORK_DIR"])

        # set start date
        control.date = stime

        # add met files
        # finder functionality may need to be updated.
        control.remove_metfile(rall=True)
        metfiles = self.metfilefinder.find(stime, duration)
        for mfile in metfiles:
            logger.debug(os.path.join(mfile[0], mfile[1]))
            control.add_metfile(mfile[0], mfile[1])
        if not metfiles:
            self.update_run_status(self.JOBID, "TERMINATED. missing met files")
            sys.exit()
        # add duration of simulation
        control.add_duration(duration)
        return control

    def set_qva_levels(self):
        # every 5000 ft (FL50 chunks)
        # Approx 1.5 km.
        levlist_fl = range(50, 650, 50)
        levlist = [FL2meters(x) for x in levlist_fl]
        rlist = []
        plev = "SFC"
        for lev in levlist_fl:
            nlev = "FL{}".format(lev)
            rlist.append("{} to {}".format(plev, nlev))
            plev = nlev
        logger.debug(rlist)
        return levlist, rlist

    def set_levels(self, cgrid):
        levlist, rlist = self.set_qva_levels()
        cgrid.set_levels(levlist)

    def additional_control_setup(self, control, stage=0):
        # for dispersion control file
        # if setting levels here then need to use the set_levels
        # function and also adjust the labeling for the levels.

        # round to nearest latitude longitude point.
        # this is for better matching with gridded volcat data.
        lat = self.inp["latitude"]
        lon = self.inp["longitude"]
        vent = self.inp["bottom"]
        height = self.inp["top"]
        emission = self.inp["emissionHours"]
        rate = self.inp["rate"]
        area = self.inp["area"]

        # add location of eruption
        # with uniform vertical line source
        control.remove_locations()
        control.add_location((lat, lon), vent, rate=rate, area=area)
        control.add_location((lat, lon), height, rate=0.0001, area=area)
        # rename cdump file
        control.concgrids[0].outdir = self.inp["WORK_DIR"]
        control.concgrids[0].outfile = self.filelocator.get_cdump_filename(stage)
        logger.info("{}".format(control.concgrids[0].outfile))
        # center concentration grid at the volcano
        # control.concgrids[0].centerlat = np.floor(lat)
        # control.concgrids[0].centerlon = np.floor(lon)
        control.concgrids[0].centerlat = int(lat)
        control.concgrids[0].centerlon = int(lon)
        # set a global grid
        control.concgrids[0].latspan = 180.0
        control.concgrids[0].lonspan = 360.0
        # set the levels
        self.set_levels(control.concgrids[0])

        # add emission duration
        rate_test = 0
        if not self.so2:
            for spec in control.species:
                spec.duration = emission
                rate_test += spec.rate
        # for so2, the sulfate emission has 0 emissions.
        else:
            spec = control.species[0]
            spec.duration = emission
            spec2 = control.species[1]
            spec2.duration = 0
            rate_test = spec.rate + spec2.rate

        if np.abs(1 - rate_test) > 0.01:
            logger.warning("Non unit emision rate {}".format(rate_test))

        # TO DO
        # set sample start to occur on the hour.
        stime = self.inp["start_date"]
        # if start past the half hour mark, set sample start at next hour.
        if stime.minute > 30:
            stime = stime + datetime.timedelta(hours=1)
        sample_start = stime.strftime("%y %m %d %H 00")
        # logger.debug("Setting sample start {}".format(sample_start))
        control.concgrids[0].sample_start = sample_start

        control.concgrids[0].interval = (self.inp["samplingIntervalHours"], 0)
        control.concgrids[0].sampletype = -1

    def compose_control(self, stage, rtype):
        control = self.setup_basic_control(stage, rtype=rtype)
        self.additional_control_setup(control, stage)
        control.write()

    def setup_setup(self, stage):
        duration = self.inp["durationOfSimulation"]
        # TO DO - may calculate particle number based on MER.
        setup = hcontrol.NameList(
            fname=self.default_setup, working_directory=self.inp["DATA_DIR"]
        )
        # read default
        setup.read(case_sensitive=False)
        print(setup.wdir, setup.fname)
        print(setup.nlist.keys())
        newname = self.filelocator.get_setup_filename(stage)
        # rename file
        setup.rename(name=newname, working_directory=self.inp["WORK_DIR"])
        # rename pardump file
        pardumpname = self.filelocator.get_pardump_filename(stage)
        setup.add("poutf", pardumpname)

        # setup.add("numpar", "2000")
        # setup.add("numpar", "20000")
        setup.add("numpar", "10000")
        # setup.add("maxpar", "500000")

        # base frequency of pardump output on length of simulation.
        if duration < 6:
            setup.add("ncycl", "1")
            setup.add("ndump", "1")
        elif duration < 12:
            setup.add("ncycl", "2")
            setup.add("ndump", "2")
        else:
            setup.add("ndump", "3")
            setup.add("ncycl", "3")

        keys = list(setup.nlist.keys())
        if "ichem" in keys:
            if int(setup.nlist["ichem"]) == 2:
                make_chemrate(self.inp["WORK_DIR"])
                logger.warning("creating chemrate file for SO2")
                self.so2 = True
        # write file
        return setup

    def compose_setup(self, stage):
        setup = self.setup_setup(stage)
        setup.write(verbose=False)

    def debug_message(self):
        # debug messages
        logger.debug("HYSPLIT_DIR     = {}".format(self.inp["HYSPLIT_DIR"]))
        # logger.debug("MAP_DIR         = {}".format(self.inp["MAP_DIR"]))
        logger.debug("WORK_DIR        = {}".format(self.inp["WORK_DIR"]))
        logger.debug("CONVERT_EXE     = {}".format(self.inp["CONVERT_EXE"]))
        logger.debug("GHOSTSCRIPT_EXE = {}".format(self.inp["GHOSTSCRIPT_EXE"]))
        logger.debug("PYTHON_EXE      = {}".format(self.inp["PYTHON_EXE"]))

    def run_model(self):
        # make control and setup files
        self.compose_control(stage=0, rtype="dispersion")
        self.compose_setup(stage=0)
        # start run and wait for it to finish..
        c = [os.path.join(self.inp["HYSPLIT_DIR"], "exec", "hycs_std"), str(self.JOBID)]
        logger.info("Running {} with job id {}".format("hycs_std", c[1]))
        Helper.execute(c)

    def doit(self):
        """
        Main work flow.
        """
        self.debug_message()
        os.chdir(self.inp["WORK_DIR"])
        if not os.path.exists("ASCDATA.CFG"):
            shutil.copyfile(
                self.inp["HYSPLIT_DIR"] + "/bdyfiles/ASCDATA.CFG", "ASCDATA.CFG"
            )
        logger.info("Please wait for further information....")
        logger.info("Model submitted on {}".format(datetime.datetime.now()))
        # if files are already there then do not run the model.
        if not self.after_run_check(update=False):
            import sys
            #sys.exit()
            # self.preprocessing()
            self.run_model()
            # self.postprocessing()
            redraw = False
        else:
            logger.info("REDRAW for run {} {}".format(self.JOBID, self.inp["WORK_DIR"]))
            redraw = True
        self.write_cxra()
        if self.after_run_check(update=True):
            self.create_plots(redraw)
        #     # create zip with cdump and pardump files etc.
        #     status = self.create_zipped_up_file(
        #         self.filelocator.get_zipped_filename(tag=""),
        #         self.filelocator.get_all_ashbase_filenames(),
        #     )
        # create zip with netcdf files for awips.
        #     logger.debug('CREATE AWIPS zipped files')
        #     status = self.create_zipped_up_file(
        #         self.filelocator.get_zipped_filename(tag="awips2_"),
        #         self.filelocator.get_awips_filenames(),
        #     )
        self.update_run_status(self.JOBID, "COMPLETED")
        self.cleanup()

    def after_run_check(self, stage=0, update=False):
        # Check for the tdump/cdump file
        rval = True
        fn = self.filelocator.get_cdump_filename(stage=stage)
        if not os.path.exists(fn):
            rval = False
            if update:
                logger.error(
                    "******************************************************************************"
                )
                logger.error(
                    "The model has crashed. Check the HYSPLIT Message file for further information."
                )
                logger.error(
                    "******************************************************************************"
                )
                self.handle_crash(stage=0)
        if not rval:
            logger.info("NOT found cdump file " + fn)
        else:
            logger.info("found cdump file " + fn)
        return rval

    def get_conc_multiplier(self):
        """
        factor to convert concentration in unit mass /m3 to mg/m3.
        Model must have been run with 1 unit mass/h released.
        """
        height = (self.inp["top"] - self.inp["bottom"]) / 1000.0
        # HT2unit outputs in grams. Then multiply by 1e3 to get mg
        conc_multiplier = 1e3 * HT2unit(height) * self.get_ash_reduction()
        return conc_multiplier

    def get_ash_reduction(self):
        eflag = float(self.inp["eflag"])
        M63 = 0.01  # fraction of  fine ash
        conc_multiplier = 10 ** (-1 * eflag)
        return M63 * conc_multiplier

    def create_plots(self, redraw=False, stage=0):
        """
        Plot creation
        """
        self.create_maptext()
        fn = self.filelocator.get_cdump_filename(stage=stage)
        self.create_massloading_plot(fn, stage)
        # DO not want maptext included in parxplot or concentration plot
        # output. only in the massloading output.
        Helper.remove(os.path.join(self.inp["WORK_DIR"], "MAPTEXT.CFG"))
        # don't need to redraw particle plots.
        if not redraw:
            self.create_parxplot(stage)
        self.create_concentration_plot(fn, stage)
        self.create_concentration_montage(stage=stage)
        # create the maptext file again for inclusion in zip file.
        self.create_maptext()
        Helper.move("MAPTEXT.CFG", self.filelocator.get_maptext_filename_for_zip())
        if self.awips:
            logger.debug("creating netcdf files for awips")
            awipsfiles = self.make_awips_netcdf()

    def cleanup(self):
        psfiles = glob.glob("*.ps")
        for psf in psfiles:
            Helper.remove(psf)
        return -1

    def create_parxplot(self, stage=0):
        # 11/13/2020 do not create kml for pardump.
        # files too large and not useful.
        fn = self.filelocator.get_pardump_filename(stage=stage)
        outputname = self.filelocator.get_parxplot_filename(ptype="ps")
        maparg = "-j" + os.path.join(self.inp["MAP_DIR"], "arlmap")
        # gisoption = self.inp['gisOption']
        gisoption = 0
        c = [
            os.path.join(self.inp["HYSPLIT_DIR"], "exec", "parxplot"),
            "-i" + fn,
            "-k1",
            "-o" + outputname,
            "-p{}".format(self.JOBID),
            "-a{}".format(gisoption),
            maparg,
        ]
        logger.debug(" ".join(c))
        Helper.execute(c)
        outputname = outputname.replace("ps", self.JOBID)
        if not os.path.exists(outputname):
            logger.debug(outputname)
            logger.error(
                "******************************************************************************"
            )
            logger.error(
                "The model was not able to create a particle plot for job {}.".format(
                    self.JOBID
                )
            )
            logger.error(
                "******************************************************************************"
            )
            self.update_run_status(self.JOBID, "PARXPLOT_FAILED")
            # logger.info(datetime.datetime.now())
        else:
            # no longer create kmz file from pardump.
            # fnames = \
            #     [self.filelocator.get_basic_kml_filename(program='parxplot')]
            # fnames = ['HYSPLITpart_ps.kml']
            # logger.debug('generating kmz from ' + ' '.join(fnames))
            # self.generate_kmz(fnames,
            #                  self.filelocator.get_kmz_filename(stage='part'))
            # create animated gif
            self.convert_ps_to_image(
                outputname,
                self.filelocator.get_parxplot_filename(ptype="gif"),
                trim=False,
                resolution=200,
            )
            # create pdf
            # self.convert_ps_to_pdf(
            #    outputname,
            #    self.filelocator.get_parxplot_filename(ptype="pdf"),
            # )

    def create_massloading_plot(self, fn, stage=0):
        # outputname = self.filelocator.get_massloading_ps_filename(ptype="ps")
        outputname = self.filelocator.get_massloading_filename(stage=stage, ptype="ps")
        # create dispersion png plot using python concplot.
        # Units in g/m2

        clrs = ConcplotColors()
        # want in grams. output in mg.
        conc_multiplier = self.get_conc_multiplier() / 1000.0
        contours = "+".join(
            [
                "50:High:{}".format(clrs.get("red")),
                "10:High:{}".format(clrs.get("magenta")),
                "5.0:Medium:{}".format(clrs.get("orange")),
                "2.0:Medium:{}".format(clrs.get("yellow")),
                "1.0:Low:{}".format(clrs.get("green")),
                "0.1::{}".format(clrs.get("tan")),
                "0.01::{}".format(clrs.get("grey")),
            ]
        )
        logger.info(
            "Creating column mass loading graphics for job {}.".format(self.JOBID)
        )
        gisopt = self.inp["gisOption"]

        mapopt = self.inp["mapBackground"]
        if mapopt == "terrain":
            maparg = "--street-map=0"
        elif mapopt == "toner":
            maparg = "--street-map=1"
        else:
            maparg = "-j" + os.path.join(self.inp["MAP_DIR"], "arlmap")

        fortran_concplot = True
        if fortran_concplot:
            # c = [self.inp['PYTHON_EXE'],
            #     os.path.join(self.inp['HYSPLIT_DIR'], 'exec', 'concplot.py'),
            c = [
                os.path.join(self.inp["HYSPLIT_DIR"], "exec", "concplot"),
                "-i" + fn,
                "-c4",
                "-v" + contours,
                "-e4",  # mass loading
                "-d4",  # over all levels
                "-ug",  # units of m2 rather than m3.
                #'+n',
                # this will change otput name of kml file.
                "-p" + str(self.JOBID),
                "-x{:2.2e}".format(conc_multiplier),
                "-o" + outputname,
                "-s0",  # sum species
                "".join(["-j", self.inp["HYSPLIT_DIR"], "/graphics/arlmap"]),
                #'-m{}'.format(self.inp['mapProjection']),
                #'-z{}'.format(self.inp['zoomFactor']),
                #'-l1',
                # maparg,
                "-a{}".format(gisopt),
                #'-g0:{}'.format(self.inp['spatialPlotRadius']),
                "-k1",
            ]
        # else
        #    c = [self.inp['PYTHON_EXE'],
        #         os.path.join(self.inp['HYSPLIT_DIR'], 'exec', 'concplot.py'),

        logger.debug("Executing concplot " + " ".join(c))
        Helper.execute(c)
        outputname = outputname.replace("ps", self.JOBID)
        if not os.path.exists(outputname):
            logger.error(
                "******************************************************************************"
            )
            logger.debug(outputname)
            logger.error(
                "The model was not able to create a graphic file for job {}.".format(
                    self.JOBID
                )
            )
            logger.error(
                "******************************************************************************"
            )
            self.update_run_status(self.JOBID, "GRAPHICS_FAILED")
            # sys.exit(3)
        else:
            # make the kmz file.
            # if the -p option used need to use the get_gelabel_filename method.
            gelist = self.make_gelabel(self.filelocator.get_gelabel_filename)
            # self.make_gelabel(self.filelocator.get_basic_gelabel_filename)
            # fnames = [self.filelocator.get_basic_kml_filename()]
            fnames = ["HYSPLIT_ps.kml"]
            fnames = [self.filelocator.get_basic_kml_filename()]
            logger.info("KML FILE NAME {}".format(fnames[0]))
            fnames.extend(gelist)
            logger.info("GELIST FILE NAME {}".format(fnames[1]))
            self.generate_kmz(
                fnames, self.filelocator.get_kmz_filename(stage="massload")
            )
            # NOT USED. convert ps to gif animation
            # convert to image and package in pdf with concentration and
            # particles.
            self.convert_ps_to_image(
                outputname,
                self.filelocator.get_massloading_filename(stage=stage, ptype="gif"),
                trim=False,
                resolution=200,
            )
            # convert ps to pdf
            # self.convert_ps_to_pdf(
            #    outputname,
            #    self.filelocator.get_massloading_filename(stage=stage, ptype="pdf"),
            # )
            # remove the kml and gelabel files
            for fn in fnames:
                Helper.remove(fn)

    def make_gelabel(self, filenamer):
        # the get_basic_gelabel_filename method should be passed if the -p
        # option was not used and files have form GELABEL_??_ps.txt.
        # otherwise the get_gelabel_filename method should be passed.

        # run gelabel executable.
        # this assumes there is a gelabel.txt file created by concplot.
        # gelabel executable creates ps files of the legend from the txt file with
        c = [
            os.path.join(self.inp["HYSPLIT_DIR"], "exec", "gelabel"),
        ]
        # check to see if -p option is needed.
        if filenamer(frame=1, ptype="txt").split(".")[0][-2:] != "ps":
            c.append("-p{}".format(self.JOBID))
        logger.debug(" ".join(c))
        Helper.execute(c)

        iii = 1
        # the ps files are then converted to gif files.
        gelabel = filenamer(frame=iii, ptype="ps")
        logger.debug(gelabel + " NAME ")
        gelist = []
        while os.path.isfile(gelabel):
            gelabelgif = filenamer(frame=iii, ptype="gif")
            gelist.append(gelabelgif)
            logger.debug(gelabel + " NAME " + gelabelgif)
            self.convert_ps_to_image(gelabel, gelabelgif, resolution=80)
            iii += 1
            gelabel = filenamer(frame=iii, ptype="ps")
            if iii > 100:
                break
        return gelist

    def create_concentration_plot(self, fn, stage):
        """
        currently generates kml file which overwrites the massloading kml file.
        but does not generate a kmz file.
        """

        # outputname = self.filelocator.get_concplot_ps_filename(ptype="ps")
        outputname = self.filelocator.get_concplot_filename(ptype="ps")
        # create dispersion png plot using python concplot.
        conc_multiplier = self.get_conc_multiplier()
        clrs = ConcplotColors()
        contours = "+".join(
            [
                "1000:High:{}".format(clrs.get("purple")),
                "100:High:{}".format(clrs.get("red")),
                "10:High:{}".format(clrs.get("magenta")),
                "5.0:Medium:{}".format(clrs.get("orange")),
                "2.0:Medium:{}".format(clrs.get("yellow")),
                "0.2:Low:{}".format(clrs.get("green")),
                "0.02:None:{}".format(clrs.get("tan")),
                "0.002:None:{}".format(clrs.get("grey")),
            ]
        )
        logger.info("Creating concentration graphics for job {}.".format(self.JOBID))
        gisopt = 3  # self.inp['gisOption']

        mapopt = self.inp["mapBackground"]
        if mapopt == "terrain":
            maparg = "--street-map=0"
        elif mapopt == "toner":
            maparg = "--street-map=1"
        else:
            maparg = "-j" + os.path.join(self.inp["MAP_DIR"], "arlmap")

        fortran_concplot = True
        if fortran_concplot:
            # c = [self.inp['PYTHON_EXE'],
            #     os.path.join(self.inp['HYSPLIT_DIR'], 'exec', 'concplot.py'),
            c = [
                os.path.join(self.inp["HYSPLIT_DIR"], "exec", "concplot"),
                "-i" + fn,
                "-c4",
                "-v" + contours,
                "-umg",  # output in mg
                #'+n',
                "-p" + str(self.JOBID),
                "-x{:2.2e}".format(conc_multiplier),
                "-o" + outputname,
                "-s0",  # sum species
                "".join(["-j", self.inp["HYSPLIT_DIR"], "/graphics/arlmap"]),
                #'-m{}'.format(self.inp['mapProjection']),
                #'-z{}'.format(self.inp['zoomFactor']),
                #'-l1',
                # maparg,
                "-a{}".format(gisopt),
                #'-g0:{}'.format(self.inp['spatialPlotRadius']),
                "-k1",
            ]
        # else
        #    c = [self.inp['PYTHON_EXE'],
        #         os.path.join(self.inp['HYSPLIT_DIR'], 'exec', 'concplot.py'),
        logger.debug("Executing concplot " + " ".join(c))
        Helper.execute(c)
        outputname = outputname.replace("ps", self.JOBID)
        if not os.path.exists(outputname):
            logger.error(
                "******************************************************************************"
            )
            logger.error(
                "The model was not able to create a concentration graphic file for job {}.".format(
                    self.JOBID
                )
            )
            logger.error(
                "******************************************************************************"
            )
            self.update_run_status(self.JOBID, "CONCENTRATION GRAPHICS_FAILED")
            logger.info(datetime.datetime.now())
            # sys.exit(3)
        else:
            # self.convert_ps_to_pdf(
            #    outputname,
            #    self.filelocator.get_concplot_filename(stage=stage, ptype="pdf"),
            # )
            self.convert_ps_to_image(
                outputname,
                self.filelocator.get_concplot_filename(stage=stage, ptype="gif"),
                resolution=200,
            )

    def create_montage_page(self, flin, flout, stage, iii, jjj):
        """
        flin : function which generates filename
        flout : function which generates filename
        jjj : frame number of output
        iii : frame number of first input file
        stage : input for flin and flout
        """
        done = False
        outputname = flout(stage, frame=jjj)
        # self.filelocator.get_concentration_montage_filename(stage,frame=jjj)
        jlist, levlist = self.set_qva_levels()
        # rlist = []
        montage = self.inp["CONVERT_EXE"].replace("convert", "montage")
        c = [montage]
        for lev in levlist:
            fname = flin(stage=stage, frame=iii)
            # self.filelocator.get_concplot_filename(stage=stage,frame=iii)
            if not os.path.exists(fname):
                done = True
                logger.debug("file does not exist {}".format(fname))
                break
            c.append("-label '{}'".format(lev))
            c.append("-pointsize 20")
            c.append("-trim {}".format(fname))
            iii += 1
        if done == True:
            return None, done, iii
        c.extend(["-geometry 200x200", "-tile 2x6", outputname])
        # logger.debug("Create montage: " + " ".join(c))
        # Sonny - not sure why I need shell=True?
        # subprocess.Popen(' '.join(c), shell=True)
        Helper.execute_with_shell(c)
        return outputname, done, iii

    def create_concentration_montage(self, stage):
        flout = self.filelocator.get_concentration_montage_filename
        flin = self.filelocator.get_concplot_filename
        file_addlist = [self.filelocator.get_massloading_filename]
        file_addlist.append(self.filelocator.get_parxplot_filename)
        stagein = [stage, stage, stage]
        self.create_montage_pdf(stagein, flin, flout, file_addlist)

    def create_montage_pdf(self, stage, flin, flout, file_addlist):
        """
        stage : list of stages.
        flin : function which gennerates file names
        flout : function which gennerates file names
        file_addlist: list of function which generate file names
        """
        # TO DO - a lot of loops and very clunky.
        # Assumes 4 levels. To change need to change 'tile'
        logger.debug("Creating pdf montage")
        done = False
        montage_list = []
        iii = 0
        jjj = 0
        while not done:
            outputname, done, iii = self.create_montage_page(
                flin, flout, stage[0], iii, jjj
            )
            # add the other files to the list for the pdf file.
            nnn = 1
            for fl in file_addlist:
                newname = fl(stage=stage[nnn], frame=jjj, ptype="gif")
                if os.path.isfile(newname):
                    montage_list.append(newname)
                else:
                    logger.debug("file not found {}".format(newname))
                nnn += 1

            if outputname:
                montage_list.append(outputname)
            if jjj >= 100:
                done = True
            jjj += 1
        if montage_list:
            c = [self.inp["CONVERT_EXE"]]
            c.extend(montage_list)
            c.append(self.filelocator.get_totalpdf_filename(stage=stage[0]))
            logger.debug("MONTAGE PDF {}".format(" ".join(c)))
            Helper.execute_with_shell(c)

    def convert_ps_to_image(self, ps_filename, gif_filename, resolution=0, trim=True):

        logger.debug("Creating images from ps {}".format(ps_filename))
        if resolution == 0:
            resolution = self.inp["graphicsResolution"]

        if not os.path.exists(ps_filename):
            logger.warn(
                "Postscript file {} does not exist. Image file {} will not be created".format(
                    ps_filename, gif_filename
                )
            )
            return
        if trim:
            c = [
                self.inp["CONVERT_EXE"],
                "+adjoin",
                "-",
                "-trim",
                "+repage",
                "GIF:{}".format(gif_filename),
            ]
        else:
            c = [
                self.inp["CONVERT_EXE"],
                "+adjoin",
                "-",
                "+repage",
                "GIF:{}".format(gif_filename),
            ]
        p1 = subprocess.Popen(
            [
                self.inp["GHOSTSCRIPT_EXE"],
                "-r{}".format(resolution),
                "-dTextAlphaBits=4",
                "-dGraphicsAlphaBits=4",
                "-dNOPAUSE",
                "-dSAFER",
                "-sDEVICE=pnmraw",
                "-q",
                "-sOutputFile=-",
                ps_filename,
                "-c",
                "quit",
            ],
            stdout=subprocess.PIPE,
            stderr=sys.stderr,
        )
        p2 = subprocess.Popen(
            c,
            stdin=p1.stdout,
            stdout=subprocess.PIPE,
            stderr=sys.stderr,
        )
        p1.stdout.close()  # allow p1 to receive a SIGPIPE if p2 exists.
        stdoutdata, stderrdata = p2.communicate()

    def convert_ps_to_pdf(self, ps_filename, pdf_filename):
        if not os.path.exists(ps_filename):
            logger.warn(
                "Postscript file {} does not exist. PDF file {} will not be created".format(
                    ps_filename, pdf_filename
                )
            )
            return

        cproc = [
            self.inp["GHOSTSCRIPT_EXE"],
            "-q",
            "-dBATCH",
            "-dSAFER",
            "-dMaxBitmap=500000000",
            "-dNOPAUSE",
            "-dAlignToPixels=0",
            "-dEPSCrop",
            "-sDEVICE=pdfwrite",
            "-sOutputFile={}".format(pdf_filename),
            "-f{}".format(ps_filename),
        ]
        logger.info("Creating PDF " + " ".join(cproc))
        Helper.execute(cproc)

    def get_maptext_info(self):
        maptexthash = {}
        rstr = "HYSPLIT dispersion run. "
        maptexthash["run_description"] = rstr

        metfiles = self.metfilefinder.find(
            self.inp["start_date"], self.inp["durationOfSimulation"]
        )
        if metfiles:
            metfiles = list(zip(*metfiles))[1]
        else:
            metfiles = "unknown"

        maptexthash["infoc"] = ",".join(metfiles)

        return maptexthash

    def create_maptext(self):
        # TO DO - write meteorology
        # vertical distribution.
        # MER and m63.
        descripline = self.maptexthash["run_description"]

        now = datetime.datetime.utcnow().strftime("%m/%d/%Y %H:%M UTC")

        emission_rate = self.get_conc_multiplier()  # gives rate in mg/hour
        emission_rate = "{:2.1e} kg/h".format(emission_rate / 1e6)

        name = self.inp["VolcanoName"]

        emission = float(self.inp["emissionHours"])
        ehours = int(emission)
        eminutes = round((emission - int(emission)) * 60)

        m63 = self.get_ash_reduction()

        lat = self.inp["latitude"]
        lon = self.inp["longitude"]
        alt = self.inp["bottom"]
        top = self.inp["top"]
        sp8 = "        "
        estart = self.inp["start_date"].strftime("%Y %m %d %H:%M UTC")

        # metid = self.inp["meteorologicalData"]
        # start_lon
        # start_lat
        # data[0] #alert type
        # vaac
        # hgts[0]
        # data[7]  #Alert detection time
        jobidline = "Job ID: {}   Job Name: {}   Job Completion: {}\n".format(
            self.JOBID, self.inp["jobname"], now
        )
        volcanoline = "Volcano: {}{} lat: {} lon: {}{}Hgt: {} to {} m\n".format(
            name, sp8, lat, lon, sp8, alt, top
        )
        poll = "Pollutant: Ash \n"
        ## TO DO calculate release quantity
        release_a = "Start: {}  Duration: {} h {} min\n".format(
            estart, ehours, eminutes
        )
        release_b = "Release: {}    m63: {:1g}\n".format(emission_rate, m63)
        info_a = "Vertical distribution: {}{}  GSD: {}{}  Particles: {}\n".format(
            "uniform", sp8, "default", sp8, "20,000"
        )
        info_b = "Meteorology: {}\n".format(self.inp["meteorologicalData"])
        info_c = self.maptexthash["infoc"]
        owner = "Run by: {}\n".format(self.inp["owner"])

        with open(self.inp["WORK_DIR"] + "MAPTEXT.CFG", "w") as fid:
            fid.write("  \n")
            fid.write("  \n")
            fid.write(jobidline)
            fid.write("  \n")
            fid.write(descripline)
            fid.write("  \n")
            fid.write("  \n")
            fid.write(volcanoline)
            fid.write(poll)
            fid.write(release_a)
            fid.write(release_b)
            fid.write(info_a)
            fid.write(info_b)
            fid.write(info_c)
            fid.write(owner)
            # current = dtime.utcnow()
            # fid.write('Job Start: '+current.strftime('%B %d, %Y')+' at '+current.strftime('%H:%M:%S')+' UTC \n')
            fid.write("\n")
            fid.close()

    # def generate_shape_file(self):
    # no shape files generated now.

    def generate_kmz(self, kml_filenames, kmz_filename):
        for f in kml_filenames:
            if not os.path.exists(f):
                logger.warn(
                    "KML file {} does not exist. Google Earth file {} will not be created".format(
                        f, kmz_filename
                    )
                )
                return

        files = [
            os.path.join(self.inp["HYSPLIT_DIR"], "guicode", "noaa_google.gif"),
            os.path.join(self.inp["HYSPLIT_DIR"], "guicode", "logocon.gif"),
            os.path.join(self.inp["HYSPLIT_DIR"], "graphics", "blueball.png"),
            os.path.join(self.inp["HYSPLIT_DIR"], "graphics", "greenball.png"),
            os.path.join(self.inp["HYSPLIT_DIR"], "graphics", "redball.png"),
        ]
        files += kml_filenames

        with zipfile.ZipFile(
            kmz_filename, "w", compresslevel=self.inp["zip_compression_level"]
        ) as z:
            for f in files:
                if os.path.exists(f):
                    bn = os.path.basename(f)
                    z.write(f, arcname=bn)

    def create_zipped_up_file(self, filename, files):
        if not files:
            return False
        logger.debug("{} files to be zipped {}".format(filename, "\n".join(files)))
        with zipfile.ZipFile(
            filename, "w", compresslevel=self.inp["zip_compression_level"]
        ) as z:
            for f in files:
                if os.path.exists(f):
                    z.write(f)
        return True

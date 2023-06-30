"""
classes
EnsembleAshRun
"""


# from abc import ABC, abstractmethod
import datetime
import logging
import os
import time

from monetio.models import hysplit

# import hysplit
from utilhysplit import metfiles

# import ashapp_plotting
from utilhysplit.evaluation import web_ensemble_plots as wep
from utilhysplit.runhandler import ProcessList

from ashbase import AshRun
from ashruninterface import ModelRunCollection, ModelRunInterface
from cdump2xml import HysplitKml
from utilvolc.runhelper import Helper

logger = logging.getLogger(__name__)


# def print_usage():
#    print(
#        """\
#    )

# test using python ash_run.py run -888
# can change test by editing runhelper module

# CHANGE LOG
# 2022 Dec 8  (AMC) changed name of ensemble_tools to ashapp_plotting
# 2023 Feb 27 (AMC) ashapp_plotting moved to web_ensemble_plots
# 2023 Feb 27 (AMC) added a 0.1 level to the ensenble relative frequency of exceedance maps
# 2023 Feb 27 (AMC) fixed the create_plots function
# 2023 Feb 27 (AMC) no longer use conprob utility so remove conprob_thresholds attribute


#class ConProbThresholds:
#    """ To be used with conprob utility """
#
#    def __init__(self, mass_converter):
#        # find probability of exceedence for the following levels.
#        self.conc_levels = [0.2, 2, 10]  # in mg
#        # this is how much 1 unit mass is in grams.
#        self.mass_converter = mass_converter
#
#    def get_unit_levels(self):
#        # this sets the levels in unit mass.
#        unit_levels = [x / self.mass_converter / 1e3 for x in self.conc_levels]
#        return unit_levels
#
#    def __str__(self):
#        unit_levels = self.get_unit_levels()
#        str_unit_levels = ["{:1.1e}".format(x) for x in unit_levels]
#        return ":".join(str_unit_levels)

#    def __repr__(self):
#        unit_levels = self.get_unit_levels()
#        str_unit_levels = ["{:1.1e}".format(x) for x in unit_levels]
#        return ":".join(str_unit_levels)

class EnsembleAshRun(ModelRunCollection):

    


    def run_model(self):
        processhandler = ProcessList()
        # redirect stdout and stderr
        processhandler.pipe_stdout()
        processhandler.pipe_stderr()
        # start all the runs
        for stage, suffix in enumerate(self.ens_suffix_list):
            logger.info("Working on {}".format(suffix))
            self.metfilefinder.set_ens_member("." + suffix)
            self.compose_control(stage + 1, rtype="dispersion")
            self.compose_setup(stage + 1)
            # start run and wait for it to finish..
            run_suffix = self.filelocator.get_control_suffix(stage + 1)
            cproc = [
                os.path.join(self.inp["HYSPLIT_DIR"], "exec", "hycs_std"),
                run_suffix,
            ]
            logger.info("Running {} with job id {}".format("hycs_std", cproc[1]))
            processhandler.startnew(cproc, self.inp["WORK_DIR"], descrip=suffix)
            # wait 5 seconds between run starts to avoid
            # runs trying to access ASCDATA.CFG at the same time.
            time.sleep(5)
        # wait for runs to finish
        done = False
        seconds_to_wait = 30
        total_time = 0
        # 60 minutes.
        max_time = 60 * 60
        # max_time = 0.5*60
        while not done:
            num_proces = processhandler.checkprocs()
            if num_proces == 0:
                done = True
            time.sleep(seconds_to_wait)
            total_time += seconds_to_wait
            if total_time > max_time:
                processhandler.checkprocs()
                processhandler.killall()
                logger.warning("HYSPLIT run Timed out")
                done = True

class EnsembleAshRunOLD(AshRun):
    def __init__(self, JOBID):
        super().__init__(JOBID)
        self.ens_suffix_list = None
        #self.conprob_thresholds = None
        self.number_of_members = 0
        self.awips = True

    def plot_massload(self):
        """
        plots ensemble mean
        """
        if self.cxra.size <= 1:
            logger.info("plot_massload cxra is empty")
            return False
        enslist = self.cxra.ens.values
        # level = self.cxra.z.values
        vlist = [self.inp["longitude"], self.inp["latitude"]]
        flin = self.filelocator.get_massloading_filename(
            stage=0, frame=999, ptype="png"
        )
        flin = flin.replace("999", "zzz")
        # flin = flin.replace('gif','pdf')
        logger.debug("Massloading FILENAME {}".format(flin))
        # massload_plot plots the ensemble mean.
        fignamelist = wep.massload_plot(self.cxra, enslist, vlist=vlist, name=flin)
        # list of figure names generated.
        return fignamelist

    def make_kml(self):
        if self.cxra.size <= 1:
            logger.info("make_kml WARNING: cxra is empty")
            return False
        # kname = 'test_{}.kml'.format(self.JOBID)
        kname = self.filelocator.get_kmz_filename(stage="massload")
        kname = kname.replace("kmz", "kml")
        mxra = hysplit.hysp_massload(self.cxra)
        # convert to g/m2
        mxra = mxra.isel(source=0).mean(dim="ens") / 1000.0
        attrs = self.cxra.attrs
        levels = wep.set_levels(mxra)
        h2xml = HysplitKml(
            levels=levels,
            sourcehash=self.inp,
            units="g/m2",
            jobid=self.JOBID,
            legend_label="Ensemble Mean Mass Loading (GEFS)",
        )
        logger.debug("creating KML file {}".format(kname))
        h2xml.create(mxra, attrs, kname)
        efiles = [
            os.path.join(self.inp["HYSPLIT_DIR"], "guicode", "noaa_google.gif"),
            os.path.join(self.inp["HYSPLIT_DIR"], "guicode", "logocon.gif"),
            os.path.join(self.inp["HYSPLIT_DIR"], "graphics", "blueball.png"),
            os.path.join(self.inp["HYSPLIT_DIR"], "graphics", "greenball.png"),
            os.path.join(self.inp["HYSPLIT_DIR"], "graphics", "redball.png"),
        ]
        comp = self.inp["zip_compression_level"]
        h2xml.create_kmz(comp, efiles=efiles)
        return True

    def plot_ATL(self):
        """
        plots probability of exceedances.
        """
        # Make sure that ATL plots are not all empty.
        # if maximum value below threshold then adjust
        # threshold so it is 1/10 the max value.
        # some time periods may be empty as the adjustment is
        # applied to the entire array.
        adjust = 10

        # must be completed after self.cxra is filled.
        if self.cxra.size <= 1:
            logger.info("plot_ATL cxra is empty")
            return False
        enslist = self.cxra.ens.values
        thresh = 0.2
        # location of volcano
        vlist = [self.inp["longitude"], self.inp["latitude"]]
        flin = self.filelocator.get_concentration_montage_filename(stage=0, frame=999)
        flin = flin.replace("999", "zzz")
        flin = flin.replace("gif", "png")
        logger.debug("NEW FILENAME{}".format(flin))
        clevels = [0.1, 5, 20, 40, 60, 80, 95]
        title = "HYSPLIT ensemble relative frequency exceeding {:0.2f}mg/m3".format(
            thresh
        )
        title += "\n GEFS {} members".format(len(enslist))

        fignamelist = wep.ATLtimeloop(
            self.cxra,
            enslist,
            thresh,
            vlist,
            name=flin,
            norm=True,
            clevels=clevels,
            title=title,
            adjust=adjust,
        )
        return fignamelist

    def get_cdump_xra(self):
        blist = []
        if self.cxra.ndim == 0:
            # ------------------------------------------------------------------------
            def make_tuple(inval):
                source_tag = "Line to {:1.0f} km".format(self.inp["top"] / 1000.0)
                suffix = inval[1]
                iii = inval[0] + 1
                cdumpname = "{}.{:03d}".format(
                    self.filelocator.get_cdump_base(stage=iii), iii
                )
                met_tag = suffix
                logger.info("adding to netcdf file :{} {}".format(met_tag, cdumpname))
                return (cdumpname, source_tag, met_tag)

            # ------------------------------------------------------------------------

            blist = [make_tuple(x) for x in enumerate(self.ens_suffix_list)]
            century = 100 * (int(self.inp["start_date"].year / 100))
            cdumpxra = hysplit.combine_dataset(blist, century=century)
            if cdumpxra.size <= 1:
                logger.debug("ENSEMBLE xra is empty")
            else:
                logger.debug("ENSEMBLE xra is full")
            self.cxra = cdumpxra
        return self.cxra

    def add_inputs(self, inp):
        logger.info("adding ensemble inputs")
        super().add_inputs(inp)
        if inp["meteorologicalData"].lower() == "gefs":
            self.ens_suffix_list = metfiles.gefs_suffix_list()
        self.number_of_members = len(self.ens_suffix_list)
        #self.conprob_thresholds = ConProbThresholds(self.get_conc_multiplier())
        self.maptexthash = self.get_maptext_info()

    def get_maptext_info(self):
        maptexthash = {}
        rstr = "HYSPLIT ensemble mean."
        maptexthash["run_description"] = rstr
        maptexthash["infoc"] = ""
        return maptexthash

    def cleanup(self):
        stage = 1
        # for suffix in self.ens_suffix_generator:
        #    run_suffix = self.filelocator.get_control_suffix(stage)

    def run_model(self):
        processhandler = ProcessList()
        # redirect stdout and stderr
        processhandler.pipe_stdout()
        processhandler.pipe_stderr()
        # start all the runs
        for stage, suffix in enumerate(self.ens_suffix_list):
            logger.info("Working on {}".format(suffix))
            self.metfilefinder.set_ens_member("." + suffix)
            self.compose_control(stage + 1, rtype="dispersion")
            self.compose_setup(stage + 1)
            # start run and wait for it to finish..
            run_suffix = self.filelocator.get_control_suffix(stage + 1)
            cproc = [
                os.path.join(self.inp["HYSPLIT_DIR"], "exec", "hycs_std"),
                run_suffix,
            ]
            logger.info("Running {} with job id {}".format("hycs_std", cproc[1]))
            processhandler.startnew(cproc, self.inp["WORK_DIR"], descrip=suffix)
            # wait 5 seconds between run starts to avoid
            # runs trying to access ASCDATA.CFG at the same time.
            time.sleep(5)
        # wait for runs to finish
        done = False
        seconds_to_wait = 30
        total_time = 0
        # 60 minutes.
        max_time = 60 * 60
        # max_time = 0.5*60
        while not done:
            num_proces = processhandler.checkprocs()
            if num_proces == 0:
                done = True
            time.sleep(seconds_to_wait)
            total_time += seconds_to_wait
            if total_time > max_time:
                processhandler.checkprocs()
                processhandler.killall()
                logger.warning("HYSPLIT run Timed out")
                done = True

    def file_not_found_error(self, fln, message=False):
        if not os.path.exists(fln):
            if message:
                logger.error(
                    "******************************************************************************"
                )
                logger.error(
                    "The model was not able to create ensemble {} file \
                              for job {}.".format(
                        fln, self.JOBID
                    )
                )
                logger.error(
                    "******************************************************************************"
                )
            rval = False
        else:
            logger.debug("cdump file found {}".format(fln))
            rval = True
        return rval

    def after_run_check(self, update=True):
        """
        Returns
        True if all cdump files are found.
        False if any of the cdump files are not found.

        If update true then will  update run status to FAILED if not all files are found.
        """
        rval = True
        # fnlist = []
        # for stage in range(1,len(self.ens_suffix_list)+1):
        #    fnlist.append(self.filelocator.get_cdump_filename(stage=stage))
        fnlist = [
            self.filelocator.get_cdump_filename(stage=x)
            for x in range(1, len(self.ens_suffix_list) + 1)
        ]
        rlist = [self.file_not_found_error(fn, update) for fn in fnlist]
        if not all(rlist):
            rval = False
            if update:
                self.update_run_status(self.JOBID, "HYSPLIT FAILED")
                logger.info(datetime.datetime.now())
        return rval

    # TO DO fix this
    def create_plots(self, redraw=False, stage=0):
        # ensemble mean of massloading to be on main display

        # make the awips2 files.
        # this will also load the data from cdump files into an xarray.
        # logger.info("creating netcdf files for awips")
        # awipsfiles = self.make_awips_netcdf()
        # return -1
        # create probability of exceedence plot.
        logger.info("creating prob of exceedance plot")
        # self.create_prob_plot()
        # logger.debug("RUNNING plot_ATL")
        ATL_fignamelist = self.plot_ATL()
        # creates kml and kmz  using cdump2xml module.

        logger.info("creating kml files")
        self.make_kml()
        # create massloading plots.
        logger.info("creating mass loading plots")
        mass_fignamelist = self.plot_massload()
        return -1
        # create parxplot for one ensemble member
        # stage would give ensemble member to use.
        self.maptexthash[
            "run_description"
        ] = "Particle Positions for 1 ensemble\
                                          member"
        self.create_maptext()
        # particle plots do not need to be re-drawn when unit mass changed.
        # if not redraw:
        #   self.create_parxplot(stage=1)

        flist = []
        iii = 0
        if len(mass_fignamelist) == len(ATL_fignamelist):
            for val in zip(mass_fignamelist, ATL_fignamelist):
                parfilename = self.filelocator.get_parxplot_filename(
                    stage=0, frame=iii, ptype="gif"
                )
                for fn in [val[0], val[1], parfilename]:
                    if os.path.exists(fn):
                        flist.append(fn)
                    else:
                        logger.warn("file {} not found".format(fn))
                iii += 1
        else:
            logger.warning(
                "Mass loading figures and ATL figures have different lengths"
            )
            flist = mass_fignamelist
        self.create_montage_pdf(flist)

        # NO longer create ensemble mean concentrations.
        # instead give probability of exceedances.
        # create montage of ensemble mean concentration.
        # Helper.remove('MAPTEXT.CFG')
        # self.create_concentration_plot(fn,stage='cmean')
        # this has been taken over by the plot_ATL function.
        # self.create_ensemble_montage()

        self.maptexthash["run_description"] = "Ensemble Run"
        self.create_maptext()
        Helper.move("MAPTEXT.CFG", self.filelocator.get_maptext_filename_for_zip())

    def create_montage_pdf(self, montage_list):
        c = [self.inp["CONVERT_EXE"]]
        c.extend(montage_list)
        c.append(self.filelocator.get_totalpdf_filename(stage=0))
        logger.info("Creating montage pdf {}".format(" ".join(c)))
        Helper.execute_with_shell(c)

    def set_qva_levels(self):
        levlist, rlist = super().set_qva_levels()
        rlist = ["Number of members\n above 0.2 mg/m3 \n" + x for x in rlist]
        return levlist, rlist

# -----------------------------------------------------------------------------
import datetime
import logging
import os
import time

import metfiles as metfile
import hysplit
from runhelper import Helper
from runhandler import ProcessList
from ashbase import AshRun
import ensemble_tools
from cdump2xml import HysplitKml

logger = logging.getLogger(__name__)

def print_usage():
    print(
        """\
USAGE: ashinverse.py JOBID

The following environment variables must be set prior to calling this script:
    RUN_API_KEY         - secret key to access Locusts web APIs.
    RUN_URL             - URL to the Locusts web application."""
    )

# Base class is AshRun

#TODO
# output needs to be more frequent than every 3 hours.
# output levels should probably be in km rather than FL.
# bottom elevation should be rounded to nearest km.

def inverse_get_suffix_list(inp):
    """
    inp: dictionay generated in

    outputs
    ---------
    suffixhash : dictionary
    key is the suffix. value is a dictionary with information about the run.
    including start date, bottom and top elevation.
    """
    suffixhash = {}
    vres = inp['inv_vertical_resolution']
    dt = datetime.timedelta(hours=inp['timeres'])
    sdate = inp['start_date']
    #edate = inp['start_date'] + datetime.timedelta(hours=inp["durationOfSimulation"])
    edate = inp['start_date'] + datetime.timedelta(hours=inp["emissionHours"])
    ndate = sdate
    done=False
    iii=0
    nnn=0
    # time loop
    while not done:
          # vertical resolution loop.
          vdone = False
          jjj=0
          bottom = inp['bottom']
          while not vdone:
              inhash = {}
              inhash['sdate'] = ndate
              inhash['edate'] = ndate+dt
              inhash['bottom'] = bottom
              inhash['top'] = bottom + vres
              suffixhash['{:03d}'.format(nnn)] = inhash.copy()
              bottom += vres
              if bottom > inp['top']: vdone=True
              jjj+=1
              nnn+=1
              print(nnn, jjj, 'bottom, top', bottom, bottom+vres)
              if jjj>50: 
                return suffixhash
          print(iii, 'dates', ndate, edate)
          iii+=1
          ndate = ndate + dt
          if ndate >= edate: done=True
    return suffixhash 

class InverseAshRun(AshRun):
    def __init__(self, JOBID):
        super().__init__(JOBID)
        self.invs_suffix_hash = {}
        self.number_of_members = 0
        self.awips = True

    def plot_massload(self):
        if self.cxra.size <= 1:
            logger.info("plot_massload cxra is empty")
            return False
        enslist = self.cxra.ens.values
        #level = self.cxra.z.values
        vlist = [self.inp["longitude"], self.inp["latitude"]]
        flin = self.filelocator.get_massloading_filename(
            stage=0, frame=999, ptype="png"
        )
        flin = flin.replace("999", "zzz")
        # flin = flin.replace('gif','pdf')
        logger.debug("Massloading FILENAME{}".format(flin))
        fignamelist = ensemble_tools.massload_plot(self.cxra, enslist, vlist=vlist, name=flin)
        # list of figure names generated.
        return fignamelist

    def make_kml(self):
        return -1

    def plot_ATL(self):
        """
        plots probability of exceedances.
        """
        return -1

    def additional_control_setup(self, control, stage=0):
        super().additional_control_setup(control,stage=stage)
        # one hour average output every hour.
        control.concgrids[0].sampletype = 1
        control.concgrids[0].interval = (1,0)

    def get_cdump_xra(self):
        blist = []

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

        suffix_list = self.invs_suffix_hash.keys()
        blist = [make_tuple(x) for x in enumerate(suffix_list)]
        century = 100 * (int(self.inp["start_date"].year / 100))
        cdumpxra = hysplit.combine_dataset(blist, century=century)
        if cdumpxra.size <= 1:
            logger.debug("ENSEMBLE xra is empty")
        else:
            logger.debug("ENSEMBLE xra is full")
        return cdumpxra

    def add_inputs(self, inp):
        logger.info("adding inverse inputs")
        self.invs_suffix_hash = inverse_get_suffix_list(inp)
        #for shash in self.invs_suffix_hash.keys():
        #    print(shash, self.invs_suffix_hash[shash])
        super().add_inputs(inp)
        #if inp["meteorologicalData"].lower() == "gefs":
        # need to get list of suffix for the inverse modeling.
        self.number_of_members = len(self.invs_suffix_hash)
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
        stage = 1
        processhandler = ProcessList()
        # redirect stdout and stderr
        processhandler.pipe_stdout()
        processhandler.pipe_stderr()
        # create a run for each suffix in this list.
        for suffix in self.invs_suffix_hash.keys():
            logger.debug("Working on {}".format(suffix))
            # not using the ensemble meteorology at this time.
            #self.metfilefinder.set_ens_member("." + suffix)
            shash = self.invs_suffix_hash[suffix]
            # Changes for each inverse modeling run.
            self.inp['bottom'] = shash['bottom']
            self.inp['top'] = shash['top']
            self.inp['emissionHours'] = self.inp['timeres']
            self.inp['start_date'] = shash['sdate']

            self.compose_control(stage, rtype="dispersion")
            self.compose_setup(stage)
            # start run and wait for it to finish..
            run_suffix = self.filelocator.get_control_suffix(stage)
            cproc = [os.path.join(self.inp["HYSPLIT_DIR"], "exec", "hycs_std"), run_suffix]
            logger.info("Running {} with job id {}".format("hycs_std", cproc[1]))
            processhandler.startnew(cproc, self.inp["WORK_DIR"], descrip=suffix)
            stage += 1
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
            rval= True
        return rval

    def after_run_check(self, update=True):
        """
        Returns
        True if all cdump files are found.
        False if any of the cdump files are not found.

        If update true then will  update run status to FAILED if not all files are found.
        """
        rval = True
        #fnlist = []
        #for stage in range(1,len(self.invs_suffix_list)+1):
        #    fnlist.append(self.filelocator.get_cdump_filename(stage=stage))
        suffix_list = list(self.invs_suffix_hash.keys())
        fnlist = [self.filelocator.get_cdump_filename(stage=x) for  x in\
                                         range(1,len(suffix_list)+1)]
        rlist = [self.file_not_found_error(fn, update) for fn in fnlist]
        if not all(rlist):
            rval = False
            if update:
                self.update_run_status(self.JOBID, "HYSPLIT FAILED")
                logger.info(datetime.datetime.now())
        return rval

    def create_plots(self, redraw=False, stage=0):
        # ensemble mean of massloading to be on main display

        # make the awips2 files.
        # this will also load the data from cdump files into an xarray.
        logger.debug("creating netcdf files")
        awipsfiles = self.make_awips_netcdf()

        # create probability of exceedence plot.
        # self.create_prob_plot()
        logger.debug("RUNNING plot_ATL")
        ATL_fignamelist = self.plot_ATL()
        # creates kml and kmz  using cdump2xml module.
        self.make_kml()
        # create massloading plots.
        mass_fignamelist = self.plot_massload()

        # create parxplot for one ensemble member
        # stage would give ensemble member to use.
        self.maptexthash[
            "run_description"
        ] = "Particle Positions for 1 ensemble\
                                          member"
        self.create_maptext()
        # particle plots do not need to be re-drawn when unit mass changed.
        if not redraw: self.create_parxplot(stage=1)
 
        flist = []
        iii=0
        if len(mass_fignamelist) == len(ATL_fignamelist):
            for val in zip(mass_fignamelist, ATL_fignamelist): 
                parfilename = self.filelocator.get_parxplot_filename(stage=0, frame=iii,ptype='gif')
                for fn in [val[0], val[1], parfilename]:
                    if os.path.exists(fn):
                        flist.append(fn)
                    else:
                        logger.warn('file {} not found'.format(fn))
                iii+=1
        else:
            logger.warning('Mass loading figures and ATL figures have different lengths')
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
        

    def create_montage_pdf(self,montage_list):
        c = [self.inp["CONVERT_EXE"]]
        c.extend(montage_list)
        c.append(self.filelocator.get_totalpdf_filename(stage=0))
        logger.info('Creating montage pdf {}'.format(" ".join(c))) 
        Helper.execute_with_shell(c)

    def set_levels_A(self):
        levlist, rlist = super().set_levels_A()
        rlist = ["Number of members\n above 0.2 mg/m3 \n" + x for x in rlist]
        return levlist, rlist

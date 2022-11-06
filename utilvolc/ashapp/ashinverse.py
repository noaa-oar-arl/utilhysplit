# -----------------------------------------------------------------------------
import datetime
import logging
import os
import time
import sys

import utilvolc.ashapp.metfiles as metfile
from monetio.models import hysplit
from utilvolc.ashapp.runhelper import Helper
from utilvolc.ashapp.runhandler import ProcessList
from utilvolc.ashapp.ashbase import AshRun
from utilvolc.ashapp import ensemble_tools
from utilvolc.ashapp.cdump2xml import HysplitKml

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

# -------LEVELS
# output levels should probably be in km rather than FL.
# If these are ONLY used for inverse modeling with the mass loading.
# we may not need so many levels.
# If height information is used or these are also used for forecasting then
# height levels are needed.
# TO DO bottom elevation should be rounded to nearest km.

# Should suffix be 
# vid_monthday_ht

# For an active volcano, unit mass runs can be done in an advance and then used
# when there is a satellite retrieval.
# This can end up being quite a bit of data though, esp. if use GEFS.
# If GMM then could only store Gaussians or only pardump output.


def get_suffix(suffix_type,dtfmt,nnn,ndate,bottom,center=None):
    if suffix_type=='int':
          suffix = '{:03d}'.format(nnn)
    elif suffix_type=='date':
          str1 = ndate.strftime(dtfmt)
          str2 = str(int(bottom))
          suffix = '{}_{}'.format(str1,str2)
    else:
          suffix = '{:03d}'.format(nnn)
    if center:
          latstr = '{:.3f}'.format(center[0]).replace('.','p')
          lonstr = '{:.3f}'.format(center[0]).replace('.','p')
          suffix += '_' + latstr + '_' + lonstr
    return suffix


def inverse_get_center_list(inp):
    import numpy as np
    # center point (lat, lon)
    center = (inp['latitude'],inp['longitude'])
    # area to cover. (meters squared)
    area = inp['area']
    # number of points to each side of center.
    num = inp['hnum']
    if num == 0:
       return center
    # here should use a square area.
    centerlist = []
    dlat = (area ** 0.5) / 111.0e3
    dlon = dlat * np.cos(center[0]*np.pi/180.0)
    lat0 = center[0]-num*dlat
    latm = center[0]+(num+0.75)*dlat
    latlist = np.arange(lat0,latm,dlat)  
    lon0 = center[1]-num*dlon
    lonm = center[1]+(num+0.75)*dlon
    lonlist = np.arange(lon0,lonm,dlon) 
    latlon = np.meshgrid(latlist, lonlist) 
    latlon = zip(latlon[0].flatten(),latlon[1].flatten())
    return list(latlon)

def inverse_get_suffix_list(inp, suffix_type='date',dtfmt="%m%d%H"):
    """
    inp: dictionay generated in
    suffix_type : str
          'int' suffix will be an integer 001, 002, 003....
          'date' suffix will be of form date_height. with date
                 determined by dtfmt.
    dtfmt : str. determines format of date when suffix_type is 'date'

    outputs
    ---------
    suffixhash : dictionary
    key is the suffix. value is a dictionary with information about the run.
    including start date, bottom and top elevation.
    """
    suffixhash = {}
    vres = inp['inv_vertical_resolution']
    dt = datetime.timedelta(hours=inp['timeres'])
    if dt.seconds%3600 > 1e-8:
       dtfmt = "%m%d%Hh%M" 
    sdate = inp['start_date']
    #edate = inp['start_date'] + datetime.timedelta(hours=inp["durationOfSimulation"])
    edate = inp['start_date'] + datetime.timedelta(hours=inp["emissionHours"])
    ndate = sdate
    done=False
    nnn=0
    if 'hnum' in inp.keys():
        latlonlist = inverse_get_center_list(inp)
    else:
        latlonlist = [1]
    # latlon loop.
    for center in latlonlist:
        iii=0
    # time loop
        if len(latlonlist) > 1: center_suffix=center
        else: center_suffix=None

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
                  suffix = get_suffix(suffix_type,dtfmt,nnn,ndate,bottom,center_suffix)
                  suffixhash[suffix] = inhash.copy()
                  bottom += vres
                  if bottom > inp['top']: vdone=True
                  jjj+=1
                  nnn+=1
                  # limit of no more than 50 heights.
                  if jjj>500: 
                     return suffixhash
              iii+=1
              ndate = ndate + dt
              if ndate >= edate: done=True
    return suffixhash 

class InverseAshRun(AshRun):
    """
    EmissionHours for the inverse run is set equal to timeres.
    """


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

    def get_conc_multiplier(self):
        """
        factor to convert concentration in unit mass /m3 to mg/m3.
        Model must have been run with 1 unit mass/h released.
        for inverse do not mutiply.
        """
        return 1

    def additional_control_setup(self, control, stage=0):
        super().additional_control_setup(control,stage=stage)
        # one hour average output every hour.
        control.concgrids[0].sampletype = -1
        control.concgrids[0].interval = (1,0)

    def get_cdump_xra(self):
        blist = []

        def make_tuple(inval):
            source_tag = "Line to {:1.0f} km".format(self.inp["top"] / 1000.0)
            suffix = inval[1]
            iii = inval[0] + 1
            #iii = inval[0] 
            cdumpname = self.filelocator.get_cdump_filename(stage=suffix)
            #cdumpname = "{}.{:03d}".format(
                #self.filelocator.get_cdump_base(stage=sxuffix), iii
            #    self.filelocator.get_cdump_base(stage=suffix), iii
            #)
            met_tag = suffix
            logger.info("adding to netcdf file :{} {}".format(met_tag, cdumpname))
            return (cdumpname, source_tag, met_tag)

        suffix_list = self.invs_suffix_hash.keys()
        blist = [make_tuple(x) for x in enumerate(suffix_list)]
        century = 100 * (int(self.inp["start_date"].year / 100))
        cdumpxra = hysplit.combine_dataset(blist, century=century,sample_time_stamp='start')
        print('attrs', cdumpxra.attrs['time description'])
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
        if inp["meteorologicalData"].lower() == "gefs":
           logger.info("ens member {}".format(inp['gefsmember']))
           self.metfilefinder.set_ens_member('.' + inp['gefsmember'])
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
        import sys
        Helper.execute('pwd')
        maxprocess=40
        #stage = 1
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
            # EmissionHours is set to timeres of the inverse modeling system.
            self.inp['emissionHours'] = self.inp['timeres']
            self.inp['start_date'] = shash['sdate']

            self.compose_control(suffix, rtype="dispersion")
            self.compose_setup(suffix)
            # start run and wait for it to finish..
            run_suffix = self.filelocator.get_control_suffix(suffix)
            cproc = [os.path.join(self.inp["HYSPLIT_DIR"], "exec", "hycs_std"), run_suffix]
            logger.info("Running {} with job id {}".format("hycs_std", cproc[1]))
            num_proces = processhandler.checkprocs()
            print('----')
            Helper.execute('pwd')
            print('----')
            print(self.inp["WORK_DIR"])
            Helper.execute(['ls', 'CONTROL.{}'.format(cproc[1])])
            while num_proces > maxprocess:
               time.sleep(10)
               num_proces = processhandler.checkprocs()
               logger.info('Waiting for process to finish before starting new {}'.format(num_proces))  
            processhandler.startnew(cproc, wdir=self.inp["WORK_DIR"], descrip=suffix)
            #stage += 1
            # wait 2 seconds between run starts to avoid
            # runs trying to access ASCDATA.CFG at the same time.
            time.sleep(2)
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
        #fnlist = [self.filelocator.get_cdump_filename(stage=x) for  x in\
        #                                 range(1,len(suffix_list)+1)]
        fnlist = [self.filelocator.get_cdump_filename(stage=x) for  x in\
                                         suffix_list]
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
        self.make_netcdf()

    def make_netcdf(self):
        import cdump2netcdf
        # convert to mg /m3
        fname = "xrfile.{}.nc".format(self.JOBID)
        mult = 1
        if os.path.isfile(fname):
            logger.info("netcdf file exists.".format(fname))
            #cxra = xr.open_dataset(fname)
            #cxra = cxra.__xarray_dataarray_variable__
            #pmult = cxra.attrs["mult"]
            #cxra = mult * cxra / pmult
            # remove netcdf file so new one can be written.
            #Helper.remove(os.path.join(self.inp["WORK_DIR"], fname))
        else:
            logger.info("netcdf file does not exist. Creating {}".format(fname))
            cxra = self.get_cdump_xra()
            cxra = cxra.assign_attrs({"mult": mult})
            logger.info("writing nc file {}".format(fname))
            cxra = cxra.assign_attrs(self.inp2attr())
            logger.debug(cxra.attrs['time description'])
            cxra.to_netcdf(fname,encoding={'zlib':True,'complevel':9})
            self.cxra = cxra

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

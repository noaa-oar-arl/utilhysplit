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
from utilvolc.runhelper import ConcplotColors, Helper, JobFileNameComposer
from ashapp.ashruninterface import MainRunInterface
#from ashapp import  utils
#from ashapp.ashnetcdf import HYSPLITAshNetcdf

#from utilvolc.volcMER import HT2unit

logger = logging.getLogger(__name__)



class MainDispersion(MainRunInterface):

    def __init__(self, inp, JOBID):

        self.JOBID = str(JOBID)  # string

        self.ilist = ['MAP_DIR','WORK_DIR','CONVERT_EXE',
                      'GHOSTSCRIPT_EXE','PYTHON_EXE']

        self.inp = inp  # dictionary from JobSetUP
        self.modelrun = inp
        self.apistr = None
        self.urlstr = None
        self.headerstr = None

        self.filelocator = None
        self.maptexthash = {}
        self.awips = True
       
        self._modeloutput = DispersionOutput(inp,[]) 
        self._modelgraphics = DispersionOutput(inp,[]) 

        # read cdump file into this xarray.
        #self._cxra = xr.DataArray()

        self.so2 = False
        utils.setup_logger()

    @property
    def inp(self):
        return self._inp

    @inp.setter
    def inp(self,inp):
        complete = True
        ilist = ['MAP_DIR','WORK_DIR','CONVERT_EXE','GHOSTSCRIPT_EXE','PYTHON_EXE']
        for iii in self.ilist:
            if iii not in inp.keys(): 
               logger.warning('Input does not contain {}'.format(iii))
               complete=False
        if complete: logger.info('Input contains all fields')
        self._inp = inp

    @property
    def modelrun(self):
        return self._modelrun

    @modelrun.setter
    def modelrun(self, inp):
        inp['jobid'] = self.JOBID
        self._modelrun = DispersionRun(inp)  

    @property
    def modeloutput(self):
        return self._modeloutput

    @property
    def modelgraphics(self):
        return self._modelgraphics


    def add_api_info(self, apistr, urlstr, headerstr):
        self.apistr = apistr
        self.urlstr = urlstr
        self.headerstr = headerstr

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

    def debug_message(self):
        # debug messages
        logger.debug("HYSPLIT_DIR     = {}".format(self.inp["HYSPLIT_DIR"]))
        logger.debug("MAP_DIR         = {}".format(self.inp["MAP_DIR"]))
        logger.debug("WORK_DIR        = {}".format(self.inp["WORK_DIR"]))
        logger.debug("CONVERT_EXE     = {}".format(self.inp["CONVERT_EXE"]))
        logger.debug("GHOSTSCRIPT_EXE = {}".format(self.inp["GHOSTSCRIPT_EXE"]))
        logger.debug("PYTHON_EXE      = {}".format(self.inp["PYTHON_EXE"]))

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
        #if not self.after_run_check(update=False):
            # self.preprocessing()
        self.modelrun.run_model()

        self.modeloutput.inputlist = self.modelrun.filelist
        self.modeloutput.postprocess()

        if self.modeloutput.check():
           self.modelgraphics.inputlist = self.modeloutput.outputlist
           self.modelgraphics.postprocess()

        self.update_run_status(self.JOBID, "COMPLETED")
        self.cleanup()

    def after_run_check(self, stage=0, update=False):
        # Check for the tdump/cdump file
        rval = True
        fnlist = [x for x in self.modelrun.filelist if 'cdump' in x]
        for fn in fnlist: 
            if not os.path.exists(fn):
                rval = False
                logger.info("NOT found cdump file " + fn)
            else:
                logger.info("found cdump file " + fn[0])
         
        if update and not rval:
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

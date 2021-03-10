#!/opt/Tools/anaconda3/envs/hysplit/bin/python
# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# runhelper.py -
#
# Helper class contains functions for executing commands.
# JobSetUp class setups the dictionary which contains information for ash runs.
# Job and JobFileNameComposer class create filenames.
#
# 18 APR 2020 (SYZ) - Initial.
# 15 Jun 2020 (AMC) - Adapted from locusts.py
# -----------------------------------------------------------------------------

# from abc import ABC, abstractmethod
import datetime
import glob
import logging
import math
import os
import pathlib

# import pytz
import shutil
import subprocess
import sys

from utilhysplit.hcontrol import NameList


logger = logging.getLogger(__name__)


class ConcplotColors:
    def __init__(self):
        colorhash = {}
        colorhash["yellow"] = "242236051"
        colorhash["orange"] = "235137052"
        colorhash["red"] = "161024014"
        colorhash["blue"] = "070051242"
        colorhash["green"] = "147219121"
        colorhash["magenta"] = "194056143"
        colorhash["purple"] = "107023156"
        colorhash["cyan"] = "075201199"
        colorhash["grey"] = "150150150"
        colorhash["tan"] = "163145131"
        self.colorhash = colorhash

    def get(self, color):
        return self.colorhash[color]


class Helper:
    def execute_with_shell(cmd, **kwargs):
        """
        cmd : string
        """
        p = subprocess.Popen(
            " ".join(cmd), shell=True, stdout=sys.stdout, stderr=sys.stderr
        )
        stdoutdata, stderrdata = p.communicate()
        if stdoutdata is not None:
            logger.info(stdoutdata)
        if stderrdata is not None:
            logger.error(stderrdata)

    def execute(cmd, **kwargs):
        """
        cmd : string
        """
        p = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)
        stdoutdata, stderrdata = p.communicate()
        if stdoutdata is not None:
            logger.info(stdoutdata)
        if stderrdata is not None:
            logger.error(stderrdata)

    def remove(f):
        """
        f : list of strings or string.
        remove file or files in list. 
        """
        if isinstance(f, list):
            for g in f:
                if os.path.exists(g):
                    os.remove(g)
        else:
            if os.path.exists(f):
                os.remove(f)

    def move(a, b):
        """
        a : string.
        b : string.
        move file a to b. 
        """
        if os.path.exists(a):
            shutil.move(a, b)

    def copy(a, b):
        """
        a : string.
        b : string.
        move file a to b. 
        """
        if os.path.exists(a):
            shutil.copy(a, b)

    def move_or_create(a, b):
        if os.path.exists(a):
            shutil.move(a, b)
        else:
            pathlib.Path(b).touch()

def make_inputs_from_file(wdir, config_file='ash_config.txt'):
    jobsetup = JobSetUp()
    config = NameList(fname=config_file,working_directory=wdir)
    config.read()
    temp = list(map(int,config.nlist['start_date'].split(":")))
    print(temp)
    config.nlist['start_date'] = datetime.datetime(temp[0],temp[1],temp[2], temp[3])
    # convert values to floats where possible.
    for key in config.nlist.keys():
        try:
           val = float(config.nlist[key])
        except:
           val = config.nlist[key]
        # get rid of any white spaces in strings
        config.nlist[key] = val
    jobsetup.inp = config.nlist
    jobsetup.add_plotting_options(config.nlist)
    return jobsetup

class JobSetUp:
    def __init__(self):
        self.inp = {}

    def parse_inputs(self, a):
        self.add_run_params(a)
        self.add_directories(a)
        self.add_plotting_options(a)
        return self.inp

    def add_input(self, inp, astr, default):
        if astr in inp.keys():
            self.inp[astr] = inp[astr]  # shapefiles
        else:
            self.inp[astr] = default

    def not_used(self, inp):
        poll_list = inp["pollutants"]
        # for  poll in poll_list:
        # poll['inCloudRemovalRate']
        # poll['belowCloudRemovalRate']
        # poll['emissionRate']
        # poll['depositVelocity']
        # poll['name']
        self.inp["wetflag"] = inp["usingWetDeposition"]
        self.inp["forwardflag"] = inp["forwardCalculation"]
        # 0 for large down to 3 for small.
        # self.inp['eflag'] = inp['eruptionSize']
        self.inp["polygon"] = None  # input polygon points to start from.

    def add_inverse_params(self):
         # time resolution for each inverse modeling run.
         # in hours. default 1 hour.
         self.inp['timeres'] = 1
         # vertical resolution for each inverse modeling run.
         # in km. default 2 km.
         self.inp['vertical_resolution'] = '1'
         self.inp['inv_vertical_resolution'] = 1000

    def add_run_params(self, inp):
        self.inp["owner"] = inp["owner"]  # TODO fix value
        self.inp["VolcanoName"] = inp["volcano"]["name"]
        self.inp["VolcanoLocation"] = inp["volcano"]["location"]
        self.inp["VolcanoName"] = inp["volcano"]["name"]
        self.inp["latitude"] = float(inp["sourceLatitude"])
        self.inp["longitude"] = float(inp["sourceLongitude"])
        self.inp["bottom"] = float(inp["volcano"]["height"])
        self.inp["top"] = float(inp["plumeHeight"])
        self.inp["durationOfSimulation"] = inp["simulationDurationInHours"]
        self.inp["emissionHours"] = inp["emissionHours"]
        # Must be in UTC.
        self.inp["start_date"] = datetime.datetime(
            inp["startDate"][0],
            inp["startDate"][1],
            inp["startDate"][2],
            inp["startTime"][0],
            inp["startTime"][1],
        )
        self.inp["samplingIntervalHours"] = inp["samplingIntervalInHours"]
        self.inp["meteorologicalData"] = inp["meteorologicalData"]
        self.inp["polygon"] = None  # input polygon points to start from.
        # runflag can be 'trajectory' or 'dispersion'
        if inp["runType"] == 0:
            self.inp["runflag"] = "trajectory"
        else:
            self.inp["runflag"] = "dispersion"
        self.inp["jobname"] = inp["runIdentifier"]
        # eflag should be a float. Large values give smaller plumes.
        # concentration value is multiplied by 10^(-1 * eflag)
        # This is because currently VAAC is used to using an 'ash reduction'
        # of 1 2 or 3 with 3 reducing the ash the most.
        # However large number to reduce plume is somewhat counter-intuitive.
        self.inp["eflag"] = float(inp["eruptionSize"])
        # below are currently always the same.
        self.inp["source_type"] = "uniform"

    def ensure_trailing_slash(self, dir):
        if dir[-1] != os.sep:
            return dir + os.sep
        return dir

    def add_directories(self, inp):
        self.inp["HYSPLIT_DIR"] = inp["readyProperties"]["directory"]["hysplit"]
        self.inp["MAP_DIR"] = inp["readyProperties"]["directory"]["map"]
        self.inp["WORK_DIR"] = self.ensure_trailing_slash(inp["workingDirectory"])
        self.inp["DATA_DIR"] = inp["dataDirectory"]

        ## TO DO CONTROL.default kept with repository in scripts directory.
        ## What will be the full path?
        self.inp["FILES_DIR"] = "./"

        self.inp["CONVERT_EXE"] = inp["readyProperties"]["executable"]["convert"]
        self.inp["GHOSTSCRIPT_EXE"] = inp["readyProperties"]["executable"][
            "ghostscript"
        ]
        self.inp["PYTHON_EXE"] = inp["readyProperties"]["executable"]["python"]

        self.inp["forecastDirectory"] = inp["readyProperties"]["directory"]["forecast"]
        self.inp["archivesDirectory"] = inp["readyProperties"]["directory"]["archives"]

    def add_plotting_options(self, inp):
        self.add_input(inp, "gisOption", 3)
        self.add_input(inp, "zoomFactor", 50)
        self.add_input(inp, "generatingPostscript", True)
        self.inp["generatingPDF"] = True
        self.inp["mapBackground"] = "arlmap"
        self.inp["mapProjection"] = 0
        self.inp["spatialPlotRadius"] = 500.0  # km
        self.inp["generatingPDF"] = True
        self.inp["graphicsResolution"] = 200
        self.inp["zip_compression_level"] = 3


    def make_test_inputs(self):
        #vname = "Kilauea"
        vname = "bezy"
        self.inp["owner"] = "A. Person"
        self.inp["top"] = 20000
        self.inp["durationOfSimulation"] = 36
        testdate = datetime.datetime.now() - datetime.timedelta(hours=24)
        # testdate = datetime.datetime(2020,10,10,11)
        testminutes = 15
        self.inp["start_date"] = datetime.datetime(
            testdate.year, testdate.month, testdate.day, testdate.hour, testminutes
        )
        if vname == 'inverse':
           self.inp["durationOfSimulation"] = 12
           self.inp["top"] = 10000
           testminutes = 0
           self.inp["start_date"] = datetime.datetime(
               testdate.year, testdate.month, testdate.day, testdate.hour, testminutes
           )
           vname='douglas'
           
        self.inp["emissionHours"] = 4
        self.inp["meteorologicalData"] = "GFS0p25"
        self.inp["EruptionSize"] = 0
        if vname.lower()=='bezy':
            # bezy data starts at 10/21 at 20:40
            #      ends at 10/22 at 21:10
            self.inp["meteorologicalData"] = "GFS0p25"
            self.inp["VolcanoName"] = "Bezymianny"
            self.inp['start_date'] = datetime.datetime(2020,10,21,19)
            self.inp["durationOfSimulation"] = 36
            self.inp["emissionHours"] = 24
            self.inp["top"] = 15000
            self.inp["bottom"] = 9455
            self.inp["latitude"] = 55.978
            self.inp["longitude"] = 160.587
              
        if vname.lower()=='raikoke':
            self.inp["meteorologicalData"] = "GFS0p25"
            self.inp["VolcanoName"] = "Raikoke"
            self.inp['start_date'] = datetime.datetime(2019,6,21,18)
            self.inp["durationOfSimulation"] = 24
            self.inp["emissionHours"] = 12
            self.inp["top"] = 15000
            self.inp["bottom"] = 1000
            self.inp["latitude"] = 48.292
            self.inp["longitude"] = 153.25
              
        if vname.lower() == "douglas":
            self.inp["meteorologicalData"] = "NAMHAK"
            self.inp["VolcanoName"] = "Douglas"
            self.inp["latitude"] = 58.855
            self.inp["longitude"] = -153.54
            self.inp["bottom"] = 7021
        if vname.lower() == "kilauea":
            self.inp["meteorologicalData"] = "NAMHHI"
            self.inp["VolcanoName"] = "Kilauea"
            self.inp["latitude"] = 19.421
            self.inp["longitude"] = -155.28
            self.inp["bottom"] = 4009
        if vname == "Reventador":
            self.inp["VolcanoName"] = "Reventador"
            self.inp["latitude"] = -0.08
            self.inp["longitude"] = -77.66
            self.inp["bottom"] = 3562
        if vname == "Veni":
            self.inp["VolcanoName"] = "Veniaminof"
            self.inp["latitude"] = 56.17
            self.inp["longitude"] = -159.38
            self.inp["bottom"] = 2507
        self.inp["samplingIntervalHours"] = 3
        self.inp["eflag"] = 0
        self.add_test_directories()
        self.add_plotting_options(inp={})
        self.inp["source_type"] = "uniform"
        self.inp["jobname"] = "ashtest"
        self.inp["runflag"] = "dispersion"
        return self.inp


class Job:
    def __init__(self, JOBID, jobName):
        self.JOBID = JOBID
        self.name = jobName

    def __str__(self):
        return "{}_{}".format(self.name, self.JOBID)

    def __repr__(self):
        return "{}_{}".format(self.name, self.JOBID)


class JobFileNameComposer:
    # amc : changed control and setup file naming.

    def __init__(self, workDir, jobId, jobname):
        self.workDirectory = workDir
        self.job = Job(jobId, jobname)

    # conprob always outputs files  with these names
    def get_original_conprob_filenames(self, stage=0):
        rlist = ["cmean", "cmax00", "cmax10", "cmax01", "cvarn", "cnumb", "ccoev"]
        rlist.extend(
            ["prob05", "prob10", "prob25", "prob50", "prob75", "prob90", "prob95"]
        )
        return rlist

    def get_conprob_filenames(self, stage=0):
        return [
            "{0!s}.{1!s}_{2:03d}".format(fn, self.job, stage)
            for fn in self.get_original_conprob_filenames(stage)
        ]

    # def get_exceedance_filename(self):
    #    tag = 'exceedance'
    #    ptype = 'ps'
    #    return '{}_{!s}.{}'.format(self.job, tag, ptype)

    def get_setup_filename(self, stage=0):
        if stage > 0:
            # return '{}_SETUP.{}.txt'.format(self.job, stage)
            return "SETUP.{}_{}".format(self.job, stage)
        return "SETUP.{}".format(self.job.JOBID)

    def get_control_suffix(self, stage=0):
        fname = self.get_control_filename(stage=stage)
        return fname.split(".")[-1]

    def get_control_filename(self, stage=0):
        if stage > 0:
            # return '{}_CONTROL.{}.txt'.format(self.job, stage)
            return "CONTROL.{}_{}".format(self.job, stage)
        return "CONTROL.{}".format(self.job.JOBID)

    def get_message_filename(self, stage=0):
        if stage > 0:
            # return '{}_MESSAGE.{}.txt'.format(self.job, stage)
            return "MESSAGE.{}_{}".format(self.job, stage)
        return "MESSAGE.{}".format(self.job.JOBID)

    def get_tdump_filename(self, stage=0):
        if stage > 0:
            return "{}_tdump.{}".format(self.job, stage)
        return "tdump.{}".format(self.job)

    def get_cdump_base(self, stage=0):
        if stage > 0:
            fname = "{0!s}_cdump.{1:03d}".format(self.job, stage)
        else:
            fname = self.get_cdump_filename().split(".")[0]
        return fname.split(".")[0]

    def get_cdump_filename(self, stage=0):
        if stage > 0:
            return "{0!s}_cdump.{1:03d}".format(self.job, stage)
        return "cdump.{}".format(self.job)

    def get_pardump_filename(self, stage=0):
        if stage > 0:
            return "{}_pardump.{}".format(self.job, stage)
        return "pardump.{}".format(self.job)

    def get_awips_filename(self, stage=0):
        if stage > 0:
            return "{}_awips2.{}".format(self.job, stage)
        return "awips2.{}".format(self.job)

    def get_tdump_unmodified_filename(self, stage=0):
        if stage > 0:
            return "{}_tdump.{}.full".format(self.job, stage)
        return "{}_tdump.full".format(self.job)

    def get_gis_status_filename(self, stage=0):
        if stage > 0:
            return "{}_gis.{}.txt".format(self.job, stage)
        return "{}_gis.txt".format(self.job)

    def get_gistmp_filename(self, stage=0):
        if stage > 0:
            return "{}_gistmp.{}.txt".format(self.job, stage)
        return "{}_gistmp.txt".format(self.job)

    def get_concplot_ps_filename(self, stage=0, ptype="ps"):
        if stage > 0:
            return "{}_conc_{}.{}".format(self.job, stage, ptype)
        return "{}_conc.{}".format(self.job, ptype)

    def get_summary_pdf_filename(self):
        return "{}_allplots.pdf".format(self.job)

    def basic_filename(self, stage, tag="conc", ptype="gif"):
        return "{}_{}_{}.{}".format(self.job, tag, stage, ptype)

    def filename_with_frame(self, stage, frame=0, tag="conc", ptype="gif"):
        return "{}_{}_{}-{:01d}.{}".format(self.job, tag, stage, frame, ptype)

    def filename_with_frame2(self, stage, frame=0, tag="conc", ptype="gif"):
        return "{}_{}_{}{:04d}.{}".format(self.job, tag, stage, frame, ptype)

    def get_trajplot_kml_filename(self, stage=0, frame=None):
        tag = "traj"
        ptype = "kml"
        return "{}_{}_{}_{:02d}.{}".format(self.job, tag, stage, frame, ptype)

    def get_trajplot_filename(self, stage=0, frame=None, ptype="gif"):
        tag = "traj"
        if frame is None:
            return self.basic_filename(stage, tag, ptype)
        else:
            return self.filename_with_frame2(stage, frame, tag, ptype)

    def get_concplot_filename(self, stage=0, frame=None, ptype="gif"):
        tag = "conc"
        if frame is None:
            return self.basic_filename(stage, tag, ptype)
        else:
            return self.filename_with_frame(stage, frame, tag, ptype)

    def get_exceedance_filename(self, stage=0, frame=None, ptype="gif"):
        tag = "exceedance"
        if frame is None:
            return self.basic_filename(stage, tag, ptype)
        else:
            return self.filename_with_frame(stage, frame, tag, ptype)

    def get_massloading_filename(self, stage=0, frame=None, ptype="gif"):
        tag = "massload"
        if frame is None:
            return self.basic_filename(stage, tag, ptype)
        else:
            return self.filename_with_frame(stage, frame, tag, ptype)

    def get_parxplot_filename(self, stage=0, frame=None, ptype="gif"):
        tag = "part"
        if frame is None:
            return self.basic_filename(stage, tag, ptype)
        else:
            return self.filename_with_frame(stage, frame, tag, ptype)

    def get_concentration_montage_filename(self, stage, frame, ptype="gif"):
        tag = "montage"
        if frame is None:
            return self.basic_filename(stage, tag, ptype)
        else:
            return self.filename_with_frame(stage, frame, tag, ptype)

    def get_totalpdf_filename(self, stage):
        tag = "allplots"
        ptype = "pdf"
        return self.basic_filename(stage, tag, ptype)

    def get_zipped_filename(self, tag=""):
        return "{}JOBID{}.zip".format(tag, self.job)

    def get_gis_zip_filename(self, stage=0):
        if stage > 0:
            return "{}_gis.{}.zip".format(self.job, stage)
        return "{}_gis.zip".format(self.job)

    def get_basic_kml_filename(self, program="concplot"):
        if program == "concplot":
            ktype = ""
        elif program == "parxplot":
            ktype = "part"
        else:
            ktype = ""
        return "HYSPLIT{}_{}.kml".format(ktype, self.job.JOBID)

    # def get_kml_filename(self, stage=0, frame=1):
    # if stage > 0:
    #    return '{}_HYSPLIT_{:02d}.{}.kml'.format(self.job, frame, stage)
    # return '{}_HYSPLIT_{:02d}.kml'.format(self.job, frame)

    def get_gelabel_filename(self, frame=1, ptype="ps"):
        if ptype == "txt":
            return "GELABEL_{:02d}_{}.{}".format(frame, self.job.JOBID, ptype)
        else:
            return "GELABEL_{:02d}_{}.{}".format(frame, self.job.JOBID, ptype)

    def get_basic_gelabel_filename(self, frame=1, ptype="ps"):
        if ptype == "txt":
            return "GELABEL_{:02d}_ps.{}".format(frame, ptype)
        else:
            return "GELABEL_{:02d}_ps.{}".format(frame, ptype)

    def get_kmz_filename(self, stage=0):
        return "{}_HYSPLIT_{}.kmz".format(self.job, stage)

    def get_maptext_filename(self):
        return "MAPTEXT.{}".format(self.job)

    def get_maptext_filename_for_zip(self):
        return "{}_MAPTEXT.txt".format(self.job)

    def get_awips_filenames(self, stage=0):
        files = glob.glob("awips2.{}*nc".format(self.job))
        return files

    def get_all_ashbase_filenames(self, stage=0):
        # Note the progress file and the run summary file are generated
        # by the Application layer written in Java.
        files = [
            "{}_progress.txt".format(self.job),
            "{}_run_setup_summary.txt".format(self.job),
        ]

        # cdump files
        files += glob.glob("{}_cdump*".format(self.job))
        files += glob.glob("cdump_{}".format(self.job))
        # netcdf files
        files += glob.glob("xrfile.{}.nc".format(self.job.JOBID))
        # pardump file
        files += glob.glob("{}_pardump*".format(self.job))
        files += glob.glob("pardump_{}".format(self.job))

        ## TEXT files which need .txt extension added.
        # MESSAGE file
        txtfiles = glob.glob("MESSAGE*{}".format(self.job.JOBID))
        # CONTROL file
        txtfiles += glob.glob("CONTROL*{}".format(self.job.JOBID))
        # SETUP file
        txtfiles += glob.glob("SETUP*{}".format(self.job.JOBID))
        ## TEXT files which need .txt extension added.
        # tdump files
        txtfiles += glob.glob("{}_tdump*".format(self.job))
        txtfiles += glob.glob("tdump_{}".format(self.job))

        # add .txt extension onto these files.
        # use copy instead of move so that if it is re-run it still works.
        [Helper.copy(x, "{}.txt".format(x)) for x in txtfiles]
        files += ["{}.txt".format(x) for x in txtfiles]

        # maptext file
        files += glob.glob("{}_MAPTEXT.txt".format(self.job))
        # KMZ files
        files += glob.glob("{}_*.kmz".format(self.job))
        # PDF files
        files += glob.glob("{}_allplots*.pdf".format(self.job))

        # png files
        files += glob.glob("{}_traj*.png".format(self.job))

        # get rid of any None
        files = [x for x in files if x]

        return files

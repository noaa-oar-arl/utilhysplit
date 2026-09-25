# ----- ------------------------------------------------------------------------
# Air Resources Laboratory
#
# runhelper.py -
#
# Classes
# Helper class contains functions for executing commands.
# JobSetUp class setups the dictionary which contains information for ash runs.
# Job and JobFileNameComposer class create filenames.
# ConcplotColors class

# Functions
# make_inputs_from_file . returns an instance of JobSetUp class.
# make_dir
# list_dirs

import datetime
import glob
import logging
import os
import pathlib

# import pytz
import shutil
import subprocess
import sys

# 18 APR 2020 (SYZ) - Initial.
# 15 Jun 2020 (AMC) - Adapted from locusts.py
# 07 Jun 2023 (AMC) - move ConcplotColors to utilhysplit.plotutils directory
# -----------------------------------------------------------------------------
from utilhysplit.hcontrol import NameList

logger = logging.getLogger(__name__)


def name2job(control_name):
    temp = control_name.splite('/')
    temp = temp[-1]
    temp = temp.split('.')
    suffix = temp[-1]
    temp = suffix.split('_')
    jobname = temp[0]
    jobid = temp[1]
    if len(temp)>=3:
       ensemble_id = temp[2]
    else:
       ensemble_id = None
    return Job(jobid,jobname,ensemble_id) 


class Job:
    def __init__(self, JOBID, jobname, ensemble_id=None):
        self.JOBID = JOBID  # created for each HYSPLIT run
        self.name = jobname  # indicates what event the run belongs to.
        self.ensemble_id=ensemble_id #indicates ensemble member

    def __str__(self):
        if isinstance(self.ensemble_id,(str,float,int)):
            return "{}_{}_{}".format(self.name,self.JOBID,self.ensemble_id)
        else:
            return "{}_{}".format(self.name, self.JOBID)

    def __repr__(self):
        if isinstance(ensemble_id,(str,float,int)):
            return "{}_{}".format(self.name, self.JOBID)
            return "{}_{}".format(self.name, self.JOBID)


def read_suffix(self, fnamestr):
    temp = fnamestr.split(".")
    suffix = temp[-1]
    temp2 = suffix.split("_")
    jobname = temp2[0]
    jobid = temp2[1]
    if len(temp2) == 3:
        stage = temp2[2]
    return jobname, jobid, stage


class JobFileNameComposer:
    def __init__(self, workDir, jobId, jobname):
        """
        jobname : identifies a particular event.
        jobId   : identifies a particular hysplit run
        stage   : identifies a particular part of a hysplit run.

        Examples:

        Eruption of Popocatepetl.
        All runs have jobname of Popocatepetl.
        jobid is set for each 'run'

        An ensemble run with the GEFS would have 31 stages (members).
        The stage would identify the GEFS member used.

        run for inverse modeling, the stages would identify the different unit source runs.
        if GEFS is used then stage identifies both the GEFS member and unit source run.
        The unit source runs are usually identified by start-time MMDDHH and start height in meters.
        e.g. 102201_11880.

        """
        self.workDirectory = workDir
        self.job = Job(jobId, jobname)
        self._stagelist = []
        self._basenamehash = self.set_defaults()

    def set_defaults(self):
        basenamehash = {}
        basenamehash["tdump"] = "tdump"
        basenamehash["cdump"] = "cdump"
        basenamehash["pardump"] = "pardump"
        basenamehash["message"] = "MESSAGE"
        basenamehash["control"] = "CONTROL"
        basenamehash["setup"] = "SETUP"
        basenamehash["concplot"] = "concplot"
        basenamehash["parxplot"] = "part"
        basenamehash["awips"] = "awips2"
        basenamehash["massload"] = "massload"
        basenamehash["trajplot"] = "traj"
        basenamehash["allpdf"] = "allplots"
        return basenamehash

    @property
    def basenamehash(self):
        return self._basenamehash

    @basenamehash.setter
    def basenamehash(self, inp):
        if isinstance(inp, dict):
            self._basenamehash.update(inp)

    @property
    def stagelist(self):
        return list(set(self._stagelist))

    def make_suffix(self, stage=None):
        if stage == 0:
            stage = None

        if isinstance(stage, (int, float, str)):
            if stage not in self._stagelist:
                self._stagelist.append(stage)
            return "{}_{}".format(self.job, stage)
        else:
            return self.job

    def get_control_filename(self, stage=None):
        suffix = self.make_suffix(stage)
        prefix = self.basenamehash["control"]
        return "{}.{}".format(prefix, suffix)

    def get_setup_filename(self, stage=None):
        suffix = self.make_suffix(stage)
        prefix = self.basenamehash["setup"]
        return "{}.{}".format(prefix, suffix)

    def get_pardump_filename(self, stage=None):
        suffix = self.make_suffix(stage)
        prefix = self.basenamehash["pardump"]
        return "{}.{}".format(prefix, suffix)

    def get_cdump_filename(self, stage=None):
        suffix = self.make_suffix(stage)
        prefix = self.basenamehash["cdump"]
        return "{}.{}".format(prefix, suffix)

    def get_tdump_filename(self, stage=None):
        suffix = self.make_suffix(stage)
        prefix = self.basenamehash["tdump"]
        return "{}.{}".format(prefix, suffix)

    def get_message_filename(self, stage=None):
        suffix = self.make_suffix(stage)
        prefix = self.basenamehash["message"]
        return "{}.{}".format(prefix, suffix)

    def get_awips_filename(self, stage=None):
        suffix = self.make_suffix(stage)
        prefix = self.basenamehash["awips"]
        return "{}.{}".format(prefix, suffix)

    def get_allplotspdf_filename(self, stage=None, frame=None, ptype="pdf"):
        prefix = self.make_suffix(stage)
        suffix = self.basenamehash["allpdf"]
        if not frame:
            return "{}_{}.{}".format(prefix, suffix, ptype)
        else:
            frame = self.check_frame()
            return "{}_{}_{:01d}.{}".format(prefix, suffix, frame, ptype)

    def get_concplot_filename(self, stage=None, frame=None, ptype="ps"):
        prefix = self.make_suffix(stage)
        suffix = self.basenamehash["concplot"]
        if not frame:
            return "{}_{}.{}".format(prefix, suffix, ptype)
        else:
            frame = self.check_frame()
            return "{}_{}_{:01d}.{}".format(prefix, suffix, frame, ptype)

    def get_parxplot_filename(self, stage=0, frame=None, ptype="gif"):
        prefix = self.make_suffix(stage)
        suffix = self.basenamehash["parxplot"]
        if not frame:
            return "{}_{}.{}".format(prefix, suffix, ptype)
        else:
            frame = self.check_frame(frame)
            return "{}_{}_{:01d}.{}".format(prefix, suffix, frame, ptype)

    def get_xrfile(self):
        # netcdf file that can be read by xarray module
        # all the stages are combined in this file so no need to input stage.
        suffix = self.make_suffix(stage=None)
        return "xrfile.{}.nc".format(suffix)

    def check_frame(self, frame):
        if isinstance(frame, float):
            frame = int(frame)
        if not isinstance(frame, int):
            frame = 0
            logger.warning(
                "get_massloading_filename frame input must be float or integer {}".format(
                    type(frame)
                )
            )
        return frame

    def get_massloading_filename(self, stage=0, frame=0, ptype="gif"):
        frame = self.check_frame(frame)
        base = self.basenamehash["massload"]
        prefix = self.make_suffix(stage)
        return "{}_{}_{:01d}.{}".format(prefix, base, frame, ptype)

    def get_zipped_filename(self):
        suffix = self.make_suffix(stage=None)
        return "JOBID{}.zip".format(suffix)

    def get_trajplot_filename(self, stage=0, frame=None, ptype="gif"):
        base = self.basenamehash["trajplot"]
        prefix = self.make_suffix(stage)
        if not frame:
            return "{}_{}.{}".format(prefix, base, ptype)
        else:
            frame = self.check_frame(frame)
            return "{}_{}_{:04d}.{}".format(prefix, base, frame, ptype)

    def get_kmz_filename(self, stage=0):
        base = "HYSPLIT"
        prefix = self.make_suffix(stage)
        suffix = "kmz"
        return "{}_{}.{}".format(prefix, base, suffix)

    def get_trajplot_kml_filename(self, stage=0, frame=None):
        # kml file is produced by trajplot. Name is same as the trajplot name.
        # but with kml instead of ps.
        trajplot = self.get_trajplot_filename(stage, frame, ptype="ps")
        return trajplot.replace(".ps", "_01.kml")

    def get_kml_filename(self, stage=0):
        # kml file is produced by trajplot. Name is same as the trajplot name.
        # but with kml instead of ps.
        base = "HYSPLIT"
        prefix = self.make_suffix(stage)
        suffix = "kml"
        return "{}_{}.{}".format(base, prefix, suffix)

    # def get_awips_filenames(self, stage=0):
    #    files = glob.glob("awips2.{}*nc".format(self.job))
    #    return files

    def get_gis_status_filename(self, stage=0):
        if stage > 0:
            return "{}_gis.{}.txt".format(self.job, stage)
        return "{}_gis.txt".format(self.job)

    # def get_gistmp_filename(self, stage=0):
    #    if stage > 0:
    #        return "{}_gistmp.{}.txt".format(self.job, stage)
    #    return "{}_gistmp.txt".format(self.job)

    # def get_concplot_ps_filename(self, stage=0, ptype="ps"):
    #    if stage > 0:
    #        return "{}_conc_{}.{}".format(self.job, stage, ptype)
    #    return "{}_conc.{}".format(self.job, ptype)

    # def get_summary_pdf_filename(self):
    #    return "{}_allplots.pdf".format(self.job)

    # def basic_filename(self, stage, tag="conc", ptype="gif"):
    #    return "{}_{}_{}.{}".format(self.job, tag, stage, ptype)

    # def filename_with_frame(self, stage, frame=0, tag="conc", ptype="gif"):
    #    return "{}_{}_{}-{:01d}.{}".format(self.job, tag, stage, frame, ptype)

    # def filename_with_frame2(self, stage, frame=0, tag="conc", ptype="gif"):
    #    return "{}_{}_{}{:04d}.{}".format(self.job, tag, stage, frame, ptype)

    # def get_trajplot_kml_filename(self, stage=0, frame=None):
    #    tag = "traj"
    #    ptype = "kml"
    #    return "{}_{}_{}_{:02d}.{}".format(self.job, tag, stage, frame, ptype)

    # def get_trajplot_filename(self, stage=0, frame=None, ptype="gif"):
    #    tag = "traj"
    #    if frame is None:
    #        return self.basic_filename(stage, tag, ptype)
    #    else:
    #        return self.filename_with_frame2(stage, frame, tag, ptype)

    # def get_concplot_filename(self, stage=0, frame=None, ptype="gif"):
    #    tag = "conc"
    #    if frame is None:
    #        return self.basic_filename(stage, tag, ptype)
    #    else:
    #        return self.filename_with_frame(stage, frame, tag, ptype)

    # def get_exceedance_filename(self, stage=0, frame=None, ptype="gif"):
    #    tag = "exceedance"
    #    if frame is None:
    #        return self.basic_filename(stage, tag, ptype)
    #    else:
    #        return self.filename_with_frame(stage, frame, tag, ptype)

    # def get_massloading_filename(self, stage=0, frame=None, ptype="gif"):
    #    tag = "massload"
    #    if frame is None:
    #        return self.basic_filename(stage, tag, ptype)
    #    else:
    #        return self.filename_with_frame(stage, frame, tag, ptype)

    # def get_parxplot_filename(self, stage=0, frame=None, ptype="gif"):
    #    tag = "part"
    #    if frame is None:
    #        return self.basic_filename(stage, tag, ptype)
    #    else:
    #        return self.filename_with_frame(stage, frame, tag, ptype)

    # def get_concentration_montage_filename(self, stage, frame, ptype="gif"):
    #    tag = "montage"
    #    if frame is None:
    #        return self.basic_filename(stage, tag, ptype)
    #    else:
    #        return self.filename_with_frame(stage, frame, tag, ptype)

    # def get_totalpdf_filename(self, stage):
    #    tag = "allplots"
    #    ptype = "pdf"
    #    return self.basic_filename(stage, tag, ptype)

    # def get_zipped_filename(self, tag=""):
    #    return "{}JOBID{}.zip".format(tag, self.job)

    # def get_gis_zip_filename(self, stage=0):
    #    if stage > 0:
    #        return "{}_gis.{}.zip".format(self.job, stage)
    #    return "{}_gis.zip".format(self.job)

    # def get_basic_kml_filename(self, program="concplot"):
    #    if program == "concplot":
    #        ktype = ""
    #    elif program == "parxplot":
    #        ktype = "part"
    #    else:
    #        ktype = ""
    #    return "HYSPLIT{}_{}.kml".format(ktype, self.job.JOBID)

    # def get_kml_filename(self, stage=0, frame=1):
    # if stage > 0:
    #    return '{}_HYSPLIT_{:02d}.{}.kml'.format(self.job, frame, stage)
    # return '{}_HYSPLIT_{:02d}.kml'.format(self.job, frame)

    def get_gelabel_filename(self, frame=1, ptype="ps"):
        if ptype == "txt":
            # return "GELABEL_{:02d}_{}.{}".format(frame, self.job.JOBID, ptype)
            return "GELABEL_ps.{}".format(ptype)
        else:
            return "GELABEL_{:02d}_ps.{}".format(frame, ptype)

    # def get_basic_gelabel_filename(self, frame=1, ptype="ps"):
    #    if ptype == "txt":
    #        return "GELABEL_{:02d}_ps.{}".format(frame, ptype)
    #    else:
    #        return "GELABEL_{:02d}_ps.{}".format(frame, ptype)

    # def get_kmz_filename(self, stage=0):
    #    return "{}_HYSPLIT_{}.kmz".format(self.job, stage)

    def get_maptext_filename(self):
        return "MAPTEXT.{}".format(self.job)

    def get_maptext_filename_for_zip(self):
        return "{}_MAPTEXT.txt".format(self.job)

    # def get_awips_filenames(self, stage=0):
    #    files = glob.glob("awips2.{}.nc".format(self.job))
    #    return files

    # def get_xrfile(self):
    # netcdf file that can be read by xarray module
    # all the stages are combined in this file so no need to input stage.
    # return "xrfile.{}.nc".format(self.job)

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


class AshDINameComposer(JobFileNameComposer):
    # used only in utildatainsertion right now for find_cdump_df_alt
    # for finding cdump files

    def get_cdump_filename(self, stage="EMIT_0"):
        cdumpfilename = self.get_di_filename(stage, "cdump.")
        return cdumpfilename

    def get_di_filename(self, stage, cstr):
        stage = str(stage)
        if "EMIT_" in stage:
            filename = stage.replace("EMIT_", cstr)
        elif "EMITIMES_" in stage:
            filename = stage.replace("EMITIMES_", cstr)
        elif "EMITIMES" in stage:
            filename = stage.replace("EMITIMES", cstr)
        elif "EMIT" in stage:
            filename = stage.replace("EMIT", cstr)
        else:
            filename = "{}.{}".format(cstr, stage)
        return filename

    def get_control_filename(self, stage="EMIT_0"):
        controlfilename = self.get_di_filename(stage, "CONTROL.")
        return controlfilename

    def get_setup_filename(self, stage=0):
        setupfilename = self.get_di_filename(stage, "SETUP.")
        return setupfilename

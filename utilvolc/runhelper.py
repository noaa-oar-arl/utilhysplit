# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# runhelper.py -
#
# Helper class contains functions for executing commands.
# JobSetUp class setups the dictionary which contains information for ash runs.
# Job and JobFileNameComposer class create filenames.

#  make_inputs_from_file . returns an instance of JobSetUp class.
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
            print(stdoutdata)
        if stderrdata is not None:
            logger.error(stderrdata)
            print(stderrdata)

    def execute(cmd, **kwargs):
        """
        cmd : string
        """
        #p = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)
        p = subprocess.Popen(cmd)
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


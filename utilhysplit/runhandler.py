# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import datetime
import logging
import os
import subprocess
import sys

"""
PRGMMR: Alice Crawford  ORG: ARL
PYTHON 3
This code written at the NOAA  Air Resources Laboratory
UID : r101
ABSTRACT: classes and functions for managing multiple HYSPLIT runs at the same time.
CTYPE: source code

Class ProcessList is used to manage multiple HYSPLIT runs.

function  is_process_running
          returns number of times a process is running

function  get_id
          returns pids associated with a process name
"""
logger = logging.getLogger(__name__)


class ProcessList(object):
    """initializes by getting all pids for process name.
    method startnew - starts a new process.
    method checkprocs - checks status of all processes started with startnew
    method updatelist - returns dictionary with information about processses which have finished.
    method move_from_temp - moves files from a temporary subdirectory where HYSPLIT was run to the directory above.
                           this method was added because HYSPLIT runs needed to be run in their own directory to prevent
                           them from occassionally trying to access the landuse file at the same time and failing.
    """

    def __init__(self):
        """gets all pids for the process name"""
        # self.currentpids = get_id(process_name)
        self.currentpids = []  # list of pid numbers
        # self.process_name = process_name
        self.pidhash = {}  # dictionary which pairs pid with description
        self.dirhash = {}  # dictionary which pairs pid with temporary directory
        self.count = len(self.currentpids)
        # for pid in self.currentpids:
        #    self.pidhash[pid] = descrip  +  process_name
        self.procra = []  ##array of processes which are started by the startnew method.
        self.err = (
            []
        )  ##array of tuples with following information, standard out, standard error, pid, description of pid.

        self.stdout = sys.stdout
        self.stderr = sys.stderr

        # self.stdout = subprocess.PIPE
        # self.stderr = subprocess.PIPE

    def pipe_stdout(self):
        self.stdout = subprocess.PIPE

    def pipe_stderr(self):
        self.stderr = subprocess.PIPE

    def startnew(self, callstr, wdir="./", descrip="new", nice=10):
        """ "starts a new process"""
        logger.info("CALLSTRING{}".format(callstr))
        xproc = subprocess.Popen(
            callstr,
            shell=False,
            cwd=wdir,
            #                         stdout=sys.stdout, stderr=sys.stderr)
            stdout=self.stdout,
            stderr=self.stderr,
        )
        pid = xproc.pid
        self.pidhash[pid] = descrip
        self.dirhash[pid] = wdir
        self.currentpids.append(pid)
        self.procra.append(xproc)
        logger.debug("Starting process " + str.join(" ", callstr) + " " + str(pid))
        return pid

    def killall(self):
        for proc in self.procra:
            subprocess.Popen.kill(proc)

    def move_from_temp(self, tempdir):
        """moves files from a temporary subdirectory where HYSPLIT was run to the directory above.
         this method was added because HYSPLIT runs needed to be run in their own directory to prevent
        them from occassionally trying to access the landuse file at the same time and failing.
        """
        os.chdir(tempdir)
        callstr = "mv CONTROL* ../"
        subprocess.call(callstr, shell=True)
        callstr = "mv SETUP* ../"
        subprocess.call(callstr, shell=True)
        callstr = "mv cdump* ../"
        subprocess.call(callstr, shell=True)
        callstr = "mv MESSAGE* ../"
        subprocess.call(callstr, shell=True)
        callstr = "mv WARNING* ../"
        subprocess.call(
            callstr, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        os.chdir("../")
        callstr = "rm -r " + tempdir
        subprocess.call(callstr, shell=True)

    def checkprocs(self, verbose=False, moveoutput=False):
        """check status of all procedures which were started using startnew method.
        Returns number of procedures started which are still running."""
        iii = 0
        newproclist = []
        logger.info("number procs " + str(len(self.procra)))
        for xproc in self.procra:
            pid = xproc.pid
            returncode = xproc.poll()
            if (
                returncode == 0
            ):  # this means process  has finished running successfully.
                stdout, stderr = xproc.communicate()
                if stderr:
                    logger.info("ERROR" + str(xproc.pid) + stderr)
                rstr = "Finished"
                # if 'Complete' not in stdout:
                #   print((xproc.pid, stdout))
                # Add a tuple with information about run which did not work out.
                #   self.err.append( (stdout, stderr, pid, self.pidhash[pid]), returncode )
                if moveoutput:
                    self.move_from_temp(self.dirhash[pid])
            elif returncode is None:
                ##this means process is still running.
                iii += 1
                newproclist.append(xproc)
                rstr = "Running"
            else:
                stdout, stderr = xproc.communicate()
                # Add a tuple with information about run which did not work out.
                self.err.append((stdout, stderr, pid, self.pidhash[pid], returncode))
                rstr = "?"
            logger.info("Status of " + self.pidhash[xproc.pid] + " " + rstr)
        self.procra = newproclist
        return iii

    def updatelist(
        self, descrip="new ", writelog=True, logname="currentpidlog.txt", verbose=False
    ):
        """returns dictionary with key=pids which have finished and value = description of pids.
        If writelog=True then will write a file with current processes running."""
        new = []
        finished = []
        still_running = []
        finished_hash = {}
        tempids = get_id(self.process_name)

        if tempids:
            for pid in tempids:
                if pid in self.currentpids:
                    still_running.append(
                        pid
                    )  # if tempid is in currentlist then process still running.
                    if verbose:
                        print(("PID", pid, " is still running"))
                else:
                    new.append(
                        pid
                    )  # if tempid is not in current list then this is a new process
                    if verbose:
                        print(("PID", pid, " has started"))

        if self.currentpids != []:
            for pid in self.currentpids:
                if pid not in tempids:
                    finished.append(pid)
                    if verbose:
                        print(
                            ("PID", pid, " has finished")
                        )  # if pid in currentlist is not in tempid then process finished.

        self.currentpids = still_running
        self.currentpids.extend(new)

        for pid in new:
            self.pidhash[pid] = descrip  # add new pids to dictionary

        for pid in finished:  # remove finished pids from dictionary
            # add them to finished hash and add current time to description.
            finished_hash[pid] = (
                self.pidhash.pop(pid)
                + " completed "
                + datetime.datetime.now().strftime("%Y:%m:%d %H:%M")
            )

        self.count = len(self.currentpids)

        if writelog:
            with open(logname, "w") as fid:
                fid.write(
                    datetime.datetime.now().strftime(
                        "Written at %Y / %m / %d  %H:%M \n"
                    )
                )
                fid.write("Processes currently running \n")
                for cpid in self.currentpids:
                    fid.write(str(cpid) + "  : " + self.pidhash[cpid] + "\n")

        return finished_hash


def is_process_running(process_name):
    """returns number of times the process is running.
    Uses information from ps -a"""
    if type(process_name) == str:
        process_name = [process_name]
    tmp = os.popen("ps -a").read()
    tcount = 0
    for pname in process_name:
        if pname in tmp:
            tcount += tmp.count(pname)
    return tcount


def get_id(process_name):
    """returns pids associated with a process name"""
    pidlist = []
    tmp = os.popen("ps -a").read().split("\n")
    for line in tmp:
        if process_name in line:
            plist = line.split()
            pidlist.append(int(plist[0]))
    return pidlist

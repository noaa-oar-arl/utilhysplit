#!/opt/Tools/anaconda3/envs/hysplit/bin/python
# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# -----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------

# import json
import datetime
import logging
import os
import sys
import traceback
from utils import setup_logger
# pylint: disable-msg=C0103

from utilvolc.runhelper import make_inputs_from_file
from maindispersion import MainEmitTimes
#from utilhysplit.hcontrol import NameList

logger = logging.getLogger(__name__)


def print_usage():
    print(
        """\
USAGE: DImain.py config.jobid.txt 1 5
The first input is the name of a configuration file for data insertion run.
The second input gives the 'time step'
The third input gives the number of 'time steps' to complete.

In the data insertion configfile the start_date tells the algorithm
to look for Emitimes files in which the file names indicate they are valid for date between
the start date and the start data + HoursToEnd.
Here HoursToEnd is calculated by adding the second input in hours to the start_date.



This program 

"""
    )


def change_config_name(cname, dt):
    temp = cname.split('.')
    new = '{}_{}'.format(temp[1], str(dt))
    return temp[0] + '.' + new + '.txt'

def get_jobid(cname):
    temp = cname.split('.')
    return temp[1]

if __name__ == "__main__":
    # Configure the logger so that log messages appears in the "Model Status" text box.
    # setup_logger(level=logging.DEBUG)
    setup_logger()
    #if len(sys.argv) != 3:
    #    print_usage()
    #    sys.exit(1)
    print(len(sys.argv))
    print(sys.argv) 
    configname = sys.argv[1]
    config = make_inputs_from_file('./',configname)
    dt = int(sys.argv[2])
    starthour= config.inp['start_date']
    done=False
    nloops = int(sys.argv[3])
    for nnn in range(0,nloops+1):
        newname = change_config_name(configname,dt*nnn)
        config.inp['start_date'] = starthour + datetime.timedelta(hours=nnn*dt)    
        config.inp['HoursToEnd'] = dt   
        config.write(newname)
        if 'qva' in config.inp['runflag']: config.inp['qvaflag'] = True
        else: config.inp['qvaflag'] = False
        jid = get_jobid(newname)
        crun = MainEmitTimes(config.inp,jid)
        crun.doit()
    sys.exit(0)

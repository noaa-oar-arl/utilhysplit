#!/opt/Tools/anaconda3/envs/hysplit/bin/python
# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# process_volcat.py -
#
# 01 JUN 2020 (AMC) - 
# 15 NOV 2022 (AMC) - updated how tests can be run
# -----------------------------------------------------------------------------
# To run in offline mode standard dispersion run  python ash_run.py -999
# To run in offline mode standard trajectory run  python ash_run.py -777
# To run in offline mode ensemble dispersion run  python ash_run.py -888
#
# TODO To run in offline mode ensemble trajectory run  python ash_run.py -555
#
#
# -----------------------------------------------------------------------------

#import json
import pandas as pd
import logging
import os
import sys
import traceback
#from abc import ABC, abstractmethod

import requests

from utilvolc import qva_logic
from utilvolc import volcat_event
from utilvolc import volcat_files
from utilhysplit.metfiles import gefs_suffix_list
from utilvolc.runhelper import JobSetUp,EventSetUp, make_inputs_from_file
from utils import setup_logger

logger = logging.getLogger(__name__)


def print_usage():
    print(
        """\
USAGE: ash_run.py RUNTYPE JOBID

where RUNTYPE is either 'run' or 'redraw' or 'test' (without quotations) and JOBID
is the job identification number.

If not using runtype of 'test' then
The following environment variables must be set prior to calling this script:
    VOLCANICASH_API_KEY         - secret key to access the hysplitash web APIs.
    VOLCANICASH_URL             - URL to the hysplitash web application.

If using 'test' then a configuration file is read. 
The call should be of the form
python ash_run.py test JOBID
"""
    )



def filterdf(sumdf):

    rstr = 'Filter by\n'
    rstr += '1: vaac\n'
    rstr += '2: country\n'
    rstr += '4: volcano name\n'
    rstr += '5: done '
    print(sumdf)
    print(rstr) 
    criteria = input()
    try:
        criteria=int(criteria)
    except:
        pass

    if criteria == 5:
       return True, sumdf
 
    elif criteria == 1:
       rstr = 'Choices are \n'
       vals = sumdf['vaac_region'].unique()
       rstr +=  str.join('\n',vals)
       rstr += '\nEnter 0 for all'
       print(rstr)
       crit2 = input()
       if crit2=='0': return False, sumdf
       return False, sumdf[sumdf['vaac_region']==crit2]

    elif criteria == 2:
       rstr = 'Choices are \n'
       vals = sumdf['volcano_country'].unique()
       rstr +=  str.join('\n',vals)
       rstr += '\nEnter 0 for all'
       print(rstr)
       crit2 = input()
       if crit2=='0': return False, sumdf
       return False, sumdf[sumdf['volcano_country']==crit2]
     
    elif criteria == 4:
       rstr = 'Choices are \n'
       vals = sumdf['volcano_name'].unique()
       rstr +=  str.join('\n',vals)
       rstr += '\nEnter 0 for all'
       print(rstr)
       crit2 = input()
       if crit2=='0': return False, sumdf
       return False, sumdf[sumdf['volcano_name']==crit2]

    else:
       print('Sorry that is not an option')
       return False, sumdf

def check_int_input(criteria, default):
    try:
        criteria=int(criteria)
    except Exception as eee:
        print('Sorry that is not one of the options {}'.format(criteria))
        print('option should be an integer')  
        criteria = default
    return criteria  
     
def main_menu():
    rstr = 'MAIN MENU \n'
    #rstr += '1: work on event\n'
    rstr += '1: check new\n'
    rstr += '2: check existing\n'
    rstr += '100: exit\n'
    print(rstr)
    choice = input()
    choice = check_int_input(choice,default=10)
    return choice 


# if choice 2 is chosen from main menu
def second_menu():
    rstr = 'MENU create new events\n'
    rstr += '1: download files\n'
    rstr += '2: summary only \n'
    rstr += '4: write emit imes \n'
    rstr += '5: write parallax corrected files \n'
    rstr += '6: setup inversion runs \n'
    rstr += '10: return to MAIN MENU \n'
    print(rstr)
    choice = input()
    choice = check_int_input(choice,default=10)
    return choice 


if __name__ == "__main__":
    setup_logger(level=logging.WARNING)
    #setup_logger()
    # create a test to run without inputs from web page.
    logging.getLogger().setLevel(20)
    logging.basicConfig(stream=sys.stdout)

    # this sets the directories.
    esetup = EventSetUp()

    done = False
    original_sumdf = pd.DataFrame()
    while not done:

        choice = main_menu() 

        #------------------------------------------------------------------------
        if choice == 100: 
           print('good bye')
           done=True

        #------------------------------------------------------------------------
        #------------------------------------------------------------------------

        if choice in [1,2]:
            choice2 = 0

            # get the summary df file
            if choice==1 or original_sumdf.empty:
                print('here because ', choice, original_sumdf.empty)
                work = qva_logic.WorkFlow(esetup.inp)
                print('How many hours back to check? (default = {}, max = {})'.format(24,24*7))
                hours = input()
                try:
                   hours = int(hours)
                except:
                   hours = 24
                if hours < 0: hours = -1*hours
                if hours > 24*7 : hours = 24*7
                sumdf = volcat_files.get_summary_file_df(esetup.inp['JPSS_DIR'],hours=hours)
                original_sumdf = sumdf.copy()
            # use the original sumdf file
            elif choice==2:
                sumdf = original_sumdf.copy()

            # filter as needed.
            filter_done = False
            iii=0
            while not filter_done:
                  filter_done, sumdf = filterdf(sumdf) 
                  iii+=1
                  if iii > 100: break

            download=False
            while choice2 != 10:
                 
                # after filtering the sumdf go to second menu.
                choice2 = second_menu() 

                # back to main menu
                if choice2 == 10: 
                   continue

                # create the Event classes
                # and download the files
                if choice2 == 1 or not download:
                    work.greenlist = sumdf['volcano_name'].unique()
                    work.get_volcat(hours=hours,verbose=False)
                    print('created events for ', work.ehash.keys())
                    download=True

                # write_parallax_corrected
                if choice2 == 5:
                   work.write_parallax_corrected(gridspace = 0.1)

                # write emit times files
                elif choice2 == 4:
                      print('WRITING EMITIMES')
                      for volc in work.ehash.keys():
                          print(volc)
                          eve = work.ehash[volc] 
                          eve.write_emit()
                          

                     
     
               

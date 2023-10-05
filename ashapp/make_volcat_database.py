import datetime
import glob
import logging
import pandas as pd
import os
import sys
import sqlite3

from utilvolc import volcat_files
from utilvolc import volcat
from utilvolc import volcano_names
from utils import setup_logger


logger = logging.getLogger(__name__)

if __name__ == "__main__":
    setup_logger()

    logdflist = []

    inp = {}
    inp['JPSS'] = '/pub/jpsss_upload/' 
    inp['VOLCAT_DIR'] = '/pub/ECMWF/JPSS/VOLCAT/Files/' 

    # this is for adding data when it is not in a summary file.
    vprops = volcano_names.VolcList('/hysplit-users/alicec/utilhysplit/ashapp/data/volclist.txt')

    edate = datetime.datetime.now()
    hours = 72
    sumdf = volcat_files.get_summary_file_df(inp['JPSS'],hours=hours,edate=edate)
    logdf = volcat_files.get_log_files(sumdf,inp)
    logger.info('Setting up new data {}'.format(sumdf.shape)) 
    logdflist.append(logdf)

    # Raikoke case
    edate = datetime.datetime(2019,6,21)
    hours = 48
    files = glob.glob(inp['VOLCAT_DIR']+'/Raikoke/*VOLCAT*nc')
    vname = volcat.VolcatName(files[0])
    khash = {'volcano_name':['Raikoke']}
    #record = vprops.get_record('Raikoke') 
    # transform into dictionary
    #record = record.to_dict(orient='records')[0]
    #logdf_raikoke = volcat_files.get_log_files(sumdf2,inp,**khash)
    #logdf_raikoke = volcat.flist2eventdf(files, record)
    #logger.info('Setting up Raikoke case {}'.format(sumdf2.shape))
    #logdflist.append(logdf_raikoke)


    # Mauna_Loa case
    edate = datetime.datetime(2022,11,29)
    hours = 24
    sumdf2 = volcat_files.get_summary_file_df(inp['JPSS'],hours=hours,edate=edate)
    khash = {'volcano_name':['Mauna Loa']}
    logdf2 = volcat_files.get_log_files(sumdf2,inp,**khash)
    logger.info('Setting up Mauna Loa case {}'.format(sumdf2.shape))
    logdflist.append(logdf2)

    # Popocatepetl case
    edate = datetime.datetime(2023,5,11)
    hours = 48
    sumdf3 = volcat_files.get_summary_file_df(inp['JPSS'],hours=hours,edate=edate)
    sumdf3 = sumdf3[sumdf3['volcano_name'] == 'Popocatepetl']
    khash = {'volcano_name':['Popocatepetl']}
    logdf3 = volcat_files.get_log_files(sumdf3,inp,**khash)
    logger.info('Setting up Popocatepetl case {}'.format(logdf3.shape))
    logdflist.append(logdf3)

    # Sheveluch case
    #files = glob.glob(inp['VOLCAT_DIR']+'/Sheveluch/*VOLCAT*nc')
    #vname = volcat.VolcatName(files[0])
    #print(vname.vhash)

    #record = vprops.get_record('Sheveluch') 
    #logdf4 = volcat.flist2eventdf(files, record[0])
    #df['observation_date'].unique()
    #logger.info('Setting up Sheveluch case {}'.format(logf4.shape))


    # Shishaldin case
    vname = 'Shishaldin'
    edate = datetime.datetime(2023,7,23)
    hours = 168
    sumdf5 = volcat_files.get_summary_file_df(inp['JPSS'],hours=hours,edate=edate)
    #sumdf5 = sumdf3[sumdf3['volcano_name'] == vname]
    khash = {'volcano_name':['Shishaldin']}
    logdf5 = volcat_files.get_log_files(sumdf5,inp,**khash)
    logger.info('Setting up {} case {}'.format(vname, sumdf5.shape))
    logdflist.append(logdf5)

    edate = datetime.datetime(2023,8,27)
    hours = 48
    sumdf5 = volcat_files.get_summary_file_df(inp['JPSS'],hours=hours,edate=edate)
    #sumdf5 = sumdf3[sumdf3['volcano_name'] == vname]
    khash = {'volcano_name':['Shishaldin']}
    logdf5 = volcat_files.get_log_files(sumdf5,inp,**khash)
    logger.info('Setting up {} case {}'.format(vname, sumdf5.shape))
    logdflist.append(logdf5)
 
    edate = datetime.datetime(2023,10,4)
    hours = 48
    sumdf5 = volcat_files.get_summary_file_df(inp['JPSS'],hours=hours,edate=edate)
    #sumdf5 = sumdf3[sumdf3['volcano_name'] == vname]
    khash = {'volcano_name':['Shishaldin']}
    logdf5 = volcat_files.get_log_files(sumdf5,inp,**khash)
    logger.info('Setting up {} case {}'.format(vname, sumdf5.shape))
    logdflist.append(logdf5)


    #some shishaldin files may have been misnamed as Pavlof

    # put them all together
    totaldf = pd.concat(logdflist,axis=0)
    basename = '/hysplit-users/alicec/utilhysplit/ashapp/VolcatDataBase'    
    totaldf.to_csv(basename + '.csv')     
    #connection = sqlite3.connect(basename + ".db")


    # if_exists may also be append or fail.
    # todo : setup dtype argument with types for the columns.

    # totaldf.to_sql(basename+".db", connection, if_exists="replace")
 

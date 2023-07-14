import datetime
import logging
import pandas as pd
import os
import sys

from utilvolc import volcat_files
from utils import setup_logger


logger = logging.getLogger(__name__)

if __name__ == "__main__":
    setup_logger()

    inp = {}
    inp['JPSS'] = '/pub/jpsss_upload/' 

    edate = datetime.datetime.now()
    hours = 168
    sumdf = volcat_files.get_summary_file_df(inp['JPSS'],hours=hours,edate=edate)
    logger.info('Setting up new data {}'.format(sumdf.shape)) 
    # Mauna_Loa case
    edate = datetime.datetime(2022,11,29)
    hours = 24
    sumdf2 = volcat_files.get_summary_file_df(inp['JPSS'],hours=hours,edate=edate)
    sumdf2 = sumdf2[sumdf2['volcano_name'] == 'Mauna Loa']
    logger.info('Setting up Mauna Loa case {}'.format(sumdf2.shape))

    # Popocatepetl case
    edate = datetime.datetime(2023,5,11)
    hours = 48
    sumdf3 = volcat_files.get_summary_file_df(inp['JPSS'],hours=hours,edate=edate)
    sumdf3 = sumdf3[sumdf3['volcano_name'] == 'Popocatepetl']
    logger.info('Setting up Popocatepetl case {}'.format(sumdf3.shape))

    # Sheveluch case
    #edate = datetime.datetime(2023,5,11)
    #hours = 48
    #sumdf3 = volcat_files.get_summary_file_df(inp['JPSS'],hours=hours,edate=edate)
    #sumdf3 = sumdf2[sumdf2['volcano_name'] == 'Popocatepetl']

    totaldf = pd.concat([sumdf,sumdf2,sumdf3],axis=0)

    totaldf.to_csv('/hysplit-users/alicec/utilhysplit/ashapp/VolcatDataBase.csv')     
 

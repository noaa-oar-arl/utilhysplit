import datetime
from utilvolc.ashapp import metfiles
import pytest
import os

def test_001():
    # tests for forecast data
    fdir = '/pub/forecast/'
    
    # tests for gfs
    metid='gfs'
    suffix='gfsf'
    finder = metfiles.MetFileFinder(metid)
    finder.set_forecast_directory(fdir)    
    finder.set_archives_directory('/pub/archive/') 

    # test on files the day before today.
    dtemp = datetime.datetime.now() - datetime.timedelta(hours=24)
    year = dtemp.year
    month = dtemp.month
    day = dtemp.day

    # start at 01 the day before today and run forward an hour.
    dstart = datetime.datetime(year, month, day, 1)
    duration = 1

    # look at '00' cycle
    rfiles = finder.find_forecast_cycle(dstart, duration, '00')
    rstr = dstart.strftime("%Y%m%d")
    answer = ['{}{}/hysplit.t00z.{}'.format(fdir,rstr,suffix)]
    assert rfiles == answer

    # look at other cycle
    dtemp2 = dstart - datetime.timedelta(hours=24)
    rstr = dtemp2.strftime("%Y%m%d")
    for cycle in ['06','12','18']:
        rfiles = finder.find_forecast_cycle(dstart, duration, cycle)
        answer = ['{}{}/hysplit.t{}z.{}'.format(fdir,rstr,cycle,suffix)]
        assert rfiles == answer

    # start at 01 the day before today and run backwards 6 hours.
    duration = -6
    for cycle in ['00','06','12','18']:
        rfiles = finder.find_forecast_cycle(dstart, duration, cycle)
        rstr = dtemp2.strftime("%Y%m%d")
        answer = ['{}{}/hysplit.t{}z.{}'.format(fdir,rstr,cycle,suffix)]
        rstr = dstart.strftime("%Y%m%d")
        answer.append('{}{}/hysplit.t{}z.{}'.format(fdir,rstr,cycle,suffix))
        assert rfiles == answer

    # start at 6z the day before today and run forward 48 hours
    duration=24 
    dstart = datetime.datetime(year, month, day, 6)
    for cycle in ['00','06','12','18']:
        answer = []
        rfiles = finder.find_forecast_cycle(dstart, duration, cycle)
        if cycle in ['12','18']:
            dtemp2 = dstart - datetime.timedelta(hours=24)
            rstr = dtemp2.strftime("%Y%m%d")
            answer.append('{}{}/hysplit.t{}z.{}'.format(fdir,rstr,cycle,suffix))
        rstr = dstart.strftime("%Y%m%d")
        answer.append('{}{}/hysplit.t{}z.{}'.format(fdir,rstr,cycle,suffix))
        dtemp2 = dstart + datetime.timedelta(hours=24)
        rstr = dtemp2.strftime("%Y%m%d")
       
        if os.path.isfile('{}{}/hysplit.t{}z.{}'.format(fdir,rstr,cycle,suffix)):
            answer.append('{}{}/hysplit.t{}z.{}'.format(fdir,rstr,cycle,suffix))
        print(rfiles)
        print(answer)
        #assert rfiles == rstr

test_001() 


   
 
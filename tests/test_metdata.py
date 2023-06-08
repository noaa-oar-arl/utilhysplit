import datetime
import os

import pytest
from utilhysplit import metfiles

def test_001():

    mf = metfiles.MetFileFinder(metid='GEFS')
    mf.set_ens_member('gec00')
    d1 = datetime.datetime.now()
    d1 = datetime.datetime(d1.year, d1.month, d1.day, 0,0)
    duration = 1
    test = mf.find_forecast(d1,duration)
    a1 = d1.strftime('/pub/forecast/%Y%m%d/')
    a2 = 'hysplit.t00z.gefs.gec00'
    assert test == [(a1,a2)]
  
    d1 = d1 - datetime.timedelta(hours=24) 
    d1 = datetime.datetime(d1.year, d1.month, d1.day, 6,0)
    test = mf.find_forecast(d1,duration)
    a1 = d1.strftime('/pub/forecast/%Y%m%d/')
    a2 = 'hysplit.t06z.gefs.gec00'
    assert test == [(a1,a2)]

    d1 = d1 - datetime.timedelta(hours=24) 
    d1 = datetime.datetime(d1.year, d1.month, d1.day, 12,0)
    test = mf.find_forecast(d1,duration)
    a1 = d1.strftime('/pub/forecast/%Y%m%d/')
    a2 = 'hysplit.t12z.gefs.gec00'
    assert test == [(a1,a2)]

    d1 = d1 - datetime.timedelta(hours=24) 
    d1 = datetime.datetime(d1.year, d1.month, d1.day, 18,0)
    test = mf.find_forecast(d1,duration)
    answer = d1.strftime('/pub/forecast/%Y%m%d/hysplit.t18z.gefs.gec00')
    a1 = d1.strftime('/pub/forecast/%Y%m%d/')
    a2 = 'hysplit.t18z.gefs.gec00'
    assert test == [(a1,a2)]


def test_002():
    # tests for gfs0p25 archive data
    adir = "/pub/archives/"
    metid = "gfs0p25"
    finder = metfiles.MetFileFinder(metid)
    finder.set_archives_directory(adir)
    year = 2020
    day = 31
    month = 12
    hour = 6
    dstart = datetime.datetime(year, month, day, hour)

    def subfunc(dstart, hours):
        mdir = "/pub/archives/gfs0p25/"
        mf = dstart.strftime("%Y%m%d_gfs0p25")
        answer = [(mdir, mf)]
        dstart += datetime.timedelta(hours=hours)
        mf = dstart.strftime("%Y%m%d_gfs0p25")
        answer.append((mdir, mf))
        dstart += datetime.timedelta(hours=hours)
        mf = dstart.strftime("%Y%m%d_gfs0p25")
        answer.append((mdir, mf))
        return answer

    # test forward
    duration = 48
    mff = finder.find_archive(dstart, duration)
    answer = subfunc(dstart, 24)
    # print(mff)
    # print(answer)
    assert answer == mff

    mff = finder.find_archive(dstart, duration)
    answer = subfunc(dstart, 24)
    assert answer == mff
    # test backward
    duration = -48
    mff = finder.find_archive(dstart, duration)
    answer = subfunc(dstart, -24)
    answer.reverse()
    assert answer == mff


# def test_001():
#    # tests for forecast data
#    fdir = '/pub/forecast/'

# tests for gfs
#    metid='gfs'
#    suffix='gfsf'


def test_003():
    # tests for forecast data
    fdir = "/pub/forecast/"

    # tests for gfs
    metid = "gfs"
    suffix = "gfsf"
    finder = metfiles.MetFileFinder(metid)
    finder.set_forecast_directory(fdir)
    finder.set_archives_directory("/pub/archive/")

    # test on files the day before today.
    dtemp = datetime.datetime.now() - datetime.timedelta(hours=24)
    year = dtemp.year
    month = dtemp.month
    day = dtemp.day

    # start at 01 the day before today and run forward an hour.
    dstart = datetime.datetime(year, month, day, 1)
    duration = 1

    # look at '00' cycle
    rfiles = finder.find_forecast_cycle(dstart, duration, "00")
    rstr = dstart.strftime("%Y%m%d")
    answer = ["{}{}/hysplit.t00z.{}".format(fdir, rstr, suffix)]
    assert rfiles == answer

    # look at other cycle
    dtemp2 = dstart - datetime.timedelta(hours=24)
    rstr = dtemp2.strftime("%Y%m%d")
    for cycle in ["06", "12", "18"]:
        rfiles = finder.find_forecast_cycle(dstart, duration, cycle)
        answer = ["{}{}/hysplit.t{}z.{}".format(fdir, rstr, cycle, suffix)]
        assert rfiles == answer

    # start at 01 the day before today and run backwards 6 hours.
    duration = -6
    for cycle in ["00", "06", "12", "18"]:
        rfiles = finder.find_forecast_cycle(dstart, duration, cycle)
        rstr = dtemp2.strftime("%Y%m%d")
        answer = ["{}{}/hysplit.t{}z.{}".format(fdir, rstr, cycle, suffix)]
        rstr = dstart.strftime("%Y%m%d")
        answer.append("{}{}/hysplit.t{}z.{}".format(fdir, rstr, cycle, suffix))
        assert rfiles == answer

    # start at 6z the day before today and run forward 48 hours
    duration = 24
    dstart = datetime.datetime(year, month, day, 6)
    for cycle in ["00", "06", "12", "18"]:
        answer = []
        rfiles = finder.find_forecast_cycle(dstart, duration, cycle)
        if cycle in ["12", "18"]:
            dtemp2 = dstart - datetime.timedelta(hours=24)
            rstr = dtemp2.strftime("%Y%m%d")
            answer.append("{}{}/hysplit.t{}z.{}".format(fdir, rstr, cycle, suffix))
        rstr = dstart.strftime("%Y%m%d")
        answer.append("{}{}/hysplit.t{}z.{}".format(fdir, rstr, cycle, suffix))
        dtemp2 = dstart + datetime.timedelta(hours=24)
        rstr = dtemp2.strftime("%Y%m%d")

        if os.path.isfile("{}{}/hysplit.t{}z.{}".format(fdir, rstr, cycle, suffix)):
            answer.append("{}{}/hysplit.t{}z.{}".format(fdir, rstr, cycle, suffix))
        # assert rfiles == rstr


test_001()

# metdata.py
# Functions designed to find the correct directories and
# meteorological files for a HYSPLIT run
import os
from datetime import timedelta as td
from datetime import datetime

# Find forecast met files for HYSPLIT run


def forecast(dstart, metdata):
    """Finds forecast meteorology data files
    dstart : datetime object of start date
    metdata : string (options: GFS, GFS0p25, NAM, NAMAK, NAMHI)
    """
    FCTDIR = '/pub/forecast/'
    ctimes = [0, 6, 12, 18]
    times = []
    metdirlist = []
    metfilelist = []
    if (metdata.lower() == 'gfs') or (metdata.lower() == 'gfs0p25'):
        met = metdata.lower()
    elif(metdata.lower() == 'nam'):
        met = 'namsf'
    elif(metdata.lower() == 'namak'):
        met = 'namsf.AK'
    elif(metdata.lower() == 'namhi'):
        met = 'namsf.HI'
    metnamefinal = 'No data found'
    for hrs in ctimes:
        times.append(datetime(dstart.year, dstart.month, dstart.day, hrs))
    for dtm in times:
        if (dtm <= dstart):
            metime = dtm
    metdir1 = FCTDIR + metime.strftime('%Y%m%d') + '/'
    metfilename = 'hysplit.' + metime.strftime('t%Hz') + '.' + met + 'f'
    if 'nam' in metdata.lower():
        metfilename = 'hysplit.' + metime.strftime('t%Hz') + '.' + met
    metname = metdir1 + metfilename
    if not os.path.exists(metname):
        metime = metime - td(hours=6)
        metdir1 = FCTDIR + metime.strftime('%Y%m%d') + '/'
        metfilename = 'hysplit.' + metime.strftime('t%Hz') + '.' + met + 'f'
        if 'nam' in metdata.lower():
            metfilename = 'hysplit.' + metime.strftime('t%Hz') + '.' + met
    metfilelist.append(metfilename)
    metdirlist.append(metdir1)
    return metdirlist, metfilelist

# Find met files for HYSPLIT archive runs


def archive(dstart, dend, metdata, direction):
    """
    dstart : datetime object start date
    dend : datetime object end date
    metdata : string : GFS0p5, GDAS1, NAM12, NAM, NAMAK, NAMHI
    direction: string : Forward, Back"""

    from datetime import timedelta as td
    from datetime import date

    DIR = '/pub/archives/'
    metdirlist = []
    metfilelist = []
    datelist = []
    if (metdata.lower() == 'gfs0p25') or (metdata.lower() == 'nam12'):
        metdir1 = DIR + metdata.lower() + '/'
        if (direction.lower() == 'forward'):
            while dstart <= dend:
                datelist.append(dstart)
                dstart += td(days=1)
        elif (direction.lower() == 'back'):
            while dstart >= dend-td(days=1):
                datelist.append(dstart)
                dstart -= td(days=1)
        for date in datelist:
            metfilelist.append(date.strftime('%Y%m%d') + '_' + metdata.lower())
            metdirlist.append(metdir1)

    if (metdata.lower() == 'nam') or (metdata.lower() == 'namak') or (metdata.lower() == 'namhi'):
        metdir1 = DIR + metdata.lower() + '/'
        ctimes = [0, 6, 12, 18]
        times = []
        if metdata.lower() == 'namak':
            end = '.namsa.AK'
        if metdata.lower() == 'namhi':
            end = '.namsa.HI'
        if metdata.lower() == 'nam':
            end = '.namsa'
        for hrs in ctimes:
            times.append(datetime(dstart.year, dstart.month, dstart.day, hrs))
        if (direction.lower() == 'forward'):
            for dtm in times:
                if (dtm <= dstart):
                    metime = dtm
            while metime <= dend:
                datelist.append(metime)
                metime += td(hours=6)
        if (direction.lower() == 'back'):
            for dtm in times:
                if (dtm <= dstart):
                    metime = dtm + td(hours=6)
            while metime >= dend:
                datelist.append(metime)
                metime -= td(hours=6)
        for date in datelist:
            metfilelist.append(date.strftime('%Y%m%d') + '_hysplit.' + date.strftime('t%Hz') + end)
            metdirlist.append(metdir1)

    if (metdata.lower() == 'gdas1'):
        metdir1 = DIR + metdata.lower() + '/'
        smonyr = dstart.strftime('%b%y').lower()
        emonyr = dend.strftime('%b%y').lower()
        if (smonyr == emonyr):
            sweek = (dstart.day - 1) // 7 + 1
            eweek = (dend.day - 1) // 7 + 1
            if (sweek == eweek):
                metdirlist.append(metdir1)
                metfilelist.append(metdata.lower() + '.' + smonyr + '.w' + str(sweek))
            if (sweek != eweek):
                if (direction.lower() == 'forward'):
                    nweeks = list(range(sweek, eweek + 1, 1))
                if (direction.lower() == 'back'):
                    nweeks = list(range(sweek, eweek - 1, -1))
                    y = 0
                    while y < len(nweeks):
                        metdirlist.append(metdir1)
                        metfilelist.append(metdata.lower() + '.' + smonyr + '.w' + str(nweeks[y]))
                        y += 1
        if (smonyr != emonyr):
            sweek = (dstart.day - 1) // 7 + 1
            eweek = (dend.day - 1) // 7 + 1
            ndays = abs((dend - dstart).days)
            if ndays < 14:
                metdirlist.append(metdir1)
                metfilelist.append(metdata.lower() + '.' + smonyr + '.w' + str(sweek))
                metdirlist.append(metdir1)
                metfilelist.append(metdata.lower() + '.' + emonyr + '.w' + str(eweek))
            if ndays > 14:
                dstep = ndays // 7 + 1
                tmpdate = dstart
                for x in range(0, dstep + 1):
                    if (direction.lower() == 'forward'):
                        tmpdate = tmpdate + td(days=7*x)
                    if (direction.lower() == 'back'):
                        tmpdate = tmpdate - td(days=7*x)
                    tmpweek = (tmpdate.day - 1) // 7 + 1
                    tmpmonyr = tmpdate.strftime('%b%y').lower()
                    metdirlist.append(metdir1)
                    metfilelist.append(metdata.lower() + '.' + tmpmonyr + '.w' + str(tmpweek))
    return metdirlist, metfilelist

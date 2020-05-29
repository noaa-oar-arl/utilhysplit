""" Obs Utilities """

import datetime
import sys

import numpy as np


def find_near(df, latlon, distance=100, sid='site_num', drange=None):
    """find all values in the df dataframe column sid which are within distance
    (km) of lat lon point. output dictionary with key as value in column sid
    and value tuple (latitude, longitude)

     Parameters
     ----------
     latlon : tuple or list
              (longitude, latitude)
     distance : float
               kilometers
     sid: string
          name of column
     drange: tuple or list with two datetimes
          consider rows with dates between these two dates.

     Returns
     --------
     lhash: dictionary
         key is the value in column sid and value is (latitude, longitude)
         position.
    """
    degree2km = 111
    if drange:
        df = timefilter(df.copy(), drange)
    lhash = get_lhash(df, sid)
    for key in lhash.keys:
        xd = (lhash[key][1] - latlon[1]) * degree2km * np.cos(
            latlon[1] * np.pi / 180.0)
        yd = (lhash[key][0] - latlon[0]) * degree2km
        dd = (xd**2 + yd**2)**0.5
        if dd > distance:
            lhash.pop(key, None)
    return lhash


def write_datem(df,
                obscolumn='obs',
                dname='datemfile.txt',
                sitename='1',
                info=None,
                drange=None,
                fillhours=1):
    """returns string in datem format (See NOAA ARL).
     datem format has the following columns:
     Year, Month, Day, Hour, Duration, lat, lon, Concentration (units), site
     id, height

     Parameters
     ----------
     obscolumn : string
            name of column with values to write in the Concentration column.
     dname : string
             name of the output file.
     sitename : string.
             If it is the name of a column in the dataframe then
             that column will be used to generate the site name column in the
             datem file. If is not the name of a column, then the string will
             be used as the site name.
     info : string
           will be written to the second line of the header.
     drange : list of two time stamp objects.
     Returns
     --------
     runstring: string
       string in datem format.
    """
    if drange:
        df = timefilter(df, drange)

    units = df['units'].tolist()
    units = list(set(units))
    if drange:
       sdate = drange[0]
    else:
       sdate = datetime.datetime(1900,1,1,2)
    if len(units) > 1:
        print('WARNING, more than one type of unit ', units)
    ustr = ''
    for uuu in units:
        ustr += uuu + ' '
    runstring = "Beginning date " + sdate.strftime(
        "%Y %m %d %H:%M") + " UTC ---"
    runstring += 'Information '
    if info:
        runstring += info + "\n"
    else:
        runstring += "\n"
    runstring += 'Year, Month, Day, Hour:Minute (UTC), Dur(hhmm) ,  LAT, LON, Concentration (' + \
        ustr + "), sid, height\n"
    lat = df['latitude']
    lon = df['longitude']
    cval = df[obscolumn]
    # print t2
    t1 = df['time']
    duration = ' 0100 '
    height = '20'
    if sitename in df.columns.values:
        sval = df[sitename]
    else:
        sval = [sitename] * len(cval)

    vlist = list(zip(t1, lat, lon, cval, sval))
    # sort by site and then by date.
    vlist.sort(key=lambda x: (x[4],x[0]))

    iii=0
    dt = datetime.timedelta(hours=fillhours)
    for val in vlist:
        runstring += get_runstring(val, height, duration)
 
        done = False
        ddd = iii
        t2 = val[0]
        while not done: 
           t2 +=  dt
           # fill in missing values with 9.
           print('val4', val[4], ddd)
           #print('vlist 4', vlist[ddd+1])
           if fillhours == 0:
              done = True
           elif ddd >= len(vlist)-1:
               done=True   
           # if next sid is not the same then continue
           elif val[4] != vlist[ddd+1][4]:
              done = True
           # if next time is an hour later then continue
           elif t2 == vlist[ddd+1][0]:
              done= True
           # if next time is missing then write a dummmy.
           # with value of -9
           else:
              val2 = (t2, val[1], val[2], -9, val[4]) 
              runstring += get_runstring(val2, height, duration)
           #if ddd >= len(vlist)-1:
           #    done=True   
           print(ddd, iii, len(vlist), val)
           if iii > 1e9:
              break
        iii+=1
            #runstring += val[0].strftime('%Y  %m  %d  %H%M') + duration
            #try:
        #    runstring += "{:.2f}".format(val[1]) + ' ' \
        #                 "{:.2f}".format(val[2]) + ' '
        #except RuntimeError:
        #    print('WARNING1', val[1])
        #    print(val[2])
        #    print(type(val[1]))
        #    print(type(val[2]))
        #    sys.exit()
        #if isinstance(val[4], str):
        #    runstring += "{:.3f}".format(
        #        val[3]) + ' ' + val[4] + ' ' + height + "\n"
        #else:
        #    runstring += "{:.3f}".format(val[3]) + ' ' + \
        #        "{0:d}".format(val[4]) + ' ' + height + "\n"
    with open(dname, 'w') as fid:
        fid.write(runstring)
    return runstring


def get_runstring(val, height, duration):
    runstring = val[0].strftime('%Y  %m  %d  %H%M') + duration
    try:
        runstring += "{:.2f}".format(val[1]) + ' ' \
                     "{:.2f}".format(val[2]) + ' '
    except RuntimeError:
        print('WARNING1', val[1])
        print(val[2])
        print(type(val[1]))
        print(type(val[2]))
        sys.exit()
    if isinstance(val[4], str):
        runstring += "{:.3f}".format(
            val[3]) + ' ' + val[4] + ' ' + height + "\n"
    else:
        runstring += "{:.3f}".format(val[3]) + ' ' + \
            "{0:d}".format(val[4]) + ' ' + height + "\n"
    return runstring

def dropna(df, inplace=True):
    """remove columns which have all Nans.
       TO DO: is this needed?"""
    return df.dropna(axis=1, inplace=inplace)


def get_lhash(df, idn):
    """returns a dictionary with the key as the input column value and the
        value a tuple of (lat, lon)  Useful for getting lat lon locations of
        different sites in a dataframe.
    """
    if 'latitude' in list(df.columns.values):
        dftemp = df.copy()
        pairs = zip(dftemp[idn], zip(dftemp['latitude'], dftemp['longitude']))
        pairs = list(set(pairs))
        lhash = dict(pairs)  # key is facility id and value is name.
        print(lhash)
    return lhash


def summarize(df, verbose=False):
    """prints list of columns. if verbose prints list of unique values in each
    column"""
    columns = list(df.columns.values)
    if verbose:
        for ccc in columns:
            print(ccc)
            print(df[ccc].unique())
    print('-------------------------------')
    for ccc in columns:
        print(ccc)


def latlonfilter(df, llcrnr, urcrnr):
    """
     removes rows from self.df with latitude longitude outside of the box
     described by llcrnr (lower left corner) and urcrnr (upper right corner)
     Parameters
     ----------
       llcrnr : tuple
                lower left corner. (latitude, longitude)
       urcrnr : tuple
                upper right corner (latittude, longitude)
       inplace: boolean
                if TRUE then replaces self.df attribute
       removes rows with latitude longitude outside of the box
       described by llcrnr (lower left corner) and urcrnr (upper right corner)

    """
    lat1 = llcrnr[0]
    lat2 = urcrnr[0]
    lon1 = llcrnr[1]
    lon2 = urcrnr[1]
    df = df[df['latitude'] < lat2]
    df = df[df['latitude'] > lat1]
    df = df[df['longitude'] > lon1]
    df = df[df['longitude'] < lon2]
    return df


def timefilter(df, daterange, inplace=True):
    """removes rows with dates outside of the daterange from self.df
     Parameters
     ----------
     daterange:  tuple
               (datetime, datetime)
     inplace: boolean
               if TRUE then replaces self.df attribute
    """
    df = df[df['time'] >= daterange[0]]
    df = df[df['time'] <= daterange[1]]
    return df

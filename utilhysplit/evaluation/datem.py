# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import numpy as np
import datetime
import os
from subprocess import call
from os import path
import sys
import os
import pandas as pd
from monetio.obs.obs_util import timefilter
# from ashfall_base_iceland import RunParams
# from arlhysplit.runh import date2dir


"""
NAME: datem.py
PRGRMMR: Alice Crawford ORG: NOAA ARL
ABSTRACT: Calls c2datem to extract information from HYSPLIT output files. Adds information to the c2datem output.
CTYPE: source code

FUNCTIONS
writedatem_sh : writes a shell script to run c2datem on the cdump files (HYPSLIT output binary concentration files).
                The shell script will also have lines to concatenate all the output into one text file and add extra information
                to the end of each line (infor about particle size and vertical concentration level which is originally in the file name).

writedatem : writes a datem file which tells c2datem which positions and times to extract concentrations from the cdump file for.

frame2datem(dfile, df,  header_str='Header', writeover=True,\

read_dataA : reads dt

"""


def writedatem_sh(
    cfile_list,
    mult="1e20",
    mdl="./",
    ofile_list=None,
    psizes=[1],
    zlra=[1],
    concat=True,
    concatfile="model.txt",
    add2ra=["none"],
):
    """Writes a .sh file to run the c2datem program on the HYSPLIT output files
       Writes separate file for each particle size. and vertical concentration level.
       psizes : a list of particle sizes to write (particle size indices 1..N)
       zlra   : a list of vertical concentration levels to write (indice 1..N)
       concat : concatenate all the files into one file with name specified by concatfile (model.txt is default).
                information about the particle size and vertical concentration level is added to the end of each line.
       concatfile : name of file with all the data
       add2ra : extra information to add to each line. The extra information is added using sed.
    """
    with open("datem.sh", "w") as fid:
        fid.write("#!/bin/sh \n")
        # removes any existing model.txt file
        fid.write("rm -f " + concatfile + "\n")
        fid.write("MDL=" + mdl + "\n")
        fid.write("mult=" + str(mult) + "\n")
        if ofile_list is None:
            ofile_list = []
            for cfile in cfile_list:
                ofile_list.append(cfile.replace("cdump", "model"))
        outfile_list = []
        for zl in zlra:
            fid.write("zl=" + str(zl) + "\n")
            iii = 0
            for cfile in cfile_list:
                for psz1 in psizes:
                    try:
                        psz = abs(int(psz1))
                    except BaseException:
                        psz1 = 1
                    else:  # if try does not raise an exception then this code is executed.
                        if zl == -1:
                            zstr = "zn1"
                        else:
                            zstr = ".z" + str(zl)
                        outfile = ofile_list[iii] + ".p" + str(int(psz)) + zstr + ".txt"
                        # print iii, outfile
                        outfile_list.append(outfile)
                    # Following block writes line in shell script to run c2datem
                    # the -h0 specifies to not write header lines.
                    # print('WRITING', cfile)
                    fid.write(
                        "$MDL/c2datem -n -h0 -i"
                        + cfile.strip()
                        + " -mdatemfile.txt -o"
                        + outfile
                        + "  -c$mult -z$zl"
                    )
                    # pollutant index select for multiple species
                    fid.write(" -p" + str(int(psz)))
                    fid.write("\n")
                    if concat:
                        temp = cfile.split(".")
                        if add2ra[0] != "none":
                            a2l = " " + add2ra[iii] + " " + str(psz) + " " + str(zl)
                            # add info to end of line and add to file.
                            fid.write(
                                "sed 's/$/"
                                + a2l
                                + "/' "
                                + outfile
                                + " >> "
                                + concatfile
                                + "\n"
                            )
                iii += 1
        fid.write("if [ ! -s model.txt ]\n")
        fid.write("then\n")
        fid.write("rm -f model.txt\n")
        fid.write("fi\n")

    return outfile_list


def frame2datem2(
    dfile,
    df,
    fillval=-999,
    header_str="Header",
    writeover=True,
    cnames=["date", "duration", "lat", "lon", "obs", "vals", "sid", "altitude"],
):
    """converts a pandas dataframe with columns names by cnames (date, duration, lat, lon, obs, vals, sid, altitude)
       to a text file in datem format.
       date should be a datetime object.
       duration should be a string format HHMM (TODO- make this more flexible?)
       lat - latitude, float
       lon - longitude, float
       obs - value of observation, float
       vals - modeled value, float
       sid  - station id, int or string
       altitude - float 
       extra columns can be written.

    """
    dlen = len(cnames)
    blen = 8
    iii = 0
    df = df.fillna(fillval)
    if writeover:
        with open(dfile, "w") as fid:
            fid.write(header_str + " (obs then model) " + "\n")
            fid.write(
                "Num Site Lat Lon  Yr Mo Da Hr Meas Calc "
            )
            #for iii in np.arange(blen,blen+(dlen-blen)):
            #    fid.write(cnames[iii] +  '  ')
            fid.write('\n')
    iii=0 
    with open(dfile, "a") as fid:
        for index, row in df.iterrows():
            fid.write(str(iii) + ' ')
            write_site(row[cnames[6]], fid)
            fid.write("%8.3f  %8.3f" % (row[cnames[2]], row[cnames[3]]))
            fid.write(row[cnames[0]].strftime(" " + "%Y %m %d %H%M") + " ")
            fid.write("%8.4f  %8.4f " % (row[cnames[4]], row[cnames[5]]))
            iii+=1
            fid.write('\n')

def write_vals(val, fid):
        if isinstance(row[cnames[iii]], int):
            fid.write("%12i" % (row[cnames[iii]]))
        elif isinstance(row[cnames[iii]], float):
            fid.write("%12.2f" % (row[cnames[iii]]))
        elif isinstance(row[cnames[6]], str):
            fid.write("%12s  " % (row[cnames[iii]]))

def write_site(site, fid):
        if isinstance(site,  int):
            fid.write("%12i" % site)
        elif isinstance(site, float):
            fid.write("%12i" % site)
        elif isinstance(site, str):
            fid.write("%12s  " % site)
        else:
            print("WARNING frame2datem function: not printing station id")
        fid.write(" ")

def frame2datem(
    dfile,
    df,
    fillval=-999,
    header_str="Header",
    writeover=True,
    cnames=["date", "duration", "lat", "lon", "obs", "vals", "sid", "altitude"],
):
    """converts a pandas dataframe with columns names by cnames (date, duration, lat, lon, obs, vals, sid, altitude)
       to a text file in datem format.
       date should be a datetime object.
       duration should be a string format HHMM (TODO- make this more flexible?)
       lat - latitude, float
       lon - longitude, float
       obs - value of observation, float
       vals - modeled value, float
       sid  - station id, int or string
       altitude - float 
       extra columns can be written.

    """
    dlen = len(cnames)
    blen = 8
    iii = 0
    df = df.fillna(fillval)
    if writeover:
        with open(dfile, "w") as fid:
            fid.write(header_str + " (obs then model) " + "\n")
            fid.write(
                "year mn dy shr dur(hhmm) LAT LON  ug/m2 ug/m2 site_id  height "
            )
            for iii in np.arange(blen,blen+(dlen-blen)):
                fid.write(cnames[iii] +  '  ')
            fid.write('\n')
    with open(dfile, "a") as fid:
        for index, row in df.iterrows():
            fid.write(row[cnames[0]].strftime("%Y %m %d %H%M") + " ")
            fid.write(str(row[cnames[1]]) + " ")
            fid.write("%8.3f  %8.3f" % (row[cnames[2]], row[cnames[3]]))
            fid.write("%8.4f  %8.4f " % (row[cnames[4]], row[cnames[5]]))
            if isinstance(row[cnames[6]], int):
                fid.write("%12i" % (row[cnames[6]]))
            elif isinstance(row[cnames[6]], float):
                fid.write("%12i" % (row[cnames[6]]))
            elif isinstance(row[cnames[6]], str):
                fid.write("%12s  " % (row[cnames[6]]))
            else:
                print("WARNING frame2datem function: not printing station id")
            fid.write("%7.2f " % (row[cnames[7]]))
            for iii in np.arange(blen,blen+(dlen-blen)):
                if isinstance(row[cnames[iii]], int):
                    fid.write("%12i" % (row[cnames[iii]]))
                elif isinstance(row[cnames[iii]], float):
                    fid.write("%12.2f" % (row[cnames[iii]]))
                elif isinstance(row[cnames[6]], str):
                    fid.write("%12s  " % (row[cnames[iii]]))


            fid.write('\n')
            # fid.write(str(row[cnames[7]]) + '\n')


def writedatem(dfile, stationlocs, sample_start, sample_end, stime, height=" 10"):
    """writes a dummy station datem file which has times for each station location.
       This file is used by c2datem to determine what concentration values to pull from the cdump files.
       stationlocs is a list of (lat,lon) tuples.

       If the -z option in c2datem is set to -1 then the height
       indicates which level will be used. It is the actual height in meters, not
       the index of the height lev
el.

       outputs 1 in the measurement concentration and sourceid columns.
    """
    iii = 0
    with open(dfile, "w") as fid:
        fid.write("DOE ASHFALL PROJECT\n")
        fid.write("year mn dy shr dur(hhmm) LAT LON g/m2  site_id  height \n")
        for iii in range(0, len(stationlocs)):
            sdate = sample_start
            while sdate < sample_end:
                fid.write(sdate.strftime("%Y %m %d %H%M") + " ")
                fid.write(str(stime[iii]) + "00 ")
                fid.write(
                    "{:0.3f}".format(stationlocs[iii][0])
                    + " "
                    + "{:0.3f}".format(stationlocs[iii][1])
                    + " "
                )
                fid.write("1 ")
                fid.write("1 ")
                # fid.write(height + '\n')
                fid.write(height + "\n")
                sdate += datetime.timedelta(hours=stime[iii])


def read_dataA(fname):
    """
    reads merged data file output by statmain.
    outputs a  dataframe with columns
    date, sid, lat, lon, obs, model
    """
    colra = ["Num", "sid", "lat", "lon", "year", "month", "day", "hour", "obs", "model"]
    dtp = {"year": int, "month": int, "day": int, "hour": int}
    try:
        datem = pd.read_csv(
            fname, names=colra, header=None, delimiter=r"\s+", dtype=dtp, skiprows=2
        )
    except pd.errors.ParserError:
        datem = pd.DataFrame()

    if not datem.empty:
        datem["minute"] = datem["hour"] % 100
        datem["hour"] = datem["hour"] / 100
        datem["hour"] = datem["hour"].astype(int)

        def getdate(x):
            return datetime.datetime(
                int(x["year"]),
                int(x["month"]),
                int(x["day"]),
                int(x["hour"]),
                int(x["minute"]),
            )

        try:
            datem["date"] = datem.apply(getdate, axis=1)
        except BaseException:
            print("EXCEPTION", fname)
            print(datem[0:10])
            sys.exit()
        datem.drop(["year", "month", "day", "hour", "minute"], axis=1, inplace=True)
    #print('DATEM FILE loaded:', fname)
    #print(datem[0:2])
    return datem


def read_datem_file(
    fname,
    dummy=False,
    verbose=False,
    header=None,
    colra=[
        "year",
        "month",
        "day",
        "hour",
        "duration",
        "meas_lat",
        "meas_lon",
        "vals",
        "stationid",
        "sourceid",
        "level",
        "thickness",
    ],
):
    """ Reads a datem file and returns a dataframe with colums described by colra
       colra : should have columns for 'year' 'month' 'day' 'hour'. 'hour' column should be in hhmm format.
       fname :
       zlevs :
       sdate :

       returns pandas dataframe.
    """
    dtp = {"year": int, "month": int, "day": int, "hour": int}
    try:
        datem = pd.read_csv(fname, names=colra, header=header, delimiter=r"\s+", dtype=dtp)
    except pd.errors.ParserError:
        print('WARNING: could not read ', fname)
        datem = pd.DataFrame()
    datem.columns = colra
    datem["minute"] = datem["hour"] % 100
    datem["hour"] = datem["hour"] / 100
    datem["hour"] = datem["hour"].astype(int)

    def getdate(x):
        return datetime.datetime(
            int(x["year"]),
            int(x["month"]),
            int(x["day"]),
            int(x["hour"]),
            int(x["minute"]),
        )

    datem["date"] = datem.apply(getdate, axis=1)
    datem.drop(["year", "month", "day", "hour", "minute"], axis=1, inplace=True)
    #colrb = [x for x in colra if x not in ["year","month","day","hour","minute"]]
    #colrb = ['date'].extend(colrb)
    #datem.columns = colrb 
    return datem



def write_datem(df,
                obscolumn='obs',
                dname='datemfile.txt',
                sitename='1',
                info=None,
                drange=None,
                fillhours=1,
                verbose=False):
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
          # print('val4', val[4], ddd)
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
           #print(ddd, iii, len(vlist), val)
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
        if verbose:
           print('writing ', dname)
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

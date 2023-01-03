#!/n-home/alicec/anaconda/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
#import numpy as np
#import string
import datetime
import os

#import subprocess
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr

from utilhysplit import hcontrol


"""
PRGMMR: Alice Crawford  ORG: ARL  
PYTHON 2.7
This code written at the NOAA  Air Resources Laboratory
UID: r102
CTYPE: source code
ABSTRACT: manages the xtrct_stn program and outputs.

CLASSES


"""


class VmixingRun:
    def __init__(
        self,
        fname,
        cname="CONTROL",
        cdir="./",
        pid=None,
        kbls=1,
        kblt=2,
        kmixd="",
        kmix0="",
        cameo=2,
        tkemin=None,
        verbose=True,
    ):
        self.control = hcontrol.HycsControl(fname=cname, rtype="vmixing")
        self.pid = self.control.replace("CONTROL.", "")
        self.kbls = kbls  # 1 fluxes  #2 wind/temp profile
        self.kblt = kblt  # 1 BJ #2 KC #3 TKE.
        self.kmixd = kmixd  # 0 Input #1 temperature #2 TKE. #3 RI#
        self.kmix0 = kmix0  # 150 is default. minimum mixing depth.
        self.cameo = cameo  # 0 no #1 yes #2 yes + wdir
        self.tkemin = tkemin
        self.woption = woption  # output extra file
        self.cdir = cdir

    def readcontrol(self):
        self.control.read()

    def writecontrol(self, cname=None, cdir=None):
        if not cdir:
            cdir = self.cdir
        if cname:
            self.control.rename(cname, working_directory=cdir)
        self.control.write()

    def assign_pid(self, pid):
        self.pid = pid
        self.control.rename("CONTROL." + str(pid))
        return -1

    def make_runstr(self, hdir):
        rstr = hdir
        if rstr[-1] != "/":
            rstr.append("/")
        rstr += "vmixing "
        if self.pid:
            rstr += "-p" + str(self.pid)
            rstr += "-s" + str(self.kbls)
            rstr += "-t" + str(self.kblt)
            rstr += "-a" + str(self.cameo)
            if self.kmixd:
                rstr += "-d" + str(self.kmixd)
            if self.kmix0:
                rstr += "-l" + str(self.kmix0)
            if tkemin:
                rstr += "-m" + str(self.tkemin)
            rstr += "-w" + str(self.woption)
        return rstr


def find_stability_files(tdir, suffix=None, szmin=0):
    rlist = []
    for (dl, dirnames, filenames) in os.walk(tdir):
        for filename in filenames:
            if filename[0:9] == "STABILITY":
                tag = filename.replace("STABILITY.", "")
                tag = tag.replace(".txt", "")
                if not suffix or suffix == tag:
                    if os.path.getsize(os.path.join(dl, filename)) > szmin:
                        message = "MESSAGE." + tag + ".txt"
                        rlist.append((dl, filename, message))
    return rlist


def find_vmix_files(tdir, suffix=None, szmin=0):
    rlist = []
    for (dl, dirnames, filenames) in os.walk(tdir):
        for filename in filenames:
            if filename[0:4] == "wmix":
                tag = filename.replace("wmix.", "")
                tag = tag.replace(".txt", "")
                if not suffix or suffix == tag:
                    if os.path.getsize(os.path.join(dl, filename)) > szmin:
                        message = "MESSAGE." + tag + ".txt"
                        rlist.append((dl, filename, message))
    return rlist


class Wclass:
    """
    helper class for WmixData.
    Helps ensure that data from different files are concatenated into
    an xarray correctly.
    """

    def __init__(self, warray, ens, sid):
        self.wf = warray  # xarray DataArray
        self.ens = str(ens)  # string
        self.sid = str(sid)  # string
        self.d1 = warray.date.values[0]
        self.d2 = warray.date.values[-1]

    def __str__(self):
        rstr = str(self.ens) + "\n"
        rstr += str(self.sid) + "\n"
        rstr += str(self.d1) + " " + str(self.d2) + "\n"
        return rstr

    def __lt__(self, other):
        """
        lt and eq so list of objects can be sorted.
        First sort by sourceid then ensemble then date.

        """
        if self.sid != other.sid:
            return self.sid < other.sid
        elif self.ens != other.ens:
            return self.ens < other.ens
        elif self.d1 != other.d1:
            return self.d1 < other.d1
        elif self.d2 != other.d2:
            return self.d2 < other.d2

    def __eq__(self, other):
        t1 = self.ens == other.ens
        t2 = self.sid == other.sid
        t3 = self.d1 == other.d1
        t4 = self.d2 == other.d2
        return t1 and t2 and t3 and t4

    def dates_equal(self, other):
        t3 = self.d1 == other.d1
        t4 = self.d2 == other.d2
        return t3 and t4


def plot_variance(dset, levels=None, ax=None):
    nset = dset**0.5
    cb = plot_ens_general(nset, levels, ax)
    return cb

## TODO trouble importing get_vmix_colors
#def plot_ens_general(dset, levels, ax=None):
#    if not ax:
#        ax = plt.gca()
    #from utilhysplit.evaluation.armdata import get_vmix_colors

#    lev, norm, cmap = get_vmix_colors(levels)

#    cb = ax.pcolormesh(dset.date, dset.dim_0, dset, cmap=cmap, norm=norm)
    # plt.colorbar(cb)
#    plt.tight_layout()
    # plt.show()
#    return cb


def plot_ens_min(dset, levels=None, ax=None):
    dset = dset**0.5
    nset = dset.min(dim="ens")
    plot_ens_general(nset, levels, ax)


def plot_ens_max(dset, levels=None, ax=None):
    dset = dset**0.5
    nset = dset.max(dim="ens")
    plot_ens_general(nset, levels, ax)


def plot_mean(dset, levels=None, ax=None):
    dset = dset**0.5
    nset = dset.mean(dim="ens")
    plot_ens_general(nset, levels, ax)


def plot_ens_var(dset, levels=None, ax=None):
    dset = dset**0.5
    nset = dset.var(dim="ens")
    plot_ens_general(nset, levels, ax)


class WmixData:

    """
    use add_data to add files.
    after adding all files use
    concat to create the xarray DataArray.
    """

    def __init__(self, century=2000, verbose=True):
        """fname : name of file output by xtrct_stn
        valra : list of values that are in fname
        century : fname lists year by last two digits only. century is needed to process date.
        """
        self.units = None
        self.df = pd.DataFrame()
        self.data = xr.DataArray()
        self.wlist = []  # list of Wclass objects
        self.enslist = []
        self.sidlist = []

    def concat_data(self):
        # sort list.
        self.wlist.sort()
        # remove duplicates from lists
        self.sidlist = list(set(self.sidlist))
        self.enslist = list(set(self.enslist))
        sidlist = []
        for sid in self.sidlist:
            # get list with all same site id.
            sublist = [x for x in self.wlist if x.sid == sid]
            # get list which all have same site and ensemble data.
            # concatenate along time axis
            enslist = []
            for ens in self.enslist:
                timelist = [x for x in sublist if x.ens == ens]
                timelist.sort()
                xrtime = xr.concat([x.wf for x in timelist], "date")
                enslist.append(xrtime)
            # if more than one ensemble concat along ensemble axis.
            if len(enslist) > 1:
                xrens = xr.concat(enslist, "ens")
            else:
                xrens = enslist[0]
            sidlist.append(xrens)
        # if more than one ensemble concat along ensemble axis.
        if len(sidlist) > 1:
            xrall = xr.concat(sidlist, "sid")
        else:
            xrall = sidlist[0]
        self.data = xrall
        return xrall

    def add_data(self, fname, century=2000, sid="None", ens="None"):
        df = self.readfile(fname, century)
        print("here", type(df))
        wf = self.make_wmix_df(df)
        print("here", type(wf))
        wf = wf.expand_dims("ens")
        wf["ens"] = [ens]
        wf = wf.expand_dims("sid")
        wf["sid"] = [sid]
        self.wlist.append(Wclass(wf, ens, sid))
        self.enslist.append(ens)
        self.sidlist.append(sid)
        # if self.wf.isnull():
        #    self.wf = wf
        # else:
        #    self.wf = xr.concat([self.wf,wf],"date")
        return wf

    def get_heights(self, messagefile):
        from utilhysplit import message_parse

        mp = message_parse.HysplitMessageFile(messagefile)
        levs = mp.get_levels()
        return levs

    def get_message_file(self, fname):
        fdir = fname.replace(fname.split("/")[-1], "")
        tag = fname.split("/")[-1].replace("wmix.", "")
        tag = tag.replace(".txt", "")
        message = "MESSAGE." + tag + ".txt"
        if os.path.isfile(os.path.join(fdir, message)):
            return os.path.join(fdir, message)
        else:
            return None

    def add_heights(self, df, messagefile):
        levs = self.get_heights(messagefile)
        print("heights are ", levs)
        ctemp = list(df.columns.values[0:10])
        ctemp.extend(levs)
        if len(ctemp) == len(df.columns.values):
            df.columns = ctemp
        else:
            print("add_heights: wrong length for header")
        return df

    @staticmethod
    def try_float(x):
        rline = x.split()

        def tofloat(y):
            try:
                rval = float(y)
            except:
                rval = y
            return rval

        rline = [tofloat(n) for n in rline]
        return rline

    def readfile(self, fname, century=2000):
        data_row = 3  # row which data starts on.
        fid = open(fname, "r")
        # fra = [x.split() for x in fid.readlines()]
        fra = [self.try_float(x) for x in fid.readlines()]
        data = fra[data_row:]
        headers = fra[data_row - 1]
        df = pd.DataFrame(data)
        headers.extend(df.columns.values[len(headers) :])
        df.columns = headers
        if self.get_message_file(fname):
            df = self.add_heights(df, self.get_message_file(fname))
        # else:
        #    headers.extend(df.columns.values[len(headers):])
        #    df.columns = headers
        # print(df[0:10])
        df["date"] = df.apply(lambda row: self.row2date(row, century), axis=1)
        return df

    def make_wmix_df(self, df):
        # dropra = ['PSQ','JDAY','YR','MO','DA','HR','MN','PSQ','U*','UMIX','VMIX','WMIX']
        dropra = [
            "PSQ",
            "JDAY",
            "YR",
            "MO",
            "DA",
            "HR",
            "MN",
            "PSQ",
            "U*",
            "UMIX",
            "VMIX",
        ]
        wf = df.drop(dropra, axis=1)
        wf = wf.transpose()
        print("here", wf.values[0][0], type(wf.values[0][0]))
        wf.columns = wf.loc["date"]
        wf = wf.drop("date", axis=0)
        wf = wf.astype(float)
        wx = xr.DataArray(wf)
        return wx
        # wx = xr.DataArray(wf)
        # return wx

    @staticmethod
    def row2date(row, century):
        """get date from a row"""
        vdate = datetime.datetime(
            int(row["YR"]) + century,
            int(row["MO"]),
            int(row["DA"]),
            int(row["HR"]),
            int(row["MN"]),
        )
        print(type(vdate))
        return vdate


class VmixingData:
    """
    add_data
    make_dummies (NOT FUNCTIONAL)
    readfile
    """

    def __init__(self, century=2000, verbose=True):
        """fname : name of file output by xtrct_stn
        valra : list of values that are in fname
        century : fname lists year by last two digits only. century is needed to process date.
        """
        self.units = None
        self.df = pd.DataFrame()

    def add_data(self, fname, vdir="./", century=2000, verbose=False, sid=None):
        df = self.readfile(fname, vdir, century, verbose)
        if sid:
            df["sid"] = sid
        if self.df.empty:
            self.df = df
        else:
            self.df = pd.concat([self.df, df], axis=0)
            # print(self.df)
            # import sys
            # sys.exit()
        return self.df

    def make_dummies(self, data_ra=[-999]):
        """instead of running,  write a dummy file like the one vmixing would write.
        Used for testing.
        """
        # sdate = datetime.datetime()
        # dt = datetime.timedelta(hour=1)
        # iii=1
        # with open(self.fname) as fid:
        #     fid.write(str(iii) + sdate.strftime(" %y %m %d %h"))
        #     iii+=1
        return -1

    def get_location(self, head1):
        # vmixing doesn't always print a space between lat and lon
        head1 = head1.replace("-", " -")
        temp1 = head1.split()
        lat = float(temp1[0])
        lon = float(temp1[1])
        met = temp1[2]
        return lat, lon, met

    def parse_header(self, head2, head3):
        temp2 = head2.split()
        temp3 = head3.split()
        cols = ["date"]
        units = ["utc"]
        cols.extend(temp2[6:])
        if "Total" in cols and "Cld" in cols:
            cols.remove("Total")
        units.extend(temp3)
        return cols, units

    def save(self, outname):
        self.df.to_csv(outname)

    def readfile(self, fname, vdir="./", century=2000, verbose=False):
        """Reads file and returns True if the file exists.
        returns False if file is not found"""
        df = pd.DataFrame()
        if os.path.isfile(vdir + fname):
            if verbose:
                print("Adding", vdir + fname)
            data = []
            with open(vdir + fname, "r") as fid:
                head1 = fid.readline()
                head2 = fid.readline()
                head3 = fid.readline()
                try:
                    lat, lon, met = self.get_location(head1)
                except:
                    print("problem with vmixing file ", fname, vdir)
                    print("header ", head1)
                    return df
                    # sys.exit()
                cols, units = self.parse_header(head2, head3)
                data = [process_line(x, century) for x in fid.readlines()]
                # data = map(process_line(x, century), fid.readlines())
                # for line in fid.readlines():
                # get the date for the line
                #    temp = line.split()
                #    try:
                #        vals = [line2date(line, century)]
                #    except:
                #        return False
                #    temp2 = []
                #    for val in temp[6:]:
                #        try:
                #            temp2.append(float(val))
                #        except:
                #            temp2.append(val)
                #    vals.extend(temp2)
                #    data.append(vals)
            df = pd.DataFrame.from_records(data)
            df.columns = cols
            df["latitude"] = lat
            df["longitude"] = lon
            df["met"] = met
            self.units = zip(cols, units)
        else:
            if verbose:
                print("Cannot Find ", vdir + fname)
        return df


def process_line(line, century):
    # get the date for the line
    temp = line.split()
    try:
        vals = [line2date(line, century)]
    except:
        return False
    temp2 = []
    for val in temp[6:]:
        try:
            temp2.append(float(val))
        except:
            temp2.append(val)
    vals.extend(temp2)
    return vals


def line2date(line, century):
    """get date from a line in output file and return datetime object"""
    temp = line.strip().split()
    year = int(temp[1]) + century
    month = int(temp[2])
    day = int(temp[3])
    hour = int(temp[4])
    minute = int(temp[5])
    vdate = datetime.datetime(year, month, day, hour, minute)
    return vdate

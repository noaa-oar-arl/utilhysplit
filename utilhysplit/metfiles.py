# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import datetime
from os import path

"""
NAME: svhy.py
PRGMMR: Alice Crawford  ORG: ARL
This code written at the NOAA  Air Resources Laboratory
ABSTRACT: choosing met files for HYSPLIT control file
"""


def getmetfiles(strfmt, sdate, runtime):
    """
    INPUTS
    strfmt : string
    sate : datetime object
    runtime : integer
    RETURNS :
    list of tuples (directory, filename)
    """
    mfiles = MetFiles(strfmt)
    return mfiles.get_files(sdate, runtime)


class MetFiles:
    """
    Class for finding metfiles used in HYSPLIT CONTROL file.
    """

    def __init__(self, strfmt, hours=-1, verbose=False):
        """
        INPUTS
        strfmt : string
        hours : integer
        verbose : boolean

        Attributes
        mdt : datetime timedelta object.
              how much time one met file spans in hours.
        verbose : boolean.
        strfmt : string

        """
        self.verbose = verbose
        self.strfmt = strfmt
        if hours >= 0:
            self.mdt = datetime.timedelta(hours=hours)
        else:
            self.mdt = -1

    def set_mdt(self, hours):
        """
        set self.mdt.
        can be used to over-ride algorithm in find_mdt method.
        """
        self.mdt = datetime.timedelta(hours=hours)

    def get_files(self, sdate, runtime):
        """
        sdate : datetime object. start date.
        runtime : integer. hours of runtime.
        RETURNS :
        list of tuples (directory, filename)
        """
        nlist = self.make_file_list(sdate, runtime)
        return self.process(nlist)

    @staticmethod
    def sub_handle_hour(sdate):
        """
        INPUTS
        sdate : datetime.datetime object
        returns
        date with 0 hour.
        """
        year = sdate.year
        month = sdate.month
        day = sdate.day
        #hour = sdate.hour
        testdate = datetime.datetime(year, month, day, 0)
        return testdate

    def handle_hour(self, sdate):
        """
        Handles met files which span less than one day.
        """
        # testdate is sdate with hour=0
        testdate = self.sub_handle_hour(sdate)
        done = False
        imax = 100
        iii = 0
        while not done:
            # add self.mdt to the testdate
            newdate = testdate + self.mdt
            # if not greater than the input date
            # then keep adding.
            if newdate < sdate:
                testdate = newdate
            # else return this date.
            else:
                done = True
            iii += 1
            if iii > imax:
                done = True
        return testdate

    def find_mdt(self, testdate):
        """
        # finds time spacing between met files by
        # seeing which spacing produces a new file name.
        # this can be problematic if there are missing files.
        # testdate = datetime.datetime(2010,1,1)
        """
        if "%H" in self.strfmt:
            mdtlist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
            testdate = self.sub_handle_hour(testdate)
        else:
            mdtlist = [1, 24, 24 * 7, 24 * 31, 24 * 356]

        file1 = testdate.strftime(self.strfmt)
        done = False
        iii = 0
        while not done:
            dttt = datetime.timedelta(hours=mdtlist[iii])
            dtt2 = testdate + dttt
            file2 = dtt2.strftime(self.strfmt)
            if file2 != file1 and path.isfile(file2):
                done = True
            iii += 1
            if iii >= len(mdtlist):
                done = True
        ##AMC ADDED TEMP
        #dttt=datetime.timedelta(hours=24)
        return dttt

    def parse_week(self, edate):
        """
        helps handle weekly files.
        """
        # used if week is in the strfmt (mostly for gdas1)
        self.mdt = datetime.timedelta(hours=7 * 24)
        temp = edate.strftime(self.strfmt)
        day = int(edate.strftime("%d"))
        # week 1 is
        if day < 8:
            temp = temp.replace("week", "w1")
        # week 2 is 8 through 14
        elif day < 15:
            temp = temp.replace("week", "w2")
        elif day < 22:
            temp = temp.replace("week", "w3")
        # week 4 is 22 through  28
        elif day < 29:
            temp = temp.replace("week", "w4")
        else:
            temp = temp.replace("week", "w5")
        return temp

    def make_file_list(self, sdate, runtime):
        """
        INPUTS
        sdate : datetime.datetime ojbect
        runtime : int (hours of runtime)
        """
        nlist = []
        sdate = sdate.replace(tzinfo=None)
        if not isinstance(self.mdt, datetime.timedelta):
            self.mdt = self.find_mdt(sdate)
        # handle backwards runs. by switching sdate and edate
        if runtime < 0:
            runtime = abs(runtime)
            end_date = sdate
            sdate = end_date - datetime.timedelta(hours=runtime)
        else:
            end_date = sdate + datetime.timedelta(hours=runtime)
        done = False
        # self.verbose=True
        if "%H" in self.strfmt:
            sdate = self.handle_hour(sdate)
        edate = sdate
        if self.verbose:
            print("GETMET", sdate, edate, end_date, runtime, self.mdt)
        zzz = 0
        while not done:
            if "week" in self.strfmt:
                temp = self.parse_week(edate)
            else:
                temp = edate.strftime(self.strfmt)
            edate = edate + self.mdt
            if not path.isfile(temp):
                temp = temp.lower()
            if not path.isfile(temp):
                print("WARNING", temp, " meteorological file does not exist")
                temp = self.find_alt(temp)
                print("REPLACE with", temp)
            if temp != "None":
                if temp not in nlist:
                    nlist.append(temp)
            if edate > end_date:
                done = True
            if zzz > 50:
                done = True
            zzz += 1
        return nlist

    @staticmethod
    def find_alt2(temp):
        """
        Finds alternative file
        """
        temp1 = temp.split(".")
        date = temp1[2]
        rval = "/ready_archives/nams/" + date + "_hysplit.t00z.namsa"
        if not path.isfile(rval):
            rval = "None"
        return rval

    @staticmethod
    def find_alt(temp):
        """
        Finds alternative file
        """
        alist = ["t03z", "t09z", "t12z", "t15z", "t18z"]
        blist = ["t03z", "t09z", "t12z", "t15z", "t18z"]
        # alist = ["t09z"]
        bbb = "t03z"
        rval = "None"
        for bbb in blist:
            for aaa in alist:
                if aaa in temp:
                    rval = temp.replace(aaa, bbb)
            if path.isfile(rval):
                return rval
        return "None"

    @staticmethod
    def process(nlist):
        """
        # convert full path to
        # list of directories and list of filenames
        # and then zips the lists to return list of tuples.
        """
        mfiles = []
        mdirlist = []
        for temp in nlist:
            siii = [x for x, char in enumerate(temp) if char == "/"]
            siii = siii[-1]
            fname = temp[siii + 1 :]
            mdir = temp[0 : siii + 1]
            mfiles.append(fname)
            mdirlist.append(mdir)
        return list(zip(mdirlist, mfiles))

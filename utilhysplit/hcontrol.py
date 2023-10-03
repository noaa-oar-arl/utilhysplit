# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import datetime
import logging
import os
from os import path

"""
PGRMMR: Alice Crawford ORG: ARL/CICS
PYTHON 3
ABSTRACT: classes and functions for creating HYSPLIT control and setup files.

   CLASSES
   HycsControl: class for reading / writing a HYSPLIT dispersion run  control
                file
   Helper classes for HycsControl class
           ControlLoc: release location for  CONTROL file.
           Species: class representing pollutant properties as defined in
                    CONTROL file
           ConcGrid: class representing concentration grid as defined in
                     CONTROL file
   NameList: class for writing SETUP.CFG file

FUNCTIONS
   writelanduse - writes ASCDATA.CFG file.
"""

# 2023 Jul 14 (AMC) added type hinting to ConcGrid class
# 2023 Sep 26 (AMC) added kwargs to NameList write method so can write without header.
#                   changed gem to a kwarg.
# 2023 Sep 26 (AMC) added a line_ending property to NameList class so can write withouth ',' at line ends.

logger = logging.getLogger(__name__)


def writeover(name, overwrite, query, verbose=False):
    """
    checks if file already exits.
    Returns 1 if file should be written.
    Returns -1 if file should not be written.

    Inputs
    name : str
           filename
    overwrite: boolean
           if True then will 1.
    query : boolean
           if True then will ask for user input before overwriting a file
    verbose : boolean
           if True will write messages

    Outputs
    rval : int
         Value is 1 if file should be written.
         Value is -1 if file should not be written.
    """
    rval = 1
    if path.isfile(name):
        if verbose:
            print("file already exists " + name)
        # fexists = True
        if query:
            istr = " Press y to overwrite file \n"
            istr += " Press any other key to continue without overwriting "
            answer = input(istr)
            if answer.strip().lower() != "y":
                overwrite = False
            else:
                overwrite = True
        if overwrite:
            if verbose:
                print("overwriting existing file")
            rval = 1
        else:
            if verbose:
                print("Continuing without overwriting file")
            rval = -1
    return rval


def writelanduse(
    landusedir,
    working_directory="./",
    overwrite=True,
    query=False,
    verbose=False,
):
    """writes an ASCDATA.CFG file in the outdir. The landusedir must
       be the name of the directory where the landuse files are located.
    Parameters
    ----------
    landusedir : string
    working_directory : string
    overwrite : boolean
    query : boolean
    verbose : boolean

    Returns
    ----------
    rval : int
    """
    rval = writeover(
        working_directory + "ASCDATA.CFG", overwrite, query, verbose=verbose
    )
    if rval == 1:
        with open(path.join(working_directory, "ASCDATA.CFG"), "w") as fid:
            fid.write("-90.0  -180.0 \n")
            fid.write("1.0    1.0    \n")
            fid.write("180    360    \n")
            fid.write("2 \n")
            fid.write("0.2 \n")
            fid.write(landusedir + " \n")
            if not path.isdir(landusedir) and verbose:
                print("writelanduse function WARNING:")
                print("landuse directory does not exist", landusedir)
    return rval


class ConcGrid:
    """concentration grid as defined by 10 lines in the HYSPLIT concentration
      CONTROL file.
    Methods
    -------
    __init__
    get_nlev
    set_annotate
    describe
    __str__
    typestrg
    definition

    Attributes
    ----------
    centerlat
    centerlon
    latdiff
    londiff
    latspan
    lonspan
    outdir
    outfile
    nlev
    sample_start
    sample_stop
    sampletype
    interval
    annotate : boolean
    """

    def __init__(
        self,
        name,
        levels=None,
        centerlat: float = 0.0,
        centerlon: float = 0.0,
        latdiff: float = -1.0,
        londiff: float = -1.0,
        latspan: float = 90.0,
        lonspan: float = 360.0,
        outdir: str = "./",
        outfile: str = "cdump",
        nlev: int = -1,
        sample_start: str = "00 00 00 00 00",
        sample_stop: str = "00 00 00 00 00",
        sampletype: int = 0,
        interval: tuple = (-1, -1),
    ):
        """
        Parameters
        ----------
        name : string
        levels : list of floats/ints
        center_lat : float
        center_lon : float
        interval is
        sample type : integer (0 is average)

        Return
        -------
        None
        """
        self.name = name
        if levels is None:
            self.levels = []
        else:
            self.levels = levels
        self.centerlat = centerlat
        self.centerlon = centerlon
        self.latdiff = latdiff
        self.londiff = londiff
        self.latspan = latspan
        self.lonspan = lonspan
        self.outdir = outdir
        self.outfile = outfile
        self.nlev = nlev
        # string (could be changed to datetime)
        self.sample_start = sample_start
        self.sample_stop = sample_stop  # string (could be changed to datetime)
        self.sampletype = sampletype
        self.interval = interval
        self.get_nlev()
        self.annotate = False


    @property
    def outdir(self):
        return self._outdir

    @outdir.setter
    def outdir(self,outdir):
        if outdir[-1] != "/":
           outdir += "/"
        self._outdir = outdir 

    def copy(self):
        return ConcGrid(
            self.name,
            self.levels,
            self.centerlat,
            self.centerlon,
            self.latdiff,
            self.londiff,
            self.latspan,
            self.lonspan,
            self.outdir,
            self.outfile,
            self.nlev,
            self.sample_start,
            self.sample_stop,
            self.sampletype,
            self.interval,
        )

    def set_annotate(self, on=True):
        """
        write
        """
        self.annotate = on

    def get_nlev(self):
        """
        computes self.nlev from self.levels
        """
        # if self.nlev == -1:
        self.nlev = len(self.levels)

    def set_levels(self, levels):
        self.levels = levels
        self.get_nlev()

    def __str__(self):
        """string method will output ten lines suitable for inserting into a
        HYSPLIT control file"""
        pnotes = self.annotate
        note = ""
        if pnotes:
            note = "  #Concentration Grid Center (latitude longitude)"
        returnstr = (
            "{:.4f}".format(self.centerlat)
            + " "
            + "{:.4f}".format(self.centerlon)
            + note
            + "\n"
        )
        if pnotes:
            note = "  #Concentration grid spacing (degrees latitude longitude)"
        returnstr += str(self.latdiff) + " " + str(self.londiff) + note + "\n"
        if pnotes:
            note = "  #Concentration grid span (degrees latitude longitude)"
        returnstr += (
            "{:.2f}".format(self.latspan)
            + " "
            + "{:.2f}".format(self.lonspan)
            + note
            + "\n"
        )
        if pnotes:
            note = "  #Directory to write concentration output file"

        if self.outdir[-1] != "/":
            self.outdir += "/"
        returnstr += self.outdir + note + "\n"
        if pnotes:
            note = "  #Filename for trajectory output"
        returnstr += self.outfile + note + "\n"
        if pnotes:
            note = "  #Number of vertical levels for concentration grid"
        returnstr += str(self.nlev) + note + "\n"
        if pnotes:
            note = "  #List of vertical levels for concentration grid"
        for lev in self.levels:
            returnstr += str(lev) + " "
        returnstr += note + "\n"
        if pnotes:
            note = "  #Sampling start time of concentration grid"
        returnstr += self.sample_start + note + "\n"
        if pnotes:
            note = "  #Sampling stop time of concentration grid"
        returnstr += self.sample_stop + note + "\n"
        if pnotes:
            note = "  # " + self.typestr()
        returnstr += "{:02.0f}".format(self.sampletype) + " "
        returnstr += "{:02.0f}".format(self.interval[0]) + " "
        returnstr += "{:02.0f}".format(self.interval[1]) + " "
        returnstr += note + "\n"
        return returnstr

    def describe(self):
        """describe prints out a description of what the lines in the
        control file mean"""
        returnstr = (
            "Center of Lat and Lon: "
            + str(self.centerlat)
            + " "
            + str(self.centerlon)
            + "\n"
        )
        returnstr += (
            "Spacing (deg) Lat, Lon: "
            + str(self.latdiff)
            + " "
            + str(self.londiff)
            + "\n"
        )
        returnstr += (
            "Span (deg Lat, Lon: " + str(self.latspan) + " " + str(self.lonspan) + "\n"
        )
        returnstr += "Output grid directory: " + self.outdir + "\n"
        returnstr += "Output grid file name: " + self.outfile + "\n"
        returnstr += "Num of vertical levels: " + str(self.nlev) + "\n"
        returnstr += "Height of levels (M Agl) : "
        for lev in self.levels:
            returnstr += str(lev) + " "
        returnstr += "\n"
        returnstr += "Sampling start (yy mm dd hh min) : " + self.sample_start + "\n"
        returnstr += "Sampling stop (yy mm dd hh min) : " + self.sample_stop + "\n"
        # returnstr +=  self.typstr() + ' ' + str(self.sampletype) + ' '
        returnstr += (
            "Interval (hh min) "
            + str(self.interval[0])
            + " "
            + str(self.interval[1])
            + "\n"
        )
        return returnstr

    def typestr(self):
        """returns a string describing what kind of sampling interval is used"""
        # print(self.interval[0], self.interval[1])
        tmstr = str(self.interval[0]).zfill(2) + ":" + str(self.interval[1]).zfill(2)
        if self.sampletype == 0:
            returnstr = "Average over  " + tmstr + " with output every " + tmstr
        elif self.sampletype == 1:
            returnstr = "Snapshot every " + tmstr
        elif self.sampletype == 2:
            returnstr = "Maximum every " + tmstr
        elif self.sampletype < 0:
            returnstr = (
                "Average over "
                + str(abs(self.sampletype))
                + " hours with output every "
                + tmstr
            )
        return returnstr

    def definition(self, lines):
        """
        Parameters
        -----------
        lines : string
        input list of 10 lines of the control file which define a concentration
        grid

        Return
        ------
        boolean
        """
        ret = True
        temp = lines[0].split()
        try:
            self.centerlat = float(temp[0])
        except TypeError:
            print("warning: center latitude not a float", temp[0])
            ret = False
        try:
            self.centerlon = float(temp[1])
        except TypeError:
            print("warning: center longitude not a float", temp[1])
            ret = False

        temp = lines[1].split()
        try:
            self.latdiff = float(temp[0])
        except TypeError:
            print("warning: spacing of latitude not a float", temp[0])
            ret = False
        try:
            self.londiff = float(temp[1])
        except TypeError:
            print("warning: spacing of longitude not a float", temp[1])
            ret = False

        temp = lines[2].split()
        try:
            self.latspan = float(temp[0])
        except TypeError:
            print("warning: span of latitude not a float", temp[0])
            ret = False
        try:
            self.lonspan = float(temp[1])
        except TypeError:
            print("warning: span of longitude not a float", temp[1])
            ret = False

        self.outdir = lines[3].strip()
        self.outfile = lines[4].strip()

        try:
            self.nlev = int(lines[5])
        except TypeError:
            print("warning: number of levels not an integer", lines[5])
            self.nlev = 0
            ret = False

        temp = lines[6].split()
        for lev in temp:
            try:
                lev = float(lev)
            except TypeError:
                print("warning: level not a float", lev)
                lev = -1
            self.levels.append(lev)

        temp = lines[7].strip()
        self.sample_start = temp
        temp = lines[8].strip()
        self.sample_stop = temp

        temp = lines[9].strip().split()
        try:
            self.sampletype = int(temp[0])
        except TypeError:
            print("warning: sample type is not an integer", temp[0])
        try:
            self.interval = (int(temp[1]), int(temp[2]))
        except TypeError:
            print("interval not integers", temp[1], temp[2])
        return ret


class Species:
    """Class which contains information to define a species or pollutant
    in a HYSPLIT control file.
    Methods
    -------
    __init__   : initialize attributes
    definition : input 3 lines from control file for defining a pollutant
    define_dep : input 5 lines from control file for defining deposition

    """

    total = 0

    @staticmethod
    def status():
        """total number of species objects"""
        return Species.total

    def __init__(
        self,
        name,
        psize=0,
        rate="1",
        duration=-1,
        density=2.5,
        shape=1,
        date="00 00 00 00 00",
        wetdepstr="0.0 0.0 0.0",
        vel="0.0 0.0 0.0 0.0 0.0",
        decay="0.0",
        resuspension="0.0",
    ):

        self.name = name
        self.rate = rate
        self.duration = duration
        self.psize = psize
        self.density = density
        self.shape = shape
        self.date = date
        # self.wetdep = wetdep
        self.wetdepstr = wetdepstr
        self.vel = vel
        self.decay = decay
        self.resuspension = resuspension
        self.datestr = "00 00 00 00 00"
        Species.total += 1

    def copy(self):
        return Species(
            self.name,
            self.psize,
            self.rate,
            self.duration,
            self.density,
            self.shape,
            self.date,
            self.wetdepstr,
            self.vel,
            self.decay,
            self.resuspension,
        )

    def definition(self, lines):
        """input 3 lines from HYSPLIT CONTROL file which define a
        pollutant/species.

        This will overwrite
        self.rate
        self.date
        self.duration
        self.datestr
        """
        try:
            self.rate = float(lines[0])
        except BaseException:
            print("warning: rate is not a float", lines[0])
            return False
        try:
            self.duration = float(lines[1])
        except BaseException:
            print("warning: duration is not a float", lines[1])
            return False
        if lines[2].strip()[0:2] == "00":
            self.date = lines[2].strip()
            self.datestr = self.date
        else:
            try:
                self.date = datetime.datetime.strptime(
                    lines[2].strip(), "%y %m %d %H %M"
                )
                self.datestr = self.date.strftime("%y %M %D %H")
            except BaseException:
                print("warning: date not valid", lines[2])
                self.date = lines[2].strip()
                self.datestr = lines[2].strip()

        return True

    def define_dep(self, lines):
        """input list of 5 lines in CONTROL file that define deposition for
        pollutant
        This will overwrite
        self.psize
        self.density
        self.shape
        self.vel
        self.wetdepstr
        self.decay
        self.resuspension
        """
        temp = lines[0].strip().split()
        try:
            self.psize = float(temp[0])
        except BaseException:
            print("warning: diameter not a float ", temp[0])
        try:
            self.density = float(temp[1])
        except BaseException:
            print("warning: density not a float ", temp[1])
        try:
            self.shape = float(temp[2])
        except BaseException:
            print("warning: shape not a float ", temp[2])
        # To do - read these in as floats
        self.vel = lines[1].strip()
        self.wetdepstr = lines[2].strip()
        # self.wetdep = 1
        self.decay = lines[3].strip()
        self.resuspension = lines[4].strip()
        return -1

    def strpollutant(self, annotate=False):
        """Prints out three lines which define a species/pollutant in HYSPLIT
        control file"""
        note = ""
        spc = " " * 20
        if annotate:
            note = spc + "#Species identifier"
        returnval = self.name + note + "\n"
        if annotate:
            note = spc + "#Rate of emission"
        returnval += str(self.rate) + note + "\n"
        if annotate:
            note = spc + "#Duration of emission"
        returnval += "%0.2f" % self.duration + note + "\n"
        if annotate:
            note = spc + "#Start date of emission"
        returnval += self.datestr + note + "\n"
        return returnval

    def add_wetdep(self, wstr):
        """add wet deposition line
        wstr : string
        """
        self.wetdepstr = wstr
        # self.wetdep = 1

    def strdep(self, annotate=True):
        """Prints out five lines which define deposition
        and gravitational settling for species/pollutant
        in HYSPLIT control file"""
        note = ""
        spc = " " * 20
        if annotate:
            if self.shape < 0:
                shapstr = "Ganser Formulation"
            if self.shape >= 0:
                shapstr = "Stokes Formulation"
            shapstr = ""
            note = spc + "#Particle diameter(um)   Density (g/cc),  Shape " + shapstr
        returnval = (
            "%05.1f" % self.psize
            + " "
            + "%02.1f" % self.density
            + " "
            + "%04.1f" % self.shape
            + note
            + "\n"
        )
        if annotate:
            note = spc + "#Dry Deposition for gas or using resistance method "
        returnval += self.vel + note + "\n"
        if annotate:
            note = spc + "#Wet deposition parameters"
        # if self.wetdep == 1:
        if self.wetdepstr == "":
            #    returnval += "0.0 4.0E+04 5.0E-06" + note + "\n"
            returnval += "0.0 0.0 0.0" + note + "\n"
        else:
            returnval += self.wetdepstr + note + "\n"
        # else:
        #    returnval += "0.0 0.0 0.0" + note + "\n"
        if annotate:
            note = spc + "#radioactive decay parameters"
        # line for radioactive decay half life
        returnval += str(self.decay) + note + "\n"
        if annotate:
            note = spc + "#resuspension from deposit"
        returnval += self.resuspension + note + "\n"  # line for resuspension factor

        return returnval


class NameList:
    """class which represents HYSPLIT SETUP.CFG file,
    This class can also be used to write GENPARM.CFG file for hycs_gem.
    In write method set gem=True"""

    def __init__(self, fname="SETUP.CFG", working_directory="./"):
        self.fname = fname
        self.nlist = {}  # dictionary of lines in the file.
        self.descrip = {}
        self._load_descrip()
        self.order = None
        if working_directory[-1] != "/":
            working_directory += "/"
        self.wdir = working_directory
        self.line_ending = ','

    @property
    def keys(self):
        keys = self.nlist.keys()
        return keys

    @property
    def line_ending(self):
        return self._line_ending

    @line_ending.setter
    def line_ending(self,line_ending):
        self._line_ending = line_ending


    def print_help(self, order=None, sep=":"):
        rstr = ""
        if not order:
            order = self.descrip.keys()
        for key in order:
            rstr += key.ljust(10) + sep
            rstr += self.descrip[key]
            rstr += "\n"
        return rstr

    def add_n(self, nlist):
        """
        add a whole dictionary.
        """
        self.nlist = nlist

    def add(self, name, value):
        """
        add one line
        """
        if name.lower() in ["poutf", "efile"]:
            if value[0] != "'":
                value = "'" + value
            if value[-1] != "'":
                value = value + "'"

        self.nlist[name.lower()] = value

    def rename(self, name, working_directory=""):
        """
        create new name and/or working directory for file
        """
        self.fname = name
        if working_directory:
            self.wdir = working_directory

    def _load_descrip(self):
        """creates dictionary with description of namelist parameters"""
        self.descrip["ichem"] = (
            "Chemistry conversion modules.\n" + "0:none, 1:matrix , 2:convert, 3:dust"
        )
        self.descrip["qcycle"] = "Cycling of emission hours"
        self.descrip["delt"] = (
            "integration time step\n" + " (0=autoset, >0= constant ,<0=minimum)"
        )
        self.descrip["kmixd"] = (
            "Mixed layer obtained from \n" + " 0:input, 1:temperature, 2: TKE",
        )
        self.descrip["kmix0"] = "mixing depth. 250 minimum"
        self.descrip["kzmis"] = (
            "Vertical mixing profile."
            + " 0:No adjustments."
            + " 1: vertical diffusivity in PBL single"
            + " average value"
        )
        self.descrip["kbls"] = (
            "Stability computed by"
            "(1) Heat and momentum fluxes," + " 2: Wind and temperature profiles"
        )
        self.descrip["kblt"] = (
            "Flag to set vertical turbulence computational"
            + "method. 1:Beljaars/Holtslag"
            + "(2):Kanthar/Clayson "
            + " 3:TKE field 4:Velocity Variances"
        )
        self.descrip["initd"] = "defines particle or puff mode"

    def summary(self):
        """prints summmary.
        Currently only prints description of INITD.
        """
        if "initd" in list(self.nlist.keys()):
            test = int(self.nlist["initd"])
            if test == 0:
                print("3D particle horizontal and vertical")
            elif test == 1:
                print("Gaussian horizontal and Top-Hat vertical puff")
            elif test == 2:
                print("Top-Hat horizontal and vertical puff")
            elif test == 3:
                print("Gaussian horizontal puff and vertical particle distribution")
            elif test == 4:
                print("Top-Hat horizontal puff and vertical particle distribution")
            else:
                print("3D particle horizontal and vertical")

    def set_dust(self):
        """
        some parameters often used for dust emissions
        """
        self.nlist["ichem"] = "3"
        self.nlist["qcycle"] = "3"

    def read(self, case_sensitive=True):
        """
        read existing SETUP.CFG file.
        """
        pname = os.path.join(self.wdir, self.fname)
        if os.path.isfile(pname):
            with open(os.path.join(self.wdir, self.fname), "r") as fid:
                content = fid.readlines()
        else:
            logger.warning("File cannot be opened {}".format(pname))
            return False
        for line in content:
            if "=" in line:
                temp = line.strip().split("=")
                key = temp[0].strip()
                if not case_sensitive:
                    key = key.lower()
                self.nlist[key] = temp[1].strip(",").strip()
        return True


    def __str__(self):
        # 9/28/2023 remove key.lower() so write with same case as input.
        order = self.order
        rstr = ""
        if order is None:
            order = list(self.nlist.keys())
        for key in order:
            kstr = True
            try:
                rstr += key + "=" + self.nlist[key] + "{}\n".format(self.line_ending)
            except BaseException:
                print("WARNING: " + str(key) + " " + str(self.nlist[key]) + " not str")
                kstr = False
            if not kstr:
                rstr += str(key) + "=" + str(self.nlist[key]) + "{}\n".format(self.line_ending)
        return rstr

    def write(self, order=None, verbose=False, overwrite=True, query=False, **kwargs):
        """
        order : list
        gem   : boolean
        verbose : boolean
        overwrite : boolean
        query : boolean

        keyword arguments
        if gem=True then will write &GENPARM at beginning of file rather than &SETUP
        """

        if 'gem' in kwargs.keys():
            gem = kwargs['gem']
        else:
            gem= False

        if 'noheader' in kwargs.keys():
            noheader = kwargs['noheader']
        else:
            noheader= False


        rval = writeover(self.wdir + self.fname, overwrite, query, verbose=verbose)
        if rval == -1:
            return rval

        if order is None:
            self.order = list(self.nlist.keys())
        else:
            self.order = order
        #  order = []
        if verbose:
            print("WRITING SETUP FILE", self.wdir + self.fname)
        with open(path.join(self.wdir, self.fname), "w") as fid:
            if gem:
                fid.write("&GEMPARM \n")
            elif noheader:
                pass
            else:
                fid.write("&SETUP \n")
            # if order == []:
            #    order = list(self.nlist.keys())
            rstr = str(self)
            fid.write(rstr)
            # for key in order:
            #    kstr = True
            #    try:
            #        fid.write(key.lower() + "=" + self.nlist[key] + ",\n")
            #    except BaseException:
            #        print(
            #            "WARNING: " + str(key) + " " + str(self.nlist[key]) + " not str"
            #        )
            #        kstr = False
            #    if not kstr:
            #        fid.write(str(key) + "=" + str(self.nlist[key]) + ",\n")
            fid.write("/ \n")
        return rval


class ControlLoc:
    """
    Release location in HYSPLIT CONTROL file
    """

    total = 0

    @staticmethod
    def status():
        """number of ControlLoc objects"""
        return ControlLoc.total

    def __init__(self, line=False, latlon=(-1, -1), alt=10.0, rate=False, area=False):
        """Can either input a string (line from HYSPLIT CONTROL file) or can enter
        latlon = tuple (default(-1,-1))
        altitude= real (default (10.0))
        rate
        area"""

        if line:
            self.definition(line)
        else:
            self.latlon = latlon
            self.alt = alt
            self.rate = rate
            self.area = area
        ControlLoc.total += 1

    def copy(self):
        return ControlLoc(False, self.latlon, self.alt, self.rate, self.area)

    def definition(self, line):
        """
        line : string
        takes line from CONTROL file and converts it to
        latitude, longitude, altitude, rate, area attributes.
        """
        temp = line.split()
        try:
            self.lat = float(temp[0])
        except BaseException:
            self.lat = -1
        try:
            self.lon = float(temp[1])
        except BaseException:
            self.lon = -1
        try:
            self.alt = float(temp[2])
        except BaseException:
            self.alt = 10.0
        try:
            self.rate = float(temp[3])
        except BaseException:
            self.rate = False
        try:
            self.area = float(temp[4])
        except BaseException:
            self.area = False
        self.latlon = (self.lat, self.lon)

    def __str__(self):
        """
        Returns string suitable for writing to CONTROL file.
        """
        spc = " "
        returnstr = "{:.4f}".format(self.latlon[0])
        returnstr += spc
        returnstr += "{:.4f}".format(self.latlon[1])
        returnstr += spc
        returnstr += "{:.1f}".format(self.alt)
        if self.rate and self.rate != -999:
            returnstr += spc
            returnstr += "{:.0f}".format(self.rate)
            # only write area if a rate was written.
            if self.area and self.area != -999:
                returnstr += spc
                returnstr += "{:.2E}".format(self.area)
        return returnstr


class HycsControl:
    """
       class which represents the HYSPLIT
       control file and all the information in it
    INPUTS
        fname : str : name of control file
        working_directory : str : directory where CONTROL file is
        rtype : string : dispersion, trajectory, vmixing

    Attributes
        fname   : str
        wdir    : str
        species   : list of objects in the Species class
        concgrids : list of objects in the ConcGrid class
        locs      : list of objects in ControlLoc class
        metfiles: list of str
        metdir  : list of str
        nlocs   : int
        num_grids : int
        num_sp  : int
        rtype   : str
        outfile : str - name of output file name
        outdir  : str
        run_duration : int
        vertical_motion : int
        ztop : int
    """

    def __init__(self, fname="CONTROL", working_directory="./", rtype="dispersion"):
        self.fname = fname
        if working_directory[-1] != "/":
            working_directory += "/"
        self.wdir = working_directory
        self.species = []  # list of objects in Species class
        self.concgrids = []  # list of object in ConcGrid class
        self.locs = []  # list of ControlLoc class objects
        self.metfiles = []  # list of str
        self.metdirs = []  # list of str
        self.nlocs = 0  # number of locations
        self.num_grids = 0  # number of concentration grids.
        self.num_sp = 0  # number of pollutants / species
        self.num_met = 0  # number of met files
        self.rtype = rtype  # dispersion or trajectory run or vmixing.
        self.outfile = "tdump"  # output file name for trajectory
        self.outdir = "./"
        self.run_duration = 1
        self.vertical_motion = 1
        self.ztop = 10000
        self.date = None  # start date of simulation

    def rename(self, name, working_directory="./"):
        """create new filename and working directory for the CONTROL file"""
        self.fname = name
        if working_directory[-1] != "/":
            working_directory += "/"
        self.wdir = working_directory

    def add_sdate(self, sdate):
        """add or overwrite the simulation start date"""
        self.date = sdate

    def remove_species(self):
        """set the species array to empty"""
        self.species = []
        self.num_sp = 0

    def add_species(self, species):
        """add new species.
        species : Species class.
        """
        self.num_sp += 1
        self.species.append(species)

    def add_cgrid(self, cgrid):
        """add new concentration grid.
        cgrid : ConcGrid class.
        """
        self.num_grids += 1
        self.concgrids.append(cgrid)

    def copy(self):
        return -1

    def add_dummy_location(self):
        newloc = self.locs[0].copy()
        self.locs.append(newloc)
        self.nlocs += 1

    def add_location(self, latlon=(0, 0), alt=10.0, rate=False, area=False):
        """add new emission location
        latlon : tuple of floats
        atl    : float
        rate   :
        area   :
        """
        cstr = str(ControlLoc(line=False, latlon=latlon, alt=alt, rate=rate, area=area))
        self.nlocs += 1
        self.locs.append(
            ControlLoc(line=False, latlon=latlon, alt=alt, rate=rate, area=area)
        )

    def print_locations(self):
        for loc in self.locs:
            print(loc)

    def add_location_str(self, locstr):
        """add new emission location
        locstr : string
        string which represents location line.
        """
        self.nlocs += 1
        self.locs.append(ControlLoc(line=locstr))

    def remove_locations(self, num=-99):
        """
        remove emission locations.
        num : integer
        default is to remove all locations.
        otherwise remove location with indice num.
        """
        if num == -99:
            self.nlocs = 0
            self.locs = []
        else:
            self.nlocs -= 1
            self.locs.pop(num)

    def add_ztop(self, ztop):
        """
        set the model top.
        ztop : integer
        """
        self.ztop = ztop

    def add_vmotion(self, vmotion):
        """
        set (or overwrite) the vertical motion method.
        vmotion : integer
        """
        self.vertical_motion = vmotion

    def add_metfile(self, metdir, metfile):
        """
        add an additional meteorological file
        metdir :  string
        metfile : string
        """
        # metdirectory needs to end with /
        if metdir[-1] != "/":
            metdir += "/"
        self.num_met += 1
        self.metfiles.append(metfile)
        self.metdirs.append(metdir)

    def remove_metfile(self, num=0, rall=False):
        """removes metfile and directory in posiiton num of the list.
        or removes all met files if rall=True"""
        if rall:
            self.num_met = 0
            self.metfiles = []
            self.metdirs = []
        else:
            self.metfiles.pop(num)
            self.metdirs.pop(num)
            self.num_met += -1

    def add_duration(self, duration):
        """will replace the duration if already exists"""
        self.run_duration = duration

    def write(
        self, annotate=False, metgrid=False, verbose=False, overwrite=True, query=False
    ):
        """writes CONTROL file to text file
        self.wdir + self.fname
        metgrid option will write a 1 before the number of met files.
        overwrite - if False then will not overwrite an existing file
        query - if True will ask before overwriting an existing file.
        """
        note = ""
        sp28 = " " * 28

        rval = writeover(self.wdir + self.fname, overwrite, query, verbose=verbose)
        if rval == -1:
            return rval

        if self.outdir[-1] != "/":
            self.outdir += "/"

        with open(path.join(self.wdir, self.fname), "w") as fid:
            fid.write(self.date.strftime("%y %m %d %H %M"))
            if annotate:
                note = " " * 18 + "#Start date of simulation"
            fid.write(note + "\n")
            if annotate:
                note = " " * 28 + "#Number of source locations"
            fid.write(str(self.nlocs) + note + "\n")
            iii = 0
            if annotate:
                note = " " * 15 + "#Lat Lon Altitude"
            for source in self.locs:
                fid.write(str(source))
                if iii > 0:
                    note = ""
                fid.write(note)
                iii += 1
                fid.write("\n")

            if annotate:
                note = sp28 + "#Duration of run"
            fid.write(str(int(self.run_duration)) + note + "\n")
            if annotate:
                note = sp28 + "#Vertical Motion"
            fid.write(str(self.vertical_motion) + note + "\n")
            if annotate:
                note = sp28 + "#Top of Model Domain"
            fid.write(str(self.ztop) + note + "\n")
            if annotate:
                note = sp28 + "#Number of Meteorological Data Files"
            if metgrid:
                fid.write("1 ")
            fid.write(str(self.num_met) + note + "\n")
            iii = 0
            for met in self.metfiles:
                if annotate:
                    note = "  #Meteorological Data Directory"
                if iii > 0:
                    note = ""
                fid.write(self.metdirs[iii])
                fid.write(note + "\n")
                if annotate:
                    note = "  #Meteorological Data Filename"
                if iii > 0:
                    note = ""
                fid.write(met)
                fid.write(note + "\n")
                iii += 1

            # done writing if using for vmixing.
            if self.rtype == "vmixing":
                return False

            # done writing if using for trajectory.
            if self.rtype == "trajectory":
                if self.outdir[-1] != "/":
                    self.outdir += "/"
                fid.write(self.outdir.strip() + "\n")
                fid.write(self.outfile.strip())
                return False

            if annotate:
                note = sp28 + "#Number of Pollutant Species"
            fid.write(str(self.num_sp) + note + "\n")
            iii = 0
            for sp in self.species:
                if iii == 0 and annotate:
                    fid.write(sp.strpollutant(annotate=True))
                else:
                    fid.write(sp.strpollutant(annotate=False))
                iii += 1
            fid.write(str(self.num_grids) + "\n")
            for cg in self.concgrids:
                if annotate:
                    cg.set_annotate()
                fid.write(str(cg))

            if annotate:
                note = sp28 + "#Number of Pollutant Species"
            fid.write(str(self.num_sp) + note + "\n")
            iii = 0
            for sp in self.species:
                if iii == 0:
                    fid.write(sp.strdep(annotate=annotate))
                else:
                    fid.write(sp.strdep(annotate=False))
                iii += 1

        return False

    def summary(self):
        """prints out summary of what is in CONTROL file"""
        print("CONTROL FILE")
        print("release start date", self.date)
        print("number of release locations", self.nlocs)
        print("run time", self.run_duration)
        print("Num of met grids ", self.num_met)
        print("Num of species ", self.num_sp)
        return True

    # def readlocs(self):
    # """
    # reads lines specifying locations.
    # """
    #    with open(self.fname, "r") as fid:
    #        for line in fid:
    #            #temp = line.strip().split(' ')
    #            #latlon = (temp[0], temp[1])
    #            self.locs.append(line.strip())
    #            self.nlocs += 1

    def parse_num_met(self, line):
        """
        line : str
              line from control file which specifies
              number of met files. Can sometimes have
              two numbers in it.
        """
        temp = line.split()
        num1 = int(temp[0])
        try:
            num2 = int(temp[1])
        except:
            num2 = 1
        return num2 * num1

    def read(self, verbose=False):
        """
        Read in control file.
        """
        with open(self.wdir + self.fname, "r") as fid:
            contentorig = fid.readlines()
            content = []
            for ln in contentorig:
                content.append(ln.split("#")[0])
            try:
                self.date = datetime.datetime.strptime(
                    content[0].strip(), "%y %m %d %H"
                )
            except BaseException:
                self.date = datetime.datetime.strptime(
                    content[0].strip(), "%y %m %d %H %M"
                )
            self.nlocs = int(content[1].strip())
            zzz = 2
            for iii in range(zzz, zzz + self.nlocs):
                temploc = content[iii].strip()
                self.locs.append(ControlLoc(line=temploc))
            zzz += self.nlocs
            self.run_duration = content[zzz].strip()
            self.vertical_motion = content[zzz + 1].strip()
            self.ztop = content[zzz + 2].strip()

            num_met = content[zzz + 3].strip()
            self.num_met = self.parse_num_met(num_met)
            # self.num_met = int(content[zzz + 3].strip())

            zzz = zzz + 4
            for iii in range(zzz, zzz + 2 * self.num_met, 2):
                self.metdirs.append(content[iii].strip())
                self.metfiles.append(content[iii + 1].strip())
            zzz = zzz + 2 * self.num_met
            # if it is a trajectory control file then just
            # two more lines
            if self.rtype == "vmixing":
                return "Vmixing"
            if self.rtype == "trajectory":
                self.outdir = content[zzz].strip()
                self.outfile = content[zzz + 1].strip()
                return "Traj"
            # this is end of trajectory file

            self.num_sp = int(content[zzz])
            zzz += 1
            for iii in range(zzz, zzz + 4 * self.num_sp, 4):
                lines = []
                spname = content[iii].strip()
                lines.append(content[iii + 1])
                lines.append(content[iii + 2])
                lines.append(content[iii + 3])
                sptemp = Species(spname)
                if sptemp.definition(lines):
                    self.species.append(sptemp)
            zzz += 4 * self.num_sp
            self.num_grids = int(content[zzz].strip())
            self.concgrids = []
            for iii in range(zzz, zzz + 10 * self.num_grids, 10):
                lines = []
                spname = content[iii].strip()
                for kkn in range(1, 11):
                    lines.append(content[iii + kkn])
                sptemp = ConcGrid(spname)
                if sptemp.definition(lines):
                    self.concgrids.append(sptemp)
            zzz += 10 * self.num_grids
            zzz += 1

            temp = int(content[zzz].strip())
            if temp != self.num_sp:
                print(
                    "warning: number of species for deposition",
                    " not equal to number of species",
                )
            for nnn, iii in enumerate(range(zzz, zzz + 5 * self.num_sp, 5)):
                lines = []
                for kkn in range(1, 6):
                    lines.append(content[iii + kkn])
                self.species[nnn].define_dep(lines)

            if verbose:
                print("---------------------------")
                print("CONTROL FILE")
                print("release start date", self.date)
                print("release locations", self.locs)
                print("run time", self.run_duration)
                print("vertical motion", self.vertical_motion)
                print("Top of model domain", self.ztop)
                print("Num of met grids ", self.num_met)
                print("Met directories ", self.metdirs)
                print("Met files ", self.metfiles)
                print("Num of species ", self.num_sp)
                for kkn, sp in enumerate(self.species):
                    print("-----Species ", str(kkn), "---------")
                    print(sp.strpollutant())
                    print(sp.strdep())
                    print("--------------")
                for kkn, grid in enumerate(self.concgrids):
                    print("-----Concentration Grid ", str(kkn), "---------")
                    print(grid)
                    print("--------------")
                print("---------------------------")

        return True


def roundtime(dto):
    """rounds input datetime to day at 00 H"""
    return datetime.datetime(dto.year, dto.month, dto.day, 0, 0)

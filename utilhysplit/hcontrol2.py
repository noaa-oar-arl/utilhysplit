# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import datetime
import logging
from os import path

import numpy as np

# import species_module
from utilhysplit.control_objects import HysplitDate, Latitude, Longitude, writeover

"""
PGRMMR: Alice Crawford ORG: ARL
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
# 2023 Dec 05 (AMC) replaced exception in parse_num_met function with test for length of array.
# 2024 Mar 20 (AMC) NameList class added remove method and modified add_n method
# 2024 Mar 25 (AMC) Fixed bug in Species datestr specification in definition method.
# 2024 Apr 14 (AMC) Updates to ConcGrid class

logger = logging.getLogger(__name__)


class PollProperties:
    """
    Pollutant Properties
    """

    def __init__(
        self,
        hysplit_name: str = "POLL",
        species: str = "generic",
        pt_diameter: float = 0.0,
        pt_density: float = 0.0,
        pt_shape: float = 0.0,
        dry_depo_rate: float = 0.0,
        mol_wgt: float = 0.0,
        aratio: float = 0.0,
        dratio: float = 0.0,
        effhenry: float = 0.0,
        henry_const: float = 0.0,
        wet_depo_in_cloud: float = 0.0,
        wet_depo_below_cloud: float = 0.0,
        half_life_days: float = 0.0,
        resuspension: float = 0.0,
    ):
        self.hysplit_name = hysplit_name  # e.g., Cs-137
        self.species = species  # e.g., Caesium

        # first line
        self.pt_diameter = pt_diameter
        self.pt_density = pt_density
        self.pt_shape = pt_shape

        # second line
        self.dry_depo_rate = dry_depo_rate
        self.mol_wgt = mol_wgt
        self.aratio = aratio
        self.dratio = dratio
        self.effhenry = effhenry

        # third line
        self.wet_depo_in_cloud = wet_depo_in_cloud
        self.wet_depo_below_cloud = wet_depo_below_cloud
        self.henry_const = henry_const

        # fourth line
        self.half_life_days = half_life_days

        # fifth line. should always be 0.
        self.resuspension = resuspension

    # def __str__(self):
    #   return f"Species(species='{self.__species}', hysplit_name='{self.hysplit_name}')"

    # @property
    # def species(self) -> str:
    #    return self.__species

    # @species.setter
    # def species(self, val: str) -> None:
    #   self.__species = val

    @property
    def primary_hysplit_name(self) -> str:
        return self.hysplit_name

    # TODO - not sure what these are
    # @property
    # def phases(self) -> tuple:
    #   return (self, )

    # TODO - not sure what these are
    # @property
    # def partitioning_ratios(self) -> tuple:
    #   return (1.0, )

    def get_hysplit_control_fragment(self) -> str:
        # particle diameter (um), density (g/cc), shape
        s = f"{self.pt_diameter:.5E} {self.pt_density:.5E} {self.pt_shape:.5E}\n"
        # dry deposition velocity (m/s), pollutant molecuar weight,
        # surface reactivity ratio, diffusivity ratio, effective Henry constant\n'
        s += f"{self.dry_depo_rate:.5E} {self.mol_wgt:0.1F} {self.aratio:0.1F} {self.dratio:0.1F} {self.effhenry:0.1F}\n"
        # actual Henry constant, in-cloud, below-cloud
        s += f"{self.henry_const:.5E} {self.wet_depo_in_cloud:.5E} {self.wet_depo_below_cloud:.5E}\n"
        # radioactive decay half-life (days)
        s += f"{self.half_life_days}\n"
        s += f"{self.resuspension}"
        return s


class PollDefinition:
    """
    Pollutant Definition
    """

    def __init__(
        self,
        hysplit_name: str = "P006",
        rate: float = 1.0,
        duration: float = 1.0,
        unit: str = "unit",
        datestr="00 00 00 00 00",
    ):
        self.hysplit_name = hysplit_name
        self.rate = rate
        self.duration = duration
        self.unit = unit
        self._datestr = HysplitDate(datestr)

    def __str__(self):
        return self.get_hysplit_control_fragment()

    @property
    def hysplit_name(self):
        return self._hysplit_name

    @hysplit_name.setter
    def hysplit_name(self, name):
        self._hysplit_name = name.strip()

    @property
    def datestr(self):
        return str(self._datestr)

    @property
    def date(self):
        return self._datestr.date

    def get_hysplit_control_fragment(self) -> str:
        # identifier
        s = f"{self.hysplit_name}\n"
        # rate
        s += f"{self.rate:.4F} # {self.unit}/h\n"
        # duration
        s += f"{self.duration:.4F}\n"
        # datestr
        s += f"{self.datestr}"
        return s


def lines2polldef(lines: list):
    """input 3 lines from HYSPLIT CONTROL file which define a"""
    if len(lines) < 4:
        logger.warning("Not enough lines for definition")
    name = lines[0]
    try:
        rate = float(lines[1])
    except BaseException:
        print("warning: rate is not a float", lines[0])
        return False
    try:
        duration = float(lines[2])
    except BaseException:
        print("warning: duration is not a float", lines[1])
        return False
    datestr = lines[3]
    return PollDefinition(
        hysplit_name=name, rate=rate, duration=duration, datestr=datestr
    )


def lines2pollprop(lines):
    """input list of 5 lines in CONTROL file that define deposition for"""

    # first line
    temp = lines[0].strip().split()
    try:
        psize = float(temp[0])
    except BaseException:
        print("warning: diameter not a float ", temp[0])
    try:
        density = float(temp[1])
    except BaseException:
        print("warning: density not a float ", temp[1])
    try:
        shape = float(temp[2])
    except BaseException:
        print("warning: shape not a float ", temp[2])

    # second line
    vel = lines[1].strip()
    vel = vel.split()
    try:
        depvel = float(vel[0])
    except ValueError as ex:
        logger.warning(ex)
        depvel = 0.0
    aratio = float(vel[1])
    dratio = float(vel[2])
    effhenry = float(vel[3])

    # second line
    wstr = lines[2].strip()
    wstr = wstr.split()
    try:
        henry = float(wstr[0])
    except ValueError as ex:
        logger.warning(ex)
        henry = 0.0
    try:
        incloud = float(wstr[1])
    except ValueError as ex:
        logger.warning(ex)
        incloud = 4e-5
    try:
        belowcloud = float(wstr[2])
    except ValueError as ex:
        logger.warning(ex)
        belowcloud = 4e-5

    # self.wetdep = 1
    half_life_days = float(lines[3].strip())
    resuspension = lines[4].strip()
    return PollProperties(
        pt_diameter=psize,
        pt_density=density,
        pt_shape=shape,
        dry_depo_rate=depvel,
        aratio=aratio,
        dratio=dratio,
        effhenry=effhenry,
        henry_const=henry,
        wet_depo_in_cloud=incloud,
        wet_depo_below_cloud=belowcloud,
        half_life_days=half_life_days,
        resuspension=resuspension,
    )


class HysplitSpecies:
    """Class which contains information to define a species or pollutant
    in a HYSPLIT control file.
    """

    total = 0

    @staticmethod
    def status():
        """total number of species objects"""
        return HysplitSpecies.total

    def __init__(
        self,
        definition: PollDefinition,
        properties: PollProperties,
    ):
        self.definition = definition
        self.properties = properties
        HysplitSpecies.total += 1

    @property
    def datestr(self):
        return self.definition.datestr

    @datestr.setter
    def datestr(self, date):
        self.definition.datestr = date

    @property
    def rate(self):
        return self.definition.rate

    @property
    def duration(self):
        return self.definition.duration

    def copy(self):
        return HysplitSpecies(
            self.definition,
            self.properties,
        )

    def strpollutant(self, annotate=False):
        return str(self.definition)

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
        return self.properties.get_hysplit_control_fragment()


class MultipleSpecies:
    def __init__(self, rate: float = 1.0, unit="unit"):
        self.rate = rate  # total rate
        self._species = []  # List of Species objects
        self.partitioning_ratios = {}  # dictionary of floats

    @property
    def specieslist(self):
        return self._species

    @property
    def number_species(self):
        return len(self._species)

    @property
    def rate(self):
        return self._rate

    @rate.setter
    def rate(self, rate: float):
        self._rate = rate

    def add_species(self, species: HysplitSpecies, partitioning_ratio: float):
        self._species.append(species)
        self.partitioning_ratios[species.definition.hysplit_name] = partitioning_ratio

    def change_ratio(self, hysplit_name: str, ratio: float):
        self.partitioning_ratios[hysplit_name] = ratio

    def check_total(self):
        total_rate = 0
        for sp in self._species:
            total_rate += sp.definition.rate
        return total_rate == self.rate

    def calculate_from_species(self):
        total_rate = 0
        unitlist = []
        for sp in self._species:
            unitlist.append(sp.definition.unit)
            total_rate += sp.definition.rate
        self.rate = total_rate
        if total_rate == 0:
            return total_rate
        for sp in self._species:
            self.partitioning_ratios[sp.definition.hysplit_name] = (
                sp.definition.rate / total_rate
            )
        return total_rate

    def apply_rate(self, rate: float):
        self.rate = rate
        for sp in self._species:
            newrate = self.rate * self.partitioning_ratios[sp.definition.hysplit_name]
            sp.definition.rate = newrate


class ConcGrid:
    """
    concentration grid as defined by 10 lines in the HYSPLIT concentration
    """

    def __init__(
        self,
        name: str = "concgrid",
        levels: list = [500],
        centerlat: float = 0.0,
        centerlon: float = 0.0,
        latdiff: float = 0.25,
        londiff: float = 0.25,
        # global grid is default
        latspan: float = 180.0,
        lonspan: float = 360.0,
        outdir: str = "./",
        outfile: str = "cdump",
        # Python 3.10+ syntax commented out
        # sample_start: datetime.datetime | str = "00 00 00 00 00",
        # sample_stop: datetime.datetime | str = "00 00 00 00 00",
        sample_start = "00 00 00 00 00",  # Changed for Python 3.7 compatibility
        sample_stop = "00 00 00 00 00",   # Changed for Python 3.7 compatibility
        sampletype: int = -1,  # <0 or 0,1,2
        # Python 3.10+ syntax commented out
        # interval: tuple[int, int] = (1, 0),
        interval = (1, 0),     # Changed for Python 3.7 compatibility
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
        self.errlist = []
        self.centerlat = Latitude(centerlat, format_specifier="2.4F")
        self.centerlon = Longitude(centerlon, format_specifier="3.4F")
        self.latdiff = latdiff
        self.londiff = londiff
        self.latspan = latspan
        self.lonspan = lonspan
        self.outdir = outdir
        self.outfile = outfile
        # string (could be changed to datetime)
        self.sample_start = HysplitDate(sample_start)
        self.sample_stop = HysplitDate(sample_stop)
        self.sampletype = sampletype
        self.interval = interval
        self.annotate = False

    @property
    def outdir(self):
        return self._outdir

    @outdir.setter
    def outdir(self, outdir):
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
            self.sample_start,
            self.sample_stop,
            self.sampletype,
            self.interval,
        )

    def add_comment(self, line: str, comment: str):
        if not self.annotate:
            return f"{line}\n"
        csc = 40  # comment start column
        clen = len(line)
        if clen >= csc:
            line += "  "
        else:
            for x in range(clen, csc):
                line += " "
        line += comment
        line += "\n"
        return line

    @property
    def nlev(self):
        return len(self.levels)

    def __str__(self):
        """string method will output ten lines suitable for inserting into a
        HYSPLIT control file"""

        # pnotes = self.annotate
        # note = ""
        # if pnotes:
        pnotes = False
        note = "#Concentration Grid Center (latitude longitude)"
        returnstr = self.add_comment(f"{self.centerlat} {self.centerlon}", note)
        note = "#Concentration grid spacing (degrees latitude longitude)"
        returnstr += self.add_comment(f"{self.latdiff} {self.londiff}", note)
        note = "#Concentration grid span (degrees latitude longitude)"
        returnstr += self.add_comment(f"{self.latspan:.2F} {self.lonspan:.2F}", note)
        note = "#Directory to write concentration output file"
        returnstr += self.add_comment(f"{self.outdir}", note)
        note = "#Filename for output"
        returnstr += self.add_comment(f"{self.outfile}", note)
        note = "#Number of vertical levels for concentration grid"
        returnstr += self.add_comment(f"{self.nlev}", note)
        note = "#List of vertical levels for concentration grid"
        levelstr = " ".join([str(x) for x in self.levels])
        returnstr += self.add_comment(levelstr, note)
        note = "#Sampling start time of concentration grid"
        returnstr += self.add_comment(str(self.sample_start), note)
        note = "#Sampling stop time of concentration grid"
        returnstr += self.add_comment(str(self.sample_stop), note)
        note = "#" + self.typestr()
        astr = "{:02.0f}".format(self.sampletype) + " "
        astr += "{:02.0f}".format(self.interval[0]) + " "
        astr += "{:02.0f}".format(self.interval[1]) + " "
        returnstr += self.add_comment(astr, note)
        return returnstr

    @property
    def sampletype(self):
        return self._sampletype

    @sampletype.setter
    def sampletype(self, value: int):
        """
        if less than 0 then indicates averaging time in hours.
        if >= 0 then can have values of
        0 (average over interval)
        1 (snapshot)
        2 (maximum value in interval)

        Values greater than 2 are invalid and the default of 0 is implemented.
        """
        if value < 0:
            self._sampletype = value
        elif value > 2:
            logger.warning(f"Invalid sample type {value}")
            self._sampletype = 0

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


def lines2concgrid(lines: list, name="ConcentrationGrid") -> ConcGrid:
    """
    Parameters
    -----------
    lines : string
    input list of 10 lines of the control file which define a concentration
    grid

    Return
    ------
    ConcGrid or None
    """
    if len(lines) < 10:
        logger.warning("Concentration Grid needs 10 lines")
        return None
    ret = "something"
    levels = []
    temp = lines[0].split()

    try:
        centerlat = float(temp[0])
    except ValueError as ex:
        logger.warning(f"centerlat: {ex}")
        ret = None
    try:
        centerlon = float(temp[1])
    except ValueError as ex:
        logger.warning(f"centerlon: {ex}")
        ret = None

    temp = lines[1].split()
    try:
        latdiff = float(temp[0])
    except ValueError as ex:
        logger.warning(f"latitude spacing: {ex}")
        ret = None
    try:
        londiff = float(temp[1])
    except ValueError as ex:
        logger.warning(f"longitude spacing: {ex}")
        ret = None

    temp = lines[2].split()
    try:
        latspan = float(temp[0])
    except ValueError as ex:
        logger.warning(f"span of latitude: {ex}")
        ret = None
    try:
        lonspan = float(temp[1])
    except ValueError as ex:
        logger.warning(f"span of longitude: {ex}")
        ret = None

    outdir = lines[3].strip()
    outfile = lines[4].strip()

    # don't use nlev
    # try:
    #    nlev = int(lines[5])
    # except ValueError as ex:
    #    logger.warning(f"number of level: {ex}")
    #    nlev = 0
    #    ret = None

    temp = lines[6].split()
    for lev in temp:
        try:
            lev = float(lev)
        except ValueError as ex:
            logger.warning(f"level: {ex}")
            lev = -1
            ret = None
        levels.append(lev)

    temp = lines[7].strip()
    sample_start = temp
    temp = lines[8].strip()
    sample_stop = temp

    temp = lines[9].strip().split()
    try:
        sampletype = int(temp[0])
    except ValueError as ex:
        logger.warning(f"sample type: {ex}")
        ret = None
    try:
        interval = (int(temp[1]), int(temp[2]))
    except ValueError as ex:
        logger.warning(f"interval: {ex}")
        ret = None
    if ret is None:
        return ret
    else:
        return ConcGrid(
            name,
            levels=levels,
            centerlat=centerlat,
            centerlon=centerlon,
            latdiff=latdiff,
            londiff=londiff,
            latspan=latspan,
            lonspan=lonspan,
            outdir=outdir,
            outfile=outfile,
            sample_start=sample_start,
            sample_stop=sample_stop,
            sampletype=sampletype,
            interval=interval,
        )


def line2control_loc(line: str):
    """
    line : string
    takes line from CONTROL file and converts it to
    latitude, longitude, altitude, rate, area attributes.
    """
    rate = None
    area = None
    temp = line.split()
    try:
        latitude = float(temp[0])
    except BaseException:
        latitude = -1
    try:
        longitude = float(temp[1])
    except BaseException:
        longitude = -1
    try:
        alt = float(temp[2])
    except BaseException:
        alt = 10.0

    if len(temp) > 3:
        try:
            rate = float(temp[3])
        except BaseException:
            pass

    if len(temp) > 4:
        try:
            area = float(temp[4])
        except BaseException:
            area = None

    return ControlLoc(latitude, longitude, alt, rate, area)


class ControlLoc:
    """
    Release location in HYSPLIT CONTROL file
    """

    total = 0

    @staticmethod
    def status():
        """number of ControlLoc objects"""
        return ControlLoc.total

    def __init__(
        self,
        latitude: float = 0.0,
        longitude: float = 0.0,
        alt: float = 10.0,
        # Python 3.10+ syntax commented out
        # rate: float | None = None,
        # area: float | None = None,
        rate = None,  # Changed for Python 3.7 compatibility
        area = None,  # Changed for Python 3.7 compatibility
    ):
        """Can either input a string (line from HYSPLIT CONTROL file) or can enter
        latlon = tuple (default(-1,-1))
        altitude= real (default (10.0))
        rate
        area"""

        self.latitude = Latitude(latitude)
        self.longitude = Longitude(longitude)
        self.alt = alt
        self.rate = rate
        self.area = area
        ControlLoc.total += 1

    def copy(self):
        return ControlLoc(
            float(self.latitude), float(self.longitude), self.alt, self.rate, self.area
        )

    def __str__(self):
        """
        Returns string suitable for writing to CONTROL file.
        """
        spc = " "
        returnstr = f"{self.latitude}"
        returnstr += spc
        returnstr += f"{self.longitude}"
        returnstr += spc
        returnstr += f"{self.alt:0.4F}"
        if self.rate is not None:
            returnstr += spc
            returnstr += f"{self.rate:.0f}"
            # only write area if a rate was written.
            if self.area is not None:
                returnstr += spc
                returnstr += f"{self.area:.4E}"
        return returnstr



class HycsControl:
    """
       class which represents the HYSPLIT
       control file and all the information in it

    INPUTS
        fname : str : name of control file
        working_directory : str : directory where CONTROL file is

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
        outfile : str - name of output file name
        outdir  : str
        run_duration : int
        vertical_motion : int
        ztop : int
        annotate : bool
    """

    def __init__(self, fname="CONTROL", working_directory="./"):
        self.fname = fname
        if working_directory[-1] != "/":
            working_directory += "/"
        self.wdir = working_directory
        self._species = []  # list of objects in Species class
        self.concgrids = []  # list of object in ConcGrid class
        self.locs = []  # list of ControlLoc class objects
        self.metfiles = []  # list of str
        self.metdirs = []  # list of str
        self.nlocs = 0  # number of locations
        self.num_grids = 0  # number of concentration grids.
        self.num_sp = 0  # number of pollutants / species
        self.num_met = 0  # number of met files
        self.outfile = "tdump"  # output file name for trajectory
        self.outdir = "./"  # str
        self.run_duration = 1  # integer (hours)
        self.vertical_motion = 1  # integer
        self.ztop = 10000  # float
        self.date = "00 00 00 00 00"  # start date of simulation

        self.annotate = False
        self.metgrid = 0

    def check(self, fix=True):
        ztop = self.ztop
        toplevel = 0
        rval = True
        for cgrid in self.concgrids:
            levels = cgrid.levels
            tlevel = levels[-1]
            toplevel = np.max([toplevel, tlevel])
        if ztop < toplevel:
            rval = False
            logger.warning("ztop below top concentration grid level")
            if fix:
                ztop = levels[-1]
        return rval

    @property
    def species(self):
        rval = MultipleSpecies()
        for sp in self._species:
            rval.add_species(sp, partitioning_ratio=1)
            rval.calculate_from_species()
        return rval

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
        self._species = []
        self.num_sp = 0

    def add_species(self, species: HysplitSpecies):
        """add new species.
        species : Species class.
        """
        self.num_sp += 1
        self.species.append(species)

    def add_multiple_species(self, multiple: MultipleSpecies):
        for sp in multiple.specieslist:
            self._species.append(sp)
            self.num_sp += 1

    def add_cgrid(self, cgrid: ConcGrid):
        """add new concentration grid.
        cgrid : ConcGrid class.
        """
        self.num_grids += 1
        self.concgrids.append(cgrid)

    # def copy(self):
    #    return -1

    def add_dummy_location(self):
        newloc = self.locs[0].copy()
        self.locs.append(newloc)
        self.nlocs += 1

    def add_location(
        self,
        # Python 3.10+ syntax commented out
        # latitude: float | None = None,
        # longitude: float | None = None,
        latitude = None,  # Changed for Python 3.7 compatibility
        longitude = None,  # Changed for Python 3.7 compatibility
        alt: float = 10.0,
        # Python 3.10+ syntax commented out
        # rate: float | None = None,
        # area: float | None = None,
        rate = None,  # Changed for Python 3.7 compatibility
        area = None,  # Changed for Python 3.7 compatibility
        latlon: tuple = (0, 0),
    ):
        """add new emission location
        latlon : tuple of floats
        atl    : float
        rate   :
        area   :
        """
        self.nlocs += 1
        if latitude is None:
            latitude = latlon[0]
        if longitude is None:
            longitude = latlon[1]
        self.locs.append(
            ControlLoc(
                latitude=latitude, longitude=longitude, alt=alt, rate=rate, area=area
            )
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
        self.locs.append(line2control_loc(locstr))

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

    @property
    def ztop(self):
        return self._ztop

    @ztop.setter
    def ztop(self, ztop: int):
        """
        set the model top.
        ztop : integer
        """
        self._ztop = ztop

    @property
    def vertical_motion(self):
        return self._vertical_motion

    @vertical_motion.setter
    def vertical_motion(self, vmotion: int):
        """
        0: data
        1: isobaric
        2: isentropic
        3: constant density
        4: constant sigma coordinate
        5: divergence
        6: special
        7: averaged
        8: damped
        """

        if vmotion > 8 or vmotion < 0:
            logger.warning(f"invalid valuefor vertical motion {vmotion}")
            vmotion = 0
        self._vertical_motion = vmotion

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

    @property
    def run_duration(self):
        return self._run_duration

    @run_duration.setter
    def run_duration(self, duration: int):
        """will replace the duration if already exists"""
        self._run_duration = int(duration)

    def __str__(self):
        """writes CONTROL file to text file
        self.wdir + self.fname
        """
        note = ""
        sp28 = " " * 28
        rval = ""
        annotate = self.annotate

        if annotate:
            note = " " * 18 + "#Start date of simulation"
        rval += f"{self.date}{note}\n"
        if annotate:
            note = " " * 28 + "#Number of source locations"
        rval += f"{self.nlocs}{note}\n"
        if annotate:
            note = " " * 15 + "#Lat Lon Altitude"
        for iii, source in enumerate(self.locs):
            rval += f"{source}{note}\n"
            if iii > 0:
                note = ""
        if annotate:
            note = sp28 + "#Duration of run"
        rval += f"{self.run_duration:}{note}\n"
        if annotate:
            note = sp28 + "#Vertical Motion"
        rval += f"{self.vertical_motion}{note}\n"
        if annotate:
            note = sp28 + "#Top of Model Domain"
        rval += f"{self.ztop}{note}\n"
        if annotate:
            note = sp28 + "#Number of Meteorological Data Files"
        if self.metgrid > 0:
            rval += "1 "
        rval += f"{self.num_met}{note}\n"
        for iii, met in enumerate(self.metfiles):
            if annotate:
                note = "  #Meteorological Data Directory"
            if iii > 0:
                note = ""
            rval += f"{self.metdirs[iii]}{note}\n"
            if annotate:
                note = "  #Meteorological Data Filename"
            if iii > 0:
                note = ""
            rval += f"{met}{note}\n"

        if annotate:
            note = sp28 + "#Number of Pollutant Species"
        rval += f"{self.num_sp}{note}\n"
        for iii, sp in enumerate(self._species):
            if iii == 0 and annotate:
                rval += sp.strpollutant(annotate=True)
            else:
                rval += sp.strpollutant(annotate=False)
            rval += "\n"
        rval += f"{self.num_grids}\n"
        for cg in self.concgrids:
            cg.annotate = self.annotate
            rval += str(cg)
        if annotate:
            note = sp28 + "#Number of Pollutant Species"
        rval += f"{self.num_sp}{note}\n"
        for iii, sp in enumerate(self._species):
            if iii == 0:
                rval += sp.strdep(annotate=annotate)
            else:
                rval += sp.strdep(annotate=False)
            rval += "\n"
        return rval

    def write(self, annotate: bool = False, overwrite=False):
        self.annotate = annotate
        rval = writeover(self.wdir + self.fname, overwrite)
        if rval:
            with open(path.join(self.wdir, self.fname), "wt") as fid:
                fid.write(str(self))
        return rval

    def summary(self):
        """prints out summary of what is in CONTROL file"""
        print("CONTROL FILE")
        print("release start date", self.date)
        print("number of release locations", self.nlocs)
        print("run time", self.run_duration)
        print("Num of met grids ", self.num_met)
        print("Num of species ", self.num_sp)
        return True

    def parse_num_met(self, line):
        """
        line : str
              line from control file which specifies
              number of met files. Can sometimes have
              two numbers in it.
        """
        temp = line.split()
        num1 = int(temp[0])
        num2 = 1
        if len(temp) > 1:
            num2 = int(temp[1])
        # try:
        #    num2 = int(temp[1])
        # except (ValueError,TypeError,IndexError):
        #    num2 = 1
        return num2 * num1

    @property
    def date(self):
        return self._date

    @date.setter
    def date(self, date):
        self._date = HysplitDate(date)

    # @date.setter
    # def date(self,date):
    #    if isinstance(date,str):
    #        date = date.strip()
    #        temp = date.split()
    #        if len(temp)==4:
    #           self._date = datetime.datetime.strptime(date, "%y %m %d %H")
    #        elif len(temp)==5:
    #           self._date = datetime.datetime.strptime(date, "%y %m %d %H %M")
    #    elif isinstance(date,datetime.datetime):
    #        self._date = date
    #    else:
    #        logger.warning(f'start date in incorrect format {type(date)}')

    def read(self, verbose=False):
        """
        Read in control file.
        """
        deflist = []
        proplist = []
        with open(self.wdir + self.fname, "r") as fid:
            contentorig = fid.readlines()
            content = []
            # remove comments
            for ln in contentorig:
                content.append(ln.split("#")[0])
            # first line is the start date.
            self.date = content[0]
            self.nlocs = int(content[1].strip())

            zzz = 2
            for iii in range(zzz, zzz + self.nlocs):
                temploc = content[iii].strip()
                self.locs.append(line2control_loc(temploc))

            zzz += self.nlocs
            self.run_duration = int(content[zzz].strip())
            self.vertical_motion = int(content[zzz + 1].strip())
            self.ztop = content[zzz + 2].strip()

            num_met = content[zzz + 3].strip()
            self.num_met = self.parse_num_met(num_met)
            # self.num_met = int(content[zzz + 3].strip())

            zzz = zzz + 4
            for iii in range(zzz, zzz + 2 * self.num_met, 2):
                self.metdirs.append(content[iii].strip())
                self.metfiles.append(content[iii + 1].strip())

            # Species definition.
            zzz = zzz + 2 * self.num_met
            self.num_sp = int(content[zzz])
            zzz += 1
            zzz2 = zzz + 4 * self.num_sp
            for iii in range(zzz, zzz2, 4):
                lines = []
                lines = content[iii : iii + 4]
                deflist.append(lines2polldef(lines))

            # Concentration Grids
            zzz += 4 * self.num_sp
            self.num_grids = int(content[zzz].strip())
            self.concgrids = []
            for iii in range(zzz, zzz + 10 * self.num_grids, 10):
                lines = []
                spname = content[iii].strip()
                for kkn in range(1, 11):
                    lines.append(content[iii + kkn])
                cgrid = lines2concgrid(lines, name=spname)
                if isinstance(cgrid, ConcGrid):
                    self.concgrids.append(cgrid)
                else:
                    logger.warning("Concentration Grid not added")
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
                proplist.append(lines2pollprop(lines))

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

        for val in zip(deflist, proplist):
            self._species.append(HysplitSpecies(val[0], val[1]))

        return True


class TrajControl(HycsControl):

    @property
    def outdir(self):
        return self._outdir

    @outdir.setter
    def outdir(self,outdir):
        self._outdir = outdir

    @property
    def outfile(self):
        return self._outfile 

    @outfile.setter
    def outfile(self,outfile):
        self._outfile = outfile

    def read(self, verbose=False):
        """
        Read in control file.
        """
        deflist = []
        proplist = []
        with open(self.wdir + self.fname, "r") as fid:
            contentorig = fid.readlines()
            content = []
            # remove comments
            for ln in contentorig:
                content.append(ln.split("#")[0])
            # first line is the start date.
            self.date = content[0]
            self.nlocs = int(content[1].strip())

            zzz = 2
            for iii in range(zzz, zzz + self.nlocs):
                temploc = content[iii].strip()
                self.locs.append(line2control_loc(temploc))

            zzz += self.nlocs
            self.run_duration = int(content[zzz].strip())
            self.vertical_motion = int(content[zzz + 1].strip())
            self.ztop = content[zzz + 2].strip()
            num_met = content[zzz + 3].strip()
            self.num_met = self.parse_num_met(num_met)
            # self.num_met = int(content[zzz + 3].strip())

            zzz = zzz + 4
            for iii in range(zzz, zzz + 2 * self.num_met, 2):
                self.metdirs.append(content[iii].strip())
                self.metfiles.append(content[iii + 1].strip())

    def __str__(self):
        """writes CONTROL file to text file
        self.wdir + self.fname
        """
        note = ""
        sp28 = " " * 28
        rval = ""
        annotate = self.annotate

        if annotate:
            note = " " * 18 + "#Start date of simulation"
        rval += f"{self.date}{note}\n"
        if annotate:
            note = " " * 28 + "#Number of source locations"
        rval += f"{self.nlocs}{note}\n"
        if annotate:
            note = " " * 15 + "#Lat Lon Altitude"
        for iii, source in enumerate(self.locs):
            rval += f"{source}{note}\n"
            if iii > 0:
                note = ""
        if annotate:
            note = sp28 + "#Duration of run"
        rval += f"{self.run_duration:}{note}\n"
        if annotate:
            note = sp28 + "#Vertical Motion"
        rval += f"{self.vertical_motion}{note}\n"
        if annotate:
            note = sp28 + "#Top of Model Domain"
        rval += f"{self.ztop}{note}\n"
        if annotate:
            note = sp28 + "#Number of Meteorological Data Files"
        if self.metgrid > 0:
            rval += "1 "
        rval += f"{self.num_met}{note}\n"
        for iii, met in enumerate(self.metfiles):
            if annotate:
                note = "  #Meteorological Data Directory"
            if iii > 0:
                note = ""
            rval += f"{self.metdirs[iii]}{note}\n"
            if annotate:
                note = "  #Meteorological Data Filename"
            if iii > 0:
                note = ""
            rval += f"{met}{note}\n"
         
        rval += f"{self.outdir }\n"
        rval += f"{self.outfile }\n"
        return rval


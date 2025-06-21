import datetime
import logging
from os import path

import os_helper

logger = logging.getLogger(__name__)


def writeover(name, overwrite, verbose=False):
    """
    checks if file already exits.
    Returns True if file should be written.
    Returns False if file should not be written.

    Inputs
    name : str
           filename
    overwrite: boolean
           if True then will 1.
    verbose : boolean
           if True will write messages

    Outputs
    rval : int
         Value is 1 if file should be written.
         Value is -1 if file should not be written.
    """
    rval = True
    if path.isfile(name):
        logger.info(f"File exists {name}")
        if overwrite:
            os_helper.remove_file(name)
            rval = True
        else:
            rval = False
    return rval


class Longitude:
    def __init__(
        self,
        longitude: float,
        format_specifier: str = "2.4F",
        default: float = -180.0,
    ):
        self._valid = True
        self._error = None
        self._default = default
        self.format_specifier = format_specifier
        self.longitude = longitude

    @property
    def longitude(self):
        return self._longitude

    @longitude.setter
    def longitude(self, value):
        # correct if from 0 to 360 and
        # apply cyclic boundary conditions
        if value > 180.0:
            value = value - 360.0
        elif value < -180:
            value = value + 360
        # if value still out of range then throw an error
        if value > 180.0 or value < -180.0:
            self._valid = False
            self._error = f"Longitude must be in range [-180.0, 180.0]"
            value = self._default
        self._longitude = value
        self._longitude = float(str(self))

    def __lt__(self, other):
        if self.longitude < other.longitude:
            return True
        else:
            return False

    def __eq__(self, other):
        if str(self) == str(other):
            return True
        else:
            return False

    def __str__(self):
        return f"{self._longitude:{self.format_specifier}}"

    def __repr__(self):
        return f"{self._longitude:{self.format_specifier}}"

    def __float__(self):
        return self._longitude

    def __add__(self, other):
        fsp = self.format_specifier
        default = self._default
        if isinstance(other, (float, int)):
            new = self._longitude + other
        elif isinstance(other, Longitude):
            new = self._longitude + other._longitude
            fsp1 = self.format_specifier
            fsp2 = other.format_specifier
            if fsp1 != fsp2:
                if fsp1 < fsp2:
                    fsp = fsp2

        return Longitude(new, format_specifier=fsp, default=default)

    @property
    def valid(self):
        return self._valid

    @property
    def error(self):
        return self._error


class Latitude:
    def __init__(
        self, latitude: float, format_specifier: str = "2.4F", default: float = -90.0
    ):
        self._default = default
        self.format_specifier = format_specifier
        self._valid = True
        self._error = None
        self.latitude = latitude

    @property
    def latitude(self):
        return self._latitude

    @latitude.setter
    def latitude(self, value):
        # apply cyclic boundary conditions
        if value > 90.0:
            value = 180.0 - value
        if value < -90.0:
            value = -180.0 - value

        # if it is still out of range, then throw an error.
        if value < -90.0 or value > 90.0:
            self._valid = False
            self._error = "Latitude must be in range [-90,90]"
            value = self._default
            logger.warning(f"Latitude out of range {value}")
        self._latitude = value
        self._latitude = float(str(self))

    def __lt__(self, other):
        if self.latitude < other.latitude:
            return True
        else:
            return False

    def __eq__(self, other):
        if str(self) == str(other):
            return True
        else:
            return False

    def __str__(self):
        return f"{self._latitude:{self.format_specifier}}"

    def __repr__(self):
        return f"{self._latitude:{self.format_specifier}}"

    def __float__(self):
        return self._latitude

    def __add__(self, other):
        fsp = self.format_specifier
        default = self._default
        if isinstance(other, (float, int)):
            new = self._latitude + other
        elif isinstance(other, Latitude):
            new = self._latitude + other._latitude
            fsp1 = self.format_specifier
            fsp2 = other.format_specifier
            if fsp1 != fsp2:
                if fsp1 < fsp2:
                    fsp = fsp2

        return Latitude(new, format_specifier=fsp, default=default)

    @property
    def valid(self):
        return self._valid

    @property
    def error(self):
        return self._error


class HysplitDate:
    def __init__(self, date):
        """
        date: datetime.datetime OR string

        hysplit date strings are in format
        YY mm DD HH
        or
        YY mm DD HH MM

        sometimes they are all 00 which
        """

        self._valid = True
        self.date = None
        self.datestr = date

    def __str__(self):
        return self.datestr

    @property
    def datestr(self):
        return self._datestr

    @property
    def valid(self):
        return self._valid

    @datestr.setter
    def datestr(self, date):
        """
        if invalid date is entered then use 00 00 00 00 00

        for sample_start_time
        The previously specified hours of emission start at this time.
        An entry of zero's in the field, when input is read from a file,
        will also result in the selection of the default values that will correspond
        with the starting time of the meteorological data file.
        Day and hour are treated as relative to the file start when month is set to zero.

        Note - day is only checked to be sure it is between 0 and 31. It does not check to
        see if it is valid for the month/year.

        """
        if isinstance(date, str):
            date = date.strip()
            temp = date.split()
            if len(temp) == 4:
                date += " 00"
                temp = date.split()
            if len(temp) == 5:
                year = int(temp[0])
                month = int(temp[1])
                day = int(temp[2])
                hour = int(temp[3])
                minute = int(temp[4])
                if month > 12 or month < 0:
                    logger.warning("Invalid month")
                    self._valid = False
                # if month is 0 then month and hour treated as relative to file start.
                # thus hour could supposedly be > 24.
                # day being greater than 31 seems unlikely though.
                if month != 0:
                    if day > 31 or day < 0:
                        logger.warning("Invalid day")
                        self._valid = False
                    if hour > 23 or hour < 0:
                        logger.warning("Invalid hour")
                        self._valid = False
                    if minute > 59 or minute < 0:
                        logger.warning("Invalid minute")
                        self._valid = False
                if self._valid:
                    self._datestr = date
                    if year != 0 and month != 0 and day != 0:
                        try:
                            self.date = datetime.datetime.strptime(
                                date, "%y %m %d %H %M"
                            )
                        except ValueError as eee:
                            logger.warning(f"'WWW {eee} {date}")
                            self._valid = False
                    else:
                        self.date = None
            else:
                self._valid = False
                logger.warning("Invalid datestr")
        elif isinstance(date, datetime.datetime):
            self.date = date
            self._datestr = date.strftime("%y %m %d %H %M")
            # this will remove any seconds that may be on the datetime object.
            self.date = datetime.datetime.strptime(self._datestr, "%y %m %d %H %M")

        else:
            logger.warning(f"Invalid type for date {type(date)}")
            self._valid = False
        if not self._valid:
            logger.warning("Invalid date using 00 00 00 00 00")
            self._datestr = "00 00 00 00 00"
            self.date = None

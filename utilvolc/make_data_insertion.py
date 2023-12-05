""" Use VOLCAT data to write emit-times files for input into HYSPLIT 
--------------
Functions:
--------------
find_emit
get_emit_name_df
make_filename2
make_filename
make_cdump_filename
make_emit_filename
parse_filename

--------------
Class: EmitName 

Class: InsertVolcat
     add_vname: adds volcano name to instance variables
     find_fnames: finds list of volcat files based on datetime object 
     (can also use volcano ID)
     make_match: makes a string containing event datetime, image datetime, and 
     volcano id which is used for generating emitimes filenames and area filenames
     get_area: calculates area of lat/lon grid
     make_1D: makes 1D arrays of lat, lon, ash height, ash mass, and area
     write_emit: writes emitimes files

     rate is written in units of g/h

--------------
"""

import os
import datetime
from glob import glob
from math import pi

import numpy as np
import numpy.ma as ma
import pandas as pd
import xarray as xr
from utilhysplit import emitimes

from utilvolc import volcat
from utilvolc.volcat import VolcatName
from utilvolc import get_area

def make_1D_sub(dset):
    """Makes compressed 1D arrays of latitude, longitude, ash height,
    mass emission rate, and area
    For use in writing emitimes files. See self.write_emit()
    Also used for writing files to generate trajectories from (see utiltraj.trajectory_input_csv)
    Input:
    areafile: (boolean) True if area netcdf is created
    correct_parallax: use parallax corrected lat/lon (default = True)
    decode_times: boolean (default=False) Must be TRUE if using L1L2 netcdf files
    clip: boolean(default=True) Use data clipped around feature, reduces array size
    Output:
    list of tuples with (latitude,longitude,mass loading (g/m2) ,height (km) ,area (m2))
    """
    # 2023 APR 07 (AMC) updated to use functions in volcat.  

    #import numpy.ma as ma

    mass = volcat.get_mass(dset,clip=False)
    height = volcat.get_height(dset,clip=False)
    if volcat.check_pc(dset):
        lat = volcat.get_pc_latitude(dset,clip=False)
        lon = volcat.get_pc_longitude(dset,clip=False)
    else:
        lat = mass.latitude.values
        lon = mass.longitude.values     

    area = get_area.get_area(mass)
    area = area * 1e6 # convert from km^2 to m^2
    time = mass.time
    tlist = volcat.matchvals(lat,lon,mass,height,area)
    return tlist 


class EmitName(VolcatName):
    """
    class for names of Emit times files written from VOLCAT files
    """

    def __init__(self, fname):
        super().__init__(fname)

    def make_datekeys(self):
        #self.datekeys = [1, 2, 4, 5]
        self.datekeys = [1,2,None,None]

    def parse(self,fname):
        #fname = self.fname.replace('EMIT_','')
        #fname = self.fname.replace('cdump_','')
        #print('zzzzzzzzzzzzzzzzzz', fname)
        super().parse(fname)

    def make_keylist(self):
        self.keylist = ["file descriptor"]
        self.keylist.append("event vid")
        self.keylist.append("satellite platform")
        self.keylist.append("feature_id")
        self.keylist.append("layer")
        self.keylist.append("particles")

    def make_filename(self,volcat_fname, prefix, suffix):
        """Makes unique string for various filenames from volcat filename
        Inputs:
        fname: name of volcat file - just volcat file, not full path (string)
        prefix : str
        suffix : str
        Outputs:
        match: (string) of important identifiying information"""
        # Parse filename for image datetime, event datetime, and volcano id
        vname = volcat.get_name_class(volcat_fname)
        pidlist = [prefix]
        pidlist.append(vname.image_date_str)
        for key in ['event vid','satellite platform','feature_id']:
            if key in vname.vhash.keys():
               mstr = vname.vhash[key]  
               pidlist.append(mstr)
            else:
               print('WARNING key not found {}'.format(key))
        if suffix: 
           pidlist.append(suffix)
        match = str.join('_',pidlist)  
        match = match.replace('.nc','')
        match = match.replace('.','')
        return match

def find_emit(tdir,etype='EMIT'):
    """
    tdir : string
    Finds files in directory indicated by tdir which conform
    to the naming convention specified in EmitName class
    Return:
    vnlist : list of EmitName objects

    """
    vnlist = []
    if not os.path.isdir(tdir):
        return vnlist
    for fln in os.listdir(tdir):
        if etype in fln:
            try:
                vn = EmitName(fln)
            except Exception as eee:
                #print('find_emit error {} {}'.format(fln, eee))
                continue
            vnlist.append(vn)
    return vnlist


def get_emit_name_df(tdir,etype='EMIT'):
    """
    tdir : string
    Return:
    vdf : pandas dataframe with information found in names of Emit-times files.
    """
    tlist = find_emit(tdir,etype=etype)
    vlist = [x.vhash for x in tlist]
    vdf = pd.DataFrame(vlist)
    return vdf

# species tag could be p006p001p002p003 for multiple particle sizes
def make_cdump_filename(volcat_fname, speciestag, layertag):
    suffix = "{}_{}".format(layertag, speciestag)
    enn = EmitName(None)
    return enn.make_filename(volcat_fname, prefix="CDUMP", suffix=suffix)


def make_emit_filename(volcat_fname, speciestag, layertag):
    suffix = "{}_{}".format(layertag, speciestag)
    enn = EmitName(None)
    return enn.make_filename(volcat_fname, prefix="EMIT", suffix=suffix)


def parse_filename(ename):
    ehash = {}
    temp = ename.split("_")
    ehash["imagedate"] = temp[1]
    ehash["imagetime"] = temp[2]
    ehash["vid"] = temp[3]
    ehash["eventdate"] = temp[1]
    ehash["eventime"] = temp[2]
    ehash["layertag"] = temp[4]
    ehash["speciestag"] = temp[5]


class InsertVolcat:
    """
    Methods:
    add_vname: adds volcano name
    find_fnames: finds volcat filename, returns list of possible filenames
    make_match: makes a string containing event datetime, image datetime, and
    volcano id which is used for generating emitimes filenames and area filenames
    get_area: calculates the domain area, gridded
    make_1D: makes 1D array of lat, lon, ash height, ash mass, area
    write_emit: writes emitimes files
    """

    # get_area was moved to its own function in get_area.py


    def __init__(
        self,
        wdir='./',
        vdir=None,
        date_time=None,
        stradd="",
        duration="0010",
        pollnum=1,
        pollpercents=[1],
        layer=0.0,
        centered=False,
        vname=None,
        vid=None,
        fname=None,
    ):
        """
        Class of tools for inserting volcat data into hysplit
        -------------
        Inputs:
        wdir: working directory where emit-times files are written
        vdir: volcat directory - where volcat data files are located (string)
        date_time: date and time of volcat file to use (datetime object)
        stradd: QUICK FIX - file formatting has changed for Nishi Data (string)
        duration: hour and minute (HHMM) of emission - default is '0010' (string)
        pollnum: number of particle size bins - default is 1 (integer)
        pollpercents: percentage of particles for each size bin,
        represented as values from [0] to [1] - default is [1] (list)
        vname: volcano name - default is None (string)
        vid: volcano id - default is None (string)
        fname: str,list or None. name of volcat file (string)  or list of names. default is None
        Outputs:
        --------------
        """

        #if wdir[-1] != "/":
        #    wdir += "/"
        #self.wdir = wdir
        #if vdir[-1] != "/":
        #    vdir += "/"
        self.wdir = wdir
        self.vdir = vdir
        self.date_time = date_time
        self.stradd = stradd
        self.set_pollnum(pollnum, pollpercents)
        self.duration = duration
        self.vname = vname
        self.vid = vid
        self.layer = layer
        self.centered = centered

        if isinstance(fname, list):
            self.fnamelist = fname
        elif isinstance(fname, str):
            self.fnamelist = [fname]
       # else:
       #     self.fnamelist = self.find_fnames()
        self.emit_name = None

    def set_duration(self, duration):
        self.duration = duration

    def set_pollnum(self, pollnum, pollpercents):
        self.pollnum = pollnum
        if isinstance(pollpercents, (float, int)):
            pollpercents = [pollpercents]
        totpercents = np.array(pollpercents).sum()
        if np.abs(totpercents - 1) > 0.02:
            print("warning particle percentages do not add to 1 {}".format(totpercents))
        self.pollpercents = pollpercents

    def add_vname(self, vname):
        """Adds volcano name"""
        self.vname = vname

    def find_fnames(self):
        """Determines filename for volcat file based on vdir, date_time and vid"""
        # vfiles = "*.nc"
        volcframe = volcat.get_volcat_name_df(volc_dir)
        volcframe1 = volcframe.loc[
            (volcframe["event date"] == "s{}".format(self.date_time.strftime("%Y%j")))
            & (volcframe["event time"] == self.date_time.strftime("%H%M%S"))
        ]
        volclist = volcframe1["filename"].tolist()
        if self.vid != None:
            volclist = volcframe1.loc[
                volcframe1["volcano id"] == "v{}".format(self.vid), "filename"
            ].tolist()
        return volclist

    def make_1D_fromlist(
        self,
        fnamelist,
        correct_parallax=True,
        area_file=False,
        decode_times=False,
        clip=True,
    ):
        lat = np.array([])
        lon = np.array([])
        hgt = np.array([])
        mass = np.array([])
        area = np.array([])
        alist = np.array([])
        for iii, fname in enumerate(fnamelist):
            tlist = self.make_1D(
                fname,
                area_file=area_file,
                decode_times=decode_times,
            )
            if iii != 0:
                alist.extend(tlist)
            else:
                alist = tlist
        temp = list(zip(*alist))
        lat = temp[0]
        lon = temp[1]
        mass = temp[2]
        hgt = np.array(temp[3])*1000  #convert to meters
        area = temp[4]
        return lat, lon, hgt, mass, area

    def make_1D(self,fname,area_file=True,decode_times=False):

        dset = volcat.open_dataset(
            os.path.join(self.vdir,fname),
            correct_parallax=False,
            decode_times=decode_times,
        )

        # TODO may need to modify. what if have different time values?
        tval = dset.time.values
        if isinstance(tval,(list,np.ndarray)):
           tval = tval[0]
        self.date_time = pd.to_datetime(tval)
        print('DATE IS {}'.format(self.date_time))
        tlist =  make_1D_sub(dset) 
        dset.close()
        return tlist

    def make_tags(self):
        speciestag = "par{}".format(self.pollnum)
        layertag = "nolayer"
        if self.layer > 0.0:
            layertag = "mlayer"
            if self.centered:
                layertag = "clayer"
        return speciestag, layertag

    def set_layer(self, layer):
        self.layer = layer

    def check_for_file(self):
        if os.path.isfile(os.path.join(self.wdir, self.make_emit_filename())):
            return True
        else:
            return False

    def make_emit_filename(self):
        numfiles = len(self.fnamelist)
        speciestag, layertag = self.make_tags()
        speciestag = "f{}{}".format(numfiles, speciestag)
        filename = make_emit_filename(self.fnamelist[0], speciestag, layertag)
        return filename

    @staticmethod
    def adjust_time(time,round_time):
        ## rounds time to the nearest round_time minute.
        if round_time == 1: return time
        minute = time.minute

        def adjust_minute(minute,a): 
            if minute%a > a/2: return (minute + a) - (minute+a)%a
            else : return minute - minute%a

        minute2 = adjust_minute(minute,round_time)
        diff = minute2 - minute
        time2 = time + datetime.timedelta(hours=diff/60.0)
        return time2 
 

    def write_emit(
        self,
        heat="0.00e+00",
        layer=None,
        area_file=True,
        centered=None,
        #correct_parallax=True,
        decode_times=False,
        clip=True,
        verbose=False,
        min_line_number=0,
        min_mass=0,  # todo allow for only writing areas with mass above min value.
        round_time = 10, # round time to nearest 10's of minutes.
    ):
        """Writes emitimes file from volcat data.
        Inputs are created using self.make_1D()
        Uses instance variables: date_time, duration, par,
        emitimes file written to wdir

        Inputs:
        lat: 1D array of latitude
        lon: 1D array of longitude
        hgt: 1D array of ash top height
        mass: 1D array of ash mass
        area: 1D array of area
        heat: (string) default=0.00e+00
        layer: (float) height in meters of ash layer. A value of 0. means no layer, ash inserted at observed VOLCAT height.
        centered: (boolean) center ash layer on volcat height
        correct_parallax: (boolean) use parallax corrected lat lon values
        areafile: (boolean) uses area file if True, calculates area if False
        min_line_number: integer. If number of lines below this number then
                         emit-times file not written and False is returned.

        Returns:
        True if file written
        False if file not written

        """
        # Call make_1D() array to get lat, lon, height, mass, and area arrays
        # try:
        filename = self.make_emit_filename()
        lat, lon, hgt, mass, area = self.make_1D_fromlist(
            self.fnamelist,
            #correct_parallax=correct_parallax,
            area_file=area_file,
            decode_times=decode_times,
            clip=clip,
        )
        date_time = self.adjust_time(self.date_time, round_time)
        #print('ADJUSTED DATETIME', date_time)
        if layer:
            self.layer = layer
        else:
            layer = self.layer
        if isinstance(centered, bool):
            self.centered = centered
        else:
            centered = self.centered

        if verbose:
            print("{} number of lines".format(len(mass)))
            if len(mass) > 0:
                print("{:2e} max mass".format(np.max(mass)))
                print("{:2e} min mass".format(np.min(mass)))
        if len(mass) < min_line_number:
            return False

        # do not write empty emit-times files
        if len(mass) == 0:
            if verbose:
                print("---- Emitimes file NOT written: " + self.wdir + filename)
            return False

        const = 60.0 / float(self.duration)
        records = np.shape(hgt)[0] * self.pollnum

        if self.layer > 0.0:
            records = records * 2
        f = open(os.path.join(self.wdir, filename), "w")
        f.write("YYYY MM DD HH DURATION(HHMM) #RECORDS \n")
        f.write(
            "YYYY MM DD HH MM DURATION(HHMM) LAT LON HGT(m) RATE(g/hr) AREA(m2) HEAT(w) \n"
        )
        f.write("{:%Y %m %d %H} {} {}\n".format(date_time, self.duration, records))
        h = 0
        while h < len(hgt):
            i = 0
            while i < len(self.pollpercents):
                # AMR: 8/3/2021 - added this for layer flag
                # AMR: 8/23/2021 - added centered flag
                if layer > 0.0 and centered:
                    layhalf = float(layer) / 2.0
                    f.write(
                        "{:%Y %m %d %H %M} {} {:9.4f} {:10.4f} {:8.2f} {:.2E} {:.2E} {} \n".format(
                            date_time,
                            self.duration,
                            lat[h],
                            lon[h],
                            hgt[h] - float(layhalf),
                            mass[h] * area[h] * const * self.pollpercents[i],
                            area[h],
                            heat,
                        )
                    )
                    f.write(
                        "{:%Y %m %d %H %M} {} {:9.4f} {:10.4f} {:8.2f} {:.2E} {:.2E} {} \n".format(
                            date_time,
                            self.duration,
                            lat[h],
                            lon[h],
                            hgt[h] + float(layhalf),
                            0.0,
                            0.0,
                            heat,
                        )
                    )
                elif layer > 0.0 and not centered:
                    f.write(
                        "{:%Y %m %d %H %M} {} {:9.4f} {:10.4f} {:8.2f} {:.2E} {:.2E} {} \n".format(
                            date_time,
                            self.duration,
                            lat[h],
                            lon[h],
                            hgt[h] - float(layer),
                            mass[h] * const * self.pollpercents[i],
                            area[h],
                            heat,
                        )
                    )
                    f.write(
                        "{:%Y %m %d %H %M} {} {:9.4f} {:10.4f} {:8.2f} {:.2E} {:.2E} {} \n".format(
                            date_time,
                            self.duration,
                            lat[h],
                            lon[h],
                            hgt[h],
                            0.0,
                            0.0,
                            heat,
                        )
                    )
                else:
                    f.write(
                        "{:%Y %m %d %H %M} {} {:9.4f} {:10.4f} {:8.2f} {:.2E} {:.2E} {} \n".format(
                            date_time,
                            self.duration,
                            lat[h],
                            lon[h],
                            hgt[h],
                            mass[h] * area[h] * const * self.pollpercents[i],
                            area[h],
                            heat,
                        )
                    )
                i += 1
            h += 1
        f.close()
        if verbose:
            print("Emitimes file written: " + self.wdir + filename)
        return os.path.join(self.wdir,filename)

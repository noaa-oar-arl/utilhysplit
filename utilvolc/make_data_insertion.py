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
--------------
"""

import os
from datetime import datetime
from glob import glob
from math import pi

import numpy as np
import numpy.ma as ma
import pandas as pd
import xarray as xr
from utilhysplit import emitimes

from utilvolc import volcat
from utilvolc.volcat import VolcatName


class EmitName(VolcatName):
    """
    class for names of Emit times files written from VOLCAT files
    """

    def __init__(self, fname):
        super().__init__(fname)

    def make_datekeys(self):
        #self.datekeys = [1, 2, 4, 5]
        self.datekeys = [2,3,None,None]

    def parse(self,fname):
        #fname = self.fname.replace('EMIT_','')
        #fname = self.fname.replace('cdump_','')
        #print('zzzzzzzzzzzzzzzzzz', fname)
        super().parse(fname)

    def make_keylist(self):
        self.keylist = ["file descriptor"]
        self.keylist.append("algorithm name")
        self.keylist.append("image date")  # should be image date (check)
        self.keylist.append("image time")
        self.keylist.append("volcano id")
        # may not need the event date in filename
        #self.keylist.append("event date")  # should be event date (check)
        #self.keylist.append("event time")
        self.keylist.append("satellite platform")
        self.keylist.append("feature id")
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
        vname = VolcatName(volcat_fname)
        pidlist = [prefix]
        for key in self.keylist:
            if key in vname.vhash.keys():
               if key == 'edate':
                  mstr = vname.vhash["edate"].strftime(vname.event_dtfmt)
                  pidlist.append(mstr)
               elif key == 'idate':
                  mstr = vname.vhash["idate"].strftime(vname.image_dtfmt)
                  pidlist.append(mstr)
               else: 
                  mstr = vname.vhash[key]  
                  pidlist.append(mstr)
        pidlist.append(suffix)
        match = str.join('_',pidlist)  
        return match

def find_emit(tdir):
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
        try:
            vn = EmitName(fln)
        except Exception as eee:
            #print('find_emit error {} {}'.format(fln, eee))
            continue
        vnlist.append(vn)
    return vnlist


def get_emit_name_df(tdir):
    """
    tdir : string
    Return:
    vdf : pandas dataframe with information found in names of Emit-times files.
    """
    tlist = find_emit(tdir)
    vlist = [x.vhash for x in tlist]
    vdf = pd.DataFrame(vlist)
    return vdf


# def make_filename2(volcat_fname, prefix, suffix):
#    """Makes unique string for various filenames from volcat filename
#    Inputs:
#    fname: name of volcat file - just volcat file, not full path (string)
#    Outputs:
#    match: (string) of important identifiying information"""
#    # Parse filename for image datetime, event datetime, and volcano id
#    vname = VolcatName(volcat_fname)
#    pid1 = vname.vhash["idate"].strftime(vname.image_dtfmt)
#    pid2 = vname.vhash["volcano id"]
#    match = "{}_{}".format(pid1, pid2, pid3)
#    filename = "{}_{}_{}".format(prefix, match, suffix)
#    return filename


# switch to using make_filename method in EmitName class.
#def make_filename_old(volcat_fname, prefix, suffix):
#    """Makes unique string for various filenames from volcat filename
#    Inputs:
#    fname: name of volcat file - just volcat file, not full path (string)
#    Outputs:
#    match: (string) of important identifiying information"""
#    # Parse filename for image datetime, event datetime, and volcano id
#    vname = VolcatName(volcat_fname)
#    #pid1 = vname.vhash["idate"].strftime(vname.image_dtfmt)
#    pid2 = vname.vhash["volcano id"]
#    pid3 = vname.vhash["edate"].strftime(vname.event_dtfmt)
#    pid4 = vname.vhash['satellite platform']
#    pid5 = vname.vhash['feature id']
#    match = "{}_{}_{}_{}".format(pid2, pid3, pid4, pid5)
    # parsedf = volcat_fname.split('_')
    # match = parsedf[4]+'_'+parsedf[5]+'_'+parsedf[7]+'_'+parsedf[12]+'_'+parsedf[13]
#    filename = "{}_{}_{}".format(prefix, match, suffix)
#    return filename


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


    TODO - should get_area method be replaced by function in get_area.py
    """

    def __init__(
        self,
        wdir,
        vdir,
        date_time,
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
        wdir: working directory - base directory(string)
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

        if wdir[-1] != "/":
            wdir += "/"
        self.wdir = wdir
        if vdir[-1] != "/":
            vdir += "/"
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
        else:
            self.fnamelist = find_fnames()
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

    def get_area(
        self, dset, write=False, correct_parallax=True, decode_times=False, clip=True
    ):
        """Calculates the area (km^2) of each volcat grid cell
        Converts degress to meters using a radius of 6378.137km.
        Input:
        write: boolean (default: False) Write area to file
        correct_parallax: boolean (default: True) Use parallax correct lat/lon values
        decode_times: boolean (default: False) Must be TRUE if using L1L2 netcdf files
        clip: boolean (default: True) Use clipped array around data, reduces domain
        output:
        area: xarray containing gridded area values
        """
        d2r = pi / 180.0  # convert degress to radians
        d2km = 6378.137 * d2r  # convert degree latitude to kilometers

        # if self.fnamelist:
        #    dset = volcat.open_dataset(
        #        self.vdir + fname,
        #        correct_parallax=correct_parallax,
        #        decode_times=decode_times,
        #    )
        # else:
        #    print("ERROR: Need volcat filename!")
        # Extracts ash mass array (two methods - one is smaller domain around feature)
        # Removes time dimension
        if clip == True:
            mass = volcat.get_mass(dset)[0, :, :]
            # mass = volcat.get_mass(dset)
        else:
            mass = dset.ash_mass_loading[0, :, :]
        lat = mass.latitude
        lon = mass.longitude
        latrad = lat * d2r  # Creating latitude array in radians
        coslat = np.cos(latrad) * d2km * d2km
        shape = np.shape(mass)

        # Make shifted lat and shifted lon arrays to use for calculations
        lat_shift = lat[1:, :].values
        lon_shift = lon[:, 1:].values
        # Adding row/column of nans to shifted lat/lon arrays
        to_add_lon = np.empty([shape[0]]) * np.nan
        to_add_lat = np.empty([shape[1]]) * np.nan
        # Back to xarray for calculations
        lat2 = xr.DataArray(np.vstack((lat_shift, to_add_lat)), dims=["y", "x"])
        lon2 = xr.DataArray(np.column_stack((lon_shift, to_add_lon)), dims=["y", "x"])

        # area calculation
        area = abs(lat - lat2) * abs(abs(lon) - abs(lon2)) * coslat
        area.name = "area"
        area.attrs["long_name"] = "area of each lat/lon grid box"
        area.attrs["units"] = "km^2"
        # Reformatting array attributes
        if write == True:
            directory = self.wdir + "area/"
            match = ""
            # match = self.make_match()
            if correct_parallax == True:
                areafname = "area_" + match + "_pc.nc"
            else:
                areafname = "area_" + match + ".nc"
            # print(directory + areafname)
            area.to_netcdf(directory + areafname)
        dset.close()
        return area

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

        for fname in fnamelist:
            print("NAME", fname)
            lat2, lon2, hgt2, mass2, area2 = self.make_1D(
                fname,
                correct_parallax=correct_parallax,
                area_file=area_file,
                decode_times=decode_times,
                clip=clip,
            )
            lat = np.append(lat, lat2)
            lon = np.append(lon, lon2)
            hgt = np.append(hgt, hgt2)
            mass = np.append(mass, mass2)
            area = np.append(area, area2)
        return lat, lon, hgt, mass, area

    def make_1D(
        self,
        fname,
        correct_parallax=True,
        area_file=True,
        decode_times=False,
        clip=True,
    ):
        """Makes compressed 1D arrays of latitude, longitude, ash height,
        mass emission rate, and area
        For use in writing emitimes files. See self.write_emit()
        Input:
        areafile: (boolean) True if area netcdf is created
        correct_parallax: use parallax corrected lat/lon (default = True)
        decode_times: boolean (default=False) Must be TRUE if using L1L2 netcdf files
        clip: boolean(default=True) Use data clipped around feature, reduces array size
        Output:
        lat: 1D array of latitude
        lon: 1D array of longitude
        hgt: 1D array of ash top height
        mass: 1D array of ash mass (grams)
        area: 1D array of area  (m2)
        """
        import numpy.ma as ma

        dset = volcat.open_dataset(
            self.vdir + fname,
            correct_parallax=correct_parallax,
            decode_times=decode_times,
        )

        # if self.fname:
        #    dset = volcat.open_dataset(
        #        self.vdir + self.fname,
        #        correct_parallax=correct_parallax,
        #        decode_times=decode_times,
        #    )
        # else:
        #    print("ERROR: Need volcat filename!")
        # Extracts ash mass array - Removes time dimension
        # TO DO - there is some error associated with the clipping.
        if clip == True:
            mass0 = volcat.get_mass(dset).isel(time=0)
            height0 = volcat.get_height(dset).isel(time=0)
        else:
            mass0 = dset.ash_mass_loading.isel(time=0)
            height0 = dset.ash_cloud_height.isel(time=0)
        # Making 0. in arrays nans
        mass = mass0.where(mass0 > 0.0)
        height = height0.where(height0 > 0.0)

        # Finds and area files. If no file detected, produces warning
        # areadir = self.vdir+'Area/'
        if area_file:
            areadir = self.wdir + "area/"
            # match = self.make_match()
            if correct_parallax == True:
                arealist = glob(areadir + "*_pc.nc")
            else:
                arealist = glob(areadir + "*0.nc")
            areafile = [f for f in arealist if match in f]
            if areafile:
                areafile = areafile[0]
                areaf = xr.open_dataset(areafile).area
        else:
            areaf = self.get_area(dset, clip=clip)

        # Calculating mass - rate is (mass/hr) in HYSPLIT
        # ashmass is g/m2. Multiply by m2 to get grams.
        areaf = areaf * 1000.0 * 1000  # converting area from (km^2) to (m^2)
        ashmass = mass * areaf
        latitude = mass.latitude
        longitude = mass.longitude

        # AMC added close statements.
        areaf.close()
        dset.close()

        # Creating compressed 1-d arrays (removing nans) to prepare for file writing
        hgt_nan = ma.compressed(height)  # arra with nan values to remove
        hgt = hgt_nan[~np.isnan(hgt_nan)]  # removing nan values
        mass = ma.compressed(ashmass)[~np.isnan(hgt_nan)]  # removing nan values
        lat = ma.compressed(latitude)[~np.isnan(hgt_nan)]  # removing nan values
        lon = ma.compressed(longitude)[~np.isnan(hgt_nan)]  # removing nan values
        area = ma.compressed(areaf)[~np.isnan(hgt_nan)]  # removing nan values
        hgt = hgt * 1000.0  # Moving hgt from (km) to (m)
        return lat, lon, hgt, mass, area

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

    def write_emit(
        self,
        heat="0.00e+00",
        layer=None,
        area_file=True,
        centered=None,
        correct_parallax=True,
        decode_times=False,
        clip=True,
        verbose=False,
        min_line_number=0,
        min_mass=0,  # todo allow for only writing areas with mass above min value.
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
            correct_parallax=correct_parallax,
            area_file=area_file,
            decode_times=decode_times,
            clip=clip,
        )
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
        f = open(self.wdir + filename, "w")
        f.write("YYYY MM DD HH DURATION(HHMM) #RECORDS \n")
        f.write(
            "YYYY MM DD HH MM DURATION(HHMM) LAT LON HGT(m) RATE(g/hr) AREA(m2) HEAT(w) \n"
        )
        f.write("{:%Y %m %d %H} {} {}\n".format(self.date_time, self.duration, records))
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
                            self.date_time,
                            self.duration,
                            lat[h],
                            lon[h],
                            hgt[h] - float(layhalf),
                            mass[h] * const * self.pollpercents[i],
                            area[h],
                            heat,
                        )
                    )
                    f.write(
                        "{:%Y %m %d %H %M} {} {:9.4f} {:10.4f} {:8.2f} {:.2E} {:.2E} {} \n".format(
                            self.date_time,
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
                            self.date_time,
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
                            self.date_time,
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
                            self.date_time,
                            self.duration,
                            lat[h],
                            lon[h],
                            hgt[h],
                            mass[h] * const * self.pollpercents[i],
                            area[h],
                            heat,
                        )
                    )
                i += 1
            h += 1
        f.close()
        if verbose:
            print("Emitimes file written: " + self.wdir + filename)
        return True

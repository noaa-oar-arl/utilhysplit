import sys
import re
import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import urllib.request, urllib.error, urllib.parse
import xml.etree.ElementTree as ET
import shapely.geometry as sgeo
from utilhysplit import geotools

"""
Read and process volcanic ash advisories in vAA forma

This code developed at NOAA Air Resources laboratory
PGRMMR: Alice Crawford 

Classes
class iwxxmFile: loads the iwxxm file from url
class iwxxmVAA:  parses the iwxxm file stores information in object attributes.
class WashingtonPage: gets url location of iwxxm file
class iwxxmCollection: processes collections of iwxxmVAA

2023 22 MAR  (AMC) added support for multiple polygons in obs or forecast.

"""


# ---------------------------------------------------------------
## Helper functions for iwxxmVAA class to parse the xml format.
# ---------------------------------------------------------------

def etree_to_dict(t):
    """
    This works as long as all children have unique keys.
    If they don't, then they will be overwritten.
    This occurrs when there are multiple polygons. 
    """
    if type(t) is ET.ElementTree: return etree_to_dict(t.getroot())
    ## concatenates 3 dictionaries together.
    return {
        **t.attrib,
        'attribute_text': t.text,
        **{e.tag: etree_to_dict(e) for e in t}
    }


def testel(root):
    try:
        for child in root:
            tag = child.tag
    except:
        return False
    return True

def listel(root):
    clist = []
    try:
        for child in root:
            print("tag", child.tag, "text", child.text)
            clist.append(child)
    except:
        print("cannot list for", type(root))
    return clist

def findel(root, tag, verbose=False):
    for child in root:
        if tag.lower() in child.tag.lower():
            return child
    if verbose:
        print("{} not found".format(tag))
    return None


def listfindel(root, tag):
    clist = []
    for child in root:
        if tag.lower() in child.tag.lower():
            clist.append(child)
    # print('{} not found'.format(tag))
    return clist


# ---------------------------------------------------------------
# ---------------------------------------------------------------
def create_new_collection(collection, vname, daterange=None, cname=None):
    newlist = collection.get_by_volcano(vname)
    if isinstance(daterange, list):
        blist = collection.get_by_date(daterange)
        # only use elements that appear in both lists.
        clist = [x for x in blist if x in newlist]
        newlist = clist
    name = "Collection of iwxxm files for {}".format(vname)
    newcollection = iwxxmCollection(name)
    newcollection.ilist = newlist
    return newcollection


class iwxxmCollection:
    """
    processing collection of vaas.
    """
    def __init__(self, name="Collection of iwxxm files"):
        """
        Attributes
        ----------
        name : str : descriptive
        ilist : list : list of iwxxmVAA objects

        Methods
        -------
        Methods work to either
        1. alter the self.ilist by adding or removing objects
        2. filter the list of iwxxm objects by certain attributes
           of the object like volcano name.
        3. compare one iwxxm object in the list to another. e.g. comparing
           observation from one object to forecast at same time in another object.
        4. summarize information in the ilist. e.g. list of volcano names, time series of areas.

        """
        self.name = "Collection of iwxxm files"
        self.ilist = []

    # METHODS to  add or remove objects from ilist
    def reset(self):
        self.ilist = []

    def add_objects(self, ilist):
        self.ilist.extend(ilist)
        # remove any duplicates
        self.ilist = list(set(self.ilist))
        # sort by issue date
        self.ilist.sort()

    def add_files(self, xlist):
        for xurl in xlist:
            try:
                xfile = iwxxmFile(xurl)
                vobj = iwxxmVAA(fstr=str(xfile))
            except Exception as eee:
                print("WARNING could not load {}".format(xurl))
                print(str(eee))
                print("-------------------------")
                continue
            self.ilist.append(vobj)
        # remove any duplicates
        self.ilist = list(set(self.ilist))
        # sort by issue date
        self.ilist.sort()

    # ----------------------------------------
    # METHODS to filter the ilist
    def get_by_volcano(self, vname):
        sublist = [x for x in self.ilist if x.volcano_name == vname]
        return sublist

    def get_by_date(self, daterange):
        sublist = [x for x in self.ilist if x.issueTime <= daterange[1]]
        sublist = [x for x in sublist if x.issueTime >= daterange[0]]
        return sublist

    # ----------------------------------------
    # METHODS to  summarize information in the objects
    def get_volcano_list(self):
        vnames = [x.volcano_name for x in self.ilist]
        unique_names = list(set(vnames))
        return unique_names

    def get_date_range(self):
        """
        Returns largest and smallest issue Time.
        """
        dates = [x.issueTime for x in self.ilist]
        return [np.min(dates), np.max(dates)]

    def get_issue_frequency(self, plot=True):
        dates = [x.issueTime for x in self.ilist]
        dates.sort()
        diff = [y - x for x, y in zip(dates, dates[1:])]
        times = [x.days * 24 + x.seconds / 3600.0 for x in diff]
        if plot:
            plt.plot(dates[1:], times, "-b.")
            plt.xticks(rotation=45)
            plt.tight_layout()
            ax = plt.gca()
            ax.set_ylabel("Time since last vaa issued (h)")
        return dates, times

    def time_series(self, vname, plot=True, ax=None, ptype='area'):
        if not ax:
            fig = plt.figure(1)
            ax = fig.add_subplot(1, 1, 1)
        sublist = self.get_by_volcano(vname)
        obstime = []
        obsarea = []
        obsarea2 = []
        f1time = []
        f1area = []
        f2time = []
        f2area = []
        f3time = []
        f3area = []
        for vaa in sublist:
            if ptype=='area':
                aaa = vaa.get_areas()
                ylabel = 'area (arbitrary unit)'
            elif ptype=='distance':
                aaa = vaa.get_distance()
                ylabel = 'distance (km)'
            elif ptype=='height':
                aaa = vaa.get_height()
                ylabel = 'Altitude (FL)'
            elif ptype=='bearing':
                aaa,bbb = vaa.get_bearing()
                ylabel = 'direction (degrees)'
            time = vaa.get_times()
            if time[0]:
                obstime.append(time[0])
                obsarea.append(aaa[0])
            if time[1]:
                f1time.append(time[1])
                f1area.append(aaa[1])
            if time[2] and len(aaa)>2:
                f2time.append(time[2])
                f2area.append(aaa[2])
            if time[3] and len(aaa)>3:
                f3time.append(time[3])
                f3area.append(aaa[3])
            if ptype=='bearing':
                obsarea2.append(bbb[0])
        if plot:
            sns.set()
            try:
                ax.plot(obstime, obsarea, "--k.")
            except:
                print("Failed")
                return obstime, obsarea
            if ptype=='bearing':
               ax.plot(obstime,obsarea2,"--k.")
            ax.plot(f1time, f1area, "r.")
            ax.plot(f2time, f2area, "b+")
            ax.plot(f3time, f3area, "yx")
            ax.set_ylabel(ylabel)
            plt.xticks(rotation=45)
            plt.tight_layout()
        return obstime, obsarea

    # METHODS to pick out a certain VAA

    # def match_vaa_by_time(intime, dt):
    #    if isinstance(dt,(int,float)):
    #       dt = datetime.timedelta(hours=dt)
    #    for iii, vaa in enumerate(sublist):
    #        timelist = vaa.get_times()
    #        ctime = timelist[0]
    #        #print('{} checking {}  to {}'.format(iii, intime, ctime))
    #        if ctime > intime:
    #           if ctime-intime < dt:
    #              matches.append(iii)
    #        elif intime <= ctime:

    # METHODS to compare information in different objects

    def check_forecasts(self, vname, fnum=1, plot=True):
        """
        RETURNS
        matches : dict : key is the index of the vaa with the  observation.
                         value is list of indices of the vaas with forecast at matching times.
        e.g. {8: [11]}
        means sublist[8] vaa has observation for same time as forecast in sublist[11] vaa.
        The vaas are generally ordered from most recent to oldest.
        """
        sublist = self.get_by_volcano(vname)
        matches = {}
        ### Want to find forecast times which overlap with observed times.
        for iii, vaa in enumerate(sublist):
            time = vaa.get_times()
            obstime = time[0]
            dt = datetime.timedelta(hours=0.5)
            matches[iii] = self.find_time_match(vname, obstime, fnum, dt)
        return matches

    def plot_checks(self, vname, obsiii, flist, fnum):

        # slist = ['Forecast0','Forecast1','Forecast2']
        # sname = slist[fnum]
        sublist = self.get_by_volcano(vname)
        polylist = []

        vaa0 = sublist[obsiii]
        #vaa0.plot_vaa()
        #plt.show()
        time0 = vaa0.get_times()[0]
        print("Reference time {}".format(time0))
        polylist.append(vaa0.get_poly_list()[0])
        print('Remark ------')
        print(vaa0.remarks)

        for fii in flist:
            vaaf = sublist[fii]
            #vaaf.plot_vaa()
            #plt.show()
            timef = vaaf.get_times()[fnum]
            obsf = vaaf.get_times()[0]
            print(
                "obstime {} Forecast time {} DIFF {}".format(obsf, timef, timef - obsf)
            )
            polylist.append(vaaf.get_poly_list()[fnum])
            print('Remark ------')
            print(vaaf.remarks)

        clr = ["-k", "-r", "-b", "-c"]
        for iii, poly in enumerate(polylist):
            for ppp in poly:
                if not ppp.is_empty:
                    plt.plot(*ppp.exterior.xy, clr[iii])
                else:
                    print("Empty polygon")
        # plt.plot(vlon,vlat,'^y',MarkerSize=10)

    def find_time_match(self, vname, intime, forecast, dt):
        """
        vname : str : volcano name
        intime : datetime.datetime object :
        forecast : integer (1,2,3) : forecast time to match (6,12,18h)
        dt : float or datetime.timedelta object : acceptable time difference for match.
        RETURNS
        -------
        matches : dict : key is the index of the vaa with the  observation.
                         value is list of indices of the vaas with forecast at matching times.
        """

        if isinstance(vname, str):
            sublist = self.get_by_volcano(vname)
        else:
            sublist = self.ilist
        matches = []
        if isinstance(dt, (int, float)):
            dt = datetime.timedelta(hours=dt)
        for iii, vaa in enumerate(sublist):
            timelist = vaa.get_times()
            ctime = timelist[forecast]
            if ctime > intime:
                if ctime - intime < dt:
                    matches.append(iii)
            elif intime <= ctime:
                if intime - ctime <= dt:
                    matches.append(iii)
        return matches


class WashingtonPage:
    """
    Scrape the washington vaac page for xml file locations
    """

    def __init__(self,year=None):
        #self.page = "https://www.ssd.noaa.gov/VAAC/messages.html"
        self.page = "https://www.ospo.noaa.gov/Products/atmosphere/vaac/messages.html"
        self.xmlpage = "https://www.ospo.noaa.gov/Products/atmosphere/vaac/volcanoes/xml_files/"
        self.archivepage = "https://www.ospo.noaa.gov/Products/atmosphere/vaac/"
        if (isinstance(year,int)):
           self.get_archive(year)
        self.hcontent = ""
        self.hcontent_loaded = False
        self.xlist = None

    def get_archive(self, year):
        self.page = '{}{}.html'.format(self.archivepage,year)

    def read(self):
        # read the web page
        hresponse = urllib.request.urlopen(self.page)
        hcode = hresponse.getcode()
        ##header = hcontent.getheaders()
        hcontent = hresponse.read()
        self.hcontent = hcontent.decode("utf-8")
        if hcode == 200:
            self.hcontent_loaded = True
        return hcode

    def find_xml(self):
        flist = []
        rlist = []
        for match in re.finditer("\.xml", self.hcontent):
            end = match.end()
            #print(self.hcontent[match.start():match.end()])
            substr = self.hcontent[end - 180 : end]
            submatch = re.search("FV", substr)
            xmlfile = self.xmlpage +  substr[submatch.start() :]
            rlist.append(substr[0 : submatch.start()])
            flist.append(xmlfile)
        self.xlist = list(zip(rlist, flist))

    # the web page changed in January 2022. this no longer works.
    def find_xml_old(self):
        flist = []
        rlist = []
        for match in re.finditer("\.xml", self.hcontent):
            end = match.end()
            substr = self.hcontent[end - 180 : end]
            submatch = re.search("https", substr)
            xmlfile = substr[submatch.start() :]
            rlist.append(substr[0 : submatch.start()])
            flist.append(xmlfile)
        self.xlist = list(zip(rlist, flist))

    def get_xml_list(self, vname=None):
        if vname:
            xlist = [x[1] for x in self.xlist if vname in x[0]]
        else:
            xlist = [x[1] for x in self.xlist]
        return xlist


class iwxxmFile:
    """
    read in iwxxm file
    """
    def __init__(self, url):
        self.page = url
        self.hcontent = ""
        self.hcode = self.read()

    def read(self):
        hresponse = urllib.request.urlopen(self.page)
        hcode = hresponse.getcode()
        self.hcontent = hresponse.read()
        if hcode == 200:
            self.hcontent_loaded = True
        return hcode

    def __str__(self):
        """
        Can be used as input to iwxxmVAA class.
        """
        return self.hcontent.decode("utf-8")

    def write(self, oname="test.xml"):
        if self.hcontent_loaded:
            with open(oname, "wb") as fid:
                fid.write(self.hcontent)
        else:
            print("content not loaded")


def get_poly(vhash):
    if "exterior" in vhash.keys():
        perim = vhash["exterior"]
        poly = sgeo.Polygon(perim)
    else:
        poly = sgeo.Polygon()
    return poly



def vaa2traj(vaa,inp={'WORK_DIR':'./','JOBID':'777'},dt=0):
    inp['jobname'] = vaa.volcano['latitude']
    inp['latitude'] = vaa.volcano['latitude']
    inp['longitude'] = vaa.volcano['longitude']
    inp['bottom'] = vaa.volcano['elevation']*0.3048
    top = vaa.get_height()[0]
    inp['top'] = top*30.48 + 2000
    inp['start_date'] = vaa.obs['date']-datetime.timedelta(hours=dt)
    inp['forecastDirectory']='/pub/forecast/'
    inp['archivesDirectory']='/pub/archives/'
    inp['meteorologicalData']='gfs0p25'
    inp['durationOfSimulation']=-12
    inp['HYSPLIT_DIR'] = '/hysplit-users/alicec/hdev/'
    inp['CONVERT_EXE'] = 'convert'
    inp['GHOSTSCRIPT_EXE']='ghostscript'
    inp['PYTHON_EXE']='python'    
    inp['DATA_DIR'] = inp['WORK_DIR']
    return inp


class iwxxmVAA:
    """
    class for parsing the VAAs in iwxxm format.
    Needs to be tested / modified when more than one polygon.
    utilizes xml.etree.ElementTree 
    
    """

    def __init__(self, fname=None, fstr=None, verbose=False):
        """
        fname : str : name of a file
        fstr  : str : contents of an iwxxm vaa
         
        reads in information and converts it into class 
        attributes listed below. attributes are generally
        strings or dictionaries.

        __hash__, __eq__, and __lt__ are defined so lists of 
        these objects can be sorted and compared.

        """

        self.strfmt = "%Y-%m-%dT%H:%M:%SZ"
        if fname:
            self.source = fname
            self.tree = ET.parse(fname)
            self.root = self.tree.getroot()
        elif fstr:
            self.source = fstr
            self.root = ET.fromstring(fstr)
        else:
            raise Exception("input strings are empty")
        self.vaaroot = self.get_vaaroot()

        # str
        self.bulletinIdentifier = self.get_bulletin_id()

        # str  This is used in the __eq__ and __hash__ functions
        # to uniquely identify the object.
        self.__bid = self.bulletinIdentifier

        # str
        self.issueTime = self.get_issue_time()

        # str
        self.issuingVolcanicAshAdvisoryCentre = self.get_vaac()

        # these are all dictionaries.
        self.volcano = self.get_volcano(verbose=verbose)
        self.details = self.get_details()
        self.obs = self.get_obs(verbose=verbose)
        self.forecast = self.get_all_forecasts(verbose=verbose)

        # str
        self.remarks = self.get_remarks()

        # str
        self.volcano_name = self.volcano["name"]
        self.volcano_id = self.volcano["volcanoID"]

    def __hash__(self):
        return hash(self.__bid)

    def __eq__(self, other):
        if self.__bid == other.__bid:
            return True
        else:
            return False

    def __lt__(self, other):
        # define so a list of vojb can be sorted by the issue date.
        if self.issueTime != other.issueTime:
            return self.issueTime < other.issueTime

    def get_times(self):
        dkey = "date"
        timelist = []

        if dkey in self.obs.keys():
            timelist.append(self.obs[dkey])
        else:
            timelist.append(None)
        for fcst in ["Forecast0", "Forecast1", "Forecast2"]:
            if dkey in self.forecast[fcst].keys():
                timelist.append(self.forecast[fcst][dkey])
            else:
                timelist.append(None)
        return timelist

    def get_areas(self):
        """
        The areas are computed from the polygons with latitude longitude
        vertices so the units are weird. Need to do some further processing
        to get real areas.
        """
        areas = []
        plist = self.get_poly_list()
        for iii, poly in enumerate(plist):
            area = 0
            for ppp in poly:
                if not ppp.is_empty:
                    area += ppp.area
            areas.append(area)
        return areas

    def get_poly_list(self):
        """
        returns list of lists.
        """
        keylist = ['Forecast0','Forecast1','Forecast2']
        plist = [[get_poly(x) for x in self.obs['poly']]]
        for key in keylist:
            plist.append([get_poly(x) for x in self.forecast[key]['poly']])
        return plist

    
    def get_height(self):
        def FL2m(fhash):
            value = fhash['value']
            #if fhash['uom']=='FL': value=value*30.48
            return value

        def htfunc(phash):
            astr = 'ashCloudExtent'
            if astr in phash.keys():
                hlist = [FL2m(phash[astr])]
            else:
                hlist = [np.nan]
            return hlist
        if self.obs['numpoly']>0:
            polygons = self.obs['poly']
            plist = []
            for poly in polygons:
                plist.append(htfunc(poly))
            # use the maximum observed polygon height.
            hlist = [np.max(plist)]
        else:
            hlist = [np.nan]

        alist = ['Forecast0','Forecast1','Forecast2']
        for aaa in alist:
            if not aaa in self.forecast.keys(): continue
            else: forecast = self.forecast[aaa]
            if forecast['numpoly']<1: 
               hlist.append(np.nan)
               continue
            polygons = forecast['poly']
            plist = []
            for poly in polygons:
               plist.append(htfunc(poly))
            hlist.append([np.max(plist)][0])
        # use the maximum observed polygon height.


            #if astr in self.forecast[aaa].keys():
            #   hlist.append(FL2m(self.forecast[aaa][astr]))
            #else:
            #   hlist.append(np.nan)
        return hlist


    def get_bearing(self):
        plist = self.get_poly_list()
        pnt = sgeo.Point(self.volcano['longitude'],self.volcano['latitude'])
        mlist = []
        alist = []
        mindiff=10
        for poly in plist:
            # if multiple polygons just use first one.
            if isinstance(poly,list):
               if poly:
                  poly = poly[0]
               else:
                  mlist.append(np.nan)
                  alist.append(np.nan)
                  continue
                   
            bearing = 0
            aaa = 400
            iii=0
            for coord in poly.exterior.coords:
                pnt2 = sgeo.Point(coord[0],coord[1])
                bbb = geotools.bearing(pnt,pnt2)
                diff = geotools.distance(pnt,pnt2)
                #bearing += geotools.bearing(pnt,pnt2)
                #iii+=1
                if bbb>bearing: bearing=bbb
                if bbb<aaa: 
                   if diff > mindiff:
                       aaa=bbb
                iii=1
            if iii>0:
                mlist.append(bearing/iii)  
                alist.append(aaa/iii)  
            else:
                mlist.append(np.nan)
                alist.append(np.nan)
        return mlist, alist

    def get_distance(self,verbose=False):
        plist = self.get_poly_list()
        pnt = sgeo.Point(self.volcano['longitude'],self.volcano['latitude'])
        mlist = []
        for poly in plist:
            # if multiple just use first one)
            if isinstance(poly,list): poly=poly[0]
            maxdiff = 0
            for coord in poly.exterior.coords:
                pnt2 = sgeo.Point(coord[0],coord[1])
                diff = geotools.distance(pnt,pnt2)
                if verbose: print(pnt, pnt2, diff)
                if diff>maxdiff:maxdiff=diff
            mlist.append(maxdiff)  
            if verbose: print('-----')
        return mlist

    def plot_vaa(self, ax=None, polylist=None, legend=True):
        """ """
        if not ax:
            fig = plt.figure(1)
            ax = fig.add_subplot(1, 1, 1)
        vlat = self.volcano["latitude"]
        vlon = self.volcano["longitude"]

        times = self.get_times()
        tstr = []
        for ttt in times:
            tstr.append(str(ttt))

        clr = ["-k", "-r", "-b", "-c"]
        plist = self.get_poly_list()
        if not polylist: polylist = np.arange(0,len(plist))
        for iii, poly in enumerate(self.get_poly_list()):
            if iii not in polylist: continue
            for jjj, ppp in enumerate(poly):
                if not ppp.is_empty:
                    if jjj==0: lbl = tstr[iii]
                    else: lbl=None
                    ax.plot(*ppp.exterior.xy, clr[iii], label=lbl)
        ax.plot(vlon, vlat, "^y", markersize=10)
        handles, labels = ax.get_legend_handles_labels()
        if legend: ax.legend(handles, labels)

    def get_vaaroot(self):
        vaa = self.root[0][0]
        return vaa

    def __str__(self):
        spc = "\n"
        rstr = "Issue Time (issueTime): {}".format(self.issueTime)
        rstr += spc
        rstr += "VAAC name: {}".format(self.issuingVolcanicAshAdvisoryCentre)
        rstr += spc
        rstr += "Bulletin ID: {}".format(self.bulletinIdentifier)
        rstr += spc

        for key in self.volcano:
            rstr += "{} : {}".format(key, self.volcano[key])
            rstr += spc
        for key in self.details:
            rstr += "{} : {}".format(key, self.details[key])
            rstr += spc
        rstr += "Observation---------\n"
        for key in ['date','numpoly','poly']:
            if key not in self.obs: continue
            if key=='poly':
               polygons = self.obs[key]
               for poly in polygons:
                   for pkey in poly:
                       rstr += "{} : {}".format(pkey, poly[pkey])
                       rstr += spc
                   rstr += '\n'
            else:        
                rstr += "{} : {}".format(key, self.obs[key])
                rstr += spc
        rstr += "---------------------------------------\n"
        for key in self.forecast:
            print('HERE', key)
            if isinstance(self.forecast[key], dict):
              rstr += "{} -------------\n".format(key)
              f1 = self.forecast[key]
              for key2 in ['date','numpoly','poly']:
                 if key2 not in f1: continue
                 if key2 == 'poly':
                   polygons = f1[key2]
                   for poly in polygons:
                      for pkey in poly:
                         rstr += "{} : {}".format(pkey, poly[pkey])
                         rstr += spc
                      rstr += '\n'
                 else:
                   rstr += "{} : {}".format(key2, f1[key2])
                   rstr += spc
 
              rstr += '\n'
        rstr += "REMARKS: {}".format(self.remarks)

        return rstr

    def basic(self):
        root = self.tree.getroot()
        for val in root.iter():
            print("TAG", val.tag)
            print("ATTR", val.attrib)
            for item in val.attrib.items():
                print("ITEM", item)
            print("TEXT", val.text)
            print("-----------------------")

    def get_bulletin_id(self):
        for child in self.root:
            if "bulletinIdentifier" in child.tag:
                bulletinIdentifier = child.tag
                return child.text
        return None

    @staticmethod
    def parse_name(name):
        """
        split the name field into volcano name and volcano id.
        """
        temp = name.text.split()
        nlen = len(temp)
        # volcano id should be last part of name
        vid = temp[-1]
        vname = temp[0:nlen-1]
        return str.join(' ',vname), vid

    def get_volcano(self, verbose=False):
        """
        get information on volcano
        """
        # strfmt = "%Y-%m-%dT%H:%M:%SZ"
        vhash = {}
        vaa = self.vaaroot
        volcano = findel(findel(vaa, "volcano"), "EruptingVolcano")

        name = findel(volcano, "name")
        vhash['name'], vhash['volcanoID'] = self.parse_name(name)
        #vhash["name"] = name.text.split()[0]
        #vhash["volcanoID"] = name.text.split()[1]

        position = findel(findel(findel(volcano, "position"), "Point"), "pos")
        vhash["position"] = position.text
        vhash["latitude"] = float(position.text.split()[0])
        vhash["longitude"] = float(position.text.split()[1])

        edate = findel(volcano, "eruptionDate")
        try:
            vhash["eruptionDate"] = datetime.datetime.strptime(edate.text, self.strfmt)
        except:
            if verbose:
                print("READING eruption DATE FAILED", edate.text)
            vhash["eruptionDate"] = edate.text
        elevation = findel(vaa, "summitElevation")
        vhash["elevation"] = float(elevation.text)
        vhash["elevation units"] = elevation.attrib
        return vhash

    def get_details(self):
        fhash = {}
        vaa = self.vaaroot
        details = findel(vaa, "eruptionDetails")
        fhash["eruptionDetails"] = details.text
        info = findel(vaa, "informationSource")
        fhash["informationSource"] = info.text
        return fhash

    def get_all_forecasts(self, verbose=False):
        vaa = self.vaaroot
        fhash = {}
        ss1 = "forecast"
        clist = listfindel(vaa, ss1)
        if verbose:
            listel(vaa)
            print("This is the get forecast method")
            print("\n {}".format(ss1))
            for temp in clist:
                listel(temp)
        for iii, forecast in enumerate(clist):
            shash = {}
            if verbose:
                print("------------------------------------")
                print("------------------------------------")
                print("------------------------------------")
                print("WORKING ON FORECAST {}".format(iii))
            shash = self.get_forecast(forecast, verbose)
            fhash["Forecast{}".format(iii)] = shash
        return fhash


    def testdate(self,date):
        """
        currently some sort of bug in iwxxm files 
        so some dates are a month ahead.
        """
        testdate = self.issueTime
        newdate = date
        # if date is more than a week ahead of issue date then clearly not correct.
        # changing month to month before works.
        if date > testdate + datetime.timedelta(hours=24*7):
           newdate = datetime.datetime(date.year,date.month-1,date.day,date.hour,date.minute)
           print('warning: date out of range {}. changed to {}'.format(date,newdate))
        return newdate

    def get_forecast(self, forecast, verbose=False):
        fhash = {}
        fhash['numpoly']=0
        fhash['poly'] = []
        ss2 = "VolcanicAshForecastConditions"
        if verbose:
            print("\n {}".format(ss2))
            tempA = findel(forecast, ss2)
            listel(tempA)

        if not testel(findel(forecast, ss2)):
            return fhash

        ss2a = "phenomenonTime"
        ss3a = "TimeInstant"
        ss4a = "timePosition"

        if testel(findel(findel(forecast, ss2), ss2a)):
            time = findel(findel(findel(findel(forecast, ss2), ss2a), ss3a), ss4a)
            try:
                fdate = datetime.datetime.strptime(time.text, self.strfmt)
             
                fhash["date"] = self.testdate(fdate)
            except:
                print("Forecast data cannot be parsed", time.text)

            if verbose:
                print("\n {}".format(ss2a))
                temp = findel(tempA, ss2a)
                listel(temp)
                print("\n {}".format(ss3a))
                temp = findel(temp, ss3a)
                listel(temp)
                print("\n {}".format(ss4a))
                temp = findel(temp, ss4a)
                listel(temp)

        ss2b = "ashCloud"
        ss3b = "VolcanicAshCloudForecast"
        ss4b = "ashCloudExtent"
        if verbose:
            print("\n {}".format(ss2b))
            temp = findel(tempA, ss2b)
            listel(temp)
            print("\n {}".format(ss3b))
            try:
                temp = findel(temp, ss3b)
                listel(temp)
            except:
                pass
        try:
            f1 = findel(findel(findel(forecast, ss2), ss2b), ss3b)
        except:
            return fhash
        forecastlist = listfindel(findel(forecast, ss2), ss2b)
        fhash['numpoly'] = len(forecastlist)
        for child in forecastlist:
            f2 = findel(child,ss3b)
            fhash2 = self.subfunc(f2, verbose=verbose)
            fhash['poly'].append(fhash2)
        return fhash

    def get_obs(self, verbose=False):
        vaa = self.vaaroot
        fhash = {}
        fhash['numpoly']=0
        # list of polygons
        fhash['poly'] = []
        fhash2 = {}
        fhash3 = {}
        ss1 = "observation"
        ss2 = "VolcanicAshObservedOrEstimatedConditions"

        if verbose:
            print("THIS IS THE get_obs method")
            print("\n {}".format(ss1))
            temp = findel(vaa, ss1, verbose)
            listel(temp)

            print("\n {}".format(ss2))
            temp = findel(temp, ss2, verbose)
            listel(temp)

        # ---------------------------------------------------------------------------
        ss3b = "phenomenonTime"
        ss4b = "TimeInstant"
        ss5b = "timePosition"
        obs = findel(findel(findel(findel(findel(vaa, ss1), ss2), ss3b), ss4b), ss5b)
        odate = datetime.datetime.strptime(obs.text, self.strfmt)
        fhash["date"] = self.testdate(odate)

        # ------------------------------------------------------------------------------
        ss3a = "ashCloud"
        ss4a = "VolcanicAshCloudObservedOrEstimated"

        if verbose:
            print("\n {}".format(ss3a))
            tempa = findel(temp, ss3a)
            listel(tempa)
            if testel(tempa):
                print("\n {}".format(ss4a))
                tempalist = listfindel(tempa, ss4a)


        tempblist = listfindel(findel(findel(vaa, ss1), ss2), ss3a)
        fhash['numpoly'] = len(tempblist) 
        for tempb in tempblist:

            testb = testel(tempb)
            if not testb:
               continue 
            tempa = findel(tempb, ss4a)
            # tempa = findel(findel(findel(findel(vaa,ss1),ss2),ss3a),ss4a)
            # ------------------------------------------------------------------------------
            ss5a = "directionOfMotion"
            motion = findel(findel(findel(findel(findel(vaa, ss1), ss2), ss3a), ss4a), ss5a)
            fhash2[ss5a] = {"value": float(motion.text), "uom": motion.attrib["uom"]}

            # ------------------------------------------------------------------------------
            ss5a = "speedOfMotion"
            motion = findel(findel(findel(findel(findel(vaa, ss1), ss2), ss3a), ss4a), ss5a)
            if verbose:
                print(type(motion))
            fhash2[ss5a] = {"value": float(motion.text), "uom": motion.attrib["uom"]}

            # ------------------------------------------------------------------------------
            fhash3 = self.subfunc(tempa, verbose)
            fhash['poly'].append({**fhash3,**fhash2}) 
            #print('zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz')
            #print(fhash2)
            #print('zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz')
 
        return  fhash

    def subfunc(self, vaa, verbose=False):
        fhash = {}
        ss5a = "ashCloudExtent"
        ss6a = "AirspaceVolume"

        if verbose:
            print("\n {}".format(ss5a))
            tempb = findel(vaa, ss5a)
            listel(tempb)
            print("\n {}".format(ss6a))
            tempb = findel(tempb, ss6a)
            listel(tempb)

        extent = findel(findel(vaa, ss5a), ss6a)
        # extent = findel(findel(findel(findel(findel(findel(vaa, ss1),ss2),ss3a),ss4a),ss5a),ss6a)
        uplimit = findel(extent, "upperLimit")
        uplimitref = findel(extent, "upperLimitReference")
        lowlimit = findel(extent, "lowerLimit")
        lowlimitref = findel(extent, "lowerLimitReference")
        projection = findel(extent, "horizontalProjection")
        if verbose:
            listel(projection)
        fhash["ashCloudExtent"] = {
            "value": float(uplimit.text),
            "uom": uplimit.attrib["uom"],
        }

        # ------------------------------------------------------------------------------

        ss7a = "horizontalProjection"
        ss8a = "Surface"
        ss9a = "patches"
        ss10a = "PolygonPatch"
        ss11a = "exterior"
        ss12a = "LinearRing"
        ss13a = "posList"

        if verbose:
            print("---------------------------------------")
            print("\n {}".format(ss7a))
            tempc = findel(extent, ss7a)
            listel(tempc)
            print("\n {}".format(ss8a))
            tempd = findel(tempc, ss8a)
            listel(tempd)
            print("\n {}".format(ss9a))
            tempe = findel(tempd, ss9a)
            listel(tempe)
            print("\n {}".format(ss10a))
            tempf = findel(tempe, ss10a)
            listel(tempf)
            print("\n {}".format(ss11a))
            tempg = findel(tempf, ss11a)
            listel(tempg)
            print("\n {}".format(ss12a))
            temph = findel(tempg, ss12a)
            listel(temph)
            # print('\n {}'.format(ss13a))
            # tempi = findel(tempg,ss13a)
            # listel(tempi)

        projection = findel(findel(findel(findel(extent, ss7a), ss8a), ss9a), ss10a)
        exterior = findel(projection, ss11a)
        lring = findel(findel(exterior, ss12a), ss13a)
        # if verbose: print('--------------')
        # if verbose: listel(lring)
        fhash["exterior"] = self.process_lring(lring.text.split())

        return fhash

    def process_lring(self, plist):
        """
        exterior is a list of [lat,lon,lat,lon,lat,lon,lat,lon]
        Need to make it into [(lat,lon),(lat,lon),(lat,lon)]
        """
        tlist = []
        for iii, elem in enumerate(plist):
            if iii % 2 == 0:
                tlist.append((float(plist[iii + 1]), float(plist[iii])))
        return tlist

    def get_remarks(self):
        vaa = self.vaaroot
        ss1 = "remarks"
        remarks = findel(vaa, ss1)
        return remarks.text

    def get_issue_time(self):
        minfo = self.root[0]
        vaa = minfo[0]
        time = vaa[0][0][0]
        try:
            itime = datetime.datetime.strptime(time.text, self.strfmt)
        except:
            print("WARNING: cannot parse issue time for {}".format(self.source))
            itime = time.txt
        return itime

    def get_vaac(self):
        minfo = self.root[0]
        vaa = minfo[0]
        vaac = vaa[1][0][0][0][2]
        return vaac.text

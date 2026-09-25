from utilvolc.helperinterface import FileNameInterface
import datetime
from utilvolc import volcano_information
import sys

class VolcatName25(FileNameInterface):
    """
    04/28/2025 works with new data format.
    """

    def __init__(self, fname):
        # if full directory path is input then just get the filename
        self.vhash = {}
        self.vhash["url"]= fname
        if isinstance(fname, str):
            if "/" in fname:
                temp = fname.split("/")
                self.fname = temp[-1]
            else:
                self.fname = fname
        print('FNAME', self.fname)
        self.vhash['filename'] = self.fname
        self.date = None
        self.image_date = None
        self.event_date = None
        self.dtfmt = "%Y%m%d%H%M%S"

        #self.make_keylist()
        #self.make_datekeys()

        # parse only if a string is given.
        if isinstance(fname, str):
            self.parse(self.fname)



    def find_name(self,ref_file):
        vid = self.vhash['event vid']
        vs = volcano_information.VolcanoSearch(ref_file)
        vname = vs.search_vid(vid)
        if vname: vname = vname[0]
        else: vname=vid
        self.vhash['volcano name'] = vname
           
    def make_filename(self):
        pass

    def parse(self,fname):
        """
        """
        temp = fname.split("_")
        for ttt in temp:
            if 'vid' in ttt:
                self.vhash['event vid'] = ttt.split('vid')[-1]
            elif ttt.startswith('s'):
                self.vhash['sdate'] = datetime.datetime.strptime(ttt[1:-1], self.dtfmt)
            elif ttt.startswith('e'):
                self.vhash['edate'] = datetime.datetime.strptime(ttt[1:-1], self.dtfmt)
            elif ttt.startswith('c'): 
                ttt = ttt.replace('.nc','')
                self.vhash['cdate'] = datetime.datetime.strptime(ttt[1:-1], self.dtfmt)
            elif ttt.startswith('v'): 
                self.vhash['version'] = ttt  
            elif ttt.startswith('g'): 
                self.vhash['g'] = ttt  

class VolcatName(FileNameInterface):
    """
    12/18/2020 works with 'new' data format.
    parse the volcat name to get information.
    attributes:
    self.fname name of file
    self.date date associated with file
    self.vhash is a dictionary which contains info
    gleaned from the naming convention.

    methods:
    compare: returns what is different between two file names.
    """

    def __init__(self, fname, original_name=None):
        # if full directory path is input then just get the filename
        self.fname = fname
        if isinstance(fname, str):
            if "/" in fname:
                temp = fname.split("/")
                self.fname = temp[-1]
        self.vhash = {}
        self.date = None
        self.image_date = None
        self.event_date = None
        self.image_dtfmt = "s%Y%j_%H%M%S"
        self.event_dtfmt = "b%Y%j_%H%M%S"

        self.make_keylist()
        self.make_datekeys()

        self.pc_corrected = False
        # parse only if a string is given.
        if isinstance(fname, str):
            self.parse(self.fname)
        if isinstance(original_name,str):
           self.vhash['filename']=original_name
        else:
           self.vhash["filename"] = fname

    def make_datekeys(self):
        self.datekeys = [3, 4, 10, 11]

    def make_keylist(self):
        self.keylist = ["algorithm name"]
        self.keylist.append("satellite platform")
        self.keylist.append("event scanning strategy")
        self.keylist.append("observation_date")  # should be image date (check)
        self.keylist.append("image time")
        self.keylist.append("feature_id")
        self.keylist.append("event vid")
        self.keylist.append("description")
        self.keylist.append("WMO satellite id")
        self.keylist.append("image scanning strategy")
        self.keylist.append("event_date")  # should be event date (check)
        self.keylist.append("event_time")
        self.keylist.append("original_feature_id")

    def __lt__(self, other):
        """
        sort by
        volcano id first.
        event date
        image date
        feature id if it exists.
        """
        if self.vhash["event vid"] < other.vhash["event vid"]:
            return True
        if "fid" in self.vhash.keys() and "fid" in other.vhash.keys():
            if self.vhash["fid"] < other.vhash["fid"]:
                return True
        if self.event_date < other.event_date:
            return True
        if self.image_date < other.image_date:
            return True
        sortlist = [
            "feature id",
            "image scanning strategy",
            "WMO satellite id",
            "description",
            "event scanning strategy",
            "satellite platform",
            "algorithm name",
        ]
        for key in sortlist:
            if key in other.vhash.keys() and key in self.vhash.keys():
                if self.vhash[key] < other.vhash[key]:
                    return True

    def compare(self, other):
        """
        other is another VolcatName object.
        Returns
        dictionary of information which is different.
        values is a  tuple of (other value, self value).
        """
        diffhash = {}
        for key in self.keylist:
            if key in other.vhash.keys() and key in self.vhash.keys():
                if other.vhash[key] != self.vhash[key]:
                    diffhash[key] = (other.vhash[key], self.vhash[key])
        return diffhash

    def __str__(self):
        # 2023 14 Jan (amc) make sure keys are in the dictionary.
        keys = self.vhash.keys()
        keylist = [x for x in self.keylist if x in keys]
        val = [str(self.vhash[x]) for x in keylist]
        return str.join("_", val)

    @staticmethod
    def split_name(fname):
        # if full_disk in filename replace with fulldisk because _ is used as separator
        fname = fname.replace("Full_Disk", "FullDisk")
        fname = fname.replace("FULL_DISK", "FullDisk")
        temp = fname.split("_")
        return temp

    def parse(self, fname):
        temp = self.split_name(fname)

        if "pc" in temp[-1]:
            self.pc_corrected = True
        jjj = 0
        for iii, key in enumerate(self.keylist):
            val = temp[jjj]
            # nishinoshima files have a g00? code before the volcano id.
            if key == "fid":
                if val[0] == "g":
                    self.vhash[key] = val
                else:
                    continue
            self.vhash[key] = val
            jjj += 1

        # Image date marks date of the data collection
        dk = self.datekeys
        if isinstance(dk[0], int) and isinstance(dk[1], int):
            dstr = "{}_{}".format(
                self.vhash[self.keylist[dk[0]]], self.vhash[self.keylist[dk[1]]]
            )
            self.image_date = datetime.datetime.strptime(dstr, self.image_dtfmt)

        # Event date is start of event
        if isinstance(dk[2], int) and isinstance(dk[3], int):
            dstr = "{}_{}".format(
                self.vhash[self.keylist[dk[2]]], self.vhash[self.keylist[dk[3]]]
            )
            self.event_date = datetime.datetime.strptime(dstr, self.event_dtfmt)
            self.vhash[self.keylist[dk[3]]] = self.vhash[self.keylist[dk[3]]].replace(
                ".nc", ""
            )

        # this is the date associated with the data
        self.vhash["observation_date"] = self.image_date
        # this date may be the same as the image date or earlier
        self.vhash["event date"] = self.event_date
        self.date = self.image_date
        return self.vhash

    @property
    def image_date_str(self):
        return self.image_date.strftime(self.image_dtfmt)

    def make_filename(self):
        """
        To do: returns filename given some inputs.
        """
        return -1


class VolcatNameA(VolcatName):
    # for the Bezymianny data and some older data the first feature id is not there.

    def make_datekeys(self):
        self.datekeys = [3, 4, 9, 10]

    def make_keylist(self):
        self.keylist = ["algorithm name"]
        self.keylist.append("satellite platform")
        self.keylist.append("event scanning strategy")
        self.keylist.append("observation_date")  # should be image date (check)
        self.keylist.append("image time")
        self.keylist.append("event vid")
        self.keylist.append("description")
        self.keylist.append("WMO satellite id")
        self.keylist.append("image scanning strategy")
        self.keylist.append("event_date")  # should be event date (check)
        self.keylist.append("event_time")
        self.keylist.append("feature_id")

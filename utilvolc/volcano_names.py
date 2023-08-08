import os
import pandas as pd


def fix_volc_name(volcname):
    """Fixes the volcano name if a comma, or space appear in the name"""
    if "," in volcname:
        s = volcname.find(",")
        tmp = volcname[:s]
        tmp2 = volcname[s + 2 :]
        volcname = tmp2 + "_" + tmp
    if " " in volcname:
        volcname = volcname.replace(" ", "_")
    return volcname

#def get_latlon(data_dir):
#    """Read csv file containing volcano name, latitude, longitude
#    Inputs:
#    data_dir: directory where Volcanoes.csv is located
#    Outputs:
#    volcdf: pandas dataframe of volcanos, latitude, longitude
#    """
#    import pandas as pd
#    volcdf = pd.read_csv(data_dir + "Volcanoes.csv", sep=",")
#    return volcdf

class VolcanoName:

    def __init__(self,name,vfilename='None'):
        self._name = name
        self.vlist = vfilename

    @property
    def name(self):
        return self._name

    @setter
    def name(self,name):
        rval = None
        if self.vlist:
           rval = self.vlist.find(name)
           if isinstance(rval,list):
              rval = rval[0]
        if not rval:
           rval = name
        return rval   

    @property
    def vfile(self):
        return self.vfile

    @setter
    def vlist(self,vfilename):
        if os.path.isfile(vfilename):
           return VolcList(vfile) 
        else:
           return None

class VolcList:
    def __init__(self,vfile):
        import pandas as pd
        self.df = pd.read_csv(vfile,
                  header=0,
                  index_col=False,
                  names=['volcano_name',
                         'volcano_region','volcano_lat',
                         'volcano_lon','volcano_elevation','type'])
        self.empty = pd.DataFrame()
        self.vlist = self.df.volcano_name.unique()

    def __bool__(self):
        if self.df.empty: return False
        else: return True

    def find_by(self,match_hash):
        """
        match_hash : dictionary
        example {'volcano_lat':45.0, 'volcano_lon':-175.2}
        """
        match = self.df.copy()
        for key in match_hash:
            match = match[match[key] = match_hash[key]]
        return match

    def find_exact(self,vname):
        vname = fix_volc_name(vname)
        match = [x for x in self.vlist if vname.lower() == x.lower()]
        if match:
            return self.df[self.df.volcano_name==match[0]]
        else:
            return self.empty

    def find_close(self,vname):
        vname = fix_volc_name(vname)
        possible = [x for x in self.vlist if vname.lower() in x.lower()]
        return self.df[self.df['volcano_name'].isin(possible)]

    def find_start(self,vname):
        vname = fix_volc_name(vname)
        lnn = len(vname)
        possible = [x for x in self.vlist if vname.lower() in x[0:lnn].lower()]
        return self.df[self.df['volcano_name'].isin(possible)]


    def find(self,vname):
        vname = fix_volc_name(vname)
        guess = self.find_exact(vname)
        if  guess.empty:
           guess = self.find_start(vname)
        if guess.empty:
           guess = self.find_close(vname)
        return guess    

    def get_record(self,vname):
        vname = fix_volc_name(vname)
        df = self.find(vname)

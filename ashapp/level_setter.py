
def get_qva_levels():
    ls = LevelSetter(bottom=50, top=650, dz=50, unit='FL')
    return ls.levlist


def get_levelsetter(inp):

    if inp['qvaflag']:
       levelsetter = LevelSetter(unit='FL')
    else:
        dz=None
        top=None
        bottom=None
        unit='m'
        if 'dzconcgrid' in inp.keys():
            dz = inp['dzconcgrid']
        if 'top' in inp.keys():
            top = inp['top']
        if 'bottom' in inp.keys():
            bottom = inp['bottom']
        levelsetter = LevelSetter(bottom,top,dz,unit)
    return levelsetter


def FL2meters(flight_level):
    meters = flight_level * 100 / 3.28084
    return int(meters)
    # return int(np.ceil(flight_level / 10.0) * 10)

class LevelSetter:

    def __init__(self, 
                 bottom=None, 
                 top=None, 
                 dz=None,
               unit='FL'):
        self._levlist = []
        self._descriptions = []
        self._unit = unit

        if unit=='FL':
           if not isinstance(bottom,(float,int)):
              bottom=50
           if not isinstance(top,(float,int)):
              top = 650 
           if not isinstance(dz,(float,int)):
              dz = 50 
           print('HERE HERE', bottom,top,dz)
           self._levlist, self._descriptions = set_qva_levels(bottom,top,dz)

        if unit=='m':
           if not isinstance(bottom,(float,int)):
              bottom=1000
           if not isinstance(top,(float,int)):
              top =  30000
           if not isinstance(dz,(float,int)):
              dz =  1000
           self._levlist, self._descriptions = set_meter_levels(bottom,top,dz)

    @property
    def unit(self):
        return self._unit
    
    @property
    def levlist(self):
        return self._levlist
 
    @property
    def descriptions(self):
        return self._descriptions


def set_qva_levels(bottom=50,top=650,dz=50):
    # every 5000 ft (FL50 chunks)
    # Approx 1.5 km.
    levlist_fl = list(range(bottom, top,dz))
    levlist = [FL2meters(x) for x in levlist_fl]
    rlist = []
    plev = "SFC"
    for lev in levlist_fl:
        nlev = "FL{}".format(lev)
        rlist.append("{} to {}".format(plev, nlev))
        plev = nlev
    return levlist, rlist

def set_meter_levels(bottom=1000, top=30000, dz=1000):
    levlist = list(range(int(bottom), int(top)+int(dz), int(dz)))
    plev = "0m"
    rlist = []
    for lev in levlist:
        nlev = "{}m".format(lev)
        rlist.append("{} to {}".format(plev, nlev))
        plev = nlev
    return levlist, rlist
    

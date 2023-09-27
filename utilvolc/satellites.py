

class SatelliteLocations:

    def __init__(self):
        shash = {} 
        shash['himawari'] = 140.7
        shash['goes17'] = -137.2
        self.shash = shash 



def viewangle(satx,vloc):
    """ 
    Estimate viewing angle for geostationary satellite.

    satx : longitude of satellite
    vloc : location of interest, such as volcano vent tuple (latitude, longitude)
    """
    yyy = np.abs(vloc[0])
    vxx = vloc[1]
    x1 = np.abs(satx-vxx)
    x2 = 360-x1
    xxx = np.min([x1,x2])
    hhh = (xxx*xxx+yyy*yyy)**0.5
    return hhh
   



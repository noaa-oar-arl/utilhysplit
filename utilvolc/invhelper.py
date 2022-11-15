# -----------------------------------------------------------------------------
import datetime
import numpy as np


def get_suffix(suffix_type, dtfmt, nnn, ndate, bottom, center=None):
    if suffix_type == "int":
        suffix = "{:03d}".format(nnn)
    elif suffix_type == "date":
        str1 = ndate.strftime(dtfmt)
        str2 = str(int(bottom))
        suffix = "{}_{}".format(str1, str2)
    else:
        suffix = "{:03d}".format(nnn)
    if center:
        latstr = "{:.3f}".format(center[0]).replace(".", "p")
        lonstr = "{:.3f}".format(center[0]).replace(".", "p")
        suffix += "_" + latstr + "_" + lonstr
    return suffix


def inverse_get_center_list(inp):

    # center point (lat, lon)
    center = (inp["latitude"], inp["longitude"])
    # area to cover. (meters squared)
    area = inp["area"]
    # number of points to each side of center.
    num = inp["hnum"]
    if num == 0:
        return center
    # here should use a square area.
    # centerlist = []
    dlat = (area ** 0.5) / 111.0e3
    dlon = dlat * np.cos(center[0] * np.pi / 180.0)
    lat0 = center[0] - num * dlat
    latm = center[0] + (num + 0.75) * dlat
    latlist = np.arange(lat0, latm, dlat)
    lon0 = center[1] - num * dlon
    lonm = center[1] + (num + 0.75) * dlon
    lonlist = np.arange(lon0, lonm, dlon)
    latlon = np.meshgrid(latlist, lonlist)
    latlon = zip(latlon[0].flatten(), latlon[1].flatten())
    return list(latlon)


def inverse_get_suffix_list(inp, suffix_type="date", dtfmt="%m%d%H"):
    """
    inp: dictionay generated in
    suffix_type : str
          'int' suffix will be an integer 001, 002, 003....
          'date' suffix will be of form date_height. with date
                 determined by dtfmt.
    dtfmt : str. determines format of date when suffix_type is 'date'

    outputs
    ---------
    suffixhash : dictionary
    key is the suffix. value is a dictionary with information about the run.
    including start date, bottom and top elevation.
    """
    suffixhash = {}
    vres = inp["inv_vertical_resolution"]
    dt = datetime.timedelta(hours=inp["timeres"])
    if dt.seconds % 3600 > 1e-8:
        dtfmt = "%m%d%Hh%M"
    sdate = inp["start_date"]
    # edate = inp['start_date'] + datetime.timedelta(hours=inp["durationOfSimulation"])
    edate = inp["start_date"] + datetime.timedelta(hours=inp["emissionHours"])
    ndate = sdate
    done = False
    nnn = 0
    if "hnum" in inp.keys():
        latlonlist = inverse_get_center_list(inp)
    else:
        latlonlist = [1]
    # latlon loop.
    for center in latlonlist:
        iii = 0
        # time loop
        if len(latlonlist) > 1:
            center_suffix = center
        else:
            center_suffix = None

        while not done:
            # vertical resolution loop.
            vdone = False
            jjj = 0
            bottom = inp["bottom"]
            while not vdone:
                inhash = {}
                inhash["sdate"] = ndate
                inhash["edate"] = ndate + dt
                inhash["bottom"] = bottom
                inhash["top"] = bottom + vres
                suffix = get_suffix(
                    suffix_type, dtfmt, nnn, ndate, bottom, center_suffix
                )
                suffixhash[suffix] = inhash.copy()
                bottom += vres
                if bottom > inp["top"]:
                    vdone = True
                jjj += 1
                nnn += 1
                # limit of no more than 50 heights.
                if jjj > 500:
                    return suffixhash
            iii += 1
            ndate = ndate + dt
            if ndate >= edate:
                done = True
    return suffixhash



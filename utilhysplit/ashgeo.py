#!/n-home/alicec/anaconda/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
#from math import *

from math import ceil, floor

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
# from netCDF4 import Dataset
#from pylab import *
import numpy as np
import numpy.ma as ma
import shapely.geometry as sgeo
#from matplotlib.path import Path
from matplotlib.collections import LineCollection
#from scipy.io import netcdf
from scipy.spatial import Delaunay
from shapely.ops import cascaded_union, polygonize, transform

# import map_projections
# from plume_models import usgstable
# from meteorology import *
# from mytools import emap
# import pandas as pd

##routines and classes which have to do with computational geometry.

##Ash could object could get info from hysplit_nc or satellite_nc or other.
##Use computational geometry for comparison rather than re-gridding and doing point-by-point comparison.


def contour2poly(qcs, plot=0, verbose=0):
    """Input a contour set object which is output by the pyplot.contour routine
    Returns a list of poly_ash objects
    Each polygon returned shows the area which has values above the contour value.
    """
    path_cntr = 0
    npoly = len(qcs.collections)
    holes = [()]
    #holelist = []
    perimeter_ilist = []
    poly_list = []
    perimeter_ra = []
    levra = []
    pathlist = []
    clr_list = ["cyan", "blue", "purple", "green", "yellow", "orange", "red", "magenta"]
    if verbose:
        print("Contour levels", qcs.levels)
        print("Number of Contour levels", npoly)
    for nnn in range(npoly):
        lev = qcs.levels[nnn]
        n_paths = len(qcs.collections[nnn].get_paths())
        if verbose:
            print("Number of paths at contour level ", lev, ":", n_paths)
        for iii in range(n_paths):
            p = qcs.collections[nnn].get_paths()[iii]
            # v = p.to_polygons()  ##some documentation claimed should use to_polygons rather than access vertices attribute directly.
            v = p.vertices
            zzz = 0
            if plot == 1:
                x = v[:, 0]
                y = v[:, 1]
                plt.plot(x, y, clr_list[zzz])
                zzz += 1
                if zzz > len(clr_list):
                    zzz = 0

            perimeter_ra.append(v)
            pathlist.append(p)
            levra.append(lev)
            holes.append([])
            path_cntr += 1
    holes[0] = []
    perimeter_ilist = list(
        range(0, path_cntr)
    )  ###Starts out assuming every path is an outer perimeter.
    levra = np.array(levra)  ##Change list to array.

    ####Each path is either an outside perimeter or a hole - an inside perimeter.
    ####This set of loops finds which contour levels are holes -
    ####These are paths contained within another path which have the same contour value.
    ####Once a path is determined to be a hole it is removed from the perimeter_ilist.
    i = 0
    for i in range(0, len(pathlist) - 1):
        for j in range(i + 1, len(pathlist)):
            if pathlist[i].contains_path(pathlist[j]) and levra[i] == levra[j]:
                holes[i].append(j)
                try:
                    perimeter_ilist.remove(j)
                except:
                    pass
            elif pathlist[j].contains_path(pathlist[i]) and levra[i] == levra[j]:
                holes[j].append(i)
                try:
                    perimeter_ilist.remove(i)
                except:
                    pass

    ####OBSOLETE This set of loops finds the actual holes####
    # for i in xrange(0,len(holes)-1):
    #     try:
    #        holes[i].remove(())
    #     except:
    #        pass
    #     for j in holes[i]:
    #       for k in holes[j]:
    #            try:
    #               holes[i].remove(k)
    #            except:
    #               pass

    ####This set of loops creates the polygon objects####
    for i in perimeter_ilist:
        hole_perimeter = []
        for j in holes[i]:
            hole_perimeter.append(perimeter_ra[j])
        poly = AshPoly(perimeter_ra[i], lvalue=levra[i], holes=hole_perimeter)
        poly_list.append(poly)
        ##Right now perimeter_ra is a numpy array and the shapely polygon needs a sequence of tuples as input
        # polyobj = sgeo.polygon([perimeter])
    if verbose:
        print("Holes", holes)
        print("Perimeter", perimeter_ilist)
    return poly_list


def suggest_contours(data, log=1, MISSING=-999, minpercent=5, binwidth=1, verbose=0):
    """uses a histogram to suggest a set of contour levels for data.
    log=1 will space levels by orders of magnitude.
    minpercent (default 5) key word sets what percent of data a bin must contain to be output as level.
    binwidth is the starting contour interval."""
    returnval_i = []
    vp = np.where(data != MISSING)
    if log == 1:
        logdata = np.log10(data[vp])
    else:
        logdata = data[vp]
    mindata = floor(np.min(logdata))
    maxdata = ceil(np.max(logdata))
    nbins = (maxdata - mindata) / binwidth
    hist, bin_edges = np.histogram(logdata, bins=nbins, range=(mindata, maxdata))
    # hist, bin_edges = np.histogram(logdata, bins=nbins  )
    bin_edges = np.array(bin_edges)
    hist = np.array(hist)
    total = np.sum(hist)
    hist = hist.astype(float) / total * 100.0
    vp2 = np.where(hist >= minpercent)
    for v in vp2[0]:
        returnval_i.append(v)
        last = v
    returnval_i.append(last + 1)
    returnval = 10 ** bin_edges[returnval_i]
    if verbose == 1:
        print("Data ranges from ", mindata, " to ", maxdata)
    return returnval


def contour_plot(x, y, z, clevs, missing=(), filled=0):
    # if missing:
    #   ma.masked_where(z == missing, z)
    clr_list = ["cyan", "blue", "purple", "green", "yellow", "orange", "red", "magenta"]
    clr_list = [
        "#99FF33",
        "yellow",
        "orange",
        "red",
        "yellow",
        "orange",
        "red",
        "magenta",
    ]
    n = len(clevs)
    clr_use = clr_list[:n]
    if filled == 0:
        cs = plt.contour(x, y, z, levels=clevs, colors=clr_use)
    else:
        cs = plt.contourf(x, y, z, levels=clevs, colors=clr_use)
    cs.cmap.set_under("grey")
    cs.cmap.set_over("purple")
    plt.colorbar(cs)
    # plt.clabel(cs,inline=1, fontsize=10)
    # plt.show()
    return cs


def cropra(dra, missing, maskthera=0, crop=1, datara=None, buffer=0):
    """Crops a 2d array to  rectangular area which contains all valid values."""
    """Set maskthearea=1 to return a masked array"""
    """if crop=0 then will mask but not crop array"""
    """if a datara is input then the datara will be cropped according to the valid area in xra"""
    if crop == 1:
        vp = np.where(dra != missing)  # Indices of valid points
        #szx = dra.shape[0]
        #szy = dra.shape[1]
        ##if the min index < 0 Python seems to set it to zero
        ##if the max index > the size of the row/column then Python seems to set it to max.
        if np.any(vp[0]):
            if np.min(vp[0]) - buffer >= 0:
                min_xi = np.min(vp[0]) - buffer
            else:
                min_xi = 0
            max_xi = np.max(vp[0]) + buffer
            if np.min(vp[1]) - buffer >= 0:
                min_yi = np.min(vp[1]) - buffer
            else:
                min_yi = 0
            max_yi = np.max(vp[1]) + buffer

            if np.any(datara):
                print("CROPPING datara", min_xi, max_xi, min_yi, max_yi)
                cropdra = datara[min_xi : max_xi + 1, min_yi : max_yi + 1]
                print(cropdra)
            else:
                print(datara)
                print("CROPPING ra", min_xi, max_xi, min_yi, max_yi)
                cropdra = dra[min_xi : max_xi + 1, min_yi : max_yi + 1]
                print(cropdra)
            if maskthera == 1:
                returnra = ma.masked_where(cropdra == missing, cropdra)
            else:
                print("input Data is cropped.")
                returnra = cropdra
        else:
            if not datara:
                returnra = dra
                print("Data could not be cropped.")
            else:
                returnra = datara
                print("Data could  not be cropped. No valid data in array.")
    # else:
    #        returnra = ma.masked_where(dra== missing,dra)
    return returnra


def histoplot(
    xval,
    bins=100,
    missing=(),
    pdf=1,
    cdf=0,
    clr="b",
    transparency=1,
    label="",
    takelog=0,
    xrange=[],
):
    xstats = {}
    """creates 1d histogram. Returns dictionary with mean and median values."""
    align = "mid"  # other values could be 'left' or 'right'
    histtype = "bar"  # other values could be 'barstacked', 'step' or 'stepfilled'
    if missing != ():
        print("MISSING", missing)
        xvalid = np.where(xval != missing)
        xval = xval[xvalid]
    if takelog == 1:
        xval = np.log10(xval)
    if xrange == []:
        plt.hist(
            xval,
            bins=bins,
            normed=pdf,
            cumulative=cdf,
            color=clr,
            alpha=transparency,
            label=label,
        )
    else:
        plt.hist(
            xval,
            bins=bins,
            normed=pdf,
            cumulative=cdf,
            color=clr,
            alpha=transparency,
            label=label,
            range=xrange,
        )
    # plt.show()
    xstats["mean"] = np.mean(xval)
    xstats["median"] = np.median(xval)
    return xstats


#####CLASSES##########################################################################################


class cloud_poly(object):
    """A collection of ash_poly objects which constitutes a cloud.
    A list of ash_poly objects.
    A dictionary key=lvalue of ash_poly object. Value is list of indices.
    A dictionary key=altitude of ash_poly object. Value is list of indices."""

    def __init__(self, value_units="", altitude_units=""):
        self.ash_poly = []
        self.values = {}
        self.value_units = value_units
        self.altitudes = {}
        self.altitude_units = altitude_units
        self.num_polys = 0
        self.area = 0

    def add_ash(self, apoly):
        num = len(apoly)
        self.ash_poly += apoly
        apoly_i = np.array(list(range(num))) + self.num_polys
        for i in apoly_i:
            value = str(self.ash_poly[i].lvalue)
            print("DEBUG", i, value)
            if value in self.values:
                self.values[value].append(i)
            else:
                self.values[value] = i
            value = str(self.ash_poly[i].altitude)
            if value != "":
                if value in self.altitudes:
                    self.values[value].append(i)
                else:
                    self.values[value] = i
        self.num_polys += num

    def __str__(self):
        values = list(self.values.keys())
        alts = list(self.altitudes.keys())
        values = ", ".join(values)
        alts = ", ".join(alts)
        outstr = (
            "Cloud contains polygons with values of "
            + values
            + " "
            + self.value_units
            + "\n"
        )
        outstr += (
            "Cloud contains polygons with altitudes of "
            + alts
            + " "
            + self.altitude_units
            + "\n"
        )
        outstr += "Cloud contains " + str(self.num_polys) + " polygons \n"
        return outstr

    def get_ash(self, lvalue="", altitude=""):
        """Returns a list of ash polygons with the specified value and altitude"""
        if lvalue != "":
            lvalue = str(lvalue)
            if lvalue in self.values:
                li = self.values[lvalue]
            else:
                print(
                    "Warning: get_ash , lvalue not valid. ",
                    lvalue,
                    " not in ",
                    list(self.values.keys()),
                )
        if alitutude != "":
            altitude = str(altitude)
            if alititude in self.altitudes:
                ai = self.values[altitude]
            else:
                print(
                    "Warning: get_ash , altitude not valid. ",
                    altitude,
                    " not in ",
                    list(self.altitudes.keys()),
                )

        i = set(li) & set(ai)
        return self.ash_poly[i]


#########################################################################################################################
class AshPoly(object):
    """A polygon which depicts a region of ash"""

    """Can be generated from the countour2poly function"""

    def __init__(
        self, perimeter, lvalue=(), hvalue=(), holes=[],  altitude=0, altunit=""
    ):
        self.holes = (
            []
        )  ##List containing sets of vertices which describe any holes in polygon.
        self.hvalue = hvalue  ##Upper limit of contoured area.
        self.lvalue = lvalue  ##lower limit of contoured area. If no upper limit then polygon represents area with value above lower limit.
        self.altitude = altitude  ##
        self.nholes = len(holes)
        self.altunit = altunit
        # self.poly = self.shapely_obj ##Shapely module polygon object.
        ##Make sure that perimeter is stored as a list of tuples.
        temp = list(zip(*perimeter))
        self.perimeter = list(
            zip(temp[0], temp[1])
        )  ##set of vertices which describes outer perimeter of polygon.
        ##Make sure that holes are stored as a list  of list of tuples.
        for h in holes:
            temp = list(zip(*h))
            self.holes.append(list(zip(temp[0], temp[1])))

        self.poly = sgeo.Polygon(self.perimeter, self.holes)
        # self.area = self.get_area()

    def __str__(self):
        outstr = "Number of holes : " + str(self.nholes) + "\n"
        outstr += "Area : " + str(self.area) + "\n"
        outstr += "Perimeter: " + str(self.perimeter) + "\n"
        return outstr

    def kml_coords(self):
        """outputs strings which can be used in a kml file to define a polygon"""
        """outputs string for outer perimeter"""
        """outputs list of strings for inner perimeters"""
        outer_str = ""
        inner_str_list = []
        if self.altitude != ():
            alt_str = str(self.altitude)
        else:
            alt_str = "0"
        for pt in self.perimeter:
            outer_str += str(pt[0]) + "," + str(pt[1]) + "," + alt_str + "\n"
        if self.nholes != 0:
            for ip in range(0, self.nholes):
                in_str = ""
                for pt in self.holes[ip]:
                    in_str += str(pt[0]) + "," + str(pt[1]) + "," + alt_str + "\n"
                inner_str_list.append(in_str)
        return outer_str, inner_str_list

    # def proj(self, proj='sinusoidal'):
    #     """returns a polygon with coordinates in projected coordinate system"""
    #     if proj == 'sinusoidal':
    #        #exterior =  map_projections.sinusoidal_proj(self.perimeter, alt=self.altitude)
    #        interiors = []
    #        for h in self.holes:
    #              interiors.append(map_projections.sinusoidal_proj(h, alt=self.altitude))
    #     poly_projected = sgeo.Polygon(exterior, interiors)
    #     return poly_projected

    # def get_area(self):
    #    """Transforms from lat-lon coords to x-y coords using a sinusoidal projection.
    #       Then uses shapely package to compute the area."""
    #    poly = self.proj(proj='sinusoidal')
    #    return poly.area


#  def get_intersection(self, other):
#      """Returns a polygon which is the intersection of self and another polygon"""
#      return self.poly.intersection(other)
#      return self.poly.difference(other)


# -----------------------------------------------------------------------------------------------------------------#


def transform2(lat, lon):
    from pyproj import Proj

    pc = Proj("+proj=longlat +datum=WGS84 +no_defs")
    x, y = pc(lon, lat)
    return [y, x, "epsg"]


def transformra(lat, lon):
    import utm

    latc = np.zeros_like(lat)
    lonc = np.zeros_like(lat)
    iii = 0
    jjj = 0
    zone = []
    for iii in np.arange(0, lat.shape[0]):
        for jjj in np.arange(0, lat.shape[1]):
            vlat = lat[iii, jjj]
            vlon = lon[iii, jjj]
            # trans = utm.from_latlon(vlat,vlon)
            trans = transform2(vlat, vlon)
            # if trans[1] == 0:
            #   print(iii, jjj, vlat, vlon, trans)
            #   import sys
            #   sys.exit(0)
            latc[iii, jjj] = trans[0]
            lonc[iii, jjj] = trans[1]
            zone.append(trans[2])
    return latc, lonc


def transformpt(lat, lon):
    import utm

    trans = []
    for val in zip(lat, lon):
        # trans.append(utm.from_latlon(val[0],val[1]))
        trans.append(transform2(val[0], val[1]))
    aaa = list(zip(*trans))
    latcoor = aaa[0]
    loncoor = aaa[1]
    print("Zone Numbers", list(set(aaa[2])))
    return latcoor, loncoor


def plotpoly(sgeo_poly):
    """xy plot of a shapely polygon"""
    x, y = sgeo_poly.exterior.xy
    # plt.plot(x,y)
    return x, y


# -------------------------------------------------------------------------------------
"""
Find concave hull of set of points.
Adapted from
gist.github.com/dwyerk/10561690
"""


def add_edge(edges, edge_points, coords, iii, jjj):
    """
    helper function for concave_hull
    """
    if (iii, jjj) in edges or (jjj, iii) in edges:
        return edges, edge_points
    else:
        edges.add((iii, jjj))
        edge_points.append(coords[[iii, jjj]])
    return edges, edge_points


def plot_delauney(ax, edge_points, mpts, hull=None):
    """
    plotting function for concave_hull output and inputs
    """
    lines = LineCollection(edge_points)
    # fig = plt.figure(1)
    # ax = fig.add_subplot(1,1,1)
    # plt.gca().add_collection(lines)

    # plots the delauney triangles
    ax.add_collection(lines)
    # plots the delauney triangles
    dpts = np.array([point.coords[0] for point in mpts])
    ax.plot(dpts[:, 0], dpts[:, 1], "k.", MarkerSize=1)
    # plots points on convex hull.
    if hull:
        x, y = plotpoly(hull)
        ax.plot(x, y, "r.")
    return dpts


def concave_hull(mpoints, alpha=1):
    """
    mpoints is a multi-point object.
    Returns
    rpoly : shapely.geometry.polygon.Polygon
    ep : list of numpy arrays.

    Note.  ep and mpoints can be fed into plot_delauney to plot the Delauney triangles.
           rpoly can be fed in as hull=rpoly.
    """
    coords = np.array([x.coords[0] for x in mpoints])
    tri = Delaunay(coords)
    edges = set()
    edge_points = []
    for ia, ib, ic in tri.vertices:
        pa = coords[ia]
        pb = coords[ib]
        pc = coords[ic]
        aaa = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        bbb = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        ccc = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)

        sss = (aaa + bbb + ccc) / 2.0
        area = np.sqrt(sss * (sss - aaa) * (sss - bbb) * (sss - ccc))
        circum_r = aaa * bbb * ccc / (4.0 * area)
        if circum_r < 1.0 / alpha:
            edges, edge_points = add_edge(edges, edge_points, coords, ia, ib)
            edges, edge_points = add_edge(edges, edge_points, coords, ib, ic)
            edges, edge_points = add_edge(edges, edge_points, coords, ic, ia)

    mmm = sgeo.MultiLineString(edge_points)
    triangles = list(polygonize(mmm))
    rpoly = cascaded_union(triangles)
    return rpoly, edge_points


def make_multi(x, y):
    """
    makes a multi-point object from list of x,y coordinates
    """
    latlon = list(zip(x, y))
    mpts = sgeo.MultiPoint(latlon)
    # poly = mpts.convex_hull
    # coords = list(poly.exterior.coords)
    return mpts


def example1():
    # creates polygon around pardump points.
    from monetio.models import pardump

    def get_thickness(df, t1, t2):
        df2 = df[df["ht"] > t1]
        df2 = df2[df2["ht"] < t2]
        return df2

    fname = "/pub/Scratch/alicec/KASATOCHI/cylindrical/e3/"
    fn = "PARDUMP.cyl.gdas1.e3"
    pd = pardump.Pardump(fname=fname + fn)
    pdict = pd.read()
    p1 = pdict["200808081400"]
    df = get_thickness(p1, 10000, 12000)
    mpts = make_multi(df["lon"], df["lat"])
    ch, ep = concave_hull(mpts, alpha=1)
    fig = plt.figure(1)
    ax = fig.add_subplot(1, 1, 1)
    plot_delauney(ax, mpts, ep, hull=ch)


class RevPar:
    def __init__(self):
        fdir = "/pub/ECMWF/JPSS/reventador/ens"
        fname = "PARDUMP.A2"
        self.tname = os.path.join(fdir, fname)

    def read_pardump(self):
        pd = pardump.Pardump(fname=self.tname)
        self.pdict = pd.read()


# proj = ccrs.Mercator(central_longitude=-175, min_latitude=40, max_latitude=70, latitude_true_scale=60)
# data_crs=ccrs.PlateCarree()
# proj.transform_point(-180, 60, data_crs)


class AshTransform:
    def __init__(self, central_longitude, central_latitude):
        self.central_longitude = central_longitude
        self.central_latitude = central_latitude

    def get_transform(self, hp1):
        hp2 = transform(self.id_func, hp1)
        return hp2

    def id_func(self, x, y, z=None):
        data_crs = ccrs.PlateCarree()
        proj = ccrs.AzimuthalEquidistant(
            central_longitude=self.central_longitude,
            central_latitude=self.central_latitude,
        )
        x, y = proj.transform_point(x, y, data_crs)
        return tuple(filter(None, [x, y, z]))


def distance(pp1, pp2):
    """
    p1 : shapely Point
    p2 : shapely Point

    x should be longitude
    y should be latitude
    """
    deg2km = 111.111  #
    aa1 = pp2.x - pp1.x  # distance in degrees
    aa2 = pp2.y - pp1.y  # distance in degrees.
    # change to meters.
    aa2 = aa2 * deg2km
    # estimate using latitude halfway between.
    aa1 = aa1 * deg2km * np.cos(np.radians(0.5 * (pp1.y + pp2.y)))
    return (aa1 ** 2 + aa2 ** 2) ** 0.5


def bearing(p1, p2):
    """
    p1 : shapely Point
    p2 : shapely Point

    x should be longitude
    y should be latitude
    """
    deg2met = 111.0  # doesn't matter.
    a1 = p2.x - p1.x  # distance in degrees
    a2 = p2.y - p1.y  # distance in degrees.
    # change to meters.
    a2 = a2 * deg2met
    # estimate using latitude halfway between.
    a1 = a1 * deg2met * np.cos(np.radians(0.5 * (p1.y + p2.y)))

    # a1 = np.cos(p1.y)*np.sin(p2.y)-np.sin(p1.y)*np.cos(p2.y)*np.cos(p2.x-p1.x)
    # a2 = np.sin(p2.x-p1.x)*np.cos(p2.y)
    angle = np.arctan2(a1, a2)
    angle = (np.degrees(angle) + 360) % 360
    return angle

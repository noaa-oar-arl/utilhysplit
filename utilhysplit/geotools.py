# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import numpy as np
from scipy.spatial import Delaunay
import shapely.geometry as sgeo
from shapely.ops import unary_union,  polygonize, transform


def plotpoly(sgeo_poly):
    """xy plot of a shapely polygon"""
    if isinstance(sgeo_poly, sgeo.multipolygon.MultiPolygon):
        for poly in sgeo_poly.geoms:
            x,y = poly.exterior.xy
            yield x,y
    else:
        x, y = sgeo_poly.exterior.xy
        # plt.plot(x,y)
        yield x, y


#def plotpoly(sgeo_poly):
#    """xy plot of a shapely polygon"""
#    x, y = sgeo_poly.exterior.xy
#    # plt.plot(x,y)
#    return x, y

def calculate_distance(lat1,lon1,lat2,lon2):
    p1 = sgeo.Point(lon1,lat1)
    p2 = sgeo.Point(lon2,lat2)
    return distance(p1,p2)

def distance(p1,p2):
    """
    p1 : shapely Point
    p2 : shapely Point

    x should be longitude
    y should be latitude

    Returns:
    approximate straitline distance in km.

    """
    deg2km = 111.111  #
    a1 = p2.x-p1.x # distance in degrees
    a2 = p2.y-p1.y # distance in degrees.
    # change to meters.
    a2 = a2 * deg2km
    # estimate using latitude halfway between.
    a1 = a1 * deg2km * np.cos(np.radians(0.5*(p1.y+p2.y))) 
    return (a1**2 + a2**2)**0.5

def bearing(p1, p2):
    """
    p1 : shapely Point
    p2 : shapely Point

    x should be longitude
    y should be latitude
    """
    deg2met = 111.0  # doesn't matter.
    a1 = p2.x-p1.x # distance in degrees
    a2 = p2.y-p1.y # distance in degrees.
    # change to meters.
    a2 = a2 * deg2met
    # estimate using latitude halfway between.
    a1 = a1 * deg2met * np.cos(np.radians(0.5*(p1.y+p2.y))) 

    #a1 = np.cos(p1.y)*np.sin(p2.y)-np.sin(p1.y)*np.cos(p2.y)*np.cos(p2.x-p1.x)
    #a2 = np.sin(p2.x-p1.x)*np.cos(p2.y)
    angle = np.arctan2(a1, a2)
    angle = (np.degrees(angle) + 360) %360
    return angle

def get_hull(z,thresh1=0.1,thresh2=1000,alpha=10):
    """
    z is a 2-d xarray data-array with coordinates latitude, longitude
    thresh1 and thresh2 are floats or ints.
    alpha : float


    """

    lon = z.longitude.values.flatten()
    lat = z.latitude.values.flatten()
    zzz = z.values.flatten()
    tlist = list(zip(lat,lon,zzz))

    # get lat lon values for values above thresh1 and below thresh2 and non nan.
    tlist = [x for x in tlist if ~np.isnan(x[2])]
    tlist = [x for x in tlist if x[2]>=thresh1]
    tlist = [x for x in tlist if x[2]<=thresh2]
    lon = [x[1] for x in tlist]
    lat = [x[0] for x in tlist]

    # create the polygons
    numpts = len(lon)
    mpts = make_multi(lon,lat)
    if numpts >= 4: 
        ch, ep = concave_hull(mpts,alpha=alpha)
    else:
        ch = mpts.convex_hull
        ep = None

    return ch, ep


def concave_hull(mpoints, alpha=1):
    """
    mpoints is a multi-point object.
    Returns
    rpoly : shapely.geometry.polygon.Polygon
    ep : list of numpy arrays.

    Note.  ep and mpoints can be fed into plot_delauney to plot the Delauney triangles.
           rpoly can be fed in as hull=rpoly. 
    """
    coords = np.array([x.coords[0] for x in mpoints.geoms]) 
    tri = Delaunay(coords)
    edges = set()
    edge_points = []
    #for ia, ib, ic in tri.vertices:
    done=False
    while not done:
      iii=0
      for ia, ib, ic in tri.simplices:
        pa = coords[ia]
        pb = coords[ib]
        pc = coords[ic]
        aaa = np.sqrt((pa[0]-pb[0])**2 + (pa[1]-pb[1])**2)
        bbb = np.sqrt((pb[0]-pc[0])**2 + (pb[1]-pc[1])**2)
        ccc = np.sqrt((pc[0]-pa[0])**2 + (pc[1]-pa[1])**2)

        sss = (aaa+bbb+ccc)/2.0
        area = np.sqrt(sss*(sss-aaa)*(sss-bbb)*(sss-ccc))
        circum_r = aaa*bbb*ccc/(4.0*area)
        if circum_r < 1.0/alpha:
           edges, edge_points = add_edge(edges, edge_points, coords, ia, ib)
           edges, edge_points = add_edge(edges, edge_points, coords, ib, ic)
           edges, edge_points = add_edge(edges, edge_points, coords, ic, ia)
           done=True
           #print('alpha', alpha,iii)
        if not done:
           alpha = alpha/2.0
        iii+=1
       
    mmm = sgeo.MultiLineString(edge_points)
    triangles = list(polygonize(mmm))
    
    #rpoly= cascaded_union(triangles)
    rpoly= unary_union(triangles)
    return rpoly, edge_points

def make_multi(x,y):
    """
    makes a multi-point object from list of x,y coordinates
    """
    from shapely.geometry import MultiPoint
    latlon = list(zip(x,y))
    mpts = MultiPoint(latlon)
    #poly = mpts.convex_hull
    #coords = list(poly.exterior.coords)
    return mpts

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



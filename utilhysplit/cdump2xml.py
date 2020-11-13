import os
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from lxml import etree
from lxml.builder import ElementMaker
import datetime
import logging
import zipfile

logger = logging.getLogger(__name__)

def color2kml(cstr):
    """
    python used RGB while KML expects BGR.
    """
    return cstr[0:2]+cstr[-2:] + cstr[4:6] + cstr[2:4]

class ColorMaker:

    def __init__(self,cname, nvals, ctype='hex', transparency='C8'):
        """
        cname : name of matplotlib colormap
        nvals : number of color values
        ctype : if 'hex' returns 8 digit hexidecimal with
                transparancy.
        transparency : str: transparency value to use in hexidecimal. 
        """
        self.transparency=transparency
        self.clist=[]
        self.get_cmap(cname,nvals)
        self.ctype = ctype
  
    def __call__(self):
        """
        Returns:
        list of nvals colors equally spaced throughout the colormap.
        and in hexidecimal format.
        """
        return self.clist

    def rgb_to_hex(self,rgb):
        """
        convert from rgb to hexidecimal.
        """
        def subfunc(val):
            rval =  hex(int(val*255)).replace('0x','').upper()
            if len(rval)==1: rval = '0'+rval
            return rval
        hval = [subfunc(x) for x in list(rgb)] 
        return '{}{}{}{}'.format(self.transparency,hval[0],hval[1],hval[2])

    def get_cmap(self,cname='viridis',nvals=10):
        cmap = matplotlib.cm.get_cmap(cname)
        cvals = cmap.N
        cspace = int(cvals/nvals)
        if self.ctype == 'hex'
            self.clist = [self.rgb_to_hex(cmap[x] for x in range(0,cvals,cpace)]
        else:
            self.clist = [cmap[x] for x in range(0,cvals,cpace)]
             
        #for iii in range(0,cvals,cspace):
        #    if ctype == 'hex':
        #        self.clist.append(self.rgb_to_hex(cmap(iii)))
        #    else:
        #        self.clist.append(cmap[iii])

def generate_kmz(filenames, kmz_filename,compression_level):
    for f in filenames:
        if not os.path.exists(f):
            logger.warn(
                "KML file {} does not exist. Google Earth file {} will not be created".format(
                    f, kmz_filename
                )
            )
            return
    logger.debug('Creating KMZ file')
    with zipfile.ZipFile(
        kmz_filename, "w", compresslevel=compression_level
    ) as z:
        for f in filenames:
            if os.path.exists(f):
                bn = os.path.basename(f)
                z.write(f,arcname=bn)


def make_legend(colors, levels, name, unit='g/m$^2$',date='',label=''):
    for val in zip(colors, levels):
        lstr = '>{} {}'.format(val[1],unit)
        plt.plot(1,1,color='#'+val[0],label=lstr,LineWidth=15)
    ax = plt.gca()
    handles,labels=ax.get_legend_handles_labels()
    plt.close()
    ht=3
    sz=10
    if len(levels)>8:ht=4 
    figlegend = plt.figure(figsize=(3,ht))
    #figlegend = plt.figure()
    ax1 = figlegend.add_subplot(1,1,1)
    ax1.legend(handles, labels, loc='center',prop={"size":sz},frameon=False)
    ax1.axis('off')
    ax1.text(0,1.1,date,transform=ax1.transAxes,size=sz)
    ax1.text(0,-0.1,label,transform=ax1.transAxes,size=sz)
    plt.savefig(name)


def add_overlay(screen, over,scr,rot,size):

    uhash = {'x':over[0],'y':over[1],'xunits':'fraction','yunits':'fraction'} 
    ovxy = etree.Element('overlayXY',uhash) 

    uhash = {'x':scr[0],'y':scr[1],'xunits':'fraction','yunits':'fraction'} 
    scxy = etree.Element('screenXY',uhash) 

    uhash = {'x':rot[0],'y':rot[1],'xunits':'fraction','yunits':'fraction'} 
    rotxy = etree.Element('rotationXY',uhash) 

    uhash = {'x':size[0],'y':size[1],'xunits':'pixels','yunits':'pixels'} 
    sz = etree.Element('size',uhash)

    screen.append(ovxy)
    screen.append(scxy)
    screen.append(rotxy)
    screen.append(sz)

def get_noaa_overlay():
    screen = etree.Element('ScreenOverlay')
    name = etree.SubElement(screen,'name')
    name.text = ('NOAA')
    snip = etree.Element('Snippet',{'maxLines':"0"})
    screen.append(snip)
    desc = etree.SubElement(screen,'description')
    desc.text = 'National Oceanic and Atmospheric Administration\
                 http:/www.noaa.gov'
    icon = etree.SubElement(screen,'Icon')
    logo = etree.SubElement(icon,'href')
    logo.text = 'noaa_google.gif'
    add_overlay(screen,["0","1"],["0.3","1"],["0","0"],["0","0"]) 
    return screen

def get_hysplit_overlay():
    screen = etree.Element('ScreenOverlay')
    name = etree.SubElement(screen,'name')
    name.text = ('HYSPLIT Information')
    desc = etree.SubElement(screen,'description')
    desc.text = 'NOAA ARL HYSPLIT Model http://www.arl.noaa.gov/HYSPLIT_info.php'
    icon = etree.SubElement(screen,'Icon')
    logo = etree.SubElement(icon,'href')
    logo.text = 'logocon.gif'
    add_overlay(screen,["1","1"],["0","1"],["0","0"],["0","0"]) 
    return screen

def add_source_style():
    style = etree.Element('Style', id='sorc')

    # IconStyle
    istyle = etree.SubElement(style, 'IconStyle')
    icolor = etree.SubElement(istyle, 'color')
    icolor.text = 'ff0000ff'
    iscale = etree.SubElement(istyle, 'scale')
    iscale.text = '0.8'
    icon = etree.SubElement(istyle,'Icon')
    href = etree.SubElement(icon,'href')
    href.text = 'icon63.png'
    hothash = {'x':"0.5","y":"0.5",'xunits':"fraction","yunits":"fraction"}
    hot = etree.SubElement(istyle, 'hotspot', hothash)    

    lab = etree.SubElement(style, 'LabelStyle')
    lcolor = etree.SubElement(lab,'color')
    lcolor.text = 'ff0000ff'
 
    line = etree.SubElement(style, 'LineStyle')
    ncolor = etree.SubElement(line,'color')
    ncolor.text = 'c8ffffff'
    nwidth = etree.SubElement(line,'width')
    nwidth.text = '2'
    return style
    
def add_point(coords):
    point = etree.Element('Point')
    ext = etree.SubElement(point,'extrude')
    ext.text = '1'
    amode = etree.SubElement(point,'altitudeMode')
    amode.text = 'clampedToGround'
    crds = etree.SubElement(point,'coordinates')
    crds.text = coords
    return point

def add_source(pname,description,coordinates):
    """
    pname : str
    description : str
    coordinates : str
    """
    folder = etree.Element('Folder')
    viz = etree.SubElement(folder,'visibility')
    # making this 0 means that contour levels do
    # not show up until they are chosen.
    viz.text = "1"
    fname = etree.SubElement(folder,'name')
    fname.text = 'Source Locations'
    pmk = etree.SubElement(folder, 'Placemark')
    pnm = etree.SubElement(pmk, 'name')
    pnm.text = pname
    descrip = etree.SubElement(pmk,'description')
    descrip.text = etree.CDATA(description)
    pmk.append(add_source_style())
    pmk.append(add_point(coordinates))
    return folder 

def add_placemark(namestr,t1,t2,stylestr):
    pmk = etree.Element('Placemark')
    name = etree.Element('name',LOC="Medium")
    name.text = namestr
    viz = etree.Element('visibility')
    viz.text = '1'
    snippet = etree.Element('Snippet',maxLines="0")
    snippet.text = ''
    style = etree.Element('styleUrl')
    style.text = stylestr
    pmk.append(name)
    pmk.append(viz)
    pmk.append(snippet)
    pmk.append(add_timespan(t1,t2))
    pmk.append(style)
    return pmk

def start_folder(iii=1,jobid=999,
             btime=datetime.datetime.now(), 
             etime=datetime.datetime.now()):
    """
    Folder up to the placemark. 
 
    iii : int
    btime : datetime object
    etime : datetime object
    Returns
    folder : etree Element
    """
    t1 = btime.strftime("%Y-%m-%dT%H:%M:00Z")
    t2 = etime.strftime("%Y-%m-%dT%H:%M:00Z")
    folder = etree.Element('Folder')
    name = etree.Element('name')
    ntext = '<pre>Mass Loading'
    validstr = '\nValid:{}UTC</pre>'.format(btime.strftime('%Y%m%d %H%M'))
    name.text = ntext + validstr
    name.text = etree.CDATA(name.text)
    folder.append(name)
    viz = etree.SubElement(folder,'visibility')
    viz.text = "1"
    op = etree.SubElement(folder,'open')
    op.text = "1"
    # TO DO add altitude range.
    description = etree.SubElement(folder,'description')
    dstr = '<pre>\nFROM {} to {}{}'.format(0,1000,validstr)
    description.text = etree.CDATA(dstr)
 
    screen = etree.Element('ScreenOverlay')
    name = etree.SubElement(screen,'name')
    name.text = ('Legend')
    snip = etree.Element('Snippet',{'maxLines':"0"})
    screen.append(snip)

    screen.append(add_timespan(btime,etime))

    desc = etree.SubElement(screen,'description')
    desc.text = 'National Oceanic and Atmospheric Administration\
                 http:/www.noaa.gov'
    icon = etree.SubElement(screen,'Icon')
    logo = etree.SubElement(icon,'href')
    logo.text = legend_namer(iii,jobid)
    #self.legend(legend_namer(iii))
    add_overlay(screen,["0","1"],["0","1"],["0","0"],["0","0"]) 
    folder.append(screen)

    return folder


def legend_namer(iii,jobid=''):
    """
    iii : int
    """
    rval = 'GELABEL_{}_{:02d}.png'.format(jobid,iii)
    return rval 

def path2str(path,ht):
    temp = ["{:0.3f},{:0.3f},{:0.2f}".format(x[0],x[1],ht) for x in path] 
    return str.join('\n',temp)

def generate_contours(cxra,levels,ht):
    qcs = plt.contour(cxra.longitude.values, 
                      cxra.latitude.values, 
                      cxra.values,
                      levels = levels)
    iii=0
    # TO DO this also yields holes and very small pieces.
    # some way to filter those out?
    for ccc in qcs.allsegs:
        for path in ccc:
            yield levels[iii], path2str(path,ht)
        iii+=1
                     
def add_timespan(t1,t2):
    tspan = etree.Element('TimeSpan')
    begin = etree.Element('begin')
    end = etree.Element('end')
    begin.text = t1.strftime("%Y-%m-%dT%H:%M:00Z")
    end.text = t2.strftime("%Y-%m-%dT%H:%M:00Z")
    tspan.append(begin)
    tspan.append(end) 
    return tspan

def add_polygon(coords):
    """
    coords : str
    """
    multi = etree.Element('MultiGeometry')
    polygon = etree.Element('Polygon')
    extrude = etree.Element('extrude')
    extrude.text = '1'
    amode = etree.Element('altitudeMode')
    amode.text = 'clampedToGround'
    boundary = etree.Element('outerBoundaryIs')
    lring = etree.Element('LinearRing')
    coord = etree.Element('coordinates')
    coord.text = coords
    polygon.append(extrude)
    polygon.append(amode)
    lring.append(coord)
    boundary.append(lring)
    polygon.append(boundary)
    multi.append(polygon)
    return multi

def add_style(sid,line_color,poly_color):
    """
    Style
    LineStyle
    color
    PolyStyle
    fill
    outline
    """
    style = etree.Element('Style',id=sid)
    line = etree.Element('LineStyle')
    lcolor = etree.Element('color')
    lcolor.text = line_color
    poly = etree.Element('PolyStyle')
    pcolor = etree.Element('color')
    pcolor.text = poly_color
    fill = etree.Element('fill')
    fill.text='1'
    outline = etree.Element('outline')
    outline.text='1'

    poly.append(pcolor)
    poly.append(fill)
    poly.append(outline)

    line.append(lcolor)

    style.append(line)
    style.append(poly)
    return style


def ppxml(xml):
    parser = etree.XMLParser(remove_blank_text=True)
    tree = etree.parse(xml,parser)
    tree.write(xml,encoding='utf-8',pretty_print=True,xml_declaration=True)


class HysplitKml:
    """
    Created specifically for the volcanic ash application.

    Attributes:
        levels : list of floats. (contour levels)
        stylehash : dictionary

    Methods
        create : creates and writes the kml file.
        write :
        set_style_dict : can be used to change colors.

    """

    # usage
    # add_root()
    # add_header()
    # add_contours_timeloop
    # write

    def __init__(self,levels,sourcehash={},
                 units='g/m2',jobid=999,
                 legend_label='Mass Loading'):
        #self.cxra = cxra      
        self.levels = levels
        self.stylehash = {}
         
        self.clist = self.get_default_colors()
        self.set_style_dict(self.clist)
 
        self.__parse_hash(sourcehash)
        self.units = units

        self.jobid = jobid
        self.legend_label = legend_label
        # used for making the kmz file.
        self.legend_file_list = []
        # set in add_root method.
        #self.root
        #self.document


    def __parse_hash(self,shash):

        tstr = 'VolcanoName'
        if tstr in shash.keys():
            self.sourcename = shash[tstr]
        else:
            self.sourcename = 'Source'

        tstr = 'start_date'
        if tstr in shash.keys():
            d1 = shash[tstr]
        else:
            d1 = datetime.datetime.now()
        self.source_time=d1.strftime("%Y-%m-%dT%H:%M:00Z")
            
        tstr = 'latitude'
        if tstr in shash.keys():
            self.latitude = shash[tstr]
        else:
            self.latitude = 45

        tstr = 'longitude'
        if tstr in shash.keys():
            self.longitude = shash[tstr]
        else:
            self.longitude = -75 

        tstr = 'top'
        if tstr in shash.keys():
            self.release_alt = shash[tstr]
        else:
            self.release_alt = 'unknown'
        tstr = 'bottom'
        if tstr in shash.keys():
            self.release_bottom = shash[tstr]
        else:
            self.release_bottom = 'unknown'

    def add_source_folder(self):
        pname = self.sourcename
        descrip = self.create_source_description()
        try:
            alt = float(self.release_alt)
        except:
            alt = 100
        coords='{},{},{}'.format(self.longitude,self.latitude,alt)
        self.document.append(add_source(pname,descrip,coords))

    def create_source_description(self):
        cstr = '<pre> Lat: {:0.3f}  LON: {:0.3f}'.format(self.latitude,self.longitude)
        cstr += '\n'
        cstr += 'Source Name {}'.format(self.sourcename)
        cstr += '\n'
        cstr += 'Released between {} and {}'.format(self.release_bottom,self.release_alt)  
        cstr += '\n'
        cstr += 'Release time {}'.format(self.source_time)
        cstr += '\n'
        cstr += '</pre>'
        return cstr

    def create(self,cra, attrs,kname):
        self.add_root()
        self.add_header()
        self.add_contours_timeloop(cra, attrs)
        logger.debug('writing {}'.format(kname))
        self.write(kname)
        self.kname = kname

    def write(self,fname):
        topstr = '<?xml version="1.0" encoding="UTF-8"?>\n'
        with open(fname,'w') as fid:
             fid.write(topstr)
             fid.write(str(self)) 

    def create_kmz(self,compression_level,efiles=None):
        zfiles = [self.kname]
        zfiles.extend(self.legend_file_list)
        if efiles: zfiles.extend(efiles)
        zname = self.kname.replace('kml','kmz')
        generate_kmz(zfiles, zname,compression_level)        

    def add_root(self):
        """
        creates beginning of document.
        up to Folder  
        """
        name1 = 'http://www.google.com/kml/ext/2.2'
        name2  = "http://www.opengis.net/kml/2.2"
        self.KML_NAMESPACE={None:name2,'gx':name1}
        # this is used for tags like <gx:TimeStamp>
        self.gxbuild = ElementMaker(namespace=name1,nsmap=self.KML_NAMESPACE)

        self.root = etree.Element('kml',nsmap=self.KML_NAMESPACE)
        self.document= etree.SubElement(self.root, 'Document')

    def add_header(self):
        self.add_intro()
        self.auto_add_styles()
        self.add_source_folder()
        self.document.append(get_hysplit_overlay()) 
        self.document.append(get_noaa_overlay()) 

    def add_intro(self):
        name = etree.SubElement(self.document,'name')
        name.text = 'NOAA HYSPLIT RESULTS'
        op = etree.SubElement(self.document,'open')
        op.text = '1'
        look = etree.SubElement(self.document,'LookAt')
        longitude = etree.SubElement(look,'longitude')
        longitude.text = '{}'.format(self.longitude)
        latitude = etree.SubElement(look,'latitude')
        latitude.text = '{}'.format(self.latitude)
        altitude = etree.SubElement(look,'altitude')
        altitude.text = '0'
        tilt = etree.SubElement(look,'tilt')
        tilt.text='0'
        lrange = etree.SubElement(look,'range')
        lrange.text='13700'
       
        ts = self.gxbuild.TimeStamp('')
        look.append(ts) 
        when = etree.SubElement(ts,'when')
        when.text = self.source_time
        look.append(ts)
       
        altmode = self.gxbuild.altitudeMode('relativeToSeaFloor')
        look.append(altmode)

    def legend(self,name='legend.png',date='',label=''):
        self.legend_file_list.append(name) 
        clist = self.clist
        clist = [x[1][2:] for x in clist]
        make_legend(clist, self.levels, name,date=date,label=label)         

    def add_style_to_root(self,tag):
        linecolor = self.stylehash[tag][0]   
        polycolor = self.stylehash[tag][1]   
        self.document.append(add_style(tag,linecolor,polycolor))


    def get_default_colors(self,cmap='plasma'):
        linecolor = 'C8000000'
        cm = ColorMaker(cmap,len(self.levels))
        alist = cm()
        clist = [(linecolor,x) for x in alist]
        return clist

    def set_style_dict(self,clist):
        # C8 sets the transparency.
        if not clist:
           clist = self.get_default_colors()
        self.clist = clist
        nmax = len(self.levels)+1
        for iii in range(1,nmax):        
            clr = (color2kml(clist[iii-1][0]),color2kml(clist[iii-1][1]))
            self.stylehash['conc{}'.format(iii)] = clr

    def auto_add_styles(self):
        """
        add a style for each contour level.
        """
        iii=1
        for lev in self.levels:
            tag = 'conc{}'.format(iii)
            self.add_style_to_root(tag)
            iii+=1

    def add_style_hash(self, tag, color1, color2):
        """
        Add a stye to the sytlehash.
        """
        self.stylehash[tag] = (color1,color2)

    def add_contours_timeloop(self, cra,attrs):
        """
        cra : xrarray DataArray. with dimensions time, x, y
              and coordiantes latitude, longitude
        attrs : dictionary
        """
        subtract=False
        if 'sample time hours' in attrs.keys():
           dtime = datetime.timedelta(hours=attrs['sample time hours'])
        else:
           dtime = datetime.timedelta(hours=1)
        if 'time description' in attrs.keys():
           if 'end' in attrs['time description']:
               subtract=True
        iii=0
        for time in cra.time.values:
            tempra = cra.sel(time=time)
            t1 = pd.to_datetime(time)
            if not subtract:
                t2 = t1 + dtime
            else:
                t2 = t1 - dtime
            self.add_contours(tempra,iii,t1,t2)
            iii+=1

    def add_contours(self, tempra, iii=1,
                    t1 = datetime.datetime.now(),
                    t2 = datetime.datetime.now()):
        levels = self.levels
        height=100
        folder = start_folder(iii,self.jobid,t1,t2)
        datestr = t1.strftime("%H:%M UTC %b %d %Y")
        datestr += ' to \n'
        datestr += t2.strftime("%H:%M UTC %b %d %Y")
        self.legend(legend_namer(iii,self.jobid),datestr,self.legend_label)
        for level, coords in generate_contours(tempra, levels,height): 
            iii = self.levels.index(level) 
            name = 'Contour Level: {} {}'.format(level,self.units)
            placemark = (add_placemark(name,t1,t2,'#conc{}'.format(iii+1)))
            placemark.append(add_polygon(coords))
            folder.append(placemark) 
        self.document.append(folder)

    def __str__(self):
        #return etree.tostring(self.root,pretty_print=True).decode() 
        return etree.tostring(self.root,pretty_print=True,encoding='unicode') 

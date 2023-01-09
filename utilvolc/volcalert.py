#volcalert.py
#Reads XML alert files
#Extracts: date, time, lat, lon, volcano name
#max height (km), max height (ft), total mass (Tg), median effective radius (um)
import lxml.etree as let
import lxml.objectify as obj
import numpy as np
from datetime import datetime

# volcalert is used to read the volcat alert files.
# There is a another volcalert.py in the autoash website.
# volcat_log


#Opening xml file
def open_xml(xmlFile):
    """Opens the xml alert file from VOLCAT"""
    with open(xmlFile) as f:
        xml = f.read()
        root = obj.fromstring(xml)    
    return root

#Extracting necessary values from alert
def get_alert_values(root):
    """Finds important alert info, returned in this order:
    1. Alert header
    2. Alert type
    3. Alert status
    4. Alert confidence
    5. Alert latitude (of radiative center)
    6. Alert longitude (of radiative center)
    7. Alert date and time (datetime object)
    8. Satellite identifier number"""

    wmo_id = root.summary.wmo_spacecraft_id.attrib.get('value').strip()
    alert_header = root.alert.alert_header.attrib.get('value').strip()
    alert_type = root.alert.alert_type.attrib.get('value')
    status = root.alert.status.attrib.get('value')
    confidence = root.alert.confidence.attrib.get('value')
    method = root.alert.method.attrib.get('value')
    lat_rc = float(root.alert.lat_rc.attrib.get('value'))
    lon_rc = float(root.alert.lon_rc.attrib.get('value'))
    date_time = root.alert.object_date_time.attrib.get('value') #Put into datetime object
    date_time = datetime.strptime(date_time, '%Y-%m-%d %H:%M:%S')
    return alert_header, alert_type, status, confidence, method, lat_rc, lon_rc, date_time, wmo_id

#Extracting array of nearby possible volcanoes
def get_nearby(root):
    """Finds information for the nearby volcanoes to the radiative center lat/lon
    Returned in this order:
    1. List of nearby volcano names
    2. List of nearby volcano latitudes
    3. List of nearby volcano longitudes
    4. List of nearby volcano distances from radiative center
    5. List of nearby volcano termal anomaly indications
    6. List of nearby volcano ID numbers"""
    values=[]
    for data in root.alert.volcanoes.getchildren():
        for volc in data.getchildren():
            values.append(volc.attrib.get('value'))
    i=0
    name=[]
    lat=[]
    lon=[]
    dist=[]
    therm=[]
    vid=[]
    #.strip() removes the whitespace in the string
    while i < len(values):
        name.append(values[i].strip())
        lat.append(values[i+1].strip())
        lon.append(values[i+2].strip())
        dist.append(values[i+3].strip())
        therm.append(values[i+4].strip())
        vid.append(values[i+5].strip())
        i+=6
    return name, lat, lon, dist, therm, vid

#Finding closest volcano (minimum distance)
#NOT NECESSARY - can easily get info from get_nearby return
def get_closest(name, lat, lon, dist, therm, vid):
    """Extracts information for the closest volcano to the radiative center.
    Uses first volcano in "nearby volcanoes" list. Values are returned in this order:
    1. Closest volcano name
    2. Closest volcano latitude
    3. Closest volcano longitude
    4. Closest volcano thermal indicator
    5. Closest volcano ID number"""
    closest = 1 #Using first volcano in list
    closest_name = name[closest-1]
    closest_lat = lat[closest-1]
    closest_lon = lon[closest-1]
    closest_therm = therm[closest-1]
    closest_vid = vid[closest-1] 
    return closest_name, closest_lat, closest_lon, closest_therm, closest_vid

#Finding VAAC region information
def get_vaac(root):
    """Returns the VAAC region of the alert"""
    vaac = root.alert.vaac_region.attrib.get('value')
    return vaac
    
#Extract plume height
#Must be an ash alert for this info!
def get_height(root):
    """Extracts the ash plume height - only valid for ash alerts!
    Values returned are in this order:
    1. Maximum ash height (kilometers)
    2. Maximum ash height (feet)
    3. Tropopause height (kilometers)
    4. Tropopause height (feet)"""
    hgt_km = root.alert.max_height.attrib.get('value')
    hgt_ft = root.alert.max_height_feet.attrib.get('value')
    tropohgt_km = root.alert.tropo_height.attrib.get('value')
    tropohgt_ft = root.alert.tropo_height_feet.attrib.get('value')
    return hgt_km, hgt_ft, tropohgt_km, tropohgt_ft

#Extract ash mass loading
def get_mass(root):
    """Extracts the total ash mass loading - only valid for ash alerts!
    Values returned are in this order:
    1. Total ash mass
    2. Total ash mass unit"""
    mass = root.alert.total_mass.attrib.get('value')
    mass_unit = root.alert.total_mass.attrib.get('units')
    return mass, mass_unit

#Extract ash effective radius
def get_radius(root):
    """Extracts the median ash effective radius - only valid for ash alerts!
    Values returned are in this order:
    1. Median ash effective radus
    2. Median ash effective radius unit"""
    effrad = root.alert.effective_radius.attrib.get('value')
    effrad_unit = root.alert.effective_radius.attrib.get('units')
    return effrad, effrad_unit

#Extract ash cloud total area
def get_area(root):
    """Extracts the ash cloud total area
    Values returned are in this order"
    1. Total ash area
    2. Total ash area unit"""
    area = root.alert.total_area.attrib.get('value')
    area_unit = root.alert.total_area.attrib.get('units')
    return area, area_unit

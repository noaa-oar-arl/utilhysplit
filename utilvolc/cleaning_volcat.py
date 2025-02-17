import os
import glob
import pandas as pd
import numpy as np

from utilvolc.volcat import flist2eventdf

import utilvolc.iwxxmVAA as ixa
from utilvolc import volcano_information
from utilvolc import volcano_names


# help remove files that are needed anymore

def inventory():
    # looks to see what VOLCAT files are currently stored.
    dhash = {}
    shash = {}
    tdir = '/pub/ECMWF/JPSS/VOLCAT/Files/'
    for dirpath,dirnames,filenames in os.walk(tdir):
        for ddir in dirnames:
            if dirpath == tdir:
               vfiles = glob.glob(tdir + ddir + '/*VOLCAT*nc')
               sz = 0
               for vf in vfiles:
                   sz += os.path.getsize(os.path.join(tdir,ddir,vf))
               print(ddir, sz/1e6)
               dhash[ddir] = flist2eventdf(vfiles, {"VOLCANO_NAME": ddir})
               shash[ddir] = sz 
    return dhash, shash


def check_vaa(years, vname):
    yhash = {}
    for year in years:
        vaac = ixa.WashingtonPage(year=year)
        vaac.read()
        vaac.find_xml()
        flist = vaac.get_xml_list(vname=vname)
        print(year, len(flist))
        yhash[year] = flist
    return yhash 


def name2vaac(name):
    nhash = {}
    nhash['Turrialba'] = 'TURR'
    nhash['Santa_Maria'] = 'SANTA MARIA'
    nhash['Concepcion'] = 'CONC'
    nhash['Conchaguita'] = 'CONCH'
    if name in nhash.keys():
       return nhash[name]
    else: 
       return name.upper()[0:4]


def washington_vaac_names():
    # in archives since 2021
    ilist = []
    ilist.append('Atitlan')
    ilist.append('')
    ilist.append('Bezymianny')
    ilist.append('Chikurachki')
    ilist.append('Concepcion')
    ilist.append('Cotopaxi')
    ilist.append('Fuego')
    #ilist.append('Fukutoku-Okanoba')
    ilist.append('Fukutoku-Oka-no-Ba')
    ilist.append('Karymsky')
    ilist.append('Katmai')
    ilist.append('Kilauea')
    #ilist.append('Kliuchevskoi')
    ilist.append('Klyuchevskoy')
    ilist.append('Mauna_Loa')
    ilist.append('Momotombo')
    ilist.append('Pacaya')
    ilist.append('Pagan')
    ilist.append('Pavlof')
    ilist.append('Poas')
    ilist.append('Popocatepetl')
    ilist.append('Quilotoa')
    ilist.append('Reventador')
    ilist.append('Rincon_de_la_Vieja')
    ilist.append('Ruiz, Nevado del')
    ilist.append('Nevado_del_Ruiz')
    ilist.append('San Cristobal')
    ilist.append('Sangay')
    ilist.append('San_Miguel')
    ilist.append('Santa_Maria')
    ilist.append('Semisopochnoi')
    ilist.append('Sheveluch')
    # this is from 2021 and shown as Soufriere hills
    # however is Soufriere St. Vincent in VAA
    # does not have xml files.
    # https://www.ospo.noaa.gov/products/atmosphere/vaac/2021.html#SOUF
    ilist.append('Soufriere_St._Vincent')
    ilist.append('Shishaldin')
    ilist.append('Telica')
    ilist.append('Turrialba')
    ilist.append('Tungurahua')
    ilist.append('Wolf')
    ilist.append('Volcan_Azul')
    return ilist


def keep_names():
    ilist = []
    ilist.append('Ruang')
    ilist.append('Raikoke')
    ilist.append('Cleveland')
    ilist.append('Veniaminof')
    ilist.append('Etna')
    ilist.append('Karymsky')
    ilist.append('La_Palma')
    ilist.append('Colima')
    ilist.append('Great_Sitkin')
    ilist.append('Nishinoshima')
    ilist.append('SheveluchA')
    ilist.append('Gareloi')
    return ilist


class Summary:

    def __init__(self):
        dhash, shash = inventory()
        self.dhash = dhash
        self.shash = shash

    def list_by_name(self):
        return sorted(self.shash.items())

    def list_by_size(self):
        return sorted(self.shash.items(), key=lambda x: x[1])

    @property
    def vnames(self):
        return list(self.dhash.keys())


    def reject(self):
        ilist =  washington_vaac_names()
        klist =  keep_names()
        keys = list(self.shash.keys())
        newlist = [x for x in keys if x not in ilist]
        newlist = [x for x in newlist if x not in klist]
        newlist = [x for x in newlist if x[0] != '3']
        with open('reject.sh','w') as fid:
            for new in sorted(newlist):
                fid.write('# {:20s}  {:5d} {:10.1f}\n'.format(new, len(self.dhash[new]), self.shash[new]/1e6)) 
                fid.write('rm -rf {}\n'.format(new))
        return newlist 

    def check(self):
        keys = list(self.shash.keys())
        for key in sorted(keys):
            vaaname = name2vaac(key)
            print(key, vaaname, len(self.dhash[key]), self.shash[key]/1e6)
            try: 
                odates = self.dhash[key].observation_date.values
            except:
                print('No files found')
                continue
            years = [pd.to_datetime(x).year for x in odates]
            years = list(set(years))
            check_vaa(years,vaaname)
            print('-----\n')       


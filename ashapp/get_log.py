import glob
import os
import pandas as pd
import re
import sys
import subprocess
from utilvolc import runhelper
from bs4 import BeautifulSoup
from utilvolc import volcat_name, volcano_information

def greenlist(name):
      """
      Check if the name is in the greenlist.
      Case-insensitive comparison by converting input to lowercase.
      """
      greenlist = ['popocatepetl', 'fuego', 'colima',
                  'pacaya',  
                  'turrialba', 'nevado del ruiz', 'reventador', 'cotopaxi', 
                  'tungurahua', 'sangay']
      return name.lower() in greenlist


def lookup_washington(key):
    dhash = {}
    dhash['342090'] = 'Fuego'
    dhash['341040'] = 'Colima'
    dhash['321010'] = 'BAKER'
    dhash['322040'] = 'SAND MOUNTAIN FIELD'
    dhash['341040'] = 'COLIMA'
    dhash['341090'] = 'POPOCATEPETL'
    dhash['341097'] = 'LA GLORIA'
    dhash['342030'] = 'SANTA MARIA'
    dhash['342090'] = 'FUEGO'
    dhash['342110'] = 'PACAYA'
    dhash['342200'] = 'None'
    dhash['344020'] = 'SAN CRISTOBAL'
    dhash['344140'] = 'None'
    dhash['345070'] = 'TURRIALBA'
    dhash['351020'] = 'NEVADO DEL RUIZ'
    dhash['351060'] = 'PURACE'
    dhash['352010'] = 'REVENTADOR'
    dhash['352050'] = 'COTOPAXI'
    dhash['352060'] = 'QUILOTOA'
    dhash['352080'] = 'TUNGURAHUA'
    dhash['352090'] = 'SANGAY'
    dhash['360170'] = 'ST. CATHERINE'
     

    if key in dhash.keys():
       return dhash[key]
    else:
       return key

def generate_links(iname):
    with open(iname, 'r') as fid:
          soup = BeautifulSoup(fid.read(), 'html.parser')
          # Find all links in the HTML
          for link in soup.find_all('a'):
              href = link.get('href')
              next_sibling = link.next_sibling
              szb = None
              if next_sibling:
                  match = re.search(r'\((\d+)\s+bytes\)', next_sibling)
                  if match:
                       szb = int(match.group(1))
              yield href, szb

if __name__ == "__main__":
   vaac='Washington'
   helper = runhelper.Helper
   helper.remove('{}.1/index.html'.format(vaac))
   helper.remove('index.html')
  
   indexfile = f'{vaac}.1'

   print(indexfile)
   # Corrected FTP URL - ensure we're listing the directory, not retrieving a file
   ftp_url = 'ftp://ftp.ssec.wisc.edu/pub/volcat/events/{}/'.format(vaac)
   try:
       result = subprocess.run(["wget", "-P", "./", ftp_url, "-O", indexfile], 
                              capture_output=True, text=True)
       # wget returns non-zero for FTP directory listings sometimes, so don't use check=True
       print(f"wget output: {result.stdout}")
       print(f"wget errors: {result.stderr}")
   except subprocess.CalledProcessError as e:
       print(f"Error running wget: {e.stderr}")
       sys.exit(1)
   except FileNotFoundError:
       print("wget command not found. Please install wget.")
       sys.exit(1)
   vlist = []
   notlist = []
   for href, sz in generate_links(indexfile):
       print(href)  
       vid = href.split('/')[-2]
       print(vid)
       ifile = f'index.{vid}' 
       result = subprocess.run(["wget", "-P", "./", href+'/netcdf/', '-O', ifile], 
                capture_output=True, text=True)
       for nref,sz in generate_links(ifile):
           #vname = nref.split('/')[-1] 
           vname = volcat_name.VolcatName25(nref)
           vname.vhash['size (bytes)'] = sz
           name = lookup_washington(vname.vhash['event vid'])
           vname.vhash['vname'] = name
           if greenlist(name):
              subdir = f'/pub/ECMWF/JPSS/VOLCAT/Files/{name}/'
              runhelper.make_dir(subdir, newdir=None, verbose=True)
              if not os.path.isfile(subdir + nref):
              # Download the file if it does not exist
                  download = subprocess.run(["wget", "-P", subdir, nref], 
                    capture_output=True, text=True)
                  print(f'download {nref} to {subdir}')
              vname.vhash['downloaded'] = True  # Fixed: Changed vhash to vname.vhash
           else:
               notlist.append(nref)
               vname.vhash['downloaded'] = False  # Fixed: Changed vhash{} to vname.vhash[] and syntax error
           vlist.append(vname.vhash)

   for name in notlist:
      print(f'Not Downloading {name} as it is not on the greenlist.')  # Fixed: Corrected message text
 
   df = pd.DataFrame(vlist) 
   print(df) 
   logfile = f'/pub/ECMWF/JPSS/VOLCAT/logs/{vaac}_log.csv'
   if not os.path.isfile(logfile):  # Fixed: Removed duplicate 'not'
      print(f'Creating new log file: {logfile}')
      df.to_csv(logfile, index=False)
   else: 
      olddf = pd.read_csv(logfile)
      df = pd.concat([olddf, df], ignore_index=True)
      df = df.drop_duplicates(subset=['url'], keep='last')  # Note: Ensure 'url' column exists
      df = df.sort_values(by='sdate', ascending=False)  # Note: Ensure 'sdate' column exists
      df.to_csv(logfile, index=False)  # Fixed: Used logfile variable instead of hardcoded path with typo

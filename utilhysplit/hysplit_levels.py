import numpy as np

def metlvl(plev):

  sfcp = 1013

  # meters per hPa by 100 hPa intervals (100 to 900)
  zph1 = [17.98,14.73,13.09,11.98,11.15,10.52,10.04,9.75,9.88 ]

  # meters per hPa by 10 hPa intervals (10 to 90)
  zph2 = [31.37,27.02,24.59,22.92,21.65,20.66,19.83,19.13,18.51]

  if plev>100:
     iii = int(plev/100.0)-1
     if iii>8: iii=8
     delz = zph1[iii]
  else:
     iii = int(plev/10.0)-1
     if iii>8: iii=8
     delz = zph2[iii]

  zlev = (sfcp-plev)*delz
  return zlev

def est():


  # meters per hPa by 100 hPa intervals (100 to 900)
  zph1 = [17.98,14.73,13.09,11.98,11.15,10.52,10.04,9.75,9.88 ]

  # meters per hPa by 10 hPa intervals (10 to 90)
  zph2 = [31.37,27.02,24.59,22.92,21.65,20.66,19.83,19.13,18.51]

  return zph1, zph2


def hypsometric(p2,p1,tbar):
    rdry = 287.04 #J/Kg-K
    grav = 9.8 #m/s2
    delz = (p2-p1)*rdry*tbar/grav
    return delz

def nasa(p):
    p = p/10 

    #if p>100:
    #    t = (1/5.256)*np.log10(p/101.29)
    #else p<100:
    h = (1.73 - np.log(p/22.65))/0.000157

    return h


def nasalvl(h):
    # h height in meters
    if h<11000:
        T = 15.04 - 0.00649*h
        p = 101.29*((T+273.1)/288.08)**5.256
    elif h<25000:
        T = -56.46
        p = 22.65*np.exp(1.73-0.000157*h)
    else:
        T = -131.21 + 0.00299*h
        p = 2.488*((T+273.1)/216.6)**-11.388
        
    pres = p*10 #convert to hPa
    return pres

#volcMER.py
#Calculates mass eruption rates in a variety of ways, eruptive volume, eruptive mass

def mastinMER(H, DRE = 2500):
    """Calculates Mass Eruption Rate (MER) using Mastin et al. 2009 equation.
    Input is height (km), dense rock equivalent (DRE) (kg/m^3). 
    dV = volumetric flow rate (m^3/s)
    Output is MER (kg/s). Optional input is DRE"""
    p=1/0.241
    dV=(H/2.00)**p
    MER=dV*DRE
    return MER 

def sparksMER(H, DRE = 2500):
    """Calculates Mass Eruption Rate (MER) using 
    Sparks equation from Volcanic Plumes (1997).
    Input is height (km). dense rock equivalent (DRE) (kg/m^3).
    dV = volumetric flow rate (m^3/s)
    Output is MER (kg/s). Optional input is DRE"""
    p=1/0.259
    dV=(H/1.67)**p
    MER=dV*DRE
    return MER

def MER2unit(MER, M63=0.1):
    """Solving for unit mass - comparable to HYSPLIT output
    Inputs: MER in kg/s, mass fraction of fine ash (M63) Fine ash is < 63um according to Mastin et al. 2009. 
    The mass fraction of fine ash is highly variable and poorly known.
    Output: unit mass in grams (unit_mass), unit mass of fine ash in grams (mass63). 
    Assume model output is one unit mass (per hour - based on MER conversion)."""
    unit_mass = MER * 3600 * 1000
    mass63 = unit_mass * M63
    return unit_mass, mass63

def HT2unit(HT, M63=0.1, verbose=True):
    """Using Mastin equation to calculate Mass Eruption Rate from plume height (km)
    Inputs: Height (HT) float (in kilometers)
    Output: unit mass in grams. 
    Assume model output is one unit mass per hour. """
    MER = mastinMER(HT)
    unit_mass = MER2unit(MER,M63=M63)
    if verbose: print('HEIGHT %0.1f km,  MER %0.3e kg/s , M63 %0.2f , unit mass=%0.3e g/hr.' %(HT, MER, M63, unit_mass))
    return unit_mass

def mastinEVEM(HT, DRE = 2500):
    """Calculating erupted volume based on plume height.
    Inputs: Plume height (HT) in km
    DRE = dense-rock equivalent (kg/m^3)
    Outputs: Eruptive Volume (km^3 DRE) and Eruptive Mass (kg)
    No time component to this??
    Equation from Mastin et al. 2009: H = 25.9 + 6.64log(V) --> log base 10"""
    p = (HT - 25.9)/6.64
    EV = 10**p
    EM = EV * DRE * (1000**3)
    return EV, EM

def sizedist():
    """
    As of 5/20/2020 operational volcanic ash runs
    utilize four particle sizes.
    The mass is distributed as in the list mfrac
    with the smallest particles receiving the least mass.
    """
    mfrac = (0.01,0.07,0.25,0.67)
    return mfrac 

def volume(dlat,dlon,dh,rlat=0):
    """
    
    dlat, dlon : float. should be in degrees
    dh : float : should be in km
    """
    import numpy as np
    adj = np.cos(rlat * np.pi / 180.0)  
    deg2km = 111.0
    km2m = 1e3
    rval = adj * dlat * dlon * (deg2km)**2 * dh * (km2m)**3
    return rval # in meters cubed. 


def est_numpar(HT, dlat, dlon, dh, 
               rlat=0, M63=0.05, 
               conc=0.2, 
               min_num = 10,
               verbose=False):
    """
    Can be used to help estimate number of particles needed to
    simulate concentrations down to 0.2 mg/m3.
    dlat : float : degrees
    dlon : float : degrees
    dh : float : km
    rlat : float : approximate latitude
    M63 : float : mass fraction fine ash
    min_num : int : number of particles needed in grid cell to get desired
                    concentration.
    conc : float : concentration of interest (mg/m3)
    """
    # Need min_num 
    # to obtan the desired concentration.

    MER = mastinMER(HT)  #mass eruption in kg/s
    MPH = MER * 3600 * 1000 * 1000 # convert to mg/hour
    mass = MPH * M63        # apply mass fraction of fine ash.
    # how much mass is distributed to each particle size (in mg).
    masslist = list(map(lambda x: x*mass, sizedist()))
   
    # volume of concentration grid.
    vol = volume(dlat,dlon,dh,rlat)
    
    # amount of mass on one particle
    mg_per_particle = conc * vol / min_num

    # how many particles are needed to get mass for one particle size 
    numpar = list(map(lambda x: x / mg_per_particle, masslist))
    if verbose: print(numpar)
   
    # sort from smallest to largest 
    numpar.sort()
    # return largest value times number of particle sizes.
    # as HYSPLIT splits particles evenly between species.
    return int(numpar[-1] * len(sizedist()))


def calc1(HT, dlat,dlon,dh,rlat=0,N=20000,M63=0.05, verbose=True):
    """
    Can be used to help estimate minimum concentration that
    can be reliably calculated.
    """
    M63 = 0.05
    MER = mastinMER(HT)  #mass eruption in kg/s
    if verbose:
        print('Number of particles emitted per hour')
        print(str(N))
        print('Mass fraction of fine ash')
        print(str(M63))
        print('Mass eruption rate (kg/s')
        print('{:0.2e}'.format(MER))
        print('Height of plume (km)')
        print(str(HT))

    MPH = MER * 3600 * 1000 * 1000 # convert to mg/hour
    mass = MPH * M63        # apply mass fraction of fine ash.
    # how much mass is distributed to each particle size.
    masslist = list(map(lambda x: x*mass, sizedist()))

    if verbose:
        print('Fraction of mass in each particle size bin')
        print(', '.join('{:0.2f}'.format(x) for x in sizedist()))   
        print('Amount of mass in each particle size bin (mg)')
        print(', '.join('{:0.2e}'.format(x) for x in masslist))   
 
    # how much mass is on each particle
    mass_per_particle = list(map(lambda x: x / (N / len(masslist)), masslist))
    vol = volume(dlat,dlon,dh,rlat)
    # concentration resulting from 1 particle in a grid cell.
    conc = list(map(lambda x: x/vol, mass_per_particle))
    # number of particles needed to obtain 0.2 mg/m3
    np = list(map(lambda x: 0.2/x, conc))

    if verbose:
        print('Amount of mass on each particle (mg)')
        print(', '.join('{:0.2e}'.format(x) for x in mass_per_particle))   
        print('Volume of grid cell (m^3)')
        print('{:0.2e}'.format(vol))
        print('dlat x dlon x dh (km)')
        print(str(dlat) + 'x' + str(dlon) + 'x' + str(dh))

        print('Concentration from 1 particle being in grid cell (mg/m3)')
        print(', '.join('{:0.2e}'.format(x) for x in conc))   

        print('Number of particles needed to obtain 0.2mg/m3')
        print(', '.join('{:0.2f}'.format(x) for x in np))   

    return conc, np





        


     
   

 












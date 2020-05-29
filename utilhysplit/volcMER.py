# volcMER.py
# Calculates mass eruption rates in a variety of ways, eruptive volume, eruptive mass


def mastinMER(H, DRE=2500):
    """Calculates Mass Eruption Rate (MER) using Mastin et al. 2009 equation.
    Input is height (km), dense rock equivalent (DRE) (kg/m^3). 
    dV = volumetric flow rate (m^3/s)
    Output is MER (kg/s). Optional input is DRE"""
    p = 1/0.241
    dV = (H/2.00)**p
    MER = dV*DRE
    return MER


def sparksMER(H, DRE=2500):
    """Calculates Mass Eruption Rate (MER) using 
    Sparks equation from Volcanic Plumes (1997).
    Input is height (km). dense rock equivalent (DRE) (kg/m^3).
    dV = volumetric flow rate (m^3/s)
    Output is MER (kg/s). Optional input is DRE"""
    p = 1/0.259
    dV = (H/1.67)**p
    MER = dV*DRE
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


def HT2unit(HT, verbose=True):
    """Using Mastin equation to calculate Mass Eruption Rate from plume height (km)
    Inputs: Height (HT) float (in kilometers)
    Output: unit mass in grams. 
    Assume model output is one unit mass per hour. """
    MER = mastinMER(HT)
    unit_mass = MER2unit(MER)
    if verbose:
        print('HEIGHT %0.1f km,  MER %0.3e kg/s , M63 %0.2f , unit mass=%0.3e g/hr.' % (HT, MER, M63, unit_mass))
    return unit_mass


def mastinEVEM(HT, DRE=2500):
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

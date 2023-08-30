# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

from math import *

import matplotlib.pyplot as plt
import numpy as np

"""
# PROGRAM:  fallspd.py   
#   PRGMMR: Alice Crawford   ORG: R/ARL       DATE:
#
# ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
#            Used to calculate terminal fall velocity, tvf, as a function of particle size.

functions for computing Stokes , Ganser and Wilson and Huang Formulations for computing tfv.

# Functions
# constants()
# wilson()
# 

"""

##02/18/2016 This file was used to calculate terminal fall velocity, tvf,  as a function of particle size
##functions for computing Stokes , Ganser and Wilson and Huang Formulations for computing tfv.
##

def constants():
    mu = 1.789e-5   #kg/m s dynamic viscosity
    grav = 9.801   #m/s2   acceleration due to gravity
    dstp = 1.2   #kg/m3  stp density
    frep = 6.53e-8   #mean free path (m at stp)
    return mu, grav, dstp, frep

def wilson(tfv, rhop, rhoa, dp, phi, verbose=True):
    """follow Wilson and Huang (1979) The influence of shape on the atmospheric settling velocity of volcanic ash particles, 
    Earth Planet. Sci. Lett., 44, 311-324.
    tfv is terminal fall velocity.
    rhop is particle density. (kg/m3)
    rhoa is air density  (kg/me)
    dp is particle diameter (m)
    phi is shape factor should be
        ##phi = (b+c) / 2a. where a = longest principle axis, b - intermediate, c=shortest
     TO DO  make sure shape parameter defined appropriately""" 

    mu, grav, dstp, frep = constants()
    condition=True
    mindiff = 0.01
    while condition: 
        Re = rhoa * tfv * dp  /  mu 
        ##F = (b+c) / 2a. where a = longest principle axis, b - intermediate, c=shortest
        F=phi
        CD = 24/Re * F**(-0.828) + 2 * sqrt(1.07 - F)
        tfv2 = sqrt( (4 * grav * dp * (rhop - rhoa))/ (3 *CD * rhoa))
        if abs(tfv2 - tfv) > mindiff:
           tfv = tfv2
        else:
           condition = False
    return tfv

def stokes(rhop, rhoa,  dp, Kshape, slip=False):
    """Stokes equation. applicable to laminar region of standard drag curve for spherical objects.
    For Reynolds number < 0.05 to 1). Corresponds to particles < 20 um.
    tfv is terminal fall velocity.
    rhop is particle density (kg/m3).
    rhoa is air density (kg/m3)
    dp is particle diameter (m)
    Kshape is shape factor
    if slip==True will implement Cunningham slip correction"""

    mu, grav, dstp, frep= constants()
    #TO DO add Cunningham slip correction?
    tfv = ((rhop - rhoa)*grav * dp**2) / (18 * mu )
    if slip:
       frea = frep * (dstp/rhoa)
       sc = 1.0 + (2.0 * frea / dp) * (1.25 + 0.4 *np.exp(-0.55 * dp / frea))
       tfv = tfv * sc  / Kshape
    return tfv


def ganser(tfv, rhop, rhoa, dp, phi, verbose=True):
    """follow Ganser (1993) A Rational approach to drag prediction for spherical and nonspherical particles, 
    Powder Technol., 77, 143-152.
    tfv is first guess from stokes terminal fall velocity.
    rhop is particle density. (kg/m3)
    rhoa is air density (kg/m3)
    dp is particle diameter (m)
    phi is shape factor should be  
    particles sphericity - ratio of the surface area of a sphere with the same volume as a particle to the actual surface area of the particle.
     """ 
    mu, grav, dstp, frep= constants()
    K1shape = 3 / (1+2 / sqrt(phi))                   #K1 is 1 for phi=1, 0.78 for phi=0.5 and 0.55 for phi=0.2
    K2shape = 10**(1.8148*((-log10(phi))**0.5742))      #K2 is 1 for phi=1, 8.1 for phi=0.5 and 30 for phi = 0.2
    condition=True
    mindiff = 0.1
    zzz =0
    while condition: 
        Re = rhoa * tfv * dp  /  mu 
    
        CD = 24/(Re * K1shape) * (1 + 0.1118*(Re*K1shape *K2shape)**0.6567) + 0.4305*K2shape / (1+ 3305/Re/K1shape/K2shape)

        tfv2 = sqrt( (4 * grav * dp * (rhop - rhoa))/ (3 *CD * rhoa))
        print('xxx', zzz, tfv2, tfv, (tfv2-tfv)/tfv*100)
        #if verbose:
           #print "Re" , Re
           #print "CD" , CD
           #print "tfv tfv2" , tfv , tfv2
        if np.abs((tfv2 - tfv)/tfv*100) < mindiff:
           condition=False
        elif zzz > 11:
           condition=False
        else:
        #if zzz <  11:
        #if abs(tfv2 - tfv) > mindiff:
           tfv = tfv2
        #else:
        #   condition = False
        zzz+=1
    return tfv2 
   
def create_curves(dpra, rhop=2500, rhoa=1.2, kshape=1 ):
    """
       rhop : float : density of particle in kg/m3.
              2500 is operational for ash at the NOAA VAACs.
       rhoa : float : density of air in kg/m3
              default value is 1.2 kg/m3 
       kshape : float : shape factor.
              default value of 1
       dpra : list of particle diameters in meters.
    """
    #m2f = 3.28
    #s2min = 1.0/60.0


    #rhop = 2500 # kg/m3 density of particle
    #rhoa = 1.2   # kg/m3 density of air
    Kshape    = 1     #shape factor

    phi =  1    #particles sphericity - ratio of the surface area of a sphere with the same volume as a particle to the actual surface area of the particle.

    #dpra   = [1e-7, 1e-6, 5e-6, 10e-6, 20e-6,  40e-6 ] #list of particle diameters in meters.
    stokes_ra = []
    ganser_ra = []
    wilson_ra = []


    for dp in dpra:
        tfv = stokes(rhop, rhoa, dp, Kshape )
        tfv2 = stokes(rhop, rhoa, dp, Kshape,slip= True)
    ##Ganser equation
        stokes_ra.append(tfv2) 
        ganser_ra.append(ganser(tfv, rhop, rhoa, dp, phi, verbose=False))
        wilson_ra.append(wilson(tfv, rhop, rhoa, dp, phi, verbose=False))

    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ylabel = 'Settling Velocity (m/s)'

    xvar = np.log10(np.array(dpra)*1e6)
    ax.plot(xvar,  np.log10(np.array(wilson_ra)), '-co', label='Wilson 2500 kg/m$^3$')
    ax.plot(xvar, np.log10(np.array(ganser_ra)), '-ko', label='Ganser 2500 kg/m$^3$')
    ax.plot(xvar, np.log10(np.array(stokes_ra)), '-bo', label='Stokes 2500 kg/m$3$')
    ax.set_xlabel('log ( Particle size ($\mu$m))')
    ax.set_ylabel('log (' + ylabel + ')')
    ax.yaxis.grid()
    ax.xaxis.grid()

    ax2.plot(np.array(dpra)*1e6, np.array(wilson_ra), '-co')
    ax2.plot(np.array(dpra)*1e6, np.array(ganser_ra), '-ko')
    ax2.plot(np.array(dpra)*1e6, np.array(stokes_ra), '-bo')
    ax2.set_xlabel('Particle size ($\mu$m)')
    ax2.set_ylabel(ylabel)

    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles, labels, loc=0, fontsize=10)

    ax2.yaxis.grid()
    ax2.xaxis.grid()

    plt.show()
    return stokes_ra, ganser_ra, wilson_ra

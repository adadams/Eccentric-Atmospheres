###################################################################
#   HAT-P-2b DATA FILE
###################################################################

import numpy as N
import astropy.constants as C
import astropy.units as U

#The published properties of the planetary system.
system_properties = {
    'name': 'HATP2b',
    #P
    'orbital period': 5.6334729 * U.d,
    #a
    'semimajor axis': 0.06878 * U.AU,
    #e
    'eccentricity': 0.50910,
    #w
    'argument of periastron': 188.09 * U.deg,
    #mp
    'planet mass': 9.09 * C.M_jup,
    #rp
    'planet radius': 1.157 * C.R_jup,
    #M
    'stellar mass': 1.36 * C.M_sun,
    #R
    'stellar radius': 1.64 * C.R_sun,
    #Teff
    'stellar temperature': 6290 * U.K
         }

#The directory for the data.
data = {'3p6': N.genfromtxt('data/planet/HATP2b/lewis_3p6.txt', comments='#', names=True),
        '4p5': N.genfromtxt('data/planet/HATP2b/lewis_4p5.txt', comments='#', names=True),
        '8p0': N.genfromtxt('data/planet/HATP2b/lewis_8p0.txt', comments='#', names=True)}

###################################################################
#   HD 80606 b DATA FILE
###################################################################

import numpy as N
import astropy.constants as C
import astropy.units as U

#The published properties of the planetary system.
system_properties = {
    'name': 'HD80606b',
    #P
    'orbital period': 111.4367 * U.d,
    #a
    'semimajor axis': 0.4564 * U.AU,
    #e
    'eccentricity': 0.93366,
    #w
    'argument of periastron': 300.8 * U.deg,
    #mp
    'planet mass': 3.94 * C.M_jup,
    #rp
    'planet radius': 0.98 * C.R_jup,
    #M
    'stellar mass': 0.97 * C.M_sun,
    #R
    'stellar radius': 0.978 * C.R_sun,
    #Teff
    'stellar temperature': 5645 * U.K
         }

#The directory for the data.
data = {'4p5': N.genfromtxt('data/planet/HD80606b/wit_4p5_alt.txt', comments='#', names=True),
        '8p0': N.genfromtxt('data/planet/HD80606b/wit_8p0_alt.txt', comments='#', names=True)}

for band in data:
    data[band]['t'] = data[band]['t'] / 24.
    data[band]['flux'] = data[band]['flux'] / 1.e6 + 1

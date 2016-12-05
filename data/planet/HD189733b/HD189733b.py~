###################################################################
#   HD 209458b DATA FILE
###################################################################

import numpy as N
import astropy.constants as C
import astropy.units as U

#The published properties of the planetary system.
system_properties = {
    'name': 'HD209458b',
    #P
    'orbital period': 3.52474859 * U.d,
    #a
    'semimajor axis': 0.04747 * U.AU,
    #e
    'eccentricity': 0.,
    #w
    'argument of periastron': 90 * U.deg,
    #mp
    'planet mass': 3.94 * C.M_jup,
    #rp
    'planet radius': 1.320 * C.R_jup,
    #M
    'stellar mass': 0.97 * C.M_sun,
    #R
    'stellar radius': 1.125 * C.R_sun,
    #Teff
    'stellar temperature': 6092 * U.K
         }

#The directory for the data.
data = {'4p5': N.genfromtxt('data/planet/HD209458b/HD209458_4p5.txt', comments='#', names=True)}

for band in data:
    data[band]['t'] = data[band]['t'] % 1
    data[band]['t'] = N.where(data[band]['t'] > 0.5, data[band]['t'] - 1, data[band]['t'])
    data[band]['t'] = data[band]['t'] * system_properties['orbital period']

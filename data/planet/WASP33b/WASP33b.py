###################################################################
#   WASP-33b DATA FILE
###################################################################

import numpy as N
import astropy.constants as C
import astropy.units as U

#The published properties of the planetary system.
system_properties = {
    'name': 'WASP-33 b',
    #P
    'orbital period': 1.2198669 * U.d,
    #a
    'semimajor axis': 0.02555 * U.AU,
    #e
    'eccentricity': 0.,
    #w
    'argument of periastron': 90. * U.deg,
    #mp
    'planet mass': 3.28 * C.M_jup,
    #rp
    'planet radius': 1.679 * C.R_jup,
    #M
    'stellar mass': 1.495 * C.M_sun,
    #R
    'stellar radius': 1.444 * C.R_sun,
    #Teff
    'stellar temperature': 7430 * U.K
         }

#The directory for the data.
data = {'3p6': N.genfromtxt('data/planet/WASP33b/wasp33b_3p6.csv', comments='#', delimiter=',', names=True, dtype=(float, float, '|S1')),
        '4p5': N.genfromtxt('data/planet/WASP33b/wasp33b_4p5.csv', comments='#', delimiter=',', names=True, dtype=(float, float, '|S1'))}

e = system_properties['eccentricity']
df = 90*U.deg - system_properties['argument of periastron']
P = system_properties['orbital period']

#The data should span -0.5 to 0.5 an orbital period, with t=0 defined as periastron passage. For circular orbits, we arbitrarily set the argument of periastron such that transit and periastron are simultaneous.
for band in data:
    data[band]['t'] = data[band]['t'] % 1
    data[band]['t'] = N.where(data[band]['t'] > 0.5, data[band]['t'] - 1, data[band]['t'])
    data[band]['t'] = data[band]['t'] * system_properties['orbital period']

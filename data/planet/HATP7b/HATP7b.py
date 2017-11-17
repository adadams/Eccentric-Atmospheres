###################################################################
#   HAT-P-7 b DATA FILE
###################################################################

import numpy as N
import astropy.constants as C
import astropy.units as U

#The published properties of the planetary system.
system_properties = {
    'name': 'HAT-P-7 b',
    #P
    'orbital period': 2.2047372 * U.d,
    #a
    'semimajor axis': 0.03676 * U.AU,
    #e
    'eccentricity': 0.016,
    #w
    'argument of periastron': 165. * U.deg,
    #mp
    'planet mass': 1.682 * C.M_jup,
    #rp
    'planet radius': 1.491 * C.R_jup,
    #M
    'stellar mass': 1.47 * C.M_sun,
    #R
    'stellar radius': 1.84 * C.R_sun,
    #Teff
    'stellar temperature': 6441 * U.K
         }

#The directory for the data.
data = {'3p6': N.genfromtxt('data/planet/HATP7b/hatp7b_3p6.csv', comments='#', delimiter=',', names=True, dtype=(float, float, '|S1')),
        '4p5': N.genfromtxt('data/planet/HATP7b/hatp7b_4p5.csv', comments='#', delimiter=',', names=True, dtype=(float, float, '|S1'))}

e = system_properties['eccentricity']
df = 90*U.deg - system_properties['argument of periastron']
P = system_properties['orbital period']

#We want to set the zero point for the time to be periastron, not transit. This corrects the offset.
for band in data:
    E = 2 * N.arctan(N.sqrt((1-e)/(1+e)) * N.tan(0.5*df))
    data[band]['t'] = (data[band]['t'] + 1/(2*N.pi) * (E/U.rad - e*N.sin(E))) % 1
    data[band]['t'] = N.where(data[band]['t'] > 0.5, data[band]['t'] - 1, data[band]['t']) * P/U.d

#Normalize eclipse depths.
data['3p6']['flux'] = data['3p6']['flux'] / 0.99927
data['4p5']['flux'] = data['4p5']['flux'] / 0.99851

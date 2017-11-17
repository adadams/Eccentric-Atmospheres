###################################################################
#   WASP-14b DATA FILE
###################################################################

import numpy as N
import astropy.constants as C
import astropy.units as U

#The published properties of the planetary system.
system_properties = {
    'name': 'WASP-14 b',
    #P
    'orbital period': 2.24376524 * U.d,
    #a
    'semimajor axis': 0.0371 * U.AU,
    #e
    'eccentricity': 0.0830,
    #w
    'argument of periastron': 252.67 * U.deg,
    #mp
    'planet mass': 7.59 * C.M_jup,
    #rp
    'planet radius': 1.240 * C.R_jup,
    #M
    'stellar mass': 1.211 * C.M_sun,
    #R
    'stellar radius': 1.306 * C.R_sun,
    #Teff
    'stellar temperature': 6475 * U.K
         }

#The directory for the data.
data = {'3p6': N.genfromtxt('data/planet/WASP14b/wasp14b_3p6.csv', comments='#', delimiter=',', names=True, dtype=(float, float, '|S1')),
        '4p5': N.genfromtxt('data/planet/WASP14b/wasp14b_4p5.csv', comments='#', delimiter=',', names=True, dtype=(float, float, '|S1'))}

e = system_properties['eccentricity']
df = 90*U.deg - system_properties['argument of periastron']
P = system_properties['orbital period']

#We want to set the zero point for the time to be periastron, not transit. This corrects the offset.
for band in data:
    data[band]['t'] = (data[band]['t'] + 1./(2*N.pi) * N.angle(e + N.cos(df) + 1j*(N.sqrt(1-e**2)*N.sin(df))) - e*N.sqrt(1-e**2)*N.sin(df)/(1+e*N.cos(df))) % 1
    data[band]['t'] = N.where(data[band]['t'] > 0.5, data[band]['t'] - 1, data[band]['t']) * P/U.d

data['3p6']['flux'] = data['3p6']['flux'] / 0.99910
data['4p5']['flux'] = data['4p5']['flux'] / 0.99866

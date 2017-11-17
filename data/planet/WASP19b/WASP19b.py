###################################################################
#   WASP-19b DATA FILE
###################################################################

import numpy as N
import astropy.constants as C
import astropy.units as U

#The published properties of the planetary system.
system_properties = {
    'name': 'WASP-19 b',
    #P
    'orbital period': 0.788838989 * U.d,
    #a
    'semimajor axis': 0.01634 * U.AU,
    #e
    'eccentricity': 0.002,
    #w
    'argument of periastron': 259. * U.deg,
    #mp
    'planet mass': 1.069 * C.M_jup,
    #rp
    'planet radius': 1.392 * C.R_jup,
    #M
    'stellar mass': 0.904 * C.M_sun,
    #R
    'stellar radius': 1.004 * C.R_sun,
    #Teff
    'stellar temperature': 5568 * U.K
         }

#The directory for the data.
data = {'3p6': N.genfromtxt('data/planet/WASP19b/wasp19b_3p6.csv', comments='#', delimiter=',', names=True, dtype=(float, float, '|S1')),
        '4p5': N.genfromtxt('data/planet/WASP19b/wasp19b_4p5.csv', comments='#', delimiter=',', names=True, dtype=(float, float, '|S1'))}

e = system_properties['eccentricity']
df = 90*U.deg - system_properties['argument of periastron']
P = system_properties['orbital period']

#We want to set the zero point for the time to be periastron, not transit. This corrects the offset.
for band in data:
    E = 2 * N.arctan(N.sqrt((1-e)/(1+e)) * N.tan(0.5*df))
    data[band]['t'] = (data[band]['t'] + 1/(2*N.pi) * (E/U.rad - e*N.sin(E))) % 1
    data[band]['t'] = N.where(data[band]['t'] > 0.5, data[band]['t'] - 1, data[band]['t']) * P/U.d

#Normalize eclipse depths.
data['3p6']['flux'] = data['3p6']['flux'] / 0.99711
data['4p5']['flux'] = data['4p5']['flux'] / 0.99657

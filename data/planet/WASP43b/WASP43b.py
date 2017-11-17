###################################################################
#   WASP-43b DATA FILE
###################################################################

import numpy as N
import astropy.constants as C
import astropy.units as U

#The published properties of the planetary system.
system_properties = {
    'name': 'WASP-43 b',
    #P
    'orbital period': 0.81347728 * U.d,
    #a
    'semimajor axis': 0.01526 * U.AU,
    #e
    'eccentricity': 0.,
    #w
    'argument of periastron': 0. * U.deg,
    #mp
    'planet mass': 2.034 * C.M_jup,
    #rp
    'planet radius': 1.036 * C.R_jup,
    #M
    'stellar mass': 0.717 * C.M_sun,
    #R
    'stellar radius': 0.667 * C.R_sun,
    #Teff
    'stellar temperature': 4400 * U.K
         }

#The directory for the data.
data = {'3p6': N.genfromtxt('data/planet/WASP43b/wasp43b_3p6.csv', comments='#', delimiter=',', names=True, dtype=(float, float, '|S1')),
        '4p5': N.genfromtxt('data/planet/WASP43b/wasp43b_4p5.csv', comments='#', delimiter=',', names=True, dtype=(float, float, '|S1'))}

e = system_properties['eccentricity']
df = 90*U.deg - system_properties['argument of periastron']
P = system_properties['orbital period']

#We want to set the zero point for the time to be periastron, not transit. This corrects the offset.
for band in data:
    E = 2 * N.arctan(N.sqrt((1-e)/(1+e)) * N.tan(0.5*df))
    data[band]['t'] = (data[band]['t'] + 1/(2*N.pi) * (E/U.rad - e*N.sin(E))) % 1
    data[band]['t'] = N.where(data[band]['t'] > 0.5, data[band]['t'] - 1, data[band]['t']) * P/U.d

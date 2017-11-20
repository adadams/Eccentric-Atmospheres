###################################################################
#   WASP-12b DATA FILE
###################################################################

import numpy as N
import astropy.constants as C
import astropy.units as U

#The published properties of the planetary system.
system_properties = {
    'name': 'WASP-12 b',
    #P
    'orbital period': 1.09142119 * U.d,
    #a
    'semimajor axis': 0.0234 * U.AU,
    #e
    'eccentricity': 0.0447,
    #w
    'argument of periastron': 272.7 * U.deg,
    #mp
    'planet mass': 1.43 * C.M_jup,
    #rp
    'planet radius': 1.825 * C.R_jup,
    #M
    'stellar mass': 1.280 * C.M_sun,
    #R
    'stellar radius': 1.630 * C.R_sun,
    #Teff
    'stellar temperature': 6360 * U.K
         }

#The directory for the data.
data = {'3p6': N.genfromtxt('data/planet/WASP12b/wasp12b_3p6.csv', comments='#', delimiter=',', names=True, dtype=(float, float, '|S1')),
        '4p5': N.genfromtxt('data/planet/WASP12b/wasp12b_4p5.csv', comments='#', delimiter=',', names=True, dtype=(float, float, '|S1'))}

e = system_properties['eccentricity']
df = 90*U.deg - system_properties['argument of periastron']
P = system_properties['orbital period']

#We want to set the zero point for the time to be periastron, not transit. This corrects the offset.
for band in data:
    E = 2 * N.arctan(N.sqrt((1-e)/(1+e)) * N.tan(0.5*df))
    data[band]['t'] = (data[band]['t'] + P/U.d/(2*N.pi) * (E/U.rad - e*N.sin(E))) % (P/U.d)
    data[band]['t'] = N.where(data[band]['t'] > 0.5*(P/U.d), data[band]['t'] - (P/U.d), data[band]['t'])

###################################################################
#   GJ 436 b DATA FILE
###################################################################

import numpy as N
import astropy.constants as C
import astropy.units as U

#The published properties of the planetary system.
system_properties = {
    'name': 'GJ 436 b',
    #P
    'orbital period': 2.64389803 * U.d,
    #a
    'semimajor axis': 0.0308 * U.AU,
    #e
    'eccentricity': 0.1616,
    #w
    'argument of periastron': 327.2 * U.deg,
    #mp
    'planet mass': 0.080 * C.M_jup,
    #rp
    'planet radius': 0.366 * C.R_jup,
    #M
    'stellar mass': 0.556 * C.M_sun,
    #R
    'stellar radius': 0.455 * C.R_sun,
    #Teff
    'stellar temperature': 3416 * U.K
         }

#The directory for the data.
data = {'8p0': N.genfromtxt('data/planet/GJ436b/gj436b_8p0.csv', comments='#', delimiter=',', names=True, dtype=(float, float, '|S1'))}

e = system_properties['eccentricity']
df = 90*U.deg - system_properties['argument of periastron']
P = system_properties['orbital period']

#We want to set the zero point for the time to be periastron, not transit. This corrects the offset.
for band in data:
    E = 2 * N.arctan(N.sqrt((1-e)/(1+e)) * N.tan(0.5*df))
    data[band]['t'] = (data[band]['t'] + P/U.d/(2*N.pi) * (E/U.rad - e*N.sin(E))) % (P/U.d)
    data[band]['t'] = N.where(data[band]['t'] > 0.5*(P/U.d), data[band]['t'] - (P/U.d), data[band]['t'])

data['8p0']['flux'] = data['8p0']['flux'] / 0.99982

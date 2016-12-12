###################################################################
#   HD 189733b DATA FILE
###################################################################

import numpy as N
import astropy.constants as C
import astropy.units as U

#The published properties of the planetary system.
system_properties = {
    'name': 'HD189733b',
    #P
    'orbital period': 2.218575143 * U.d,
    #a
    'semimajor axis': 0.0312 * U.AU,
    #e
    'eccentricity': 0.,
    #w
    'argument of periastron': 90 * U.deg,
    #mp
    'planet mass': 1.162 * C.M_jup,
    #rp
    'planet radius': 1.216 * C.R_jup,
    #M
    'stellar mass': 0.846 * C.M_sun,
    #R
    'stellar radius': 0.805 * C.R_sun,
    #Teff
    'stellar temperature': 4875 * U.K
         }

#The directory for the data.
data = {'3p6': N.genfromtxt('data/planet/HD189733b/3p6.csv', comments='#', delimiter=',', names=True),
        '4p5': N.genfromtxt('data/planet/HD189733b/4p5.csv', comments='#', delimiter=',', names=True)}

for band in data:
    data[band]['t'] = data[band]['t'] % 1
    data[band]['t'] = N.where(data[band]['t'] > 0.5, data[band]['t'] - 1, data[band]['t'])
    data[band]['t'] = data[band]['t'] * system_properties['orbital period']

data['3p6']['flux'] = data['3p6']['flux'] / 0.999
data['4p5']['flux'] = data['4p5']['flux'] / 0.9987

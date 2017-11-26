###################################################################
#   HD 149026b DATA FILE
###################################################################

import numpy as N
import astropy.constants as C
import astropy.units as U

#The published properties of the planetary system.
system_properties = {
    'name': 'HD 149026 b',
    #P
    'orbital period': 2.8758911 * U.d,
    #a
    'semimajor axis': 0.04288 * U.AU,
    #e
    'eccentricity': 0.,
    #w
    'argument of periastron': 90 * U.deg,
    #mp
    'planet mass': 0.368 * C.M_jup,
    #rp
    'planet radius': 0.813 * C.R_jup,
    #M
    'stellar mass': 1.345 * C.M_sun,
    #R
    'stellar radius': 1.541 * C.R_sun,
    #Teff
    'stellar temperature': 6160 * U.K
         }

#The directory for the data.
data = {'3p6': N.genfromtxt('data/planet/HD149026b/hd149026b_3p6.csv', comments='#', delimiter=',', names=True, dtype=(float, float, '|S1')),
        '4p5': N.genfromtxt('data/planet/HD149026b/hd149026b_4p5.csv', comments='#', delimiter=',', names=True, dtype=(float, float, '|S1'))}

#The data should span -0.5 to 0.5 an orbital period, with t=0 defined as periastron passage. For circular orbits, we arbitrarily set the argument of periastron such that transit and periastron are simultaneous.
for band in data:
    data[band]['t'] = data[band]['t'] % 1
    data[band]['t'] = N.where(data[band]['t'] > 0.5, data[band]['t'] - 1, data[band]['t'])
    data[band]['t'] = data[band]['t'] * system_properties['orbital period']

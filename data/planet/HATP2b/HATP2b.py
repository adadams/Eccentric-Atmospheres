###################################################################
#   HAT-P-2 b DATA FILE
###################################################################

import numpy as N
import astropy.constants as C
import astropy.units as U

#The published properties of the planetary system.
system_properties = {
    'name': 'HAT-P-2 b',
    #P
    'orbital period': 5.6334729 * U.d,
    #a
    'semimajor axis': 0.06878 * U.AU,
    #e
    'eccentricity': 0.5171,
    #w
    'argument of periastron': 185.22 * U.deg,
    #mp
    'planet mass': 9.09 * C.M_jup,
    #rp
    'planet radius': 1.157 * C.R_jup,
    #M
    'stellar mass': 1.36 * C.M_sun,
    #R
    'stellar radius': 1.64 * C.R_sun,
    #Teff
    'stellar temperature': 6290 * U.K
         }

#The directory for the data.
data = {'3p6': N.genfromtxt('data/planet/HATP2b/hatp2b_3p6.csv', comments='#', delimiter=',', names=True, dtype=(float, float, '|S1')),
        '4p5': N.genfromtxt('data/planet/HATP2b/hatp2b_4p5.csv', comments='#', delimiter=',', names=True, dtype=(float, float, '|S1')),
        '8p0': N.genfromtxt('data/planet/HATP2b/hatp2b_8p0.csv', comments='#', delimiter=',', names=True, dtype=(float, float, '|S1'))}

#The data should span -0.5 to 0.5 an orbital period, with t=0 defined as periastron passage. For circular orbits, we arbitrarily set the argument of periastron such that transit and periastron are simultaneous.
for band in data:
    data[band]['t'] = data[band]['t'] / system_properties['orbital period']
    data[band]['t'] = data[band]['t'] % 1
    data[band]['t'] = N.where(data[band]['t'] > 0.5, data[band]['t'] - 1, data[band]['t'])
    data[band]['t'] = data[band]['t'] * system_properties['orbital period']

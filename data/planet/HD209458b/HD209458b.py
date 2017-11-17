###################################################################
#   HD 209458b DATA FILE
###################################################################

import numpy as N
import astropy.constants as C
import astropy.units as U

#The published properties of the planetary system.
system_properties = {
    'name': 'HD 209458 b',
    #P
    'orbital period': 3.52474859 * U.d,
    #a
    'semimajor axis': 0.04747 * U.AU,
    #e
    'eccentricity': 0.,
    #w
    'argument of periastron': 90 * U.deg,
    #mp
    'planet mass': 0.714 * C.M_jup,
    #rp
    'planet radius': 1.380 * C.R_jup,
    #M
    'stellar mass': 0.97 * C.M_sun,
    #R
    'stellar radius': 1.125 * C.R_sun,
    #Teff
    'stellar temperature': 6092 * U.K
         }

#The directory for the data.
data = {'4p5': N.genfromtxt('data/planet/HD209458b/hd209458b_4p5.csv', comments='#', delimiter=',', names=True, dtype=(float, float, '|S1'))}

#The data should span -0.5 to 0.5 an orbital period, with t=0 defined as periastron passage. For circular orbits, we arbitrarily set the argument of periastron such that transit and periastron are simultaneous.
for band in data:
    data[band]['t'] = data[band]['t'] % 1
    data[band]['t'] = N.where(data[band]['t'] > 0.5, data[band]['t'] - 1, data[band]['t'])
    data[band]['t'] = data[band]['t'] * system_properties['orbital period']

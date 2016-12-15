###################################################################
#   WASP-19b DATA FILE
###################################################################

import numpy as N
import astropy.constants as C
import astropy.units as U

#The published properties of the planetary system.
system_properties = {
    'name': 'WASP19b',
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
data = {'3p6': {'t': None, 'flux': None}, '4p5': {'t': None, 'flux': None}}
data['3p6']['t'], data['3p6']['flux'], data['4p5']['t'], data['4p5']['flux'] = N.load('data/planet/WASP19b/wasp19_ch1&2.dat')

#Correct offset using transit ephemerides, so that t=0 corresponds to transit midpoint.
offset_JD = 2455750
transit_JD = 2455708.534626

#In this case we bin the data in 5-minute intervals.
bin_step = 5./1440.

for band in data:
    t_unbin = (data[band]['t'] + offset_JD - transit_JD - (system_properties['argument of periastron']/U.deg-90)/360.*system_properties['orbital period']/U.d) % (system_properties['orbital period']/U.d)
    f_unbin = data[band]['flux']

    t_min = N.min(t_unbin)
    t_max = N.max(t_unbin)
    time_bin = N.arange(start=t_min, stop=t_max, step=bin_step)

    t_bin = N.array([t + 0.5*bin_step for t in time_bin[:-1]])
    f_bin = N.array([N.mean(f_unbin[(t_unbin >= t) & (t_unbin <= t + bin_step)]) for t in time_bin[:-1]])
    
    #data[band]['t'] = t_bin
    data[band]['t'] = N.where(t_bin > (system_properties['orbital period']/U.d)/2, t_bin - (system_properties['orbital period']/U.d), t_bin)
    
    data[band]['flux'] = f_bin

#Normalize eclipse depths.
data['3p6']['flux'] = data['3p6']['flux'] / 0.99711
data['4p5']['flux'] = data['4p5']['flux'] / 0.99657

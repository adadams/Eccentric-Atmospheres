#########################################################
##### QUICK RUN SCRIPT FOR ONLY GRID SEARCH ROUTINE #####
#########################################################

# ##### The typical packages to import.

import astropy.constants as C
import astropy.units as U
from astropy.analytic_functions import blackbody_lambda
import numpy as N

import matplotlib
from matplotlib import pyplot as plt

import datetime

# ##### The planet class contains all the methods used for the geometry of the planet-star system.:

from planet_class import Planet

# ##### Import the modules for the thermal model (currently a blackbody model).:

from thermal import blackbody as thermal_model

# ##### Import the system properties and data of the desired planet (contained in a Python dictionary in a separate file).

from data.planet.HD189733b import HD189733b as exoplanet

# ##### Import the instrumental response data (currently for the Spitzer IRAC bands) and the routine to convert surface temperatures to observed planet-star flux ratios in a given band.

from data.bandpass.spitzer_IRAC import spitzer_IRAC as instrument
from data.bandpass.response import light_curve

# ##### Import the likelihood calculation routine.

from stats.gaussian import log_likelihood

# ##### Create a planet instance with the user-supplied system properties.

planet = Planet(exoplanet.system_properties)

print("{0} loaded as planet to model.".format(planet.name))

# ##### Specify the spatial and time resolution for the calculations, including the number of orbits to run.

planet.set_resolution(longitude_resolution = 72,
                      latitude_resolution = 36,
                      time_resolution = 200,
                      num_orbits = 3)

# ##### Specify ranges to be used for the grid in parameter space to be sampled. Make sure that the rotation period is the first parameter. For the blackbody model one also specifies the radiative timescale at 1000 K, the "nightside" (baseline) temperature of the planet, and the global Bond albedo.

prot = N.linspace(0.4, 1.5, num=12) * planet.pseudosynchronous_period()
t1000 = N.linspace(0, 50, num=26)
t1000[0] += 0.001
Tn = N.linspace(500, 2000, num=16) * U.K
albedo = N.linspace(0.0, 1.0, num=20)

parameters = [prot, t1000*U.hr, Tn, albedo]

# ##### A wrapper function using the thermal model and light curve routines to start with a set of specified parameter values and return a light curve. This function will be passed to the statistical likelihood function.

def generate_model(parameters, spectral_array):
    temperature_map = thermal_model.temperatures(planet=planet, parameters=parameters)
    model = light_curve(planet, temperature_map, spectral_array, parameters=parameters, use_tidal_distortion=False)
    return {'temp': temperature_map, 'model': model}

# Likelihood grid search routine
print("Starting likelihood grid search routine.")
for par in parameters:
    print("{0} to {1}, {2} divisions".format(N.min(par), N.max(par), len(par)))

logls = {}
for band in exoplanet.data:
    print("Band: {0} um".format(band.replace('p','.')))
    logls[band] = log_likelihood(planet = planet,
                                 data = exoplanet.data[band],
                                 spectral_array = instrument.bandpass[band],
                                 model_function = generate_model,
                                 parameters = parameters,
                                 opt = 'parallel')

# ##### Option to save the output in a numpy array file.

with open('files/{0}_{1}_grid.npy'.format(datetime.date.today(), planet.name), 'wb') as grid_file:
    N.save(grid_file, logls)

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

from data.planet.HD209458b import HD209458b as exoplanet

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

prot = N.linspace(0.1, 2.5, num=50) * planet.pseudosynchronous_period()
t1000 = N.linspace(2, 100, num=50)
#t1000[0] += 0.001
Tn = N.linspace(1145.93, 2000, num=1) * U.K
albedo = N.linspace(0.3829, 1.0, num=1)

parameters = [prot, t1000*U.hr, Tn, albedo]

#Get a list of the search lengths along each parameter dimension, as well as the axis with the greatest number of elements and the total number of search elements.
search_dims = [len(l) for l in parameters]
max_axis = N.argmax(search_dims)
num_elements = N.prod(search_dims)

#The above information allows us to determine how we should handle the search parallelization (in a broadcasting sense). Running each point in the grid in series is certainly the slowest, but running everything broadcast as one single array can throw memory errors if the number of elements is too high. Here we set an estimate of the maximum number of elements that will feasibly run at once, and break the original array into components of at most this size. We choose to break the array along the axis of greatest length (if multiple are equal, we choose the first of those).
max_blocksize = 100
num_chunks = int(N.ceil(num_elements/max_blocksize))

if num_chunks > 1:
    param_blocksplit = N.array_split(parameters[max_axis], num_chunks)
    run_blocks = []
    for block in param_blocksplit:
        run_blocks.append(parameters[:])
        run_blocks[-1][max_axis] = block

else:
    run_blocks = [parameters]

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
    logls_parts = []
    #Since we are splitting our original array into chunks, the resulting partitioned container will have its first dimension as the number of chunks. The subsequent dimensions match the dimensions of the parameter grid, which have their 1st and 2nd axes swapped via numpy's meshgrid. We swap them back here to preserve the position of the axis along which we originally split the array.
    for block in run_blocks:
        logls_parts.append(N.swapaxes(log_likelihood(planet = planet,
                                                     data = exoplanet.data[band],
                                                     spectral_array = instrument.bandpass[band],
                                                     model_function = generate_model,
                                                     parameters = block,opt = 'parallel')['logl'], 0,1))
    logls[band] = N.concatenate(logls_parts, axis=max_axis+1)

grid_dict = {'logl': logls, 'par': parameters}

# ##### Option to save the output in a numpy array file.

with open('files/{0}_{1}_grid.npy'.format(datetime.date.today(), planet.name), 'wb') as grid_file:
    N.save(grid_file, grid_dict)

##################################################
##### QUICK RUN SCRIPT FOR ONLY MCMC ROUTINE #####
##################################################

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

from data.planet.GJ436b import GJ436b as exoplanet

# ##### Import the instrumental response data (currently for the Spitzer IRAC bands) and the routine to convert surface temperatures to observed planet-star flux ratios in a given band.

from data.bandpass.spitzer_IRAC import spitzer_IRAC as instrument
from data.bandpass.response import light_curve

# ##### Import the likelihood calculation routine.

from stats.gaussian import log_likelihood
from stats.metropolis import MCMC

# ##### Create a planet instance with the user-supplied system properties.

planet = Planet(exoplanet.system_properties)

print('{0} loaded as planet to model.'.format(planet.name))

# ##### Specify the spatial and time resolution for the calculations, including the number of orbits to run.

planet.set_resolution(longitude_resolution = 18,
                      latitude_resolution = 8,
                      time_resolution = 200,
                      num_orbits = 2)

# ##### A wrapper function using the thermal model and light curve routines to start with a set of specified parameter values and return a light curve. This function will be passed to the statistical likelihood function. (Note: if modeling a planet with expected tidal distortion, consider setting "use_tidal_distortion" to True.)

def generate_model(parameters, spectral_array):
    temperature_map = thermal_model.temperatures(planet=planet, parameters=parameters)
    model = light_curve(planet, temperature_map, spectral_array, parameters=parameters, use_tidal_distortion=False)
    return {'temp': temperature_map, 'model': model}

# ### MCMC LIKELIHOOD SEARCH
# ##### First define a (log) prior function that unpacks a set of parameters.

def mcmc_prior(parameters):
    p, r, T, a = parameters
    lnp = N.where(p>10*U.hr/planet.pseudosynchronous_period(), 0, -N.inf)
    lnp = N.where(r>0, lnp, -N.inf)
    lnp = N.where(T>0, lnp, -N.inf)
    lnp = N.where(0<=a, lnp, -N.inf)
    lnp = N.where(a<=1, lnp, -N.inf)
    return lnp


starting_position = {'values': N.array([0.35, 1.3, 589, 0.06]),
                     'units': N.array([planet.pseudosynchronous_period(), U.hr, U.K, 1])}

num_walkers = 1
num_steps = 10
step_size = [0.05, 5, 100, 0.05]

print('Starting MCMC optimization.')
print('{0} walker(s) executing {1} steps per band.'.format(num_walkers, num_steps))
print('Step sizes: {0}'.format(step_size))

simultaneous = False

if simultaneous:
    mcmc_model = MCMC(planet = planet,
                      data = [exoplanet.data[band] for band in exoplanet.data],
                      spectral_array = [instrument.bandpass[band] for band in exoplanet.data],
                      model_function = generate_model)
    mcmc_model.set_prior(mcmc_prior)
    mcmc_model.set_initial_position(starting_position)
    samples = mcmc_model.run_samples(num_walkers=num_walkers, num_steps=num_steps, step_size=step_size)

else:
    mcmc_model = {}
    samples = {}
    uncertainties = {}
    for band in exoplanet.data:
        print('Band: {0} um'.format(band.replace('p', '.')))
        mcmc_model[band] = MCMC(planet = planet,
                                data = [exoplanet.data[band]],
                                spectral_array = [instrument.bandpass[band]],
                                model_function = generate_model)
        mcmc_model[band].set_prior(mcmc_prior)
        mcmc_model[band].set_initial_position(starting_position)   
        samples[band] = mcmc_model[band].run_samples(num_walkers=num_walkers, num_steps=num_steps, step_size=step_size, annealing=True)

        final_position = {'values': samples[band]['pos'][-1][0],
                          'units': starting_position['units']}
        print('Final position for {0} um: {1}'.format(band.replace('p', '.'), final_position['values']*final_position['units']))
        print('Starting MCMC chain for uncertainties.')
        print(final_position)
        mcmc_model[band].set_initial_position(final_position)
        uncertainties[band] = mcmc_model[band].run_samples(num_walkers=num_walkers, num_steps=num_steps, step_size=step_size, annealing=False)

# ##### Save MCMC outputs: the record of walker positions in parameter space, and the record of log likelihoods at those positions.

planet_name = planet.name.replace(' ','')

with open('files/{0}_{1}_mcmc.npy'.format(datetime.date.today(), planet_name), 'wb') as mcmc_file:
    N.save(mcmc_file, samples)

with open('files/{0}_{1}_mcmc_uncertainties.npy'.format(datetime.date.today(), planet_name), 'wb') as mcmc_uncert_file:
    N.save(mcmc_uncert_file, uncertainties)

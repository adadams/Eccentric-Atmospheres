##################################################
##### QUICK RUN SCRIPT TO GENERATE LIGHT CURVE  ##
##### DATA FOR A SERIES OF PARAMETER POSITIONS  ##
##################################################

# ##### Import shared packages.

import datetime
from glob import glob
from importlib import import_module

import astropy.constants as C
import astropy.units as U
import numpy as N


# ##### Import instrumental response data.

from data.bandpass.spitzer_IRAC import spitzer_IRAC as spitzer_IRAC
from data.bandpass.JHKs.TwoMASS import TwoMASS as TwoMASS
from data.bandpass.kepler import kepler as kepler
bandpasses = {**kepler.bandpass, **TwoMASS.bandpass, **spitzer_IRAC.bandpass}


# ***
# ## Planet-specific Modules

# ##### The planet class contains all the methods used for the geometry of the planet-star system.

from planet_class import Planet


# ##### Import the system properties and data of the planets (contained in Python dictionarys in a separate file), into a dictionary of planet classes. For some plots we just want to pick out a specific planet, which we specify here.

paths = [s.split('/')[-1] for s in glob('data/planet/*') if '__' not in s]
planets = {}
for path in paths:
    planets[path] = Planet(import_module('data.planet.{0}.{0}'.format(path)))

specific_name = 'HD209458b'


# ##### Import the module for the thermal model.

from thermal import blackbody as thermal_model


# ##### Import the routine to convert surface temperatures to observed planet-star flux ratios in a given band.

from data.bandpass.response import light_curve


# ##### Import the likelihood calculation routine.

from stats.gaussian import log_likelihood
from stats.metropolis import MCMC


# ***
# ## Best-Fit Model Light Curves from Existing MCMC Results

# ##### A wrapper function using the thermal model and light curve routines to start with a set of specified parameter values and return a light curve. This function will be passed to the statistical likelihood function.

def generate_model(planet, parameters, spectral_array):
    temperature_map = thermal_model.temperatures(planet=planet, parameters=parameters)
    model = light_curve(planet=planet,
                        temp_map=temperature_map,
                        spectral_array=spectral_array,
                        parameters=parameters,
                        use_tidal_distortion=True)
    return {'temp': temperature_map, 'model': model}


# ##### If there are saved MCMC outputs we want to use, import them here. Usually we have one run which converges (hopefully) to the optimum point in parameter space, as well as a run that starts at that optimum and freely explores the region around the optimum to get an estimate on the uncertainty.
# 
# ##### For the MCMC routine we feed the parameters in as a list. The parameter order is:
# 1. Rotation period
# 2. Radiative timescale
# 3. Minimum temperature
# 4. Albedo

mcmc_uncertainty_file = {}
mcmc_uncertainty = {}

for name, planet in sorted(planets.items()):
    try:
        mcmc_uncertainty_file[name] = N.load('files/{0}_mcmc_uncertainties.npy'.format(planet.name.replace(" ", "")), encoding='latin1', fix_imports=True)[()]
    except:
        continue
    
    mcmc_uncertainty[name] = {}

unique_models = {}
for name, planet in planets.items():
    if name == 'HD209458b':
        continue
    planet.set_resolution(longitude_resolution = 72,
                          latitude_resolution = 36,
                          time_resolution = 200,
                          num_orbits = 5)
    unique_models[name] = {}
    for band in planet.data:
        try:
            mcmc_uncertainty[name][band] = mcmc_uncertainty_file[name][band]['pos'].squeeze()
            product = N.product(mcmc_uncertainty[name][band], axis=1)
        except:
            continue
        index, count = N.unique(product, return_index=True, return_counts=True)[1:]
        unique_parameters = mcmc_uncertainty[name][band][index]
        print('{0} at {1}: {2}'.format(name, band, N.shape(unique_parameters)[0]))
        unique_models[name][band] = {'count': count, 'model': []}
        for i, step in enumerate(unique_parameters):
            print('Step {0}'.format(i))
            step = step * [planets[name].pseudosynchronous_period(), U.hr, U.K, U.dimensionless_unscaled]
            step_paramgrid = N.meshgrid(*step)
            unique_models[name][band]['model'].append(generate_model(planet, step_paramgrid, spitzer_IRAC.bandpass[band])['model'])

with open('files/{0}_uncertainty_lightcurves.npy'.format(datetime.date.today()), 'wb') as uncert_file:
    N.save(uncert_file, unique_models)

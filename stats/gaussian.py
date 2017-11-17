###################################################################
#   GAUSSIAN LIKELIHOOD STATISTICAL MODEL
###################################################################

import numpy as N
import astropy.units as U
from scipy import interpolate

###########################################################
# Calculate the log likelihoods given some data. 'logl' gives the full array of likelihoods over parameter space,
# 'index' gives the triplet of indices for the most favorable parameter values based on the likelihood.
########################################################### 
def log_likelihood(planet, data, spectral_array, model_function, parameters, opt='serial', quiet=False):
    #Store the masks for whether the planet is in an eclipse or transit at a given time, as well as the index corresponding to the flux value at eclipse and transit.
    #One can either flag occulted data points by calculating a window in which we expect the planet to transit/eclipse geometrically, or have a special column in the data file itself that manually marks the points in transit/eclipse. I prefer the latter for consistency, and it's not terribly time intensive.
    #flag_eclipse = N.array([planet.calculate_occultation(d*U.d)['eclipse flag'] for d in data['t']])
    flag_eclipse = (data['occultation'] == b'e')
    eclipse_index = N.array(list(range(N.shape(planet.times)[0])))[planet.calculate_occultation(planet.times)['eclipse flag']][0]
    #flag_transit = N.array([planet.calculate_occultation(d*U.d)['transit flag'] for d in data['t']])
    flag_transit = (data['occultation'] == b't')
    transit_index = N.array(list(range(N.shape(planet.times)[0])))[planet.calculate_occultation(planet.times)['transit flag']][0]
    
    #Arbitrary uncertainty in the Gaussian errors. This number is currently chosen to be roughly on order the spread in data for most light curves in units of the star-planet contrast ratio, so that the arguments in the exponent will be of order unity in most cases. In the future perhaps this will no longer be arbitrary.
    data_spread = 5e-4

    #Store the dimensions of the parameter ranges.
    search_dims = N.array([len(prange) for prange in parameters])

    #This mask contains enough of the time data to include the final periastron-periastron orbit simulated. (The mask actually cuts from third-to-final apastron to final apastron, but we'll need some extra on either side for some phase overlap. Note that this requires one to generate at least 2 orbits for the model.
    assert planet.num_orbits >= 2, \
        'Number of orbits generated must be at least 2 for likelihood calculation.'
    time_cut = (planet.times >= (planet.num_orbits-2.5)*planet.P) & (planet.times <= (planet.num_orbits-0.5)*planet.P)

    #The "serial" option uses a for loop to fill the likelihood array.
    if opt == 'serial':
        model_interps = N.zeros(N.append(search_dims, len(data['t'])))
        log_likelihoods = N.zeros(search_dims)

        for i in range(N.prod(search_dims)):
            #A single for loop takes care of looping over multiple possible dimensions, which we leave as arbitrary by using an "unravel" routine to recover the N-D index here.
            idx = N.unravel_index(i, search_dims)

            #Get the parameter values corresponding to that N-D index.
            parameter_sample = [par[j] for par, j in zip(parameters, idx)]
            print('Parameter position: {0}'.format(parameter_sample))

            #The model light curve.
            model = model_function(parameter_sample, spectral_array)['model']

            #Generate an interpolation function from the model light curve. Apply to the data to get interpolated values.
            interp_func = interpolate.interp1d(planet.times[time_cut] - (planet.P*(planet.num_orbits-2)), model[...,time_cut])
            model_interps[idx] = interp_func(data['t'])

            #The interpolation might not have enough time resolution to catch the dips during eclipse and transit. Here we swap out the interpolated values during eclipses and transits with what the values should be.
            model_interps[idx][flag_eclipse] = model_function(parameter_sample, spectral_array)['model'].T[eclipse_index]
            model_interps[idx][flag_transit] = model_function(parameter_sample, spectral_array)['model'].T[transit_index]

            #Sum up the log likelihoods over time.
            residual = 0.5*((data['flux'] - model_interps[idx])/data_spread)**2
            log_likelihoods[idx] = N.einsum('...t->...', residual)

    #The "parallel" option uses Numpy broadcasting to do all calculations at once.
    if opt == 'parallel':
        #We need meshgrids in N-D of each parameter.
        parameter_mesh = N.meshgrid(*parameters)

        #The model light curve.
        model = model_function(parameters, spectral_array)['model']

        #Generate an interpolation function from the model light curve. Apply to the data to get interpolated values.
        interp_func = interpolate.interp1d(planet.times[time_cut] - (planet.P*(planet.num_orbits-2)), model[...,time_cut])
        model_interps = interp_func(data['t'])

        #The interpolation might not have enough time resolution to catch the dips during eclipse and transit. Here we swap out the interpolated values during eclipses and transits with what the values should be.
        model_interps[...,flag_eclipse] = model[...,eclipse_index,N.newaxis]
        model_interps[...,flag_transit] = model[...,transit_index,N.newaxis]

        #Sum up the log likelihoods over time.
        residual = 0.5*((data['flux'] - model_interps)/data_spread)**2
        log_likelihoods = N.einsum('...t->...', residual)

    #Return the index in the likelihood array corresponding to the smallest (most favorable) value.
    bestfit_index = N.unravel_index(N.argmin(log_likelihoods), search_dims)

    if not quiet:
        print("Most favorable parameter values: {0}".format([par[j] for par, j in zip(parameters, bestfit_index)]))
        
    return {'logl': log_likelihoods,
            'res': residual,
            'index': bestfit_index}

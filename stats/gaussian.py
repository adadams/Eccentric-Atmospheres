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
def log_likelihood(planet_object, data, spectral_array, model_function, parameters, opt='serial', quiet=False):
    #Store the masks for whether the planet is in an eclipse or transit at a given time, as well as the index corresponding to the flux value at eclipse and transit.
    flag_eclipse = N.array([planet_object.calculate_occultation(d*U.d)['eclipse flag'] for d in data['t']])
    eclipse_index = N.array(range(N.shape(planet_object.times)[0]))[planet_object.calculate_occultation(planet_object.times)['eclipse flag']][0]
    flag_transit = N.array([planet_object.calculate_occultation(d*U.d)['transit flag'] for d in data['t']])
    transit_index = N.array(range(N.shape(planet_object.times)[0]))[planet_object.calculate_occultation(planet_object.times)['transit flag']][0]
    
    #Arbitrary uncertainty in the Gaussian errors.
    data_spread = 5e-4

    #Store the dimensions of the parameter ranges.
    search_dims = N.array([len(prange) for prange in parameters])

    #This mask contains enough of the time data to include the final periastron-periastron orbit simulated. (The mask actually cuts from third-to-final apastron to final apastron, but we'll need some extra on either side for some phase overlap.
    time_cut = (planet_object.times > (planet_object.num_orbits-2)*planet_object.P) & (planet_object.times < (planet_object.num_orbits)*planet_object.P)

    #The "serial" option uses a for loop to fill the likelihood array.
    if opt == 'serial':
        model_interps = N.zeros(N.append(search_dims, len(data['t'])))
        log_likelihoods = N.zeros(search_dims)

        for i in range(N.prod(search_dims)):
            #A single for loop takes care of looping over multiple possible dimensions, which we leave as arbitrary by using an "unravel" routine to recover the N-D index here.
            idx = N.unravel_index(i, search_dims)

            #Get the parameter values corresponding to that N-D index.
            parameter_sample = [par[j] for par, j in zip(parameters, idx)]
            print parameter_sample

            #The model light curve.
            model = model_function(parameter_sample, spectral_array)['model']

            #Generate an interpolation function from the model light curve. Apply to the data to get interpolated values.
            interp_func = interpolate.interp1d(planet_object.times - (planet_object.P*(planet_object.num_orbits)), model)
            model_interps[idx] = interp_func(data['t'])

            #The interpolation might not have enough time resolution to catch the dips during eclipse and transit. Here we swap out the interpolated values during eclipses and transits with what the values should be.
            model_interps[idx][flag_eclipse] = model_function(parameter_sample, spectral_array)['model'][eclipse_index]
            model_interps[idx][flag_transit] = model_function(parameter_sample, spectral_array)['model'][transit_index]

            #Sum up the log likelihoods over time.
            log_likelihoods[idx] = N.einsum('...t->...', 0.5*((data['flux'] - model_interps[idx])/data_spread)**2)

    #The "parallel" option uses Numpy broadcasting to do all calculations at once.
    if opt == 'parallel':
        #We need meshgrids in N-D of each parameter.
        parameter_mesh = N.meshgrid(*parameters)

        #The model light curve.
        model = model_function(parameter_mesh, spectral_array)['model']

        #Generate an interpolation function from the model light curve. Apply to the data to get interpolated values.
        interp_func = interpolate.interp1d(planet_object.times - (planet_object.P*(planet_object.num_orbits)), model)
        model_interps = interp_func(data['t'])

        #The interpolation might not have enough time resolution to catch the dips during eclipse and transit. Here we swap out the interpolated values during eclipses and transits with what the values should be.
        model_interps[...,flag_eclipse] = model_function(parameter_mesh, spectral_array)['model'][...,eclipse_index,N.newaxis]
        model_interps[...,flag_transit] = model_function(parameter_mesh, spectral_array)['model'][...,transit_index,N.newaxis]

        #Sum up the log likelihoods over time.
        log_likelihoods = N.einsum('...t->...', 0.5*((data['flux'] - model_interps)/data_spread)**2)

    #Return the index in the likelihood array corresponding to the smallest (most favorable) value.
    bestfit_index = N.unravel_index(N.argmin(log_likelihoods), search_dims)

    if not quiet:
        print "Most favorable parameter values: {0}".format([par[j] for par, j in zip(parameters, bestfit_index)])
        
    return {'logl': log_likelihoods,
            'index': bestfit_index}

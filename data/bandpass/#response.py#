###################################################################
#   SPECTRAL RESPONSE
###################################################################

# INPUT: A grid of planet luminosities (time x phi x theta x parameters).
# Gives the observed flux ratios of planet to star given an input set of spectral responses from a telescope filter.

import numpy as N
import astropy.units as U
from astropy.analytic_functions import blackbody_lambda

###########################################################
# Given a grid of luminosities, calculate the ratio of planet to stellar brightness in the provided bandpass. 
###########################################################       
def light_curve(planet_object,
                temp_map,
                spectral_array,
                prot   = N.array([10])*U.h,
                nbins=None):
    
    # The spectral array will have two columns, one for the wavelengths and one for the spectral response in units proportional to energy (i.e. electron count "weighted" by the wavelength).
    num_waves = N.shape(spectral_array)[0]
    if (nbins == None) or (nbins == num_waves):
        wavelengths = N.array(spectral_array['wavelength']) * U.um
        responses = N.array(spectral_array['weighted_spectrum']) / U.eV
    else:
        wavelengths = N.array([N.average(spectral_array['wavelength'][i*num_waves/nbins:(i+1)*num_waves/nbins]) for i in range(nbins)]) * U.um
        responses = N.array([N.average(spectral_array['weighted_spectrum'][i*num_waves/nbins:(i+1)*num_waves/nbins]) for i in range(nbins)]) / U.eV

    #Mask for the times where the planet is in transit.
    transit_flag = planet_object.calculate_occultation(planet_object.times)['transit flag']

    #Calculate the planet-star cross-sectional area ratio.
    area_ratio = (planet_object.rp/planet_object.R)**2
    
    #For each wavelength bin in the instrumental response curve, calculate the corresponding ratio of planet-star flux. FUN TIP: some Numpy methods like "einsum", while useful, will REMOVE Astropy units from your arrays, so be careful to preserve them using the .unit method to re-unit your data.
    for k, (l, E) in enumerate(zip(wavelengths, responses)):
        
        specific_luminosity = l * observable_luminosity(planet_object, temp_map, l, prot)['integrated'] * E
        if k == 0:
            planet_luminosity = N.einsum('t...->...t', specific_luminosity)
        else:
            planet_luminosity += N.einsum('t...->...t', specific_luminosity)
            
        specific_stellar_luminosity = l * blackbody_lambda(l, planet_object.Teff) * N.pi*U.sr*planet_object.R**2 * E
        if k == 0:
            stellar_luminosity = specific_stellar_luminosity
        else:
            stellar_luminosity += specific_stellar_luminosity

        observed_stellar_luminosity = stellar_luminosity * (1 - area_ratio*transit_flag)
            
        bandpass_ratio = (planet_luminosity*specific_luminosity.unit+observed_stellar_luminosity) / stellar_luminosity.to(specific_luminosity.unit)
        return bandpass_ratio

###########################################################
# Create a visibility mask for the observer, to be multiplied by a flux map to get the observable flux.
########################################################### 
def visibility(planet_object,
               prot   = N.array([10])*U.h):

    #We need the longitudes that face the observer to determine what hemisphere of the globe is visible at a given time step.
    #Additionally, calculate trig functions ahead of time so they won't be computationally expensive in the loops.
    #This requires calculating the surface area in each cell.
    cos_theta = N.cos(planet_object.theta_range)
    long_obs = ((1.5*N.pi + (2*N.pi*(planet_object.times/prot[...,N.newaxis]).decompose() - ((planet_object.w)/U.rad).decompose())) % (2*N.pi))*U.rad
    eclipse_flag = planet_object.calculate_occultation(planet_object.times)['eclipse flag']
    
    #Create a mask that incorporates both the cell area and the angle from the observer. Any grid cell more than 90 deg in longitude from the observer-facing longitude will not be visible.
    #This loop calculates the separation angle between a cell and subsolar longitude at any time, computes the cosine, and creates a mask that returns True and False depending on whether the cell would be illuminated at all.
    cos_sep_angle = N.cos(long_obs[..., N.newaxis]-planet_object.phi_range)
    cell_visible = cos_sep_angle > 0
    horizontal_visibility = N.einsum('...tv->...vt', cos_sep_angle*cell_visible) * ~eclipse_flag
    
    return N.einsum('...vt,u->tuv...', horizontal_visibility, cos_theta)

###########################################################
# Calculate an observable luminosity map and integrated observable luminosity at a given wavelength, using a Planck spectrum.
# An option is available to return the values in ratios of the stellar observable luminosity at the same wavelength.
# The stellar spectrum currently uses a simple blackbody model.
###########################################################
def observable_luminosity(planet_object,
                          temp_map,
                          wavelength,
                          prot   = N.array([10])*U.h,
                          contrast_ratio=False):
    
    #Calculates an array of specific intensities. This is the Planck function.
    intensity = blackbody_lambda(wavelength, temp_map)
    
    #Since we're calling the phi values directly they need to be in assumed radian values. Theta only goes into trig functions so it just needs some angle unit.
    phi_step = ((planet_object.phi_range[1] - planet_object.phi_range[0])/U.rad).decompose()
    theta_step = planet_object.theta_range[1] - planet_object.theta_range[0]
    
    #Function to calculate the cell areas as a function of phi and theta, given the resolution specified.
    cell_area = N.array([phi_step * (N.sin(planet_object.theta_range+theta_step) - N.sin(planet_object.theta_range))]*(planet_object.longitude_resolution+1)).T * planet_object.rp**2
    
    #Calculates the luminosity per cell at a given time and desired wavelength, and then an integrated luminosity over the visible area at the given time.
    luminosity = N.einsum('tuv...->...tuv', intensity * visibility(planet_object, prot)) * intensity.unit * cell_area
    integrated_luminosity = N.nansum(N.einsum('...tuv->tuv...', luminosity), axis=(1,2)) * luminosity.unit * U.sr
    
    #contrast_ratio sets the luminosities in ratios of the stellar luminosity at the desired wavelength.
    if contrast_ratio:
        stellar_luminosity = blackbody_lambda(wavelength, planet_object.Teff) * N.pi*U.sr*planet_object.R**2
        return {'map': luminosity/stellar_luminosity, 'integrated': integrated_luminosity/stellar_luminosity}
    else:
        return {'map': luminosity, 'integrated': integrated_luminosity}

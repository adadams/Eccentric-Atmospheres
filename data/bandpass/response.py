###################################################################
#   SPECTRAL RESPONSE
###################################################################

# INPUT: A grid of planet luminosities (time x phi x theta x parameters).
# Gives the observed flux ratios of planet to star given an input set of spectral responses from a telescope filter.

import numpy as N
import astropy.constants as C
import astropy.units as U
from astropy.analytic_functions import blackbody_lambda

###########################################################
# Given a grid of luminosities, calculate the ratio of planet to stellar brightness in the provided bandpass. 
###########################################################
def light_curve(planet,
                temp_map,
                spectral_array,
                parameters,
                use_tidal_distortion=False,
                nbins=None):
    
    prot, trad_EQ, Tn, albedo = parameters
    prot, trad_EQ, Tn, albedo = N.meshgrid(prot, trad_EQ, Tn, albedo)

    if use_tidal_distortion:
        tidal_axes = planet.tidal_deformation_axes()
        ellipticity = N.sqrt(1-(tidal_axes['y']/tidal_axes['x'])**2)
    else:
        ellipticity = 0
    
    # The spectral array will have two columns, one for the wavelengths and one for the spectral response in units proportional to energy (i.e. electron count "weighted" by the wavelength).
    num_waves = N.shape(spectral_array)[0]
    if (nbins == None) or (nbins == num_waves):
        wavelengths = N.array(spectral_array['wavelength']) * U.um
        responses = N.array(spectral_array['weighted_spectrum']) / U.eV
    else:
        wavelengths = N.array([N.average(spectral_array['wavelength'][i*num_waves/nbins:(i+1)*num_waves/nbins]) for i in range(nbins)]) * U.um
        responses = N.array([N.average(spectral_array['weighted_spectrum'][i*num_waves/nbins:(i+1)*num_waves/nbins]) for i in range(nbins)]) / U.eV

    #Mask for the times where the planet is in transit.
    transit_flag = planet.calculate_occultation(planet.times)['transit flag']

    #Calculate the planet-star cross-sectional area ratio.
    area_ratio = (planet.rp/planet.R)**2
    
    #For each wavelength bin in the instrumental response curve, calculate the corresponding ratio of planet-star flux. FUN TIP: some Numpy methods like "einsum", while useful, will REMOVE Astropy units from your arrays, so be careful to preserve them using the .unit method to re-unit your data.
    for k, (l, E) in enumerate(zip(wavelengths, responses)):
        
        specific_luminosity = l * observable_luminosity(planet, temp_map, l, prot, ellipticity)['integrated'] * E
        if k == 0:
            planet_luminosity = N.einsum('t...->...t', specific_luminosity)
        else:
            planet_luminosity += N.einsum('t...->...t', specific_luminosity)
            
        specific_stellar_luminosity = l * blackbody_lambda(l, planet.Teff) * N.pi*U.sr*planet.R**2 * E
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
def visibility(planet,
               prot   = N.array([10])*U.h):

    #We need the longitudes that face the observer to determine what hemisphere of the globe is visible at a given time step.
    #Additionally, calculate trig functions ahead of time so they won't be computationally expensive in the loops.
    #This requires calculating the surface area in each cell.
    #cos_theta = N.cos((planet.theta_range[1:] + planet.theta_range[:-1])/2)
    cos_theta = N.sin(planet.i-(planet.theta_range[1:] + planet.theta_range[:-1])/2)
    long_obs = planet.observer_longitude(planet.times, prot[...,N.newaxis])
    #long_obs = ((1.5*N.pi - (2*N.pi*(planet.times/prot[...,N.newaxis]).decompose() + ((planet.w)/U.rad).decompose())) % (2*N.pi))*U.rad
    eclipse_flag = planet.calculate_occultation(planet.times)['eclipse flag']
    
    #Create a mask that incorporates both the cell area and the angle from the observer. Any grid cell more than 90 deg in longitude from the observer-facing longitude will not be visible.
    #This loop calculates the separation angle between a cell and subsolar longitude at any time, computes the cosine, and creates a mask that returns True and False depending on whether the cell would be illuminated at all.
    cos_sep_angle = N.cos(long_obs[..., N.newaxis]-((planet.phi_range[1:]+planet.phi_range[:-1])/2))
    cell_visible = cos_sep_angle > 0
    horizontal_visibility = N.einsum('...tv->...vt', cos_sep_angle*cell_visible) * ~eclipse_flag
    
    return N.einsum('...vt,u->tuv...', horizontal_visibility, cos_theta)

###########################################################
# Calculate an observable luminosity map and integrated observable luminosity at a given wavelength, using a Planck spectrum.
# An option is available to return the values in ratios of the stellar observable luminosity at the same wavelength.
# The stellar spectrum currently uses a simple blackbody model.
###########################################################
def observable_luminosity(planet,
                          temp_map,
                          wavelength,
                          prot   = N.array([10])*U.h,
                          ellipticity = 0,
                          contrast_ratio=False):
    
    #Calculates an array of specific intensities. This is the Planck function.
    #intensity = blackbody_lambda(wavelength, temp_map)
    intensity = (2 * C.h * C.c**2) / (wavelength**5) / (N.exp((C.h * C.c)/(wavelength * C.k_B * temp_map)) - 1) / U.sr
    
    #Since we're calling the phi values directly they need to be in assumed radian values. Theta only goes into trig functions so it just needs some angle unit.
    theta_start = planet.theta_range[:-1]
    theta_end = planet.theta_range[1:]
    
    #Function to calculate the cell areas as a function of phi and theta, given the resolution specified.
    if ellipticity > 0:
        a2 = planet.rp**2 / (1-ellipticity**2)

        subsolar_longitude = planet.subsolar_longitude(planet.times, prot[...,N.newaxis])

        phi_start = (planet.phi_range[:-1] - subsolar_longitude[...,N.newaxis]).to(U.rad)
        phi_end = (planet.phi_range[1:] - subsolar_longitude[...,N.newaxis]).to(U.rad)

        #The phi array for the elliptical case has dimensions of ((parameters), time, lon), and theta is just an array of latitudes, so we need to add enough extra dimensions to theta to have them broadcast.
        theta_start_expand = theta_start.reshape((1,)*phi_start.ndim + N.shape(theta_start))
        theta_end_expand = theta_end.reshape((1,)*phi_end.ndim + N.shape(theta_end))

        #The ellipticity array has dimensions of the parameter array, so it needs 3 more for time (t), longitude (v), and latitude (u).
        ep = ellipticity.reshape(N.shape(ellipticity) + (1,)*3)

        #The phi array needs just one more dimension, for latitude (u).
        phi_start_expand = phi_start[...,N.newaxis]
        phi_end_expand = phi_end[...,N.newaxis]

        #Integrals for the surface area of a cell at a given lat/lon, based on the prolate geometry.
        def spheroid_integrand(phi, theta):
            cosp2 = N.cos(phi)**2
            sinp2 = N.sin(phi)**2
            cost2 = N.cos(theta)**2
            sint2 = N.sin(theta)**2
            zeta2 = 1 - ep**2
            ep2 = ep**2
            return N.sqrt( cost2*(1 - ep2*(2-sint2) + ep2**2 * (cost2 + sint2**2*sinp2*cosp2)) )

        def integrate_area(integrand, p1, t1, n):
            def integrate_p(integrand, t, p1, n):
                dp = 2*N.pi*U.rad / planet.longitude_resolution
                p_span = N.linspace(dp/(2*n), dp*(1-1/(2*n)), n)
                p = p1[...,N.newaxis,N.newaxis] + p_span
                integrand_p = integrand(p, t)
                return N.sum(integrand_p, axis=-1)*dp/n
            
            dt = N.pi*U.rad / planet.latitude_resolution
            t_span = N.linspace(dt/(2*n), dt*(1-1/(2*n)), n)
            t = (t1[...,N.newaxis] + t_span)[...,N.newaxis]
            integrand_t = integrate_p(integrand, t, p1, n)
            return N.sum(integrand_t, axis=-1)*dt/n
        
        solid_angle = integrate_area(spheroid_integrand, phi_start_expand, theta_start_expand, n=10)
        
        cell_area = a2 * N.einsum('...tvu->tuv...', solid_angle) * U.rad**2
        
    else:
        phi_start = planet.phi_range[:-1]
        phi_end = planet.phi_range[1:]
        phi_term = ((phi_end - phi_start)/U.rad).decompose()

        theta_term = N.sin(theta_end) - N.sin(theta_start)
        
        cell_area = (phi_term[...,N.newaxis] * theta_term).T * U.rad**2 * planet.rp**2

    phi_step = float(planet.phi_range[1]/U.rad - planet.phi_range[0]/U.rad)
    cell_area_spherical = N.array([phi_step * (N.sin(theta_end) - N.sin(theta_start))]*(planet.longitude_resolution)).T * U.rad**2 * planet.rp**2
    
    #Calculates the luminosity per cell at a given time and desired wavelength, and then an integrated luminosity over the visible area at the given time.

    if ellipticity > 0:
        luminosity = N.einsum('tuv...->...tuv', intensity[:,:-1,:-1,...] * visibility(planet, prot) * N.abs(cell_area)) * intensity.unit * cell_area.unit
        luminosity2 = N.einsum('tuv...->...tuv', intensity[:,:-1,:-1,...] * visibility(planet, prot)) * intensity.unit * N.abs(cell_area_spherical)
        
    else:
        luminosity = N.einsum('tuv...->...tuv', intensity[:,:-1,:-1,...] * visibility(planet, prot)) * intensity.unit * cell_area
        
    integrated_luminosity = N.nansum(N.einsum('...tuv->tuv...', luminosity), axis=(1,2)) * luminosity.unit
    
    #contrast_ratio sets the luminosities in ratios of the stellar luminosity at the desired wavelength.
    if contrast_ratio:
        stellar_luminosity = blackbody_lambda(wavelength, planet.Teff) * N.pi*U.sr*planet.R**2
        return {'map': luminosity/stellar_luminosity, 'integrated': integrated_luminosity/stellar_luminosity}
    else:
        return {'map': luminosity, 'integrated': integrated_luminosity}

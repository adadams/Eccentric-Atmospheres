###################################################################
#   BLACKBODY THERMAL MODEL
###################################################################

import numpy as N
import astropy.units as U

###########################################################
# Calculates the fourth power of the equilibrium temperature at a given time and lat/lon, given a rotation period and "nightside" minimum temperature.
########################################################### 
def equilibrium_temperatures4(planet,
                             prot,
                             Tn,
                             albedo):

    #Define the stellar luminosity and the instantaneous star-planet separations.
    l = planet.R**2 * planet.Teff**4
    r = planet.a * (1-planet.e**2)/(1+planet.e*N.cos(planet.anomaly(planet.times)['true']))
    
    #Alpha is the attenuation factor due to the apparent angle of the star to the horizon at a given point on the planet. Set to zero on the night side.
    sub_sol = planet.subsolar_longitude(planet.times, rotation_period = prot[...,N.newaxis])
    alpha = N.cos(planet.thetas) * N.cos(sub_sol[...,N.newaxis,N.newaxis] - planet.phis)
    alpha[alpha < 0] = 0
    visible_flux = (l/r**2) * N.einsum('...tuv->...uvt', alpha)
    Teq4 = (1 - albedo) * N.einsum('...uvt->tuv...', visible_flux)*visible_flux.unit + Tn**4

    return Teq4.to('K4')

###########################################################
# Produces a new 2D array of surface temperatures from an initial temperature array and a time step.
########################################################### 
def temperatures(planet, 
                 parameters = [N.array([5])*U.h, N.array([10])*U.h, N.array([1000])*U.K, N.array([0.2])]):

    prot, trad_EQ, Tn, albedo = parameters
    prot, trad_EQ, Tn, albedo = N.meshgrid(prot, trad_EQ, Tn, albedo)
    
    eq_T4 = equilibrium_temperatures4(planet, prot, Tn, albedo)
    
    temperature_timeseries = []
    temperature_timeseries.append(N.einsum('...uv->uv...', Tn[...,N.newaxis,N.newaxis]/U.K * N.ones_like(planet.thetas/U.deg))*U.K)

    trad = trad_EQ * ((planet.orbital_equilibrium_temperature()*U.K)**4 / eq_T4)**0.75
    timestep = planet.times[1] - planet.times[0]

    for t, time in enumerate(planet.times):
        if t > 0:
            dT = 0.25 * eq_T4[t]**0.25/trad[t] * (1 - temperature_timeseries[t-1]**4/eq_T4[t]) * timestep
            evolved_temperatures = N.where(N.abs(dT)>N.abs(temperature_timeseries[t-1]-(eq_T4[t])**0.25), eq_T4[t]**0.25, temperature_timeseries[t-1]+dT)*U.K
            temperature_timeseries.append(evolved_temperatures)
        
    return N.array(temperature_timeseries)*U.K

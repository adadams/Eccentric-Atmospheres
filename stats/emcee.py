###################################################################
#   MCMC ROUTINE (USING EMCEE PACKAGE)
###################################################################

import emcee as E
import numpy as N
import astropy.units as U
from gaussian import *

###########################################################
#
########################################################### 

class MCMC:

    def __init__(self, planet_object, data, spectral_array, model_function):
        self.planet_object = planet_object
        self.data = data
        self.spectral_array = spectral_array
        self.model_function = model_function

    def set_prior(self, prior_function):
        self.lnprior = prior_function

    def set_initial_position(self, position):
        self.init_pos = position['values']
        self.units = position['units']
        
    def lnlike(self, par_vals):
        pars = [[par]*unit for par, unit in zip(par_vals, self.units)]
        return log_likelihood(self.planet_object, self.data, self.spectral_array, self.model_function, pars)['logl'].reshape(1)[0]

    def lnprob(self, par_vals):
        lp = self.lnprior(par_vals)
        if not N.isfinite(lp):
            return -N.inf
        return lp + self.lnlike(par_vals)

    def run_samples(self, num_walkers, num_steps, step_size):
        ndim = len(self.init_pos)
        pos = [[par*(1+N.random.randn()*step_size) for par in self.init_pos] for i in range(num_walkers)]

        sampler = E.EnsembleSampler(num_walkers, ndim, self.lnprob)
        sampler.run_mcmc(pos, num_steps)

        return sampler

    

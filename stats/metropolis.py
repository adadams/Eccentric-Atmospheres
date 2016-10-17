###################################################################
#   MCMC ROUTINE (USING METROPOLIS-HASTINGS)
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

        return log_likelihood(self.planet_object, self.data, self.spectral_array, self.model_function, pars, opt='serial', quiet=True)['logl']

    def lnprob(self, par_vals):
        lp = self.lnprior(par_vals)
        lp = N.where(N.isfinite(lp), lp, -N.inf)
        return lp + self.lnlike(par_vals)

    def run_samples(self, num_walkers, num_steps, step_size):
        ndim = len(self.init_pos)

        #This array contains the position evolution of all the MCMC walkers in each parameter dimension.
        pos = N.zeros((num_steps, num_walkers, ndim))
        logls = N.zeros((num_steps, num_walkers))

        #Set the initial position of each walker to the same value.
        pos[0] = N.tile(self.init_pos, (num_walkers,1))
        logls[0] = N.array([self.lnprob(c) for c in pos[0]]).reshape(num_walkers)

        #Implement M-H algorithm to begin walking.
        for i in xrange(1, num_steps):
            print i
            #The candidate position.
            cand = N.array([[par+par0*N.random.randn()*step_size for par, par0 in zip(pos[i-1][j], pos[0][j])] for j in range(num_walkers)])
            cand_check = N.array([self.lnprior(c) for c in cand])
            while any(~N.isfinite(entry) for entry in cand_check):
                cand = N.array([[par+par0*N.random.randn()*step_size for par, par0 in zip(pos[i-1][j], pos[0][j])] for j in range(num_walkers)])
                cand_check = N.array([self.lnprior(c) for c in cand])

            #Acceptance ratio.
            cand_prob = N.array([self.lnprob(c) for c in cand])
            curr_prob = N.array([self.lnprob(c) for c in pos[i-1]])
            #ratio = self.lnprob(cand.T) - self.lnprob(pos[i-1].T)
            ratio = (curr_prob - cand_prob).reshape(num_walkers)
            for j in range(num_walkers):
                if N.log(N.random.uniform())<ratio[j]:
                    pos[i][j] = cand[j]
                    logls[i][j] = cand_prob[j]
                else:
                    pos[i][j] = pos[i-1][j]
                    logls[i][j] = curr_prob[j]
            #pos[i] = N.where(N.log(N.random.uniform(num_walkers))<ratio, cand, pos[i-1])

        return {'pos': pos, 'logl': logls}

    

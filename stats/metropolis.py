###################################################################
#   MCMC ROUTINE (USING METROPOLIS-HASTINGS)
###################################################################

import numpy as N
import astropy.units as U
from .gaussian import *

###########################################################
#
########################################################### 

class MCMC:

    def __init__(self, planet, data, spectral_array, model_function):
        self.planet = planet
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

        likelihoods = [log_likelihood(self.planet, d, s, self.model_function, pars, opt='serial', quiet=True) for d, s in zip(self.data, self.spectral_array)]
        self.likelihoods = likelihoods
        logl = N.array([l['logl'] for l in likelihoods])
        return N.sum(logl, axis=0)

    def lnprob(self, par_vals):
        lp = self.lnprior(par_vals)
        lp = N.where(N.isfinite(lp), lp, -N.inf)
        return lp + self.lnlike(par_vals)

    def run_samples(self, num_walkers, num_steps, step_size, annealing=True):
        ndim = len(self.init_pos)

        #This array contains the position evolution of all the MCMC walkers in each parameter dimension.
        pos = N.zeros((num_steps, num_walkers, ndim))
        logls = N.zeros((num_steps, num_walkers))

        #Set the initial position of each walker to the same value.
        pos[0] = N.tile(self.init_pos, (num_walkers,1))
        logls[0] = N.array([self.lnprob(c) for c in pos[0]]).reshape(num_walkers)

        #For simulated annealing we have 2 characteristic temperatures (high/low), where initially we allow many wild proposals, but gradually reduce the temperature such that the acceptances become much more conservative. We start by defining these temperature limits in terms of 2 characteristic acceptance probabilities.
        if annealing:
            P0 = 0.9
            Pf = 0.001
            T0 = -0.03 / N.log(P0)
            Tf = -0.03 / N.log(Pf)
            T = N.array([T0 * (Tf/T0 ** (i/(num_steps-1))) for i in range(num_steps)])

        #Implement M-H algorithm to begin walking.
        for i in range(1, num_steps):
            print('Step: {0}'.format(i))
            #The candidate position.
            cand = N.array([[par+N.random.randn()*step for par, step in zip(pos[i-1][j], step_size)] for j in range(num_walkers)])
            cand_check = N.array([self.lnprior(c) for c in cand])
            while any(~N.isfinite(entry) for entry in cand_check):
                cand = N.array([[par+N.random.randn()*step for par, step in zip(pos[i-1][j], step_size)] for j in range(num_walkers)])
                cand_check = N.array([self.lnprior(c) for c in cand])

            #Acceptance ratio.
            cand_prob = N.array([self.lnprob(c) for c in cand])
            curr_prob = N.array([self.lnprob(c) for c in pos[i-1]])
            if annealing:
                ratio = ((curr_prob - cand_prob)/T[i]).reshape(num_walkers)
            else:
                ratio = (curr_prob - cand_prob).reshape(num_walkers)
            for j in range(num_walkers):
                if N.log(N.random.uniform())<ratio[j]:
                    pos[i][j] = cand[j]
                    logls[i][j] = cand_prob[j]
                else:
                    pos[i][j] = pos[i-1][j]
                    logls[i][j] = curr_prob[j]

        return {'pos': pos, 'logl': logls}

    

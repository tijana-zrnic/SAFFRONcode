import numpy as np

from scipy.stats import norm
from scipy.stats import bernoulli

class rowexp_new_batch:

    def __init__(self, NUMHYP, numdraws, alt_vec, mu0, mu_alt_vec, markov_lag):
        self.numhyp = NUMHYP
        self.alt_vec = alt_vec
        self.mu0 = mu0
        self.mu_vec = mu0*np.ones(NUMHYP) + np.multiply(alt_vec, mu_alt_vec)
        self.pvec = np.zeros(NUMHYP)
        self.numdraws = numdraws
        self.markov_lag = markov_lag

        '''
        Function drawing p-values: Mixture of two Gaussians
        '''
    def gauss_two_mix(self, sigma = 1, markov_lag = 0, rndsd = 0):

        np.random.seed(rndsd)
        Z = np.zeros(self.numhyp)

        # Draw Z values according to lag
        if (markov_lag == 0):
            Z = self.mu_vec + np.random.randn(self.numhyp)*sigma # draw gaussian acc. to hypothesis, if sigma are all same
        else:
            cov_mat = np.ones([markov_lag + 1,markov_lag + 1])*0.2
            np.fill_diagonal(cov_mat,1)
            # compute samples before L + 1
            Z[0:(markov_lag + 1)] = np.random.multivariate_normal(self.mu_vec[0:(markov_lag + 1)], cov_mat)
            # compute samples after L + 1
            cov_vec = cov_mat[markov_lag,0:markov_lag]
            cov_submat = cov_mat[0:markov_lag,0:markov_lag]
            for i in range(markov_lag + 1,self.numhyp):
                mu_i = self.mu_vec[i] + np.dot(np.dot(cov_vec,np.linalg.inv(cov_submat)),(Z[i-markov_lag:i]-self.mu_vec[i-markov_lag:i]))
                var_i = sigma*sigma - np.dot(np.dot(cov_vec,np.linalg.inv(cov_submat)),np.transpose(cov_vec))
                Z[i] = mu_i + np.random.randn(1)*np.sqrt(var_i)


        # Compute p-values and save
        self.pvec = [(1 - norm.cdf(z)) for z in Z]

    def beta_draws(self, rndsd = 0):
        np.random.seed(rndsd)
        self.pvec = [(np.random.beta(0.5,5,1)*self.alt_vec[i]+np.random.uniform(0,1,1)*(1-self.alt_vec[i])) for i in range(self.numhyp)]

# SMC Resampling Frequency Experiments

Aiming to find the effect of resampling frequency on Sequential Monte Carlo (SMC) and its variations w.r.t different resampling schemes. We are now experimenting with Stochastic Vitality Model(Doucet. 2012) and Gaussian Mixture Model(Richardson and Green. 1997) with SMC and Sequential Monte Carlo Samplers (SMCS, Del Moral. 2006). 

## Gaussian Mixture Model

Bayesian estimation of gaussian mixture's parameters.


$$
f = \sum_{i=1}^r \omega_iN(\mu_i,\lambda_i^{-1})
$$


Now using the data from a $r=2$ gaussian mixture with weakly informative priors similar with Richardson and Green. 1997 and update $\omega_i,\mu_i,\lambda_i$ respectively with, 

1. update $\omega_i$ with additive normal random walk on logit scale, in our code $N(0,0.1^2)$
2. update $\mu_i$ with additive normal random walk, in our code $N(0,0.1^2)$
3. update $\lambda_i$ with multiplicative normal random walk, in our code $\log N(0,0.1^2)$

Use simulated annealing to generate the sequence of distribution to sample from.
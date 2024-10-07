# BayesMVP 0.4

* Allows to choose MH or Gibbs sampler for the coefficient's variance in HRR-type models
* Add `swapBeta()`  in `HRR_Chain.cpp`
* (maybe) Change `var_w_proposal=0.5` and `var_w0_proposal=0.5`  in `HRR_Chain.cpp`
* (maybe) Disable `updateProposalVariances()` for updating the MH proposal variance of w, w0, o and pi in `HRR_Chain.cpp`

* Should we swap Z in GA? But should work for one chain without swap!
* Should we swap Z and beta in GA in HRR.cpp?
* Must check again and again the changes, should something wrong, so that always computing huge values!!!
* It might be wrong to pass hyperparameters via xlm2 tool. Try remove it and use a vector to include all hyperparameters with a fixed order to pass into C++.

* Should we consider gamma prior for the precision parameter, which is equivalent to IG for the variance parameter? If using gamma prior for variance, its posterior is not conjugate, how to sample the posterior?

* TODO1 (done): add variance factor c1 to spike-and-slab prior
* TODO2: adapt the HRR code to the probit models
* TODO3: adapt prediction R-function for binary responses
* TODO4: check whether the calculation of elpd needs to be changed or not
* TODO5 (?): make gamma-prior alternative to inverse-gamma for random effects

# BayesMVP 0.3

* Restrict residual variances to be 1 

# BayesMVP 0.2

* Included the latent responses Z and related calculations
* Fixed a bug in function `SUR_Chain::stepWGibbs()` for extracting a submatrix
* Simplified the calculation in function `SUR_Chain::logPGamma`

# BayesMVP 0.1

* First version
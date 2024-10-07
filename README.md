# BayesMVP

<!-- badges: start --

[![CRAN](http://www.r-pkg.org/badges/version/BayesMVP)](https://cran.r-project.org/package=BayesMVP)
[![r-universe](https://zhizuio.r-universe.dev/badges/BayesMVP)](https://zhizuio.r-universe.dev/BayesMVP)
[![DOI](https://img.shields.io/badge/doi-10.32614%2FCRAN.package.BayesMVP-brightgreen)](https://doi.org/10.32614/CRAN.package.BayesMVP)

-- badges: end -->

[![R-CMD-check](https://github.com/zhizuio/BayesMVP/workflows/R-CMD-check/badge.svg)](https://github.com/zhizuio/BayesMVP/actions)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)


This R package is for high-dimensional Bayesian multivariate probit (BayesMVP) models with general variable selection and dense/sparse covariance matrix, including extended methods from [Bottolo et al. (2021)](https://doi.org/10.1111/rssc.12490), [Zhao et al. (2021)](https://doi.org/10.18637/jss.v100.i11) and [Zhao et al. (2024)](https://doi.org/10.1093/jrsssc/qlad102). 

## Installation

Install the latest development version from [GitHub](https://github.com/zhizuio/BayesMVP)

```r
#install.packages("remotes")
remotes::install_github("zhizuio/BayesMVP")
```

## Simulation 1

Load a simulated data set with structured $\Gamma$ and $G$ (see source code in file `BayesMVP/data-raw/simData.R`) and run the iMVP-MRF model. 

```r
rm(list=ls())
library(BayesMVP)
data("simData", package = "BayesMVP")
data("simDataTest", package = "BayesMVP")
attach(simData)

# here hyperparameters might be not optimal
hyperpar <- list(mrf_d=-5, mrf_e=1, a_w=2, b_w=500, a_w0=2, b_w0=2) 
set.seed(28173)
fit<-BayesMVP(Y=Y, X=X, X_0=X0,
              outFilePath="iMVP-MRF1",nIter=50000,
              nChains=3, burnin=10000,
              hyperpar = hyperpar,
              standardize = TRUE,
              covariancePrior = "IG",
              wSampler = "gibbs",
              gammaPrior = "MRF",
              betaPrior = "reGroup",
              mrfG = mrfG,
              output_CPO = TRUE,
              maxThreads = 3)
# show estimates
plotEstimator(fit, estimator = c("beta", "gamma"), 
              name.predictors = "auto", Pmax = 0.5, beta.type = "conditional")
# show MCMC diagnosis
plot(fit, estimator = "logP", type = "diagnostics")

# summarize results
gammas=getEstimator(fit, estimator = "gamma")
(accuracy_b=sum((gammas>0.5) == (simData$Gamma==1))/prod(dim(Gamma)))
(sensitivity_b=sum(gammas>0.5 & simData$Gamma==1)/sum(simData$Gamma==1))
(specificity_b=sum(gammas<0.5 & simData$Gamma==0)/sum(simData$Gamma==0))

betas=getEstimator(fit, estimator = "beta", beta.type = "conditional", Pmax = 0.5)
(beta_rmse <- sqrt( mean((betas[-c(1:4),] - simData$B)^2) ))
z_pred <- cbind(X0, X) %*% betas
(accuracy_y=sum((z_pred>0) == (simData$Y==1))/prod(dim(z_pred)))
(sensitivity_y=sum(z_pred>0 & simData$Y==1)/sum(simData$Y==1))
(specificity_y=sum(z_pred<=0 & simData$Y==0)/sum(simData$Y==0))
z_test <- cbind(simDataTest$X0, simDataTest$X) %*% betas
(accuracy_y_test=sum((z_test>0) == (simDataTest$Y==1))/prod(dim(z_test)))
(sensitivity_y_test=sum(z_test>0 & simDataTest$Y==1)/sum(simDataTest$Y==1))
(specificity_y_test=sum(z_test<=0 & simDataTest$Y==0)/sum(simDataTest$Y==0))
```

Similarly, we can run the iMVP-Ber model.

```r
rm(list=ls())
library(BayesMVP)
data("simData", package = "BayesMVP")
data("simDataTest", package = "BayesMVP")
attach(simData)

# here hyperparameters might be not optimal
hyperpar <- list(a_omega=2, b_omega=1.5, a_w=2, b_w=1.5, a_w0=2, b_w0=1.5)
set.seed(28173)
fit2 <- BayesMVP(Y=Y, X=X, X_0=X0,
                 outFilePath="iMVP-Ber1",nIter=50000,
                 nChains=3, burnin=10000,
                 hyperpar = hyperpar,
                 standardize = TRUE,
                 covariancePrior = "IG",
                 wSampler = "gibbs",
                 gammaPrior = "hierarchical",
                 betaPrior = "reGroup",
                 output_CPO = TRUE,
                 maxThreads = 3)
```

## Simulation 2

Load another simulated data set with random $\Gamma$ and independent responses (see source code in file `BayesMVP/data-raw/simData2.R`) and run the iMVP-Ber model. 

```r
rm(list=ls())
library(BayesMVP)
data("simData2", package = "BayesMVP")
data("simDataTest2", package = "BayesMVP")
attach(simData2)

# here hyperparameters might be not optimal
hyperpar <- list(a_omega=2, b_omega=1.5, a_w=2, b_w=1.5, a_w0=2, b_w0=1.5)
set.seed(28173)
fit2 <- BayesMVP(Y=Y, X=X, X_0=X0,
                 outFilePath="iMVP-Ber2",nIter=50000,
                 nChains=3, burnin=30000,
                 hyperpar = hyperpar,
                 standardize = TRUE,
                 covariancePrior = "IG",
                 wSampler = "gibbs",
                 gammaPrior = "hierarchical",
                 betaPrior = "reGroup",
                 output_CPO = TRUE,
                 maxThreads = 3)
```


# References

> Leonardo Bottolo, Marco Banterle, Sylvia Richardson, Mika Ala-Korpela, Marjo-Riitta Järvelin, Alex Lewin (2021).
> A computationally efficient Bayesian seemingly unrelated regressions model for high-dimensional quantitative trait loci discovery.
> _Journal of the Royal Statistical Society: Series C (Applied Statistics)_, 70(4):886-908. DOI: [10.1111/rssc.12490](https://doi.org/10.1111/rssc.12490).

> Zhi Zhao, Marco Banterle, Leonardo Bottolo, Sylvia Richardson, Alex Lewin, Manuela Zucknick (2021).
> BayesMVP: An R package for high-dimensional multivariate Bayesian variable and covariance selection in linear regression.
> _Journal of Statistical Software_, 100(11):1-32. DOI: [10.18637/jss.v100.i11](https://doi.org/10.18637/jss.v100.i11).

> Zhi Zhao, Marco Banterle, Alex Lewin, Manuela Zucknick (2023).
> Multivariate Bayesian structured variable selection for pharmacogenomic studies.
> _Journal of the Royal Statistical Society: Series C (Applied Statistics)_, 73(2):420-443 qlad102. DOI: [10.1093/jrsssc/qlad102](https://doi.org/10.1093/jrsssc/qlad102).

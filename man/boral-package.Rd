\name{boral-package}
\alias{boral-package}
\docType{package}
\title{Bayesian Ordination and Regression AnaLysis (boral)}
\description{boral is a package offering Bayesian model-based approaches for analyzing multivariate data in ecology. Estimation is performed using Bayesian/Markov Chain Monte Carlo (MCMC) methods via JAGS (Plummer, 2003). Three ``types" of models may be fitted: 1) With covariates and no latent variables, \code{boral} fits independent column GLMs such that the columns of \eqn{y} are assumed to be independent; 2) With no covariates, \code{boral} fits a purely latent variable model (Skrondal and Rabe-Hesketh, 2004) to perform model-based unconstrained ordination (Hui et al., 2014); 3) With covariates and latent variables, \code{boral} fits correlated column GLMs, with latent variables accounting for correlations between the columns of \eqn{y}.}

\details{
\tabular{ll}{
Package: \tab boralv2\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2014-02-18\cr
License: \tab GPL-2\cr
}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\references{
\itemize{
\item Hui et al., (2014). Model-based approaches to unconstrained ordination. Methods in Ecology and Evolution: In press.
\item Plummer, M. (2003). JAGS: A program for analysis of Bayesian graphical models using Gibbs sampling. In Proceedings of the 3rd International Workshop on Distributed Statistical Computing. March (pp. 20-22).
\item Skrondal, A., and Rabe-Hesketh, S. (2004). Generalized latent variable modeling: Multilevel, longitudinal, and structural equation models. CRC Press.
\item Yi W. et al. (2013). mvabund: statistical methods for analysing multivariate abundance data. R package version 3.8.4. } 
}

\examples{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun

## Example 1 - model with 2 latent variables, site effects, 
## 	and no environmental covariates
spider.fit.nb <- boral(y, family = "negative.binomial", num.lv = 2, 
     row.eff = TRUE, n.burnin = 10, n.iteration = 100, 
     n.thin = 1, save.model = TRUE, calc.ics = FALSE)

## Example 2 - model with latent variables, no site effects, 
## 	and environmental covariates
## This is quite similar to the manyglm() function in the mvabund package, except:
## the latent variables are permitted provide a method of accounting for
## 	correlation between species. 
spider.fit.nb2 <- boral(y, X = spider$x, family = "negative.binomial", num.lv = 2, 
     row.eff = FALSE, n.burnin = 10, n.iteration = 100, 
     n.thin = 1, save.model = TRUE, calc.ics = FALSE)
}

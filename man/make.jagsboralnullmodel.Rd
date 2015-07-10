\name{make.jagsboralnullmodel}
\alias{make.jagsboralnullmodel}
\title{Write a text file containing an boral model for use into JAGS}
\description{This function is designed to write boral models with no latent variables (so-called "null" models).}

\usage{make.jagsboralnullmodel(family, num.X = 0, row.eff = "none", n, p, 
	hypparams = c(100,20,100,100), ssvs.index = -1, model.name = NULL)}

\arguments{
\item{family}{Either a single element, or a vector of length equal to the number of columns in \eqn{y}. The former assumes all columns of \eqn{y} come from this distribution. The latter option allows for different distributions for each column of \eqn{y}. Elements can be one of "binomial" (with logit link), "poisson" (with log link), "negative.binomial" (with log link), "normal" (with identity link), "lnormal" for lognormal (with log link), "tweedie" (with log link), "exponential" (with log link), "gamma" (with log link), "beta" (with logit link), "ordinal" which performs cumulative logit regression. For the negative binomial distribution, the variance is parameterized as \eqn{Var(y) = \mu + \phi \mu^2}, where \eqn{\phi} is the column-specific overdispersion parameter. For the normal distribution, the variance is parameterized as \eqn{Var(y) = \phi}, where \eqn{\phi} is the column-specific variance. For the tweedie distribution, the variance is parameterized as \eqn{Var(y) = \phi \mu^p} where \eqn{\phi} is the column-specific dispersion parameter and \eqn{p} is a power parameter common to all columns assumed to be tweedie, with \eqn{1 < p < 2}. For the gamma distribution, the variance is parameterized as \eqn{Var(y) = \mu/\phi} where \eqn{\phi} is the column-specific rate (henceforth referred to also as dispersion parameter). For the beta distribution, the parameterization is in terms of the mean \eqn{\mu} and sample size \eqn{\phi} (henceforth referred to also as dispersion parameter), so that the two shape parameters are given by \eqn{a = \mu\phi} and \eqn{b = (1-\mu)\phi}.

All columns assumed to have ordinal responses are constrained to have the same cutoffs points, with a column-specific intercept to account for differences between the columns (please see \emph{Details} for formulation). 
}

\item{num.X}{Number of columns in the model matrix \code{X}. Defaults to 0, in which case it is assumed not model matrix is included in model.}

  \item{row.eff}{Single element indicating whether row effects are included as fixed effects ("fixed"), random effects ("random") or not included ("none") in the boral model. If random effects, they are drawn from a normal distribution with mean zero and unknown variance. Defaults to "none". } 

\item{n}{The number of rows in the response matrix \eqn{y}.}

  \item{p}{The number of columns in the response matrix \eqn{y}.}

 \item{hypparams}{Vector of 4 hyperparameters used in the set up of prior distributions. The first element is the variance for the normal priors of all (column-specific) intercepts, row effects, cutoff points for ordinal data. It also controls the variance of the normal prior for the mean, and the maximum of the uniform prior for the variance, for the random effects normal distribution if \code{row.eff = "random"}. The second element is the variance for the normal priors of all column-specific coefficients relating to the latent variables (ignored if \code{num.lv = 0}). The third element is the variance for the normal priors of all column-specific coefficients relating to the model matrix \code{X} (ignored if \code{X = NULL}). The fourth element is the maximum of the uniform prior used for the dispersion parameters i.e, for the negative binomial, normal, lognormal, tweedie, gamma, beta (the dispersion parameter here refers to the sample size) distributions. Note the common power parameter in the tweedie distribution is assumed to have uniform prior from 1 to 2.
 
Finally, if all columns of \eqn{y} are assumed to be ordinal responses, a sum-to-zero constraint is imposed on \eqn{\beta_{0j}} for model identifiability.}

 \item{ssvs.index}{Indices to be used for Stochastic Search Variable Selection (SSVS, George and McCulloch, 1993). Either a single element or a vector with length equal to the number of columns in the implied model matrix \code{X}. Each element can take values of -1 (no SSVS is performed on this covariate), 0 (SSVS is performed on individual coefficients for this covariate, or any integer exceeding 1 (SSVS is performed on collectively all coefficients on this covariate/s.) Defaults to -1, in which case no model selection is performed on the fitted model at all. 
 
 This argument is only read if \code{X.eff = TRUE}, and is necessary to establish the prior distributions used for any explanatory variables. Please see the \code{\link{boral}} help file for more information regarding the implementation of SSVS.}

 \item{model.name}{Name of the text file that the JAGS model is written to. Defaults to NULL, in which case the default of "jagsboralmodel.txt" is used.}

 }

\details{
This function is automatically executed inside \code{\link{boral}}, and therefore does not need to be run separately before fitting the boral model. It can however be run independently if one is: 1) interested in what the actual JAGS file for a particular boral model looks like, 2) wanting to modify a basic JAGS model file to construct more complex model e.g., include environmental variables. Note however \code{\link{boral}} currently does not allow the user to manually enter a script to be run. 

When running the main function \code{\link{boral}}, setting \code{save.model = TRUE} which automatically save the JAGS model file as a text file (with name based on the \code{model.name}) in the current working directory.
}

\value{
A text file is created, containing the JAGS model to be called by the boral function for entering into jags. This file is automatically deleted once boral has finished running unless \code{save.model = TRUE}.}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\seealso{
\code{\link{make.jagsboralmodel}} for writing boral model JAGS scripts with one or more latent variables.}

\examples{
\dontrun{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
n <- nrow(y); p <- ncol(y); 

## Create a boral "null" model JAGS script, where distributions alternative 
## between Poisson and negative distributions 
##   across the rows of y...cool!
make.jagsboralnullmodel(family = rep(c("poisson","negative.binomial"),n), 
     num.X = ncol(spider$x), row.eff = "fixed", n = n, p = p)

## Create a boral "null" model JAGS script, where distributions are all negative 
## 	binomial distributions and covariates will be included!
make.jagsboralnullmodel(family = rep("negative.binomial",length=n), 
     num.X = ncol(spider$x), n = n, p = p, model.name = "myawesomeordnullmodel.txt")
}
}
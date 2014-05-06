\name{make.jagsboralnullmodel}
\alias{make.jagsboralnullmodel}
\title{Write a text file containing an boral model for use into JAGS}
\description{This function is designed to write boral models with no latent variables (so-called "null" models).}

\usage{make.jagsboralnullmodel(family, X.eff, site.eff, n, p, hypparams = c(100,100))}

\arguments{
\item{family}{Either a single element, or a vector of length equal to the number of columns in \eqn{y}. The former assumes all columns of \eqn{y} come from this distribution. The latter option allows for different distributions for each column of \eqn{y}. Elements can be one of "binomial" (with logit link), "poisson" (with log link), "negative.binomial" (with log link), "normal" (with identity link), "lnormal" for lognormal (with log link), "tweedie" (with log link), "ordinal" which performs cumulative logit regression. For the negative binomial distribution, the variance is parameterized as \eqn{Var(y) = \mu + \phi \mu^2}, where \eqn{\phi} is the column-specific overdispersion parameter. For the normal distribution, the variance is parameterized as \eqn{Var(y) = \phi}, where \eqn{\phi} is the column-specific variance. For the tweedie distribution, the variance is parameterized as \eqn{Var(y) = \phi \mu^p} where \eqn{\phi} is the column-specific dispersion parameter and \eqn{p} is a power parameter common to all columns assumed to be tweedie, with \eqn{1 < p < 2}. 

All columns assumed to have ordinal responses are constrained to have the same cutoffs points, with a column-specific intercept to account for differences between the columns (please see \emph{Details} for formulation). 
}

  \item{X.eff}{A logical value indicating whether covariates are to be included in the model.}

 \item{site.eff}{A logical value indicating whether to include row effects or not.}    

\item{n}{The number of rows in the response matrix \eqn{y}.}

  \item{p}{The number of columns in the response matrix \eqn{y}.}

 \item{hypparams}{Vector of 2 hyperparameters in the Bayesian estimation. The first hyperparameter is the variance for the normal priors of all coefficients i.e., row effects, column-specific intercepts, column-specific regression coefficients corresponding to all latent variables and covariates, cutoff points for ordinal data. The second hyperparameter is the maximum of the uniform prior used for the dispersion parameters i.e, for the negative binomial, normal, and tweedie distributions. Defaults to 100 for both elements. The common power parameter in the tweedie distributon is assumed to have uniform prior from 1 to 2.
 
If all columns of \eqn{y} are assumed to be ordinal responses, a sum-to-zero constaint is imposed on \eqn{\theta_{0j}} for model identifiability. Defaults to 100 for both hyperparameters.}
}

\details{
This function is automatically executed inside \code{\link{boral}}, and therefore does not need to be run separately before fitting the boral model. It can however be run independently if one is: 1) interested in what the actual JAGS file for a particukar boral model looks like, 2) wanting to modify a basic JAGS model file to construct more complex model e.g., include environmental variables. Note however \code{\link{boral}} currently does not allow the user to manually enter a script to be run. 

When running the main function \code{\link{boral}}, setting \code{save.model = TRUE} which automatically save the JAGS model file as a text file called "jagsboralmodel.txt" in the current working directory.
}

\value{
A text file named "jagsboralmodel.txt" is created, containing the JAGS model to be called by the boral function for entering into jags. This file is automatically deleted once boral has finished running unless \code{save.model = TRUE}.}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\seealso{
\code{\link{make.jagsboralmodel}} for writing boral model JAGS scripts with one or more latent variables.}

\examples{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
n <- nrow(y); p <- ncol(y); 

## Create a boral "null" model JAGS script, where distributions alternative 
## between Poisson and negative distributions 
##   across the rows of y...cool!
make.jagsboralnullmodel(family = rep(c("poisson","negative.binomial"),n), 
     site.eff = TRUE, X.eff = TRUE, n, p)

## Create a boral "null" model JAGS script, where distributions are all negative 
## 	binomial distributions and covariates will be included!
make.jagsboralnullmodel(family = rep("negative.binomial",length=n), 
     site.eff = FALSE, X.eff = TRUE, n, p)
}
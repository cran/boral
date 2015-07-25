\name{make.jagsboralmodel}
\alias{make.jagsboralmodel}
\title{Write a text file containing an boral model for use into JAGS}
\description{This function is designed to write boral models with one or more latent variables.}

\usage{make.jagsboralmodel(family, num.X = 0, num.traits = 0, 
   which.traits = NULL, row.eff = "none", trial.size = 1, n, p, 
   hypparams = c(100,20,100,100), ssvs.index = -1, model.name = NULL)}

\arguments{
\item{family}{Either a single element, or a vector of length equal to the number of columns in \code{y}. The former assumes all columns of \code{y} come from this distribution. The latter option allows for different distributions for each column of \code{y}. Elements can be one of "binomial" (with probit link), "poisson" (with log link), "negative.binomial" (with log link), "normal" (with identity link), "lnormal" for lognormal (with log link), "tweedie" (with log link), "exponential" (with log link), "gamma" (with log link), "beta" (with logit link), "ordinal" (cumulative probit regression). 

For the negative binomial distribution, the variance is parameterized as \eqn{Var(y) = \mu + \mu^2/\phi}, where \eqn{\phi} is the column-specific dispersion parameter (often referred to as size). For the normal distribution, the variance is parameterized as \eqn{Var(y) = \phi}, where \eqn{\phi} is the column-specific variance. For the tweedie distribution, the variance is parameterized as \eqn{Var(y) = \phi \mu^p} where \eqn{\phi} is the column-specific dispersion parameter and \eqn{p} is a power parameter common to all columns assumed to be tweedie, with \eqn{1 < p < 2}. For the gamma distribution, the variance is parameterized as \eqn{Var(y) = \mu/\phi} where \eqn{\phi} is the column-specific rate (henceforth referred to also as dispersion parameter). For the beta distribution, the parameterization is in terms of the mean \eqn{\mu} and sample size \eqn{\phi} (henceforth referred to also as dispersion parameter), so that the two shape parameters are given by \eqn{a = \mu\phi} and \eqn{b = (1-\mu)\phi}.

All columns assumed to have ordinal responses are constrained to have the same cutoffs points, with a column-specific intercept to account for differences between the columns (please see \emph{Details} for formulation). 
}

\item{num.X}{Number of columns in the model matrix \code{X}. Defaults to 0, in which case it is assumed that no covariates are included in the model.}

\item{num.traits}{Number of columns in the model matrix \code{traits}. Defaults to 0, in which case it is assumed no traits are included in model.}

\item{which.traits}{A list of length equal to (number of columns in \code{X} + 1), informing which columns of \code{traits} the column-specific intercepts and each of the column-specific regression coefficients should be regressed against. The first element in the list applies to the column-specific intercept, while the remaining elements apply to the regression coefficients. Each element of \code{which.traits} is a vector indicating which traits are to be used. For example, if \code{which.traits[[2]] = c(2,3)}, then the regression coefficients corresponding to the first column in \code{X} are regressed against the second and third columns of \code{traits}. If \cr \code{which.traits[[2]] = 0}, then the regression coefficients are treated as independent. 

Defaults to \code{NULL}, in conjunction with \code{num.traits = 0}).} 

\item{row.eff}{Single element indicating whether row effects are included as fixed effects ("fixed"), random effects ("random") or not included ("none") in the boral model. If random effects, they are drawn from a normal distribution with mean zero and unknown variance. Defaults to "none". } 

 \item{trial.size}{Either equal to a single element, or a vector of length equal to the number of columns in y. If a single element, then all columns assumed to be binomially distributed will have trial size set to this. If a vector, different trial sizes are allowed in each column of y. The argument is ignored for all columns not assumed to be binomially distributed. Defaults to 1, i.e. Bernoulli distribution.}

 \item{n}{The number of rows in the response matrix \code{y}.}

  \item{p}{The number of columns in the response matrix \code{y}.}

 \item{hypparams}{Vector of four hyperparameters used in the set up of prior distributions. The first element is the variance for the normal priors of all column-specific intercepts, row effects, and cutoff points for ordinal data. It also controls the the maximum of the uniform prior for the variance of the random effects normal distribution if \code{row.eff = "random"}. The second element is the variance for the normal priors of all latent variable coefficients (ignored if \code{num.lv = 0}). The third element is the variance for the normal priors of all column-specific coefficients relating to the model matrix \code{X} (ignored if \code{X = NULL}). The fourth element is the maximum of the uniform prior used for any dispersion parameters. Note the common power parameter in the tweedie distribution is assumed to have uniform prior from 1 to 2. Also, if all columns of \code{y} are assumed to be ordinal responses, a sum-to-zero constraint is imposed on \eqn{\beta_{0j}} for model identifiability.}

 \item{ssvs.index}{Indices to be used for Stochastic Search Variable Selection (SSVS, George and McCulloch, 1993). Either a single element or a vector with length equal to the number of columns in the implied model matrix \code{X}. Each element can take values of -1 (no SSVS is performed on this covariate), 0 (SSVS is performed on individual coefficients for this covariate), or any integer exceeding 1 (SSVS is performed on collectively all coefficients on this covariate/s.) Defaults to -1, in which case no model selection is performed on the fitted model at all. 
 
 This argument is only read if \code{X.eff = TRUE}, and is necessary to establish the prior distributions used for any explanatory variables. Please see the \code{\link{boral}} help file for more information regarding the implementation of SSVS.}

 \item{model.name}{Name of the text file that the JAGS model is written to. Defaults to NULL, in which case the default of "jagsboralmodel.txt" is used.}
 
 }

\details{
This function is automatically executed inside \code{\link{boral}}, and therefore does not need to be run separately before fitting the boral model. It can however be run independently if one is: 1) interested in what the actual JAGS file for a particular boral model looks like, 2) wanting to modify a basic JAGS model file to construct more complex model e.g., include environmental variables. 

Please note that \code{\link{boral}} currently does not allow the user to manually enter a script to be run. 

When running the main function \code{\link{boral}}, setting \code{save.model = TRUE} which automatically save the JAGS model file as a text file (with name based on the \code{model.name}) in the current working directory.
}

\value{
A text file is created, containing the JAGS model to be called by the boral function for entering into jags. This file is automatically deleted once boral has finished running \code{save.model = TRUE}.}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\seealso{
\code{\link{make.jagsboralnullmodel}} for writing boral models JAGS scripts with no latent variables (so-called "null models").}

\examples{
\dontrun{
library(mvtnorm)
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
n <- nrow(y); p <- ncol(y); 

## Create a boral model JAGS script, where distributions alternative 
## between Poisson and negative binomial distributions 
##   across the rows of y...cool!
make.jagsboralmodel(family = rep(c("poisson","negative.binomial"),length=n), 
	row.eff = "fixed", num.X = 0, n = n, p = p)

## Create a boral model JAGS script, where distributions are all negative 
## 	binomial distributions and covariates will be included.
make.jagsboralmodel(family = "negative.binomial", num.X = ncol(spider$x),
	n = n, p = p)

	
## Simulate some ordinal data and create a JAGS model script
true.lv <- rbind(rmvnorm(10,mean=c(-2,-2)),rmvnorm(10,mean=c(2,2)))
## 20 rows (sites) with two latent variables 
true.lv.coefs <- rmvnorm(30,mean = rep(0,3)); 
## 30 columns (species)
true.lv.coefs[nrow(true.lv.coefs),1] <- -sum(true.lv.coefs[-nrow(true.lv.coefs),1])
## Impose a sum-to-zero constraint on the column effects
true.ordinal.cutoffs <- seq(-2,10,length=10-1)
## Cutoffs for proportional odds regression (must be in increasing order)

sim.y <- create.life(true.lv = true.lv, lv.coefs = true.lv.coefs, 
     n.burnin = 10, n.iteration = 100, n.thin = 1,
	family = "ordinal", cutoffs = true.ordinal.cutoffs) 

make.jagsboralmodel(family = "ordinal", num.X = 0, row.eff = FALSE, n=20, p=30,
	model.name = "myawesomeordmodel.txt")


## Have a look at the JAGS model file for a boral model involving traits,
## based on the ants data from mvabund.
library(mvabund)
data(antTraits)

y <- antTraits$abun
X <- as.matrix(antTraits$env)
## Include only traits 1, 2, and 5, plus an intercept
traits <- as.matrix(cbind(1,antTraits$traits[,c(1,2,5)]))
## Please see help file for boral regarding the use of which.traits
which.traits <- vector("list",ncol(X)+1)
for(i in 1:length(which.traits)) which.traits[[i]] <- 1:ncol(traits)

fit.traits <- boral(y, X = X, traits = traits, which.traits = which.traits, 
	family = "negative.binomial", num.lv = 2, do.fit = FALSE, 
	model.name = "anttraits.txt")

}

}

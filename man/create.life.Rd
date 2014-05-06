\name{create.life}
\alias{create.life}
\title{Simulate a Multivariate Response Matrix}
\description{Simulate a multivariate response matrix, given parameters such as but not necessarily all of: family, number of latent variables and related coefficients, an matrix of explanatory variables and related coefficients, row effects, cutoffs for proportional odds regression of ordinal responses.}

\usage{ 
create.life(true.lv = NULL, lv.coefs, X = NULL, X.coefs = NULL, 
     family, site.coefs = NULL, trial.size = NULL, 
     cutoffs = NULL, powerparam = NULL, manual.dim = NULL)
}

\arguments{
  \item{true.lv}{A matrix of true latent variables. With multivariate abundance data in ecology for instance, each row corresponds to the true site ordination coordinates. Defaults to \code{NULL}, in which case no latent variables are included.}

 \item{lv.coefs}{The column-specific intercepts, regression coefficients relating to the latent variables \code{true.lv}, and dispersion parameters.}

  \item{X}{An model matrix of covariates, which can be included as part of the data generation. Defaults to \code{NULL}, in which case no model matrix is used. No intercept column should be included in \code{X}.}  

  \item{X.coefs}{The coefficients relating to the model matrix \code{X}.}

\item{family}{Either a single element, or a vector of length equal to the number of columns in \eqn{y}. The former assumes all columns of \eqn{y} come from this distribution. The latter option allows for different distributions for each column of \eqn{y}. Elements can be one of "binomial" (with logit link), "poisson" (with log link), "negative.binomial" (with log link), "normal" (with identity link), "lnormal" for lognormal (with log link), "tweedie" (with log link), "ordinal" which performs cumulative logit regression. For the negative binomial distribution, the variance is parameterized as \eqn{Var(y) = \mu + \phi \mu^2}, where \eqn{\phi} is the column-specific overdispersion parameter. For the normal distribution, the variance is parameterized as \eqn{Var(y) = \phi}, where \eqn{\phi} is the column-specific variance. For the tweedie distribution, the variance is parameterized as \eqn{Var(y) = \phi \mu^p} where \eqn{\phi} is the column-specific dispersion parameter and \eqn{p} is a power parameter common to all columns assumed to be tweedie, with \eqn{1 < p < 2}. 

All columns assumed to have ordinal responses are constrained to have the same cutoffs points, with a column-specific intercept to account for differences between the columns (please see \emph{Details} for formulation). 
}

\item{site.coefs}{Row effects. Defaults to \code{NULL}, in which case it is assumed there are no row effects.} 

 \item{trial.size}{Either equal to NULL, a single element, or a vector of length equal to the number of columns in \eqn{y}. If a single element, then all columns assumed to be binomially distributed will have trial size set to this. If a vector, different trial sizes are allowed in each column of \eqn{y}.}

  \item{cutoffs}{A vector of common common cutoffs for proportional odds regression when any of \code{family} is ordinal. They should be increasing order. Defaults to \code{NULL}.}

  \item{powerparam}{A common power parameter for tweedie regression when any of \code{family} is tweedie. Defaults to \code{NULL}.}

  \item{manual.dim}{A vector of length 2, containing the number of rows (\eqn{n}) and columns (\eqn{p}) for the multivariate response matrix. This is a "backup" argument only required when \code{create.life} can not determine how many rows or columns the multivariate response matrix should be.}
}

\value{
A multivariate response matrix of dimension \eqn{n} times \eqn{p}.
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\seealso{
\code{\link{boral}} for the default function for fitting a boral (Bayesian Ordination and Regression AnaLysis) model. 
}

\examples{
## Let's simulate a multivariate response matrix of normally distributed data
## (normal data doesn't occur too often in ecology!)
library(mvtnorm)

true.lv <- rbind(rmvnorm(n=15,mean=c(1,2)),rmvnorm(n=15,mean=c(-3,-1))) 
## 30 rows (sites) with two latent variables 
lv.coefs <- cbind(matrix(runif(30*3),30,3),1)
## 30 columns (species)

X <- matrix(rnorm(30*4),30,4) 
## 4 explanatory variables
X.coefs <- matrix(rnorm(30*4),30,4)

sim.y <- create.life(true.lv, lv.coefs, X, X.coefs, family = "normal")

\dontrun{
fit.boral <- boral(sim.y, X = X, family = "normal", num.lv = 2, 
     site.eff = FALSE, save.model = FALSE)

lvsplot(fit.boral)
}

## Simulate a multivariate response matrix of ordinal data

true.lv <- rbind(rmvnorm(10,mean=c(-2,-2)),rmvnorm(10,mean=c(2,2)))
## 20 rows (sites) with two latent variables 
true.lv.coefs <- rmvnorm(30,mean = rep(0,3)); 
## 30 columns (species)
true.lv.coefs[nrow(true.lv.coefs),1] <- -sum(true.lv.coefs[-nrow(true.lv.coefs),1])
## Impose a sum-to-zero constraint on the column effects
true.ordinal.cutoffs <- seq(-2,10,length=10-1)
## Cutoffs for proportional odds regression (must be in increasing order)

sim.y <- create.life(true.lv = true.lv, lv.coefs = true.lv.coefs, 
     family = "ordinal", cutoffs = true.ordinal.cutoffs) 

\dontrun{
fit.boral <- boral(y = sim.y, family = "ordinal", num.lv = 2, site.eff = F, 
     n.thin = 5, save.model = FALSE, calc.ics = TRUE)
}
}
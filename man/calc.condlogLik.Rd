\name{calc.condlogLik}
\alias{calc.condlogLik}
\title{Conditional log-likelihood for an boral model}
\description{Calculates the conditional log-likelihood for a set of parameter estimates from an boral model, whereby the latent variables (and row effects if applicable) are treated as "fixed effects".}

\usage{
calc.condlogLik(y, X = NULL, family, trial.size = 1, lv.coefs, 
	X.coefs = NULL, row.coefs = NULL, lv, cutoffs = NULL, 
     powerparam = NULL)
}

\arguments{
  \item{y}{The response matrix the boral model was fitted to.}

  \item{X}{The model matrix used in the boral model. Defaults to \code{NULL}, in which case it is assumed no model matrix was used.}  

\item{family}{Either a single element, or a vector of length equal to the number of columns in \eqn{y}. The former assumes all columns of \eqn{y} come from this distribution. The latter option allows for different distributions for each column of \eqn{y}. Elements can be one of "binomial" (with logit link), "poisson" (with log link), "negative.binomial" (with log link), "normal" (with identity link), "lnormal" for lognormal (with log link), "tweedie" (with log link), "exponential" (with log link), "gamma" (with log link), "beta" (with logit link), "ordinal" which performs cumulative logit regression. For the negative binomial distribution, the variance is parameterized as \eqn{Var(y) = \mu + \phi \mu^2}, where \eqn{\phi} is the column-specific overdispersion parameter. For the normal distribution, the variance is parameterized as \eqn{Var(y) = \phi}, where \eqn{\phi} is the column-specific variance. For the tweedie distribution, the variance is parameterized as \eqn{Var(y) = \phi \mu^p} where \eqn{\phi} is the column-specific dispersion parameter and \eqn{p} is a power parameter common to all columns assumed to be tweedie, with \eqn{1 < p < 2}. For the gamma distribution, the variance is parameterized as \eqn{Var(y) = \mu/\phi} where \eqn{\phi} is the column-specific rate (henceforth referred to also as dispersion parameter). For the beta distribution, the parameterization is in terms of the mean \eqn{\mu} and sample size \eqn{\phi} (henceforth referred to also as dispersion parameter), so that the two shape parameters are given by \eqn{a = \mu\phi} and \eqn{b = (1-\mu)\phi}.

All columns assumed to have ordinal responses are constrained to have the same cutoffs points, with a column-specific intercept to account for differences between the columns (please see \emph{Details} for formulation). 
}

 \item{trial.size}{ This argument is ignored for all columns not assumed to be binomially distributed. Either equal to a single element, or a vector of length equal to the number of columns in y. If a single element, then all columns assumed to be binomially distributed will have trial size set to this. If a vector, different trial sizes are allowed in each column of y. Defaults to 1.}

  \item{lv.coefs}{The column-specific intercept, coefficient estimates relating to the latent variables, and dispersion parameters from the boral model.}

  \item{X.coefs}{The coefficients estimates relating to the model matrix \code{X} from the boral model. Defaults to \code{NULL}, in which it is assumed there are no covariates in the model.}

  \item{row.coefs}{ Row effect estimates for the boral model. Even if they were estimated as random effects, the conditional log-likelihood is defined conditional on these estimates i.e., they are (also) treated as ``fixed effects". Defaults to NULL, in which case it is assumed there are no row effects in the model.} 

  \item{lv}{Latent variables "estimates" from the boral model, which the conditional log-likelihood is based on. For boral models with no latent variables, please use \code{\link{calc.logLik.lv0}} to calculate the conditional log-likelihood.}

  \item{cutoffs}{Common cutoff estimates from the boral model when any of the columns of \eqn{y} are ordinal responses. Defaults to \code{NULL}.}

\item{powerparam}{Common power parameter from the boral model when any of the columns of \eqn{y} are tweedie responses. Defaults to \code{NULL}.}
}

\details{
For an \eqn{n x p} response matrix \eqn{y}, suppose we fit an boral model with \eqn{q} latent variables. If we denote the latent variables by \eqn{\bm{z}_i; i = 1,\ldots,n}, then the conditional log-likelihood is given by (with parameters where appropriate),

\deqn{ \log(f) = \sum_{i=1}^n \sum_{j=1}^p \log (f(y_{ij} | \alpha_i, \tau_k, \beta_{0j}, \bm{\theta}_j, \bm{z}_i, \bm{x}_i, \bm{\beta}_j, \phi_j)), }

where \eqn{f(y_{ij}|\cdot)} is the assumed distribution for column \eqn{j}, \eqn{\alpha_i} is the row effect, \eqn{\tau_k} are the common cutoffs for proportional odds regression, \eqn{\beta_{0j}} is the column-specific intercept, \eqn{\bm{z}'_i} are the latent variables, \eqn{\bm{\theta}_j} are the column-specific coefficients relating to the latent variables, \eqn{\bm{x}'_i} is row \eqn{i} of the model matrix,, \eqn{\bm{\beta}_j} are the column-specific coefficients relating to the model matrix of explanatory variables, \eqn{\phi_j} are-column-specific dispersion parameters.

The key difference between this and the marginal likelihood (see \code{\link{calc.marglogLik}}) is that the conditional log-likelihood treats the latent variables (and row effects if applicable) as "fixed effects", while the marginal log-likelihood treats them as "random effects" which require integrating over.
}

\value{
A list with the following components:
\item{logLik}{Value of the conditional log-likelihood.}
\item{logLik.comp}{A vector of the log-likelihood values for each row of \eqn{y}, \cr such that \code{sum(logLik.comp) = logLik}.}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\note{
The conditional DIC, WAIC, EAIC, and EBIC returned from \code{\link{get.measures}} are based on the conditional log-likelihood calculated from this function. Additionally, \code{\link{get.measures}} returns the conditional log-likelihood evaluated at all MCMC samples of a fitted boral model.
}

\seealso{
\code{\link{get.measures}} for some information criteria based on the conditional log-likelihood; \code{\link{calc.marglogLik}} for calculation of the marginal log-likelihood; \code{\link{calc.logLik.lv0}} to calculate the conditional/marginal log-likelihood for an boral model with no latent variables.
}

\examples{
\dontrun{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
n <- nrow(y); p <- ncol(y); 

## Example 1 - model with 2 latent variables, site effects, 
## 	and no environmental covariates
spider.fit.nb <- boral(y, family = "negative.binomial", num.lv = 2, 
     n.burnin = 10, n.iteration = 100, n.thin = 1,
	row.eff = "fixed", save.model = TRUE, calc.ics = FALSE)

## Extract all MCMC samples
fit.mcmc <- mcmc(spider.fit.nb$jags.model$BUGSoutput$sims.matrix, 
	start = 1, thin = spider.fit.nb$n.thin) 

## Find the posterior medians
coef.mat <- matrix(apply(fit.mcmc[,grep("all.params",colnames(fit.mcmc))],
     2,median),nrow=p)
site.coef.median <- apply(fit.mcmc[,grep("row.params", colnames(fit.mcmc))],
     2,median)
lvs.mat <- matrix(apply(fit.mcmc[,grep("lvs",colnames(fit.mcmc))],2,median),nrow=n)

## Caculate the conditional log-likelihood at the posterior median
calc.condlogLik(y, family = "negative.binomial", 
     lv.coefs =  coef.mat, row.coefs = site.coef.median, lv = lvs.mat)


## Example 2 - model with one latent variable, no site effects, 
## 	and environmental covariates
spider.fit.nb2 <- boral(y, X = spider$x, family = "negative.binomial", num.lv = 1, 
	n.burnin = 10, n.iteration = 100, n.thin = 1,
     save.model = TRUE, calc.ics = FALSE)

## Extract all MCMC samples
fit.mcmc <- mcmc(spider.fit.nb2$jags.model$BUGSoutput$sims.matrix, 
	start = 1, thin = spider.fit.nb2$n.thin) 

## Find the posterior medians
coef.mat <- matrix(apply(fit.mcmc[,grep("all.params",colnames(fit.mcmc))],
     2,median),nrow=p)
X.coef.mat <- matrix(apply(fit.mcmc[,grep("X.params",colnames(fit.mcmc))],
	2,median),nrow=p)
lvs.mat <- matrix(apply(fit.mcmc[,grep("lvs",colnames(fit.mcmc))],2,median),nrow=n)

## Caculate the log-likelihood at the posterior median
calc.condlogLik(y, X = spider$x, family = "negative.binomial", 
	lv.coefs =  coef.mat, X.coefs = X.coef.mat, lv = lvs.mat)
}
}

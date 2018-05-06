\name{calc.condlogLik}

\alias{calc.condlogLik}

\title{Conditional log-likelihood for a boral model}

\description{Calculates the conditional log-likelihood for a set of parameter estimates from an boral model, where everything is treated as "fixed effects" including latent variables, row effects, and so on.}

\usage{
calc.condlogLik(y, X = NULL, family, trial.size = 1, lv.coefs, 
	X.coefs = NULL, row.coefs = NULL, row.ids = NULL,
	offset = NULL, lv = NULL, cutoffs = NULL, powerparam = NULL)
}

\arguments{
  \item{y}{The response matrix the boral model was fitted to.}

  \item{X}{The model matrix used in the boral model. Defaults to \code{NULL}, in which case it is assumed no model matrix was used.}  

  \item{family}{Either a single element, or a vector of length equal to the number of columns in \code{y}. The former assumes all columns of \code{y} come from this distribution. The latter option allows for different distributions for each column of \code{y}. Elements can be one of "binomial" (with probit link), "poisson" (with log link), "negative.binomial" (with log link), "normal" (with identity link), "lnormal" for lognormal (with log link), "tweedie" (with log link), "exponential" (with log link), "gamma" (with log link), "beta" (with logit link), "ordinal" (cumulative probit regression). 

  Please see \code{\link{about.distributions}} for information on distributions available in boral overall.
  }

 \item{trial.size}{Either equal to a single element, or a vector of length equal to the number of columns in y. If a single element, then all columns assumed to be binomially distributed will have trial size set to this. If a vector, different trial sizes are allowed in each column of y. The argument is ignored for all columns not assumed to be binomially distributed. Defaults to 1, i.e. Bernoulli distribution.}

  \item{lv.coefs}{The column-specific intercept, coefficient estimates relating to the latent variables, and dispersion parameters from the boral model.}

  \item{X.coefs}{The coefficients estimates relating to \code{X} from the boral model. Defaults to \code{NULL}, in which it is assumed there are no covariates in the model.}

  \item{row.coefs}{Row effect estimates for the boral model. The conditional likelihood is defined conditional on these estimates i.e., they are also treated as ``fixed effects". Defaults to \code{NULL}, in which case it is assumed there are no row effects in the model.} 

 \item{row.ids}{A matrix with the number of rows equal to the number of rows in \code{y}, and the number of columns equal to the number of row effects to be included in the model. Element \eqn{(i,j)} indicates to the cluster ID of row \eqn{i} in \code{y} for random effect eqn{j}; please see the \code{\link{boral}} function for details. Defaults to \code{NULL}, so that if \code{row.coefs = NULL} then the argument is ignored, otherwise if \code{row.coefs} is supplied then \code{row.ids = matrix(1:nrow(y), ncol = 1)} i.e., a single, row effect unique to each row. An internal check is done to see \code{row.coefs} and \code{row.ids} are consistent in terms of arguments supplied.}    

 \item{offset}{A matrix with the same dimensions as the response matrix \code{y}, specifying an a-priori known component to be included in the linear predictor during fitting. Defaults to \code{NULL}.}

 \item{lv}{Latent variables "estimates" from the boral model, which the conditional likelihood is based on. Defaults to \code{NULL}, in which case it is assumed no latent variables were included in the model.}

  \item{cutoffs}{Common cutoff estimates from the boral model when any of the columns of \code{y} are ordinal responses. Defaults to \code{NULL}.}

\item{powerparam}{Common power parameter from the boral model when any of the columns of \code{y} are tweedie responses. Defaults to \code{NULL}.}
}


\details{
For an \eqn{n x p} response matrix \code{y}, suppose we fit an boral model with one or more latent variables. If we denote the latent variables by \eqn{\bm{z}_i; i = 1,\ldots,n}, then the conditional log-likelihood is given by,

\deqn{ \log(f) = \sum_{i=1}^n \sum_{j=1}^p \log \{f(y_{ij} | \bm{z}_i, \bm{\theta}_j, \beta_{0j}, \ldots)\}, }

where \eqn{f(y_{ij}|\cdot)} is the assumed distribution for column \eqn{j}, \eqn{\bm{z}_i} are the latent variables and \eqn{\bm{\theta}_j} are the coefficients relating to them, \eqn{\beta_{0j}} are column-specific intercepts, and \eqn{\ldots} denotes anything else included in the model, such as row effects, regression coefficients related \code{X} and \code{traits}, etc...

The key difference between this and the marginal likelihood (see \code{\link{calc.marglogLik}}) is that the conditional likelihood treats everything as "fixed effects" i.e., conditions on them. These include the latent variables \eqn{\bm{z}_i} and other parameters that were included in the model as random effects e.g., row effects if \code{row.eff = "random"}, regression coefficients related to \code{X} if traits were included in the model, and so on.

The conditional DIC, WAIC, EAIC, and EBIC returned from \code{\link{get.measures}} are based on the conditional likelihood calculated from this function. Additionally, \code{\link{get.measures}} returns the conditional likelihood evaluated at all MCMC samples of a fitted boral model.
}


\value{
A list with the following components:
\item{logLik}{Value of the conditional log-likelihood.}
\item{logLik.comp}{A matrix of the log-likelihood values for each element in \code{y}, \cr such that \code{sum(logLik.comp) = logLik}.}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\seealso{
\code{\link{calc.logLik.lv0}} to calculate the conditional/marginal log-likelihood for an boral model with no latent variables; \code{\link{calc.marglogLik}} for calculation of the marginal log-likelihood; 
}

\examples{
\dontrun{
## NOTE: The values below MUST NOT be used in a real application;
## they are only used here to make the examples run quick!!!
example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
     n.thin = 1)

library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
n <- nrow(y)
p <- ncol(y)

## Example 1 - model with 2 latent variables, site effects, 
## 	and no environmental covariates
spiderfit_nb <- boral(y, family = "negative.binomial", num.lv = 2, 
	row.eff = "fixed", save.model = TRUE, mcmc.control = example_mcmc_control)

## Extract all MCMC samples
fit_mcmc <- get.mcmcsamples(spiderfit_nb) 
mcmc_names <- colnames(fit_mcmc)

## Find the posterior medians
coef_mat <- matrix(apply(fit_mcmc[,grep("lv.coefs",mcmc_names)],
     2,median),nrow=p)
site_coef <- list(ID1 = apply(fit_mcmc[,grep("row.coefs.ID1", mcmc_names)],
     2,median))
lvs_mat <- matrix(apply(fit_mcmc[,grep("lvs",mcmc_names)],2,median),nrow=n)

## Calculate the conditional log-likelihood at the posterior median
calc.condlogLik(y, family = "negative.binomial", 
     lv.coefs = coef_mat, row.coefs = site_coef, lv = lvs_mat)


## Example 2 - model with no latent variables and environmental covariates
X <- scale(spider$x)
spiderfit_nb2 <- boral(y, X = X, family = "negative.binomial", num.lv = 0, 
     save.model = TRUE, mcmc.control = example_mcmc_control)

## Extract all MCMC samples
fit_mcmc <- get.mcmcsamples(spiderfit_nb2) 
mcmc_names <- colnames(fit_mcmc)

## Find the posterior medians
coef_mat <- matrix(apply(fit_mcmc[,grep("lv.coefs",mcmc_names)],
     2,median),nrow=p)
X_coef_mat <- matrix(apply(fit_mcmc[,grep("X.coefs",mcmc_names)],
	2,median),nrow=p)

## Calculate the log-likelihood at the posterior median
calc.condlogLik(y, X = X, family = "negative.binomial", 
	lv.coefs =  coef_mat, X.coefs = X_coef_mat)
}
}

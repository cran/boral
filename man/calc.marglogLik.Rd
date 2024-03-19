\name{calc.marglogLik}

\alias{calc.marglogLik}

\docType{package}

\title{Marginal log-likelihood for a fitted model}

\description{
\Sexpr[results=rd, stage=render]{lifecycle::badge("deprecated")}

Calculates the marginal log-likelihood for a set of parameter estimates from a fitted model, whereby the latent variables and random effects (if applicable) are integrated out. The integration is performed using Monte Carlo integration. WARNING: As of version 1.9, this function is no longer being maintained (and probably does not work properly, if at all)!}

\usage{
calc.marglogLik(y, X = NULL, family, trial.size = 1, lv.coefs, 
     X.coefs = NULL, row.eff = "none", row.params = NULL, 
     row.ids = NULL,offset = NULL, num.lv, lv.mc = NULL, 
     cutoffs = NULL, powerparam = NULL)
}

\arguments{
\item{y}{The response matrix that the model was fitted to.}

\item{X}{The covariate matrix used in the model. Defaults to \code{NULL}, in which case it is assumed no model matrix was used.}  

  \item{family}{Either a single element, or a vector of length equal to the number of columns in the response matrix. The former assumes all columns of the response matrix come from this distribution. The latter option allows for different distributions for each column of the response matrix. Elements can be one of "binomial" (with probit link), "poisson" (with log link), "negative.binomial" (with log link), "normal" (with identity link), "lnormal" for lognormal (with log link), "tweedie" (with log link), "exponential" (with log link), "gamma" (with log link), "beta" (with logit link), "ordinal" (cumulative probit regression). 

  Please see \code{\link{about.distributions}} for information on distributions available in boral overall.
  }

 \item{trial.size}{Either equal to a single element, or a vector of length equal to the number of columns in y. If a single element, then all columns assumed to be binomially distributed will have trial size set to this. If a vector, different trial sizes are allowed in each column of y. The argument is ignored for all columns not assumed to be binomially distributed. Defaults to 1, i.e. Bernoulli distribution.}

  \item{lv.coefs}{The response-specific intercept, coefficient estimates relating to the latent variables, and dispersion parameters from the fitted model.}

  \item{X.coefs}{The coefficients estimates relating to the covariate matrix from the fitted model. Defaults to \code{NULL}, in which it is assumed there are no covariates in the model.}

  \item{row.eff}{Single element indicating whether row effects are included as fixed effects ("fixed"), random effects ("random") or not included ("none") in the fitted model. If fixed effects, then for parameter identifiability the first row effect is set to zero, which analogous to acting as a reference level when dummy variables are used. If random effects, they are drawn from a normal distribution with mean zero and standard deviation given by \code{row.params}. Defaults to "none". } 

  \item{row.params}{Parameters corresponding to the row effect from the fitted model. If \cr \code{row.eff = "fixed"}, then these are the fixed effects and should have length equal to the number of columns in the response matrix. If \code{row.eff = "random"}, then this is standard deviation for the random effects normal distribution. If \code{row.eff = "none"}, then this argument is ignored.}

 \item{row.ids}{A matrix with the number of rows equal to the number of rows in the response matrix, and the number of columns equal to the number of row effects to be included in the model. Element \eqn{(i,j)} indicates the cluster ID of row \eqn{i} in the response matrix for random effect \eqn{j}; please see \code{\link{boral}} for details. Defaults to \code{NULL}, so that if \code{row.params = NULL} then the argument is ignored, otherwise if \code{row.params} is supplied then \cr \code{row.ids = matrix(1:nrow(y), ncol = 1)} i.e., a single, row effect unique to each row. An internal check is done to see \code{row.params} and \code{row.ids} are consistent in terms of arguments supplied.}

 \item{offset}{A matrix with the same dimensions as the response matrix, specifying an a-priori known component to be included in the linear predictor during fitting. Defaults to \code{NULL}.}

 \item{num.lv}{The number of latent variables used in the fitted model. For models with no latent variables, please use \code{\link{calc.logLik.lv0}} to calculate the log-likelihood.}

  \item{lv.mc}{A matrix used for performing the Monte Carlo integration. Defaults to \code{NULL}, in which case a matrix is generated within the function.}

  \item{cutoffs}{Common cutoff estimates from the fitted model when any of the columns of the response matrix are ordinal responses. Defaults to \code{NULL}.}

\item{powerparam}{Common power parameter from the fitted model when any of the columns of the response matrix are tweedie responses. Defaults to \code{NULL}.}
}

\details{
For an \eqn{n x p} response matrix \eqn{\bm{Y}}, suppose we fit a model with one or more latent variables. If we denote the latent variables by \eqn{\bm{u}_i; i = 1,\ldots,n}, then the marginal log-likelihood is given by

\deqn{ \log(f) = \sum_{i=1}^n \log ( \int \prod_{j=1}^p \{f(y_{ij} | \bm{u}_i, \beta_{0j}, \bm{\theta}_j, \ldots) \} f(\bm{u}_i) d\bm{u}_i), }

where \eqn{f(y_{ij}|\cdot)} is the assumed distribution for column \eqn{j}, \eqn{\beta_{0j}} are the response-specific intercepts, \eqn{\bm{\theta}_j} are the response-specific latent variable coefficients, and \eqn{\ldots} generically denotes anything else included in the model, e.g. row effects, dispersion parameters etc... The quantity \eqn{f(\bm{u}_i)} denotes the distribution of the latent variable, which is assumed to be standard multivariate Gaussian. Standard Monte Carlo integration is used for calculating the marginal likelihood. If \code{lv.mc = NULL}, the function automatically generates a matrix as \cr \code{lv.mc <- rmvnorm(1000, rep(0,num.lv))}. If there is a need to apply this function numerous times, we recommend a matrix be inserted into \code{lv.mc} to speed up computation.

The key difference between this and the conditional likelihood (see \code{calc.condlogLik}) is that the marginal likelihood treats the latent variables as "random effects" and integrates over them, whereas the conditional likelihood treats the latent variables as "fixed effects". 

Please note the function is written conditional on all regression coefficients. Therefore, if traits are included in the model, in which case the regression coefficients \eqn{\beta_{0j}, \bm{\beta}_j} become random effects instead (please see \code{\link{about.traits}}), then the calculation of the log-likelihood does NOT take this into account, i.e. does not marginalize over them! Likewise if more than two columns are ordinal responses, then the regression coefficients \eqn{\beta_{0j}} corresponding to these columns become random effects, and the calculation of the log-likelihood also does NOT take this into account, i.e. does not marginalize over them!

When a single \eqn{\alpha_i} random row effect is inclued, then the log-likelihood is calculated by integrating over this,

\deqn{ \log(f) = \sum_{i=1}^n \log ( \int \prod_{j=1}^p \{f(y_{ij} | \bm{u}_i, \beta_{0j}, \alpha_i, \ldots)\} f(\bm{u}_i) f(\alpha_i) d\alpha_i ), } 
     
where \eqn{f(\alpha_i)} is the random effects distribution with mean zero and standard deviation given by the \code{row.params}. The integration is again performed using standard Monte Carlo integration. This naturally extends to multiple random row effects structures.
}

\value{
A list with the following components:
     \item{logLik}{Value of the marginal log-likelihood.}
     \item{logLik.comp}{A vector of the log-likelihood values for each row of the response matrix, \cr such that \code{sum(logLik.comp) = logLik}.}
}

\author{
\packageAuthor{boral}

Maintainer: \packageMaintainer{boral}
}

\section{Warnings}{As of version 1.6, this function is longer updated!}


\seealso{
\code{\link{calc.condlogLik}} for calculation of the conditional log-likelihood; 
\code{\link{calc.logLik.lv0}} to calculate the conditional/marginal log-likelihood for a model with no latent variables.
}

\examples{
\dontrun{
## NOTE: The values below MUST NOT be used in a real application;
## they are only used here to make the examples run quick!!!
example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
     n.thin = 1)

testpath <- file.path(tempdir(), "jagsboralmodel.txt")


library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
n <- nrow(y)
p <- ncol(y)
    
## Example 1 - model with two latent variables, site effects, 
## 	and no environmental covariates
spiderfit_nb <- boral(y, family = "negative.binomial", 
    lv.control = list(num.lv = 2), row.eff = "fixed", save.model = TRUE, 
    mcmc.control = example_mcmc_control, model.name = testpath)

## Extract all MCMC samples
fit_mcmc <- get.mcmcsamples(spiderfit_nb) 
mcmc_names <- colnames(fit_mcmc)

## Find the posterior medians
coef_mat <- matrix(apply(fit_mcmc[,grep("lv.coefs",mcmc_names)],
    2,median),nrow=p)
site_coef <- list(ID1 = apply(fit_mcmc[,grep("row.coefs.ID1", mcmc_names)],
    2,median))
     
## Calculate the marginal log-likelihood at the posterior median
calc.marglogLik(y, family = "negative.binomial",
    lv.coefs = coef_mat, row.eff = "fixed", row.params = site_coef, 
    num.lv = 2)

	
## Example 2 - model with one latent variable, no site effects, 
## 	and environmental covariates
spiderfit_nb2 <- boral(y, X = spider$x, family = "negative.binomial", 
     lv.control = list(num.lv = 2), save.model = TRUE, 
     mcmc.control = example_mcmc_control, model.name = testpath)

## Extract all MCMC samples
fit_mcmc <- get.mcmcsamples(spiderfit_nb2) 
mcmc_names <- colnames(fit_mcmc)

## Find the posterior medians
coef_mat <- matrix(apply(fit_mcmc[,grep("lv.coefs",mcmc_names)],
    2,median),nrow=p)
X_coef_mat <- matrix(apply(fit_mcmc[,grep("X.coefs",mcmc_names)],
    2,median),nrow=p)

## Calculate the log-likelihood at the posterior median
calc.marglogLik(y, X = spider$x, family = "negative.binomial", 
    lv.coefs = coef_mat, X.coefs = X_coef_mat, num.lv = 2)	
}
}

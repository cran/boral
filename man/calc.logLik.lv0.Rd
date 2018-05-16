\name{calc.logLik.lv0}

\alias{calc.logLik.lv0}

\title{Log-likelihood for a boral model with no latent variables}

\description{Calculates the log-likelihood for a set of parameter estimates from an boral model with no latent variables. If the row effects are assumed to be random, they are integrated over using Monte Carlo integration.}

\usage{
calc.logLik.lv0(y, X = NULL, family, trial.size = 1, lv.coefs, 
	X.coefs = NULL, row.eff = "none", row.params = NULL, 
	row.ids = NULL, offset = NULL, cutoffs = NULL,
	powerparam = NULL)
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

  \item{row.eff}{Single element indicating whether row effects are included as fixed effects ("fixed"), random effects ("random") or not included ("none") in the boral model. If random effects, they are drawn from a normal distribution with mean zero and standard deviation given by \code{row.params}. Defaults to "none". } 

  \item{row.params}{Parameters corresponding to the row effect from the boral model. If \cr \code{row.eff = "fixed"}, then these are the fixed effects and should have length equal to the number of columns in \code{y}. If \code{row.eff = "random"}, then this is the standard deviation for the random effects normal distribution. If \code{row.eff = "none"}, then this argument is ignored.}

 \item{row.ids}{A matrix with the number of rows equal to the number of rows in \code{y}, and the number of columns equal to the number of row effects to be included in the model. Element \eqn{(i,j)} indicates to the cluster ID of row \eqn{i} in \code{y} for random effect eqn{j}; please see \code{\link{boral}} for details. Defaults to \code{NULL}, so that if \code{row.params = NULL} then the argument is ignored, otherwise if \code{row.params} is supplied then \cr \code{row.ids = matrix(1:nrow(y), ncol = 1)} i.e., a single, row effect unique to each row. An internal check is done to see \code{row.params} and \code{row.ids} are consistent in terms of arguments supplied.}    

 \item{offset}{A matrix with the same dimensions as the response matrix \code{y}, specifying an a-priori known component to be included in the linear predictor during fitting. Defaults to \code{NULL}.}

 \item{cutoffs}{Common cutoff estimates from the boral model when any of the columns of \code{y} are ordinal responses. Defaults to \code{NULL}.}

\item{powerparam}{Common power parameter from the boral model when any of the columns of \code{y} are tweedie responses. Defaults to \code{NULL}.}
}


\details{
For an \eqn{n x p} response matrix \code{y}, the log-likelihood for a model with no latent variables included is given by,

\deqn{ \log(f) = \sum_{i=1}^n \sum_{j=1}^p \log \{f(y_{ij} | \beta_{0j}, \alpha_i, \ldots)\}, }

where \eqn{f(y_{ij}|\cdot)} is the assumed distribution for column \eqn{j}, \eqn{\beta_{0j}} is the column-specific intercepts, \eqn{\alpha_i} is the row effect, and \eqn{\ldots} generically denotes anything else included in the model, e.g. row effects, dispersion parameters etc...

Please note the function is written conditional on all regression coefficients. Therefore, if traits are included in the model, in which case the regression coefficients \eqn{\beta_{0j}, \bm{\beta}_j} become random effects instead (please see \code{\link{about.traits}}), then the calculation of the log-likelihood does NOT take this into account, i.e. does not marginalize over them!

Likewise if more than two columns are ordinal responses, then the regression coefficients \eqn{\beta_{0j}} corresponding to these columns become random effects, and the calculation of the log-likelihood also does NOT take this into account, i.e. does not marginalize over them!

When a single \eqn{\alpha_i} random row effect is inclued, then the log-likelihood is calculated by integrating over this,

\deqn{ \log(f) = \sum_{i=1}^n \log ( \int \prod_{j=1}^p \{f(y_{ij} | \beta_{0j}, \alpha_i, \ldots)\}f(\alpha_i) d\alpha_i ), } 
     
where \eqn{f(\alpha_i)} is the random effects distribution with mean zero and standard deviation given by the \code{row.params}. The integration is performed using standard Monte Carlo integration. This naturally extends to multiple random row effects structures.
}

\value{
A list with the following components:
    \item{logLik}{Value of the log-likelihood}
	\item{logLik.comp}{A vector of the log-likelihood values for each row of \code{y}, \cr such that \code{sum(logLik.comp) = logLik}.}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}; Wade Blanchard \email{wade.blanchard@anu.edu.au}
}

\seealso{
\code{\link{calc.marglogLik}} for calculation of the log-likelihood marginalizing over one or more latent variables, and \code{\link{calc.condlogLik}} for calculation of the conditional log-likelihood for models where everything is treated as "fixed effects", including latent variables, row effects, and so on.
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

## Example 1 - NULL model with site effects only
spiderfit_nb <- boral(y, family = "negative.binomial", 
     row.eff = "fixed", save.model = TRUE, mcmc.control = example_mcmc_control)

## Extract all MCMC samples
fit_mcmc <- get.mcmcsamples(spiderfit_nb) 
mcmc_names <- colnames(fit_mcmc)

## Find the posterior medians
coef_mat <- matrix(apply(fit_mcmc[,grep("lv.coefs",mcmc_names)],
     2,median),nrow=p)
site_coef <- list(ID1 = apply(fit_mcmc[,grep("row.coefs.ID1", mcmc_names)],
     2,median))

## Calculate the log-likelihood at the posterior median
calc.logLik.lv0(y, family = "negative.binomial",
     lv.coefs =  coef_mat, row.eff = "fixed", row.params = site_coef)


## Example 2 - Model with environmental covariates and random row effects
X <- scale(spider$x)
spiderfit_nb2 <- boral(y, X = X, family = "negative.binomial", 
     row.eff = "random",save.model = TRUE, mcmc.control = example_mcmc_control)

## Extract all MCMC samples
fit_mcmc <- get.mcmcsamples(spiderfit_nb2) 
mcmc_names <- colnames(fit_mcmc)

## Find the posterior medians
coef_mat <- matrix(apply(fit_mcmc[,grep("lv.coefs",mcmc_names)],
     2,median),nrow=p)
X_coef_mat <- matrix(apply(fit_mcmc[,grep("X.coefs",mcmc_names)],
	2,median),nrow=p)
site.sigma <- list(ID1 = 
	median(fit_mcmc[,grep("row.sigma.ID1", mcmc_names)]))

	
## Calculate the log-likelihood at the posterior median
calc.logLik.lv0(y, X = spider$x, family = "negative.binomial", 
     row.eff = "random",lv.coefs =  coef_mat, X.coefs = X_coef_mat, 
     row.params = site.sigma)
}
}

\name{calc.logLik.lv0}
\alias{calc.logLik.lv0}
\title{Log-likelihood for a boral model with no latent variables}
\description{Calculates the log-likelihood for a set of parameter estimates from an boral model with no latent variables. In such a model, the conditional and marginal log-likelihood are equivalent.}

\usage{
calc.logLik.lv0(y, X = NULL, family, trial.size = NULL, lv.coefs, 
	X.coefs = NULL, row.coefs = NULL, cutoffs = NULL,
     powerparam = NULL)
}

\arguments{
\item{y}{The response matrix the boral model was fitted to.}

\item{X}{The model matrix used in the boral model. Defaults to \code{NULL}, in which case it is assumed no model matrix was used.}  

\item{family}{Either a single element, or a vector of length equal to the number of columns in \eqn{y}. The former assumes all columns of \eqn{y} come from this distribution. The latter option allows for different distributions for each column of \eqn{y}. Elements can be one of "binomial" (with logit link), "poisson" (with log link), "negative.binomial" (with log link), "normal" (with identity link), "lnormal" for lognormal (with log link), "tweedie" (with log link), "exponential" (with log link), "gamma" (with log link), "beta" (with logit link), "ordinal" which performs cumulative logit regression. For the negative binomial distribution, the variance is parameterized as \eqn{Var(y) = \mu + \phi \mu^2}, where \eqn{\phi} is the column-specific overdispersion parameter. For the normal distribution, the variance is parameterized as \eqn{Var(y) = \phi}, where \eqn{\phi} is the column-specific variance. For the tweedie distribution, the variance is parameterized as \eqn{Var(y) = \phi \mu^p} where \eqn{\phi} is the column-specific dispersion parameter and \eqn{p} is a power parameter common to all columns assumed to be tweedie, with \eqn{1 < p < 2}. For the gamma distribution, the variance is parameterized as \eqn{Var(y) = \mu/\phi} where \eqn{\phi} is the column-specific rate (henceforth referred to also as dispersion parameter). For the beta distribution, the parameterization is in terms of the mean \eqn{\mu} and sample size \eqn{\phi} (henceforth referred to also as dispersion parameter), so that the two shape parameters are given by \eqn{a = \mu\phi} and \eqn{b = (1-\mu)\phi}.

All columns assumed to have ordinal responses are constrained to have the same cutoffs points, with a column-specific intercept to account for differences between the columns (please see \emph{Details} for formulation). 
}

 \item{trial.size}{Either equal to NULL, a single element, or a vector of length equal to the number of columns in \eqn{y}. If a single element, then all columns assumed to be binomially distributed will have trial size set to this. If a vector, different trial sizes are allowed in each column of \eqn{y}.}

  \item{lv.coefs}{The column-specific intercept, coefficient estimates relating to the latent variables, and dispersion parameters from the boral model.}

  \item{X.coefs}{The coefficients estimates relating to the model matrix \code{X} from the boral model. Defaults to \code{NULL}, in which it is assumed there are no covariates in the model.}

  \item{row.coefs}{Row effects estimates for the boral model. Defaults to \code{NULL}, in which case it is assumed there are no row effects in the model.} 

  \item{cutoffs}{Common cutoff estimates from the boral model when any of the columns of \eqn{y} are ordinal responses. Defaults to \code{NULL}.}

\item{powerparam}{Common power parameter from the boral model when any of the columns of \eqn{y} are tweedie responses. Defaults to \code{NULL}.}
}

\details{
A separate function to \code{\link{calc.marglogLik}} and \code{\link{calc.condlogLik}} has been written for boral models with with no latent variables, because the marginal and conditional log-likelihoods functions for such a model are equivalent. For an \eqn{n x p} response matrix \eqn{y}, the log-likelihood is given by (with parameters where appropriate),

\deqn{ \log(f) = \sum_{i=1}^n \sum_{j=1}^p \log (f(y_{ij} | \alpha_i, \tau_k, \theta_{0j}, \bm{x}_i, \bm{\beta}_j, \phi_j)), }

where \eqn{f(y_{ij}|\cdot)} is the assumed distribution for column \eqn{j}, \eqn{\alpha_i} is the row effect, \eqn{\tau_k} are the cutoffs for proportional odds regression, \eqn{\theta_{0j}} is the column-specific intercepts, \eqn{\bm{x}'_i} is row \eqn{i} of the model matrix, \eqn{\bm{\beta}_j} are the column-specific coefficients relating to the model matrix of covariates, \eqn{\phi_j} are column-specific dispersion parameters.

The function returns log-likelihood values both on a per-row and per-column basis. The latter is useful as it returns the log-likelihood for the separate GLMs fitted to each column of \eqn{y}. 
}

\value{
A list with the following components:
     \item{logLik}{Value of the log-likelihood.}
     \item{logLik.row.comp}{A vector of the log-likelihood values for each row of \eqn{y}, \cr such that \code{sum(logLik.row.comp) = logLik}.}
     \item{logLik.col.comp}{A vector of the log-likelihood values for each column of \eqn{y}, \cr such that \code{sum(logLik.row.comp) = logLik}.}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\seealso{
\code{\link{calc.marglogLik}} for calculation of the marginal log-likelihood, and \code{\link{calc.condlogLik}} for calculation of the conditional log-likelihood for boral models with one or more latent variables.
}

\examples{
\dontrun{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
n <- nrow(y); p <- ncol(y); 

## Example 1 - NULL model with site effects and species specific intercepts
spider.fit.nb <- boral(y, family = "negative.binomial", num.lv = 0, 
     row.eff = TRUE, save.model = TRUE, calc.ics = FALSE)

## Extract all MCMC samples
fit.mcmc <- as.mcmc(spider.fit.nb$jags.model)[[1]]

## Find the posterior medians
coef.mat <- matrix(apply(fit.mcmc[,grep("all.params",colnames(fit.mcmc))],
     2,median),nrow=p)
site.coef.median <- apply(fit.mcmc[,grep("site.params", colnames(fit.mcmc))],
     2,median)

## Calculate the log-likelihood at the posterior median
calc.logLik.lv0(y, family = "negative.binomial",
     lv.coefs =  coef.mat, row.coefs = site.coef.median)


## Example 2 - Model without site effects, latent variables, 
##   but includes environmental covariates
spider.fit.nb2 <- boral(y, X = spider$x, family = "negative.binomial", num.lv = 0, 
     row.eff = FALSE, save.model = TRUE, calc.ics = FALSE)

## Extract all MCMC samples
fit.mcmc <- as.mcmc(spider.fit.nb2$jags.model)[[1]]

## Find the posterior medians
coef.mat <- matrix(apply(fit.mcmc[,grep("all.params",colnames(fit.mcmc))],
     2,median),nrow=p)
X.coef.mat <- matrix(apply(fit.mcmc[,grep("X.params",colnames(fit.mcmc))],
	2,median),nrow=p)

## Calculate the log-likelihood at the posterior median
calc.logLik.lv0(y, X = spider$x, family = "negative.binomial", 
	lv.coefs =  coef.mat, X.coefs = X.coef.mat)
}     
}

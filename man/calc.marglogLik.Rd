\name{calc.marglogLik}
\alias{calc.marglogLik}
\title{Marginal log-likelihood for an boral model}
\description{Calculates the marginal log-likelihood for a set of parameter estimates from an boral model, whereby the latent variables are integrated out. The integration is performed using Monte Carlo integration.}

\usage{
calc.marglogLik(y, X = NULL, family, trial.size = NULL, lv.coefs, 
     X.coefs = NULL, site.coefs = NULL, num.lv, 
     X.mc = NULL, cutoffs = NULL, powerparam = NULL)
}

\arguments{
  \item{y}{The response matrix that the boral model was fitted to.}

  \item{X}{The model matrix used in the boral model. Defaults to \code{NULL}, in which case it is assumed no model matrix was used.}  

\item{family}{Either a single element, or a vector of length equal to the number of columns in \eqn{y}. The former assumes all columns of \eqn{y} come from this distribution. The latter option allows for different distributions for each column of \eqn{y}. Elements can be one of "binomial" (with logit link), "poisson" (with log link), "negative.binomial" (with log link), "normal" (with identity link), "lnormal" for lognormal (with log link), "tweedie" (with log link), "ordinal" which performs cumulative logit regression. For the negative binomial distribution, the variance is parameterized as \eqn{Var(y) = \mu + \phi \mu^2}, where \eqn{\phi} is the column-specific overdispersion parameter. For the normal distribution, the variance is parameterized as \eqn{Var(y) = \phi}, where \eqn{\phi} is the column-specific variance. For the tweedie distribution, the variance is parameterized as \eqn{Var(y) = \phi \mu^p} where \eqn{\phi} is the column-specific dispersion parameter and \eqn{p} is a power parameter common to all columns assumed to be tweedie, with \eqn{1 < p < 2}. 

All columns assumed to have ordinal responses are constrained to have the same cutoffs points, with a column-specific intercept to account for differences between the columns (please see \emph{Details} for formulation). 
}

 \item{trial.size}{Either equal to NULL, a single element, or a vector of length equal to the number of columns in \eqn{y}. If a single element, then all columns assumed to be binomially distributed will have trial size set to this. If a vector, different trial sizes are allowed in each column of \eqn{y}.}

  \item{lv.coefs}{The column-specific intercept, coefficient estimates relating to the latent variables, and dispersion parameters from the boral model.}

  \item{X.coefs}{The coefficients estimates relating to the model matrix \code{X} from the boral model. Defaults to \code{NULL}, in which it is assumed there are no covariates in the model.}

  \item{site.coefs}{Row effects estimates for the boral model. Defaults to \code{NULL}, in which case it is assumed there are no row effects in the model.} 

  \item{num.lv}{The number of latent variables used in the boral model. For boral models with no latent variables, please use \code{\link{calc.logLik.lv0}} to calculate the marginal log-likelihood.}

  \item{X.mc}{A matrix used for performing the Monte Carlo integration. Defaults to \code{NULL}, in which case a matrix is generated within the function.}

  \item{cutoffs}{Common cutoff estimates from the boral model when any of the columns of \eqn{y} are ordinal responses. Defaults to \code{NULL}.}

\item{powerparam}{Common power parameter from the boral model when any of the columns of \eqn{y} are tweedie responses. Defaults to \code{NULL}.}
}

\details{
For an \eqn{n x p} response matrix \eqn{y}, suppose we fit an boral model with \eqn{q} latent varables. If we denote the latent varibles by \eqn{b_i; i = 1,\ldots,n}, then the marginal log-likelihood is given by (with parameters where appropriate)

\deqn{ \log(f) = \sum_{i=1}^n \log ( \int \prod_{j=1}^p f(y_{ij} | \alpha_i, \tau_k, \theta_{0j}, \bm{\theta}_j, \bm{z}_i, \bm{x}_i, \bm{\beta}_j, \phi_j) f(\bm{z}_i) d\bm{z}_i), }

where \eqn{f(y_{ij}|\cdot)} is the assumed distribution for column \eqn{j}, \eqn{\alpha_i} is the row effect, \eqn{\tau_k} are the cutoffs for proportional odds regression, \eqn{\theta_{0j}} is the column-specific intercepts, \eqn{\bm{z}'_i} are the latent variables, \eqn{\bm{\theta}_j} are the column-specific coefficients relating to the latent variables, \eqn{\bm{x}'_i} is row \eqn{i} of the model matrix, \eqn{\bm{\beta}_j} are the column-specific coefficients relating to the model matrix of covariates, \eqn{\phi_j} are column-specific dispersion parameters.

The quantity \eqn{f(\bm{z}_i)} denotes the distribution of the latent variable, which is assumed to be standard multivariate Gaussian for boral models.

The key difference between this and the conditional likelihood (see \code{calc.condlogLik}) is that the marginal log-likelihood treats the latent variables as "random effects" and integrates over them, whereas the conditional log-likelihood treats the latent variables as "fixed effects". 

Monte Carlo integration is used for calculating the marginal log-likelihood. If \code{X.mc = NULL}, the function automatically generates a matrix as \code{X.mc <- cbind(1, rmvnorm(5000, rep(0,num.lv)))}. If there is need to apply this function numerous times, we recommend a matrix be inserted into \code{X.mc} to speed up computation.
}

\value{
A list with the following components:
     \item{logLik}{Value of the marginal log-likelihood.}
     \item{logLik.comp}{A vector of the log-likelihood values for each row of \eqn{y}, \cr such that \code{sum(logLik.comp) = logLik}.}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\note{
The AIC at posterior median, BIC at posterior median, and compound Laplce-Metropolis estimator at posterior median, as returned from \code{\link{get.measures}}, are all based on the marginal log-likelihood calculated from this function. Additionally, \code{\link{get.more.measures}} returns even more information criteria based on the marginal log-likelihod.}

\seealso{
\code{\link{get.measures}} and \code{\link{get.more.measures}} for information criteria based on the marginal log-likelihood; \code{\link{calc.condlogLik}} for calculation of the conditional log-likelihood; \code{\link{calc.logLik.lv0}} to calculate the conditional/marginal log-likelihood for an boral model with no latent variables.
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
     site.eff = TRUE, save.model = TRUE, calc.ics = FALSE)

## Extract all MCMC samples
fit.mcmc <- as.mcmc(spider.fit.nb$jags.model)[[1]]

## Find the posterior medians
coef.mat <- matrix(apply(fit.mcmc[,grep("all.params",colnames(fit.mcmc))],
     2,median),nrow=p)
site.coef.median <- apply(fit.mcmc[,grep("site.params", colnames(fit.mcmc))],
     2,median)
     
## Caculate the marginal log-likelihood at the posterior median
calc.marglogLik(y, family = "negative.binomial",
	lv.coefs = coef.mat, site.coefs = site.coef.median, num.lv = 2, X.mc = NULL)

	
## Example 2 - model with one latent variable, no site effects, 
## 	and environmental covariates
spider.fit.nb2 <- boral(y, X = spider$x, family = "negative.binomial", 
     num.lv = 1, site.eff = FALSE, save.model = TRUE, calc.ics = FALSE)

## Extract all MCMC samples
fit.mcmc <- as.mcmc(spider.fit.nb2$jags.model)[[1]]

## Find the posterior medians
coef.mat <- matrix(apply(fit.mcmc[,grep("all.params",colnames(fit.mcmc))],
     2,median),nrow=p)
X.coef.mat <- matrix(apply(fit.mcmc[,grep("X.params",colnames(fit.mcmc))],
	2,median),nrow=p)

## Caculate the log-likelihood at the posterior median
calc.marglogLik(y, X = spider$x, family = "negative.binomial", 
	lv.coefs = coef.mat, X.coefs = X.coef.mat, num.lv = 1)	
	}
}

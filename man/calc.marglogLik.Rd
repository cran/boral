\name{calc.marglogLik}
\alias{calc.marglogLik}
\title{Marginal log-likelihood for an boral model}
\description{Calculates the marginal log-likelihood for a set of parameter estimates from an boral model, whereby the latent variables and random effects (if applicable) are integrated out. The integration is performed using Monte Carlo integration.}

\usage{
calc.marglogLik(y, X = NULL, family, trial.size = 1, lv.coefs, 
     X.coefs = NULL, row.eff = "none", row.params = NULL, row.ids = NULL,
     num.lv, lv.mc = NULL, cutoffs = NULL, powerparam = NULL)
}

\arguments{
\item{y}{The response matrix that the boral model was fitted to.}

\item{X}{The model matrix used in the boral model. Defaults to \code{NULL}, in which case it is assumed no model matrix was used.}  

  \item{family}{Either a single element, or a vector of length equal to the number of columns in \code{y}. The former assumes all columns of \code{y} come from this distribution. The latter option allows for different distributions for each column of \code{y}. Elements can be one of "binomial" (with probit link), "poisson" (with log link), "negative.binomial" (with log link), "normal" (with identity link), "lnormal" for lognormal (with log link), "tweedie" (with log link), "exponential" (with log link), "gamma" (with log link), "beta" (with logit link), "ordinal" (cumulative probit regression). 

  For the negative binomial distribution, the variance is parameterized as \eqn{Var(y) = \mu + \phi\mu^2}, where \eqn{\phi} is the column-specific dispersion parameter. For the normal distribution, the variance is parameterized as \eqn{Var(y) = \phi^2}, where \eqn{\phi} is the column-specific standard deviation. For the tweedie distribution, the variance is parameterized as \eqn{Var(y) = \phi \mu^p} where \eqn{\phi} is the column-specific dispersion parameter and \eqn{p} is a power parameter common to all columns assumed to be tweedie, with \eqn{1 < p < 2}. For the gamma distribution, the variance is parameterized as \eqn{Var(y) = \mu/\phi} where \eqn{\phi} is the column-specific rate (henceforth referred to also as dispersion parameter). For the beta distribution, the parameterization is in terms of the mean \eqn{\mu} and sample size \eqn{\phi} (henceforth referred to also as dispersion parameter), so that the two shape parameters are given by \eqn{a = \mu\phi} and \eqn{b = (1-\mu)\phi}.

 All columns assumed to have ordinal responses are constrained to have the same cutoffs points, with a column-specific intercept to account for differences between the columns (please see \emph{Details} for formulation). 
}

 \item{trial.size}{Either equal to a single element, or a vector of length equal to the number of columns in y. If a single element, then all columns assumed to be binomially distributed will have trial size set to this. If a vector, different trial sizes are allowed in each column of y. The argument is ignored for all columns not assumed to be binomially distributed. Defaults to 1, i.e. Bernoulli distribution.}

  \item{lv.coefs}{The column-specific intercept, coefficient estimates relating to the latent variables, and dispersion parameters from the boral model.}

  \item{X.coefs}{The coefficients estimates relating to the model matrix \code{X} from the boral model. Defaults to \code{NULL}, in which it is assumed there are no covariates in the model.}

  \item{row.eff}{Single element indicating whether row effects are included as fixed effects ("fixed"), random effects ("random") or not included ("none") in the boral model. If random effects, they are drawn from a normal distribution with mean zero and standard deviation given by \code{row.params}. Defaults to "none". } 

  \item{row.params}{Parameters corresponding to the row effect from the boral model. If \cr \code{row.eff = "fixed"}, then these are the fixed effects and should have length equal to the number of columns in \code{y}. If \code{row.eff = "random"}, then this is standard deviation for the random effects normal distribution. If \code{row.eff = "none"}, then this argument is ignored.}

 \item{row.ids}{A matrix with the number of rows equal to the number of rows in \code{y}, and the number of columns equal to the number of row effects to be included in the model. Element \eqn{(i,j)} indicates to the cluster ID of row \eqn{i} in \code{y} for random effect eqn{j}; please see the help file for the main boral function for details. Defaults to \code{NULL}, so that if \code{row.params = NULL} then the argument is ignored, otherwise if \code{row.params} is supplied then \code{row.ids = matrix(1:nrow(y), ncol = 1)} i.e., a single, row effect unique to each row. An internal check is done to see \code{row.params} and \code{row.ids} are consistent in terms of arguments supplied.}    

 \item{num.lv}{The number of latent variables used in the boral model. For boral models with no latent variables, please use \code{\link{calc.logLik.lv0}} to calculate the log-likelihood.}

  \item{lv.mc}{A matrix used for performing the Monte Carlo integration. Defaults to \code{NULL}, in which case a matrix is generated within the function.}

  \item{cutoffs}{Common cutoff estimates from the boral model when any of the columns of \code{y} are ordinal responses. Defaults to \code{NULL}.}

\item{powerparam}{Common power parameter from the boral model when any of the columns of \code{y} are tweedie responses. Defaults to \code{NULL}.}
}

\details{
For an \eqn{n x p} response matrix \code{y}, suppose we fit an boral model with one or more latent variables. If we denote the latent variables by \eqn{\bm{z}_i; i = 1,\ldots,n}, then the marginal log-likelihood is given by

\deqn{ \log(f) = \sum_{i=1}^n \log ( \int \prod_{j=1}^p f(y_{ij} | \bm{z}_i, \beta_{0j}, \bm{\theta}_j, \ldots) f(\bm{z}_i) d\bm{z}_i), }

where \eqn{f(y_{ij}|\cdot)} is the assumed distribution for column \eqn{j}, \eqn{\beta_{0j}} are the column-specific intercepts, \eqn{\bm{\theta}_j} are the column-specific latent variable coefficients, and \eqn{\ldots} generically denotes anything else included in the model, e.g. row effects, dispersion parameters etc... The quantity \eqn{f(\bm{z}_i)} denotes the distribution of the latent variable, which is assumed to be standard multivariate Gaussian. Standard Monte Carlo integration is used for calculating the marginal likelihood. If \code{lv.mc = NULL}, the function automatically generates a matrix as \cr \code{lv.mc <- cbind(1, rmvnorm(2000, rep(0,num.lv)))}. If there is a need to apply this function numerous times, we recommend a matrix be inserted into \code{lv.mc} to speed up computation.

The key difference between this and the conditional likelihood (see \code{calc.condlogLik}) is that the marginal likelihood treats the latent variables as "random effects" and integrates over them, whereas the conditional likelihood treats the latent variables as "fixed effects". 

Please note the function is written conditional on all regression coefficients. Therefore, if traits are included in the model, in which case the regression coefficients \eqn{\beta_{0j}, \bm{\beta}_j} become random effects instead (please see section titled ``Including species traits" in the main \code{boral} function), then the calculation of the log-likelihood does NOT take this into account, i.e. does not marginalize over them!

Likewise if more than two columns are ordinal responses, then the regression coefficients \eqn{\beta_{0j}} corresponding to these columns become random effects, and the calculation of the log-likelihood also does NOT take this into account, i.e. does not marginalize over them!

In the special case where \eqn{\alpha_i} is a random row effect, then the log-likelihood is calculated by integrating over this as well,

\deqn{\log(f) = \sum_{i=1}^n \log ( \int \prod_{j=1}^p (f(y_{ij} | \bm{z}_i, \beta_{0j}, \bm{\theta}_j, \alpha_i, \ldots)) f(\bm{z}_i) f(\alpha_i) d\bm{z}_i d\alpha_i ),
} 
     
where \eqn{f(\alpha_i)} is the random effects distribution with mean zero and standard deviation given by the \code{row.params}. The integration is performed using standard Monte Carlo integration. Beyond this however, for more complicated random row effects structures the function currently does not permit the calculation of the log-likelihood in such cases, as the integration over these random effects becomes too complicated to perform in such case. Sorry!
}

\value{
A list with the following components:
     \item{logLik}{Value of the marginal log-likelihood.}
     \item{logLik.comp}{A vector of the log-likelihood values for each row of \code{y}, \cr such that \code{sum(logLik.comp) = logLik}.}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\note{
The AIC and BIC at posterior median returned from \code{\link{get.measures}} are all based on the marginal log-likelihood calculated from this function. Additionally, \code{\link{get.more.measures}} returns even more information criteria based on the marginal log-likelihood. As mentioned in the details though, these information criteria do not take into account that traits are included in the model!
}

\seealso{
\code{\link{get.measures}} and \code{\link{get.more.measures}} for information criteria based on the marginal log-likelihood; \code{\link{calc.condlogLik}} for calculation of the conditional log-likelihood; \code{\link{calc.logLik.lv0}} to calculate the conditional/marginal log-likelihood for an boral model with no latent variables.
}

\examples{
\dontrun{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
n <- nrow(y); p <- ncol(y); 
    
## Example 1 - model with two latent variables, site effects, 
## 	and no environmental covariates
spider.fit.nb <- boral(y, family = "negative.binomial", num.lv = 2, 
     row.eff = "fixed", save.model = TRUE, calc.ics = FALSE)

## Extract all MCMC samples
fit.mcmc <- mcmc(spider.fit.nb$jags.model$BUGSoutput$sims.matrix) 

## Find the posterior medians
coef.mat <- matrix(apply(fit.mcmc[,grep("all.params",colnames(fit.mcmc))],
     2,median),nrow=p)
site.coef <- list(ID1 = apply(fit.mcmc[,grep("row.params", colnames(fit.mcmc))],
     2,median))
     
## Calculate the marginal log-likelihood at the posterior median
calc.marglogLik(y, family = "negative.binomial",
	lv.coefs = coef.mat, row.eff = "fixed", row.params = site.coef, 
	num.lv = 2)

	
## Example 2 - model with one latent variable, no site effects, 
## 	and environmental covariates
spider.fit.nb2 <- boral(y, X = spider$x, family = "negative.binomial", 
     num.lv = 2, save.model = TRUE, calc.ics = FALSE)

## Extract all MCMC samples
fit.mcmc <- mcmc(spider.fit.nb2$jags.model$BUGSoutput$sims.matrix) 

## Find the posterior medians
coef.mat <- matrix(apply(fit.mcmc[,grep("all.params",colnames(fit.mcmc))],
     2,median),nrow=p)
X.coef.mat <- matrix(apply(fit.mcmc[,grep("X.params",colnames(fit.mcmc))],
	2,median),nrow=p)

## Calculate the log-likelihood at the posterior median
calc.marglogLik(y, X = spider$x, family = "negative.binomial", 
	lv.coefs = coef.mat, X.coefs = X.coef.mat, num.lv = 2)	
}
}

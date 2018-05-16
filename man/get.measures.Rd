\name{get.measures}

\alias{get.measures}

\title{Information Criteria for boral models}

\description{Calculates some information criteria for an boral model, which could be used for model selection. WARNING: As of version 1.5, this function will no longer be updated...use at your own peril!!!}

\usage{
get.measures(y, X = NULL, family, trial.size = 1, row.eff = "none", 
	row.ids = NULL, offset = NULL, num.lv, fit.mcmc)
}

\arguments{
  \item{y}{The response matrix that the boral model was fitted to.}
  
  \item{X}{The model matrix used in the boral model. Defaults to \code{NULL}, in which case it is assumed no model matrix was used.}  
  
  \item{family}{Either a single element, or a vector of length equal to the number of columns in \code{y}. The former assumes all columns of \code{y} come from this distribution. The latter option allows for different distributions for each column of \code{y}. Elements can be one of "binomial" (with probit link), "poisson" (with log link), "negative.binomial" (with log link), "normal" (with identity link), "lnormal" for lognormal (with log link), "tweedie" (with log link), "exponential" (with log link), "gamma" (with log link), "beta" (with logit link), "ordinal" (cumulative probit regression). 

  Please see \code{\link{about.distributions}} for information on distributions available in boral overall.
  }
  
 \item{trial.size}{Either equal to a single element, or a vector of length equal to the number of columns in y. If a single element, then all columns assumed to be binomially distributed will have trial size set to this. If a vector, different trial sizes are allowed in each column of y. The argument is ignored for all columns not assumed to be binomially distributed. Defaults to 1, i.e. Bernoulli distribution.}
  
  \item{row.eff}{Single element indicating whether row effects are included as fixed effects ("fixed"), random effects ("random") or not included ("none") in the boral model. If random effects, they are drawn from a normal distribution with mean zero and unknown standard deviation. Defaults to "none". } 
  
 \item{row.ids}{A matrix with the number of rows equal to the number of rows in \code{y}, and the number of columns equal to the number of row effects to be included in the model. Element \eqn{(i,j)} indicates to the cluster ID of row \eqn{i} in \code{y} for random effect eqn{j}; please see \code{\link{boral}} for details. Defaults to \code{NULL}, so that if \code{row.eff = "none"} then the argument is ignored, otherwise if \cr \code{row.eff = "fixed"} or \code{"random"}, \cr then \code{row.ids = matrix(1:nrow(y), ncol = 1)} i.e., a single, row effect unique to each row.}    

 \item{offset}{A matrix with the same dimensions as the response matrix \code{y}, specifying an a-priori known component to be included in the linear predictor during fitting. Defaults to \code{NULL}.}

 \item{num.lv}{The number of latent variables used in the fitted boral model.}
  
  \item{fit.mcmc}{All MCMC samples for the fitted boral model. These can be extracted by fitting an boral model using \code{\link{boral}} with \code{save.model = TRUE}, and then applying \code{get.mcmcsamples(fit)}.} 
}

\details{
The following information criteria are currently calculated, when permitted: 1) Widely Applicable Information Criterion (WAIC, Watanabe, 2010) based on the conditional log-likelihood; 2) expected AIC (EAIC, Carlin and Louis, 2011); 3) expected BIC (EBIC, Carlin and Louis, 2011); 4) AIC (using the marginal likelihood) evaluated at the posterior median; 5) BIC (using the marginal likelihood) evaluated at the posterior median.

1) WAIC has been argued to be more natural and extension of AIC to the Bayesian and hierarchical modeling context (Gelman et al., 2013), and is based on the conditional log-likelihood calculated at each of the MCMC samples. 

2 & 3) EAIC and EBIC were suggested by (Carlin and Louis, 2011). Both criteria are of the form -2*mean(conditional log-likelihood) + penalty*(no. of parameters in the model), where the mean is averaged all the MCMC samples. EAIC applies a penalty of 2, while EBIC applies a penalty of \eqn{log(n)}.

4 & 5) AIC and BIC take the form -2*(marginal log-likelihood) + penalty*(no. of parameters in the model), where the log-likelihood is evaluated at the posterior median. If the parameter-wise posterior distributions are unimodal and approximately symmetric, these will produce similar results to an AIC and BIC where the log-likelihood is evaluated at the posterior mode. EAIC applies a penalty of 2, while EBIC applies a penalty of \eqn{log(n)}.

Intuitively, comparing boral models with and without latent variables (using information criteria such as those returned) amounts to testing whether the columns of the response matrix \code{y} are correlated. With multivariate abundance data for example, where \code{y} is a matrix of \eqn{n} sites and \eqn{p} species, comparing models with and without latent variables tests whether there is any evidence of correlation between species.

Please note that criteria 4 and 5 are not calculated all the time. In models where traits are included in the model (such that the regression coefficients \eqn{\beta_{0j}, \bm{\beta}_j} are random effects), or more than two columns are ordinal responses (such that the intercepts \eqn{\beta_{0j}} for these columns are random effects), then criteria 4 and 5 are will not calculated. This is because the calculation of the marginal log-likelihood in such cases currently fail to marginalize over such random effects; please see the details in \code{calc.logLik.lv0} and \code{calc.marglogLik}.
}

\value{
A list with the following components:
\item{waic}{WAIC based on the conditional log-likelihood.}
\item{eaic}{EAIC based on the mean of the conditional log-likelihood.}
\item{ebic}{EBIC based on the mean of the conditional log-likelihood.}
\item{all.cond.logLik}{The conditional log-likelihood evaluated at all MCMC samples. This is done via repeated application of \code{\link{calc.condlogLik}}.}
\item{cond.num.params}{Number of estimated parameters used in the fitted model, when all parameters are treated as "fixed" effects.}
\item{do.marglik.ics}{A boolean indicating whether marginal log-likelihood based information criteria are calculated.}

If \code{do.marglik.ics = TRUE}, then we also have:
\item{median.logLik}{The marginal log-likelihood evaluated at the posterior median.}
\item{marg.num.params}{Number of estimated parameters used in the fitted model, when all parameters are treated as "fixed" effects.}
\item{aic.median}{AIC (using the marginal log-likelihood) evaluated at the posterior median.}
\item{bic.median}{BIC (using the marginal log-likelihood) evaluated at the posterior median.}
}

\section{Warning}{
As of version 1.5, this function will no longer be updated...use at your own peril!!!

Using information criterion for variable selection should be done with extreme caution, for two reasons: 1) The implementation of these criteria are both \emph{heuristic} and experimental. 2) Deciding what model to fit for ordination purposes should be driven by the science. For example, it may be the case that a criterion suggests a model with 3 or 4 latent variables. However, if we interested in visualizing the data for ordination purposes, then models with 1 or 2 latent variables are far more appropriate. As an another example, whether or not we include row effects when ordinating multivariate abundance data depends on if we are interested in differences between sites in terms of relative species abundance (\code{row.eff = FALSE}) or in terms of species composition (\code{row.eff = "fixed"}).  

Also, the use of information criterion in the presence of variable selection using SSVS is questionable.
}

\references{
\itemize{
\item Carlin, B. P., and Louis, T. A. (2011). Bayesian methods for data analysis. CRC Press.

\item Gelman et al. (2013). Understanding predictive information criteria for Bayesian models. Statistics and Computing, 1-20.

\item Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. The Journal of Machine Learning Research, 11, 3571-3594.
}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}; Wade Blanchard \email{wade.blanchard@anu.edu.au}
}

\note{
When a boral model is fitted using \code{\link{boral}} with \code{calc.ics = TRUE}, then this function is applied and the information criteria are returned as part of the model output. 
}

\seealso{
\code{\link{get.dic}} for calculating the Deviance Information Criterion (DIC) based on the conditional log-likelihood; \code{\link{get.more.measures}} for even more information criteria.}

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

spiderfit_pois <- boral(y, family = "poisson", 
	num.lv = 2, row.eff = "random",
	mcmc.control = example_mcmc_control)

spiderfit_pois$ics ## Returns information criteria

spiderfit_nb <- boral(y, family = "negative.binomial", 
	num.lv = 2, row.eff = "random",
	mcmc.control = example_mcmc_control)

spiderfit_nb$ics ## Returns the information criteria 
}
}
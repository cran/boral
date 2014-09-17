\name{get.measures}
\alias{get.measures}
\title{Information Criteria for boral models}
\description{Calculates some information criteria for an boral model, which could be used for model selection.}

\usage{
get.measures(y, X = NULL, family, trial.size  = NULL, row.eff, num.lv, 
	fit.mcmc, more.measures = FALSE)
}

\arguments{
  \item{y}{The response matrix that the boral model was fitted to.}
  
  \item{X}{The model matrix used in the boral model. Defaults to \code{NULL}, in which case it is assumed no model matrix was used.}  
  
\item{family}{Either a single element, or a vector of length equal to the number of columns in \eqn{y}. The former assumes all columns of \eqn{y} come from this distribution. The latter option allows for different distributions for each column of \eqn{y}. Elements can be one of "binomial" (with logit link), "poisson" (with log link), "negative.binomial" (with log link), "normal" (with identity link), "lnormal" for lognormal (with log link), "tweedie" (with log link), "exponential" (with log link), "gamma" (with log link), "beta" (with logit link), "ordinal" which performs cumulative logit regression. For the negative binomial distribution, the variance is parameterized as \eqn{Var(y) = \mu + \phi \mu^2}, where \eqn{\phi} is the column-specific overdispersion parameter. For the normal distribution, the variance is parameterized as \eqn{Var(y) = \phi}, where \eqn{\phi} is the column-specific variance. For the tweedie distribution, the variance is parameterized as \eqn{Var(y) = \phi \mu^p} where \eqn{\phi} is the column-specific dispersion parameter and \eqn{p} is a power parameter common to all columns assumed to be tweedie, with \eqn{1 < p < 2}. For the gamma distribution, the variance is parameterized as \eqn{Var(y) = \mu/\phi} where \eqn{\phi} is the column-specific rate (henceforth referred to also as dispersion parameter). For the beta distribution, the parameterization is in terms of the mean \eqn{\mu} and sample size \eqn{\phi} (henceforth referred to also as dispersion parameter), so that the two shape parameters are given by \eqn{a = \mu\phi} and \eqn{b = (1-\mu)\phi}.

All columns assumed to have ordinal responses are constrained to have the same cutoffs points, with a column-specific intercept to account for differences between the columns (please see \emph{Details} for formulation). 
}
  
 \item{trial.size}{Either equal to NULL, a single element, or a vector of length equal to the number of columns in \eqn{y}. If a single element, then all columns assumed to be binomially distributed will have trial size set to this. If a vector, different trial sizes are allowed in each column of \eqn{y}.}
  
  \item{row.eff}{A logical value indicating whether to row effects were included in the model.}    
  
  \item{num.lv}{The number of latent variables used in the fitted boral model.}
  
  \item{fit.mcmc}{All MCMC samples for the fitted boral model, as obtained from JAGS. These can be extracted by fitting an boral model using \code{\link{boral}} with \code{save.model = TRUE}, and then applying \code{as.mcmc} on the \code{jags.model} component of the output.} 

  \item{more.measures}{A logical value indicating whether to run \code{\link{get.more.measures}} to obtain additional information criteria.}    
}

\details{
Currently, six information criteria has been implemented: 1) Widely Applicable Information Criterion (WAIC, Watanabe, 2010) based on the conditional log-likelihood; 2) expected AIC (EAIC, Carlin and Louis, 2011); 3) expected BIC (EBIC, Carlin and Louis, 2011); 4) AIC (using the marginal likelihood) evaluated at the posterior median; 5) BIC (using the marginal likelihood) evaluated at the posterior median; 6) Compound Laplace-Metropolis estimator of the model likelihood (Lewis and Raftery, 2011). 

1) WAIC has been argued to be more natural and extension of AIC to the Bayesian and hierarchical modeling context (Gelman et al., 2013), and is based on the conditional log-likelihood calculated at each of the MCMC samples. 

2 & 3) EAIC and EBIC were suggested by (Carlin and Louis, 2011). Both criteria are of the form -2*mean(conditional log-likelihood) + penalty*(no. of parameters in the model), where the mean is averaged all the MCMC samples. EAIC applies a penalty of 2, while EBIC applies a penalty of \eqn{log(n)}.

4 & 5) AIC and BIC take the form -2*(marginal log-likelihood) + penalty*(no. of parameters in the model), where the log-likelihood is evaluated at the posterior median. If the parameter-wise posterior distributions are unimodal and approximately symmetric, these will produce similar results to an AIC and BIC where the log-likelihood is evaluated at the posterior mode. EAIC applies a penalty of 2, while EBIC applies a penalty of \eqn{log(n)}.

6) The model likelihood is the probability of the data given a model, and both BIC and the compound Laplace-Metropolis estimator are based on asymptotic approximations to the this. However, while the first term in both criteria are of the same form, namely -2*(marginal log-likelihood), where the log-likelihood is evaluated at the posterior median, the compound Laplace-Metropolis estimator explicitly calculates the determinant of the relevant Hessian matrix (evaluated at the posterior median) is to derive its penalty. 

In our very limited experience, if information criteria are to be used for model selection between boral models, we found BIC at the posterior median and the compound Laplace-Metropolis estimator tend to perform best. WAIC, AIC, and DIC (see \code{\link{get.dic}}) tend to over select the number of latent variables. For WAIC and DIC, part of this overfitting could be due to the fact both criteria are calculated from the conditional rather than the marginal log-likelihood (see Millar, 2009).  

Intuitively, comparing boral models with and without latent variables (using information criteria such as those returned) amounts to testing whether the columns of the response matrix \eqn{y} are correlated. With multivariate abundance data for example, where \eqn{y} is a matrix of \eqn{n} sites and \eqn{p} species, comparing models with and without latent variables tests whether there is any evidence of correlation between species.

}

\value{
A list with the following components:
\item{waic}{WAIC based on the conditional log-likelihood.}
\item{eaic}{EAIC based on the mean of the conditional log-likelihood.}
\item{ebic}{EBIC based on the mean of the conditional log-likelihood.}
\item{aic.median}{AIC (using the marginal log-likelihood) evaluated at the posterior median.}
\item{bic.median}{BIC (using the marginal log-likelihood) evaluated at the posterior median.}
\item{comp.lm}{Compound Laplace-Metropolis estimator of the model likelihood, evaluated at the posterior median.}
\item{all.cond.logLik}{The conditional log-likelihood evaluated at all MCMC samples. This is done via repeated application of \code{\link{calc.condlogLik}}.}
\item{num.params}{Number of estimated parameters used in the fitted model.}
}

\section{Warning}{
Using information criterion for variable selection should be done with extreme caution, for two reasons: 1) The implementation of these criteria are both \emph{heuristic} and experimental. 2) Deciding what model to fit for ordination purposes should be driven by the science. For example, it may be the case that a criterion suggests a model with 3 or 4 latent variables. However, if we interested in visualizing the data for ordination purposes, then models with 1 or 2 latent variables are far more appropriate. As an another example, whether or not we include row effects when ordinating multivariate abundance data depends on if we are interested in differences between sites in terms of relative species abundance (\code{site.eff = FALSE}) or in terms of species composition (\code{site.eff = TRUE}).  

Also, the use of information criterion in the presence of variable selection using SSVS is questionable.
}

\references{
\itemize{
\item Carlin, B. P., and Louis, T. A. (2011). Bayesian methods for data analysis. CRC Press.
\item Gelman et al. (2013). Understanding predictive information criteria for Bayesian models. Statistics and Computing, 1-20.
\item Lewis, S. M., and Raftery, A. E. (1997). Estimating Bayes factors via posterior simulation with the Laplace-Metropolis estimator. Journal of the American Statistical Association, 92, 648-655.
\item Millar, R. B. (2009). Comparison of hierarchical Bayesian models for overdispersed count data using DIC and Bayes' factors. Biometrics, 65, 962-969.
\item Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. The Journal of Machine Learning Research, 11, 3571-3594.
}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\note{
When a boral model is fitted using \code{\link{boral}} with \code{calc.ics = TRUE}, then this function is applied with \code{more.measures = FALSE}, and the information criteria are returned as part of the model output. 
}

\seealso{
\code{\link{get.dic}} for calculating the Deviance Information Criterion (DIC) based on the conditional log-likelihood; \code{\link{get.more.measures}} for even more information criteria.}

\examples{
\dontrun{

library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
n <- nrow(y); p <- ncol(y); 
    
spider.fit.pois <- boral(y, family = "poisson", num.lv = 2, 
     row.eff = TRUE, save.model = FALSE, calc.ics = TRUE)

spider.fit.pois$ics ## Returns information criteria

spider.fit.nb <- boral(y, family = "negative.binomial", num.lv = 2, 
     row.eff = TRUE, save.model = FALSE, calc.ics = TRUE)

spider.fit.nb$ics ## Returns the information criteria 
}
}
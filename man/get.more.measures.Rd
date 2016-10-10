\name{get.more.measures}
\alias{get.more.measures}
\title{Additional Information Criteria for boral models}
\description{Calculates some information criteria beyond those from \code{\link{get.measures}} for a boral model, although this set of criteria takes much longer to compute!}

\usage{
get.more.measures(y, X = NULL, family, trial.size = 1, 
	row.eff = "none", row.ids = NULL, num.lv, fit.mcmc, 
	verbose = TRUE)
}

\arguments{
  \item{y}{The response matrix that the boral model was fitted to.}
  
  \item{X}{The model matrix used in the boral model. Defaults to \code{NULL}, in which case it is assumed no model matrix was used.}  
  
  \item{family}{Either a single element, or a vector of length equal to the number of columns in \code{y}. The former assumes all columns of \code{y} come from this distribution. The latter option allows for different distributions for each column of \code{y}. Elements can be one of "binomial" (with probit link), "poisson" (with log link), "negative.binomial" (with log link), "normal" (with identity link), "lnormal" for lognormal (with log link), "tweedie" (with log link), "exponential" (with log link), "gamma" (with log link), "beta" (with logit link), "ordinal" (cumulative probit regression). 

  For the negative binomial distribution, the variance is parameterized as \eqn{Var(y) = \mu + \phi\mu^2}, where \eqn{\phi} is the column-specific dispersion parameter. For the normal distribution, the variance is parameterized as \eqn{Var(y) = \phi^2}, where \eqn{\phi} is the column-specific standard deviation. For the tweedie distribution, the variance is parameterized as \eqn{Var(y) = \phi \mu^p} where \eqn{\phi} is the column-specific dispersion parameter and \eqn{p} is a power parameter common to all columns assumed to be tweedie, with \eqn{1 < p < 2}. For the gamma distribution, the variance is parameterized as \eqn{Var(y) = \mu/\phi} where \eqn{\phi} is the column-specific rate (henceforth referred to also as dispersion parameter). For the beta distribution, the parameterization is in terms of the mean \eqn{\mu} and sample size \eqn{\phi} (henceforth referred to also as dispersion parameter), so that the two shape parameters are given by \eqn{a = \mu\phi} and \eqn{b = (1-\mu)\phi}.

 All columns assumed to have ordinal responses are constrained to have the same cutoffs points, with a column-specific intercept to account for differences between the columns (please see \emph{Details} for formulation). 
}

 \item{trial.size}{Either equal to a single element, or a vector of length equal to the number of columns in y. If a single element, then all columns assumed to be binomially distributed will have trial size set to this. If a vector, different trial sizes are allowed in each column of y. The argument is ignored for all columns not assumed to be binomially distributed. Defaults to 1, i.e. Bernoulli distribution.}

  \item{row.eff}{Single element indicating whether row effects are included as fixed effects ("fixed"), random effects ("random") or not included ("none") in the boral model. If random effects, they are drawn from a normal distribution with mean zero and unknown standard deviation. Defaults to "none". } 

 \item{row.ids}{A matrix with the number of rows equal to the number of rows in \code{y}, and the number of columns equal to the number of row effects to be included in the model. Element \eqn{(i,j)} indicates to the cluster ID of row \eqn{i} in \code{y} for random effect eqn{j}; please see the help file for the main boral function for details. Defaults to \code{NULL}, so that if \code{row.eff = "none"} then the argument is ignored, otherwise if \code{row.eff = "fixed"} or \code{"random"}, then \code{row.ids = matrix(1:nrow(y), ncol = 1)} i.e., a single, row effect unique to each row.}    

 \item{num.lv}{The number of latent variables used in the fitted boral model.}

  \item{fit.mcmc}{All MCMC samples for the fitted boral model, as obtained from JAGS. These can be extracted by fitting an boral model using \code{\link{boral}} with \code{save.model = TRUE}, and then accessing the \code{jags.model} component of the output.} 

  \item{verbose}{If TRUE, a notice is printed every 100 samples indicating progress in calculation of the marginal log-likelihood. Defaults to \code{TRUE}.}    
}

\details{
Currently, four information criteria are calculated using this function, when permitted: 1) AIC (using the marginal likelihood) evaluated at the posterior mode; 2) BIC (using the marginal likelihood) evaluated at the posterior mode; 3) Deviance information criterion (DIC) based on the marginal log-likelihood; 4) Widely Applicable Information Criterion (WAIC, Watanabe, 2010) based on the marginal log-likelihood. Since flat priors are used in fitting boral models, then the posterior mode should be approximately equal to the maximum likelihood estimates.

All four criteria require computing the marginal log-likelihood across all MCMC samples. This takes a very long time to run, since Monte Carlo integration needs to be performed for all MCMC samples. Consequently, this function is currently not implemented as an argument in main \code{\link{boral}} fitting function, unlike \code{\link{get.measures}} which is available via the \code{calc.ics = TRUE} argument.  

Moreover, note these criteria are not calculated all the time. In models where traits are included in the model (such that the regression coefficients \eqn{\beta_{0j}, \bm{\beta}_j} are random effects), more than two columns are ordinal responses (such that the intercepts \eqn{\beta_{0j}} for these columns are random effects), or more than one random row effect is included (such that \code{row.ids} contains more than one column), then these extra information criteria are will not calculated, and the function returns nothing except a simple message. This is because the calculation of the marginal log-likelihood in such cases currently fail to marginalize over such random effects; see the details in the help files for \code{calc.logLik.lv0} and \code{calc.marglogLik}. 

The two main differences between the criteria and those returned from \code{\link{get.measures}} are:
\itemize{
\item The AIC and BIC computed here are based on the log-likelihood evaluated at the posterior mode, whereas the AIC and BIC from \code{\link{get.measures}} are evaluated at the posterior median. The posterior mode and median will be quite close to one another if the component-wise posterior distributions are unimodal and symmetric. Furthermore, given uninformative priors are used, then both will be approximate maximum likelihood estimators. 
\item The DIC and WAIC computed here are based on the marginal log-likelihood, whereas the DIC and WAIC from \code{\link{get.measures}} are based on the conditional log-likelihood. Criteria based on the two types of log-likelihood are equally valid, and to a certain extent, which one to use depends on the question being answered i.e., whether to condition on the latent variables or treat them as "random effects" (see discussions in Spiegelhalter et al. 2002, and Vaida and Blanchard, 2005). 
}
}

\value{
If calculated, then a list with the following components:
\item{marg.aic}{AIC (using on the marginal log-likelihood) evaluated at posterior mode.}
\item{marg.bic}{BIC (using on the marginal log-likelihood) evaluated at posterior mode.}
\item{marg.dic}{DIC based on the marginal log-likelihood.}
\item{marg.waic}{WAIC based on the marginal log-likelihood.}
\item{all.marg.logLik}{The marginal log-likelihood evaluated at all MCMC samples. This is done via repeated application of \code{\link{calc.marglogLik}}.}
\item{num.params}{Number of estimated parameters used in the fitted model.}
}

\section{Warning}{
Using information criterion for variable selection should be done with extreme caution, for two reasons: 1) The implementation of these criteria are both \emph{heuristic} and experimental. 2) Deciding what model to fit for ordination purposes should be driven by the science. For example, it may be the case that a criterion suggests a model with 3 or 4 latent variables. However, if we interested in visualizing the data for ordination purposes, then models with 1 or 2 latent variables are far more appropriate. As an another example, whether or not we include row effects when ordinating multivariate abundance data depends on if we are interested in differences between sites in terms of relative species abundance (\code{row.eff = FALSE}) or in terms of species composition (\code{row.eff = "fixed"}).  

Also, the use of information criterion in the presence of variable selection using SSVS is questionable.
}

\references{
\itemize{
\item Spiegelhalter et al. (2002). Bayesian measures of model complexity and fit. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 64, 583-639.

\item Vaida, F., and Blanchard, S. (2005). Conditional Akaike information for mixed-effects models. Biometrika, 92, 351-370.

\item Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. The Journal of Machine Learning Research, 11, 3571-3594.
}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\seealso{\code{\link{get.measures}} for several information criteria which take considerably less time to compute, and are automatically implemented in \code{\link{boral}} with \code{calc.ics = TRUE}.}

\examples{
\dontrun{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
n <- nrow(y); p <- ncol(y); 
    
spider.fit.nb <- boral(y, family = "negative.binomial", num.lv = 2, 
     row.eff = "fixed", save.model = TRUE, calc.ics = TRUE)

## Extract MCMC samples
fit.mcmc <- mcmc(spider.fit.nb$jags.model$BUGSoutput$sims.matrix)

## WATCH OUT! The following takes a very long time to run!
get.more.measures(y, family = "negative.binomial", 
     num.lv = spider.fit.nb$num.lv, fit.mcmc = fit.mcmc, 
     row.eff = "fixed", row.ids = spider.fit.nb$row.ids)		

## Illustrating what happens in a case where these criteria will 
## 	not be calculated.
data(antTraits)
y <- antTraits$abun
X <- as.matrix(scale(antTraits$env))
## Include only traits 1, 2, and 5
traits <- as.matrix(antTraits$traits[,c(1,2,5)])
which.traits <- vector("list",ncol(X)+1)
for(i in 1:length(which.traits)) which.traits[[i]] <- 1:ncol(traits)

fit.traits <- boral(y, X = X, traits = traits, num.lv = 2, 
	which.traits = which.traits, family = "negative.binomial", 
	calc.ics = FALSE, save.model = TRUE)
     
## Extract MCMC samples
fit.mcmc <- mcmc(fit.traits$jags.model$BUGSoutput$sims.matrix)

get.more.measures(y, X = X, family = "negative.binomial", 
     num.lv = fit.traits$num.lv, fit.mcmc = fit.mcmc)		

}
}
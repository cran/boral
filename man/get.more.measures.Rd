\name{get.more.measures}
\alias{get.more.measures}
\title{Additional Information Criteria for boral models}
\description{Calculates some information criteria beyond those from \code{\link{get.measures}} for an boral model, although this set of criteria takes much longer to compute!}

\usage{
get.more.measures(y, X = NULL, family, trial.size = NULL, num.lv, fit.mcmc, 
	site.eff, verbose = TRUE)
}

\arguments{
  \item{y}{The response matrix that the boral model was fitted to.}
  
  \item{X}{The model matrix used in the boral model. Defaults to \code{NULL}, in which case it is assumed no model matrix was used.}  
  
\item{family}{Either a single element, or a vector of length equal to the number of columns in \eqn{y}. The former assumes all columns of \eqn{y} come from this distribution. The latter option allows for different distributions for each column of \eqn{y}. Elements can be one of "binomial" (with logit link), "poisson" (with log link), "negative.binomial" (with log link), "normal" (with identity link), "lnormal" for lognormal (with log link), "tweedie" (with log link), "ordinal" which performs cumulative logit regression. For the negative binomial distribution, the variance is parameterized as \eqn{Var(y) = \mu + \phi \mu^2}, where \eqn{\phi} is the column-specific overdispersion parameter. For the normal distribution, the variance is parameterized as \eqn{Var(y) = \phi}, where \eqn{\phi} is the column-specific variance. For the tweedie distribution, the variance is parameterized as \eqn{Var(y) = \phi \mu^p} where \eqn{\phi} is the column-specific dispersion parameter and \eqn{p} is a power parameter common to all columns assumed to be tweedie, with \eqn{1 < p < 2}. 

All columns assumed to have ordinal responses are constrained to have the same cutoffs points, with a column-specific intercept to account for differences between the columns (please see \emph{Details} for formulation). 
}
  
 \item{trial.size}{Either equal to NULL, a single element, or a vector of length equal to the number of columns in \eqn{y}. If a single element, then all columns assumed to be binomially distributed will have trial size set to this. If a vector, different trial sizes are allowed in each column of \eqn{y}.}

  \item{num.lv}{The number of latent variables used in the fitted boral model.}

  \item{fit.mcmc}{All MCMC samples for the fitted boral model, as obtained from JAGS. These can be extracted by fitting an boral model using \code{\link{boral}} with \code{save.model = TRUE}, and then applying \code{as.mcmc} on the \code{jags.model} component of the output.} 

  \item{site.eff}{A logical value indicating whether to include row effects were included in the model.}    

  \item{verbose}{If \code{TRUE}, a notice is printed every 100 samples indicating progress in calculation of the marginal log-likelihood. Defaults to \code{TRUE}.}    
}

\details{
Currently, four information criteria has been implemented in this function: 1) AIC (using the marginal likelihood) evaluated at the posterior mode; 2) BIC (using the marginal likelihood) evalulated at the posterior mode; 3) Deviance information criterion (DIC) based on the marginal log-likelihood; 4) Widely Applicable Information Criterion (WAIC, Watanabe, 2010) based on the marginal log-likelihood. Since flat priors are used in fitting boral models, then the posterior mode should be approximately equal to the maximum likelihood estimates.

All four criteria require computing the marginal log-likelihood across all MCMC samples. This takes a very long time to run, since Monte Carlo integration needs to be performed for all MCMC samples. Consequently, this function is currently not implemented as an argument in main \code{\link{boral}} fitting function, unlike \code{\link{get.measures}} which is available via the \code{calc.ics = TRUE} argument.  

The two main differences between the criteria and those returned from \code{\link{get.measures}} are:
\itemize{
\item The AIC and BIC computed here are based on the log-likelihood evalulated at the posterior mode, whereas the AIC and BIC from \code{\link{get.measures}} are evalulated at the posterior median. The posterior mode and median will be quite close to one another if the component-wise posterior distributions are unimodal and symmetric. Furthermore, given uninformative priors are used, then both will be approximate maximum likelihood estimators. 
\item The DIC and WAIC computed here are based on the marginal log-likelihood, whereas the DIC and WAIC from \code{\link{get.measures}} are based on the conditional log-likelihood. Criteria based on the two types of log-likelihood are equally valid, and to a certain extent, which one to use depends on the question being answered i.e., whether to condition on the latent variables or treat them as "random effects" (see discussions in Spiegelhalter et al. 2002, and Vaida and Blanchard, 2005). Having said that, there is evidence to suggests, for models with overdispersed count data, conditional DIC/WAIC may not perform as well as than marginal DIC/WAIC for overdispered abundance data (Millar, 2009).
}

In our very limited experience, we found BIC evalulated at the posterior mode tends to be quite stable, whereas marginal DIC and WAIC tend to overfit the number of latent variables. 
}

\value{
A list with the following components:
\item{marg.aic}{AIC (using on the marginal log-likelihood) evaluated at posterior mode.}
\item{marg.bic}{BIC (using on the marginal log-likelihood) evaluated at posterior mode.}
\item{marg.dic}{DIC based on the marginal log-likelihood.}
\item{marg.waic}{WAIC based on the marginal log-likelihood.}
\item{all.marg.logLik}{The marginal log-likelihood evaluated at all MCMC samples. This is done via repeated application of \code{\link{calc.marglogLik}}.}
\item{num.params}{Number of estimated parameters used in the fitted model.}
}

\section{Warning}{
Using information criteria for variable selection should be done with extreme caution, for two reasons: 1) The implementation of these criteria are both \emph{heuristic} and experimental. 2) Deciding what model to fit for ordination purposes should be driven by the science. For example, it may be the case that criteria suggests a model with 3 or 4 latent variables. However, if we interested in visualizing the data for ordination purposes, then models with 1 or 2 latent variables are far more appropriate. As an another example, whether or not we include row (site) effects when ordinating multivariate abundance data depends on if we are interested in differences between sites in terms of abundance (\code{site.eff = FALSE}) or in terms of species composition (\code{site.eff = TRUE}).  
}

\references{
\itemize{
\item Millar, R. B. (2009). Comparison of hierarchical Bayesian models for overdispersed count data using DIC and Bayes' factors. Biometrics, 65, 962-969.
\item Spiegelhalter, D. J., Best, N. G., Carlin, B. P., & Van Der Linde, A. (2002). Bayesian measures of model complexity and fit. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 64, 583-639.
\item Vaida, F., & Blanchard, S. (2005). Conditional Akaike information for mixed-effects models. Biometrika, 92, 351-370.
\item Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. The Journal of Machine Learning Research, 11, 3571-3594.
}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\seealso{\code{\link{get.measures}} for several information criteria which take less time to compute, and are automatically implemented in \code{\link{boral}} with \code{calc.ics = TRUE}.}

\examples{
\dontrun{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
n <- nrow(y); p <- ncol(y); 
    
spider.fit.nb <- boral(y, family = "negative.binomial", num.lv = 2, 
     site.eff = TRUE, save.model = TRUE, calc.ics = TRUE)

all.mcmc <- as.mcmc(spider.fit.nb$jags.model)

## WATCH OUT! The following takes a very long time to run!
get.more.measures(y, family = "negative.binomial", 
     num.lv = 2, fit.mcmc = all.mcmc, site.eff = TRUE)		
}
}
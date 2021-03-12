\name{get.more.measures}

\alias{get.more.measures}

\docType{package}

\title{Additional Information Criteria for models}

\description{
\Sexpr[results=rd, stage=render]{lifecycle::badge("defunct")}

Calculates some information criteria beyond those from \code{\link{get.measures}} for a fitted model, although this set of criteria takes much longer to compute! WARNING: As of version 1.6, this function is no longer maintained (and probably doesn't work properly, if at all)!}

\usage{
get.more.measures(y, X = NULL, family, trial.size = 1, 
	row.eff = "none", row.ids = NULL, offset = NULL,
	num.lv, fit.mcmc, verbose = TRUE)
}

\arguments{
  \item{y}{The response matrix that the model was fitted to.}
  
  \item{X}{The covariate matrix used in the model. Defaults to \code{NULL}, in which case it is assumed no model matrix was used.}  
  
  \item{family}{Either a single element, or a vector of length equal to the number of columns in the response matrix. The former assumes all columns of the response matrix come from this distribution. The latter option allows for different distributions for each column of the response matrix. Elements can be one of "binomial" (with probit link), "poisson" (with log link), "negative.binomial" (with log link), "normal" (with identity link), "lnormal" for lognormal (with log link), "tweedie" (with log link), "exponential" (with log link), "gamma" (with log link), "beta" (with logit link), "ordinal" (cumulative probit regression), "ztpoisson" (zero truncated Poisson with log link), "ztnegative.binomial" (zero truncated negative binomial with log link). 

  Please see \code{\link{about.distributions}} for information on distributions available in boral overall.
  }

 \item{trial.size}{Either equal to a single element, or a vector of length equal to the number of columns in y. If a single element, then all columns assumed to be binomially distributed will have trial size set to this. If a vector, different trial sizes are allowed in each column of y. The argument is ignored for all columns not assumed to be binomially distributed. Defaults to 1, i.e. Bernoulli distribution.}

  \item{row.eff}{Single element indicating whether row effects are included as fixed effects ("fixed"), random effects ("random") or not included ("none") in the fitted model. If fixed effects, then for parameter identifiability the first row effect is set to zero, which analogous to acting as a reference level when dummy variables are used. If random effects, they are drawn from a normal distribution with mean zero and estimated standard deviation. Defaults to "none". } 

 \item{row.ids}{A matrix with the number of rows equal to the number of rows in the response matrix, and the number of columns equal to the number of row effects to be included in the model. Element \eqn{(i,j)} indicates the cluster ID of row \eqn{i} in the response matrix for random effect eqn{j}; please see \code{\link{boral}} for details. Defaults to \code{NULL}, so that if \code{row.eff = "none"} then the argument is ignored, otherwise if \cr \code{row.eff = "fixed"} or \code{"random"}, \cr then \code{row.ids = matrix(1:nrow(y), ncol = 1)} i.e., a single, row effect unique to each row.}    

 \item{offset}{A matrix with the same dimensions as the response matrix, specifying an a-priori known component to be included in the linear predictor during fitting. Defaults to \code{NULL}.}

 \item{num.lv}{The number of latent variables used in the fitted model.}

  \item{fit.mcmc}{All MCMC samples for the fitted model. These can be extracted by fitting a model using \code{\link{boral}} with \code{save.model = TRUE}, and then applying \code{get.mcmcsamples(fit)}.} 

  \item{verbose}{If TRUE, a notice is printed every 100 samples indicating progress in calculation of the marginal log-likelihood. Defaults to \code{TRUE}.}    
}

\details{
Currently, four information criteria are calculated using this function, when permitted: 1) AIC (using the marginal likelihood) evaluated at the posterior mode; 2) BIC (using the marginal likelihood) evaluated at the posterior mode; 3) Deviance information criterion (DIC) based on the marginal log-likelihood; 4) Widely Applicable Information Criterion (WAIC, Watanabe, 2010) based on the marginal log-likelihood. When uninformative priors are used in fitting models, then the posterior mode should be approximately equal to the maximum likelihood estimates.

All four criteria require computing the marginal log-likelihood across all MCMC samples. This takes a very long time to run, since Monte Carlo integration needs to be performed for all MCMC samples. Consequently, this function is currently not implemented as an argument in main \code{\link{boral}} fitting function, unlike \code{\link{get.measures}} which is available via the \code{calc.ics = TRUE} argument.  

Moreover, note these criteria are not calculated all the time. In models where traits are included in the model (such that the regression coefficients \eqn{\beta_{0j}, \bm{\beta}_j} are random effects), or more than two columns are ordinal responses (such that the intercepts \eqn{\beta_{0j}} for these columns are random effects), then these extra information criteria are will not calculated, and the function returns nothing except a simple message. This is because the calculation of the marginal log-likelihood in such cases currently fail to marginalize over such random effects; please see the details in \code{calc.logLik.lv0} and \code{calc.marglogLik}. 

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
As of version 1.6, this function is no longer maintained (and probably doesn't work)!

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
\packageAuthor{boral}

Maintainer: \packageMaintainer{boral}
}

\seealso{\code{\link{get.measures}} for several information criteria which take considerably less time to compute, and are automatically implemented in \code{\link{boral}} with \code{calc.ics = TRUE}.}

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

spiderfit_nb <- boral(y, family = "negative.binomial", lv.control = list(num.lv = 2),
    row.eff = "fixed", save.model = TRUE, calc.ics = TRUE,
    mcmc.control = example_mcmc_control, model.name = testpath)

## Extract MCMC samples
fit_mcmc <- get.mcmcsamples(spiderfit_nb)

## NOTE: The following takes a long time to run!
get.more.measures(y, family = "negative.binomial", 
    num.lv = spiderfit_nb$num.lv, fit.mcmc = fit_mcmc, 
    row.eff = "fixed", row.ids = spiderfit_nb$row.ids)		

     
## Illustrating what happens in a case where these criteria will 
## 	not be calculated.
data(antTraits)
y <- antTraits$abun
X <- as.matrix(scale(antTraits$env))
## Include only traits 1, 2, and 5
traits <- as.matrix(antTraits$traits[,c(1,2,5)])
example_which_traits <- vector("list",ncol(X)+1)
for(i in 1:length(example_which_traits))
    example_which_traits[[i]] <- 1:ncol(traits)

fit_traits <- boral(y, X = X, traits = traits, lv.control = list(num.lv = 2), 
    which.traits = example_which_traits, family = "negative.binomial", 
    save.model = TRUE, mcmc.control = example_mcmc_control,
    model.name = testpath)
     
## Extract MCMC samples
fit_mcmc <- get.mcmcsamples(fit_traits)

get.more.measures(y, X = X, family = "negative.binomial", 
    num.lv = fit_traits$num.lv, fit.mcmc = fit_mcmc)	
}
}

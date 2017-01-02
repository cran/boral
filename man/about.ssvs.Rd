\name{about.ssvs}

\alias{about.ssvs}

\title{Stochastic search variable selection (SSVS) in boral}

\description{This help file provides more information regarding the implementation of the stochastic search variable selection (SSVS, George and McCulloch, 1993) as implemented in the boral package.}

\details{
An alternative to using information criterion for model selection, which is applicable specifically to the Bayesian MCMC framework, is stochastic search variable selection (SSVS, George and McCulloch, 1993). In boral, SSVS is implemented on the column-specific coefficients \eqn{\bm{\beta}_j}. Basically, SSVS works by placing a spike-and-slab priors on these coefficients, such that the spike is a narrow normal distribution concentrated around zero and the spike is a normal distribution with a large variance.

\deqn{\rho(\beta) = I_{\beta = 1}\times\mathcal{N}(0,\sigma^2) + (1-I_{\beta = 1})\times \mathcal{N}(0,g*\sigma^2),}

where \eqn{\sigma^2} is determined by \code{prior.control$hypparams[3]}, \eqn{g} is determined by \code{ssvs.g}, and \eqn{I_{\beta = 1} = P(\beta = 1)} is an indicator function representing whether coefficient is included in the model. It is given a Bernoulli prior with probability of inclusion 0.5. After fitting, the posterior probability of \eqn{\beta} being included in the model is returned based on posterior mean of the indicator function \eqn{I_{\beta = 1}}. Note this is NOT the same as a \emph{p}-value seen in maximum likelihood estimation: a \emph{p}-value provides an indication of how much evidence there is against the null hypothesis of \eqn{\beta = 0}, while the posterior probability provides a measure of how likely it is for \eqn{\beta \ne 0} given the data.

In boral, SSVS can be applied at a grouped or individual coefficient level, and this is governed by \code{prior.control$ssvs.index}: 
\itemize{
\item For elements of \code{ssvs.index} equal to -1, SSVS is not applied on the corresponding covariate of the model matrix \code{X}. 

\item For elements equal to 0, SSVS is applied to each individual coefficients of the corresponding covariate in \code{X}. That is, the fitted model will return posterior probabilities for this covariate, one for each column of \code{y}.

\item For elements taking positive integers {1,2,...}, SSVS is applied to each group of coefficients of the corresponding covariate in \code{X}. That is, the fitted model will return a single posterior probability for this covariate, indicating whether this covariate should be included for all columns of \code{y}; see O'Hara and Sillanpaa (2009) and Tenan et al. (2014) among many others for an discussion of Bayesian variable selection methods.
}


Note the last application of SSVS allows multiple covariates to be selected \emph{simultaneously}. For example, suppose \code{X} consists of five columns: the first two columns are environmental covariates, while the last three correspond to quadratic terms of the two covariates as well as their interaction. If we want to "test" whether any quadratic terms are required, then we can set \cr \code{prior.control$ssvs.index = c(-1,-1,1,1,1)}, so a single posterior probability of inclusion is returned for the last three columns of \code{X}. 

Finally, note using information criterion (and possibly residual analysis) should probably not be done at the same as when SSVS is used, and it is advised to separate out their applications e.g., choose the explanatory variables first using SSVS, and then use information criterion to select the number of latent variables? Obtaining summaries such as posterior medians and HPD intervals of the coefficients from a boral model that is implementing SSVS is also (perhaps) problematic, because the posterior distribution is multi-modal. 
}


\references{
\itemize{
\item George, E. I. and McCulloch, R. E. (1993). Variable selection via Gibbs sampling. Journal of the American Statistical Association, 85, 398-409.

\item O'Hara, B., and Sillianpaa, M.J. (2009). A Review of Bayesian Variable Selection Methods: What, How and Which. Bayesian Analysis, 4, 85-118.

\item Tenan, S., et al. (2014). Bayesian model selection: The steepest mountain to climb. Ecological Modelling, 283, 62-69.
	}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\section{Warnings}{
\itemize{
\item Summaries of the coefficients such as posterior medians and HPD intervals may also be problematic when SSVS is being used, since the posterior distribution will be multi-modal. 

\item If \code{save.model = TRUE}, the raw jags model is also returned. This can be quite very memory-consuming, since it indirectly saves all the MCMC samples.
}
}

\seealso{
\code{\link{boral}} for the main boral fitting function which implemented SSVS. 
}

\examples{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
X <- scale(spider$x)
n <- nrow(y); p <- ncol(y); 

## NOTE: The two examples below and taken directly from the boral help file

example.mcmc.control <- list(n.burnin = 10, n.iteration = 100, 
     n.thin = 1)

\dontrun{
## Example 3a - Extend example 2 to demonstrate grouped covariate selection
## on the last three covariates. 
set.prior <- list(type = c("normal","normal","normal","uniform"), 
	hypparams = c(100, 20, 100, 50), ssvs.index = c(-1,-1,-1,1,2,3))
spider.fit.nb2 <- boral(y, X = X, family = "negative.binomial", 
	num.lv = 0, calc.ics = FALSE, mcmc.control = example.mcmc.control,
	prior.control = set.prior)
     
summary(spider.fit.nb2) 


## Example 3b - Extend example 2 to demonstrate individual covariate selection
## on the last three covariates. 
set.prior <- list(type = c("normal","normal","normal","uniform"), 
	hypparams = c(100, 20, 100, 50), ssvs.index = c(-1,-1,-1,0,0,0))
spider.fit.nb3 <- boral(y, X = X, family = "negative.binomial", 
	num.lv = 0, calc.ics = FALSE, mcmc.control = example.mcmc.control,
	prior.control = set.prior)
summary(spider.fit.nb3) 
}

}
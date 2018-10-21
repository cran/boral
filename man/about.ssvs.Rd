\name{about.ssvs}

\alias{about.ssvs}

\docType{package}

\title{
Stochastic search variable selection (SSVS) in boral
}

\description{This help file provides more information regarding the implementation of the stochastic search variable selection (SSVS, George and McCulloch, 1993) as implemented in the boral package.}

\details{
Stochastic search variable selection (SSVS, George and McCulloch, 1993) is a approach for model selection, which is applicable specifically to the Bayesian MCMC framework. As of boral version 1.5, SSVS is implemented in two ways.

\bold{SSVS on coefficients in \code{X}:}
SSVS is implemented on the column-specific coefficients \eqn{\bm{\beta}_j}. Basically, SSVS works by placing a spike-and-slab priors on these coefficients, such that the spike is a narrow normal distribution concentrated around zero and the spike is a normal distribution with a large variance.

\deqn{\rho(\beta) = I_{\beta = 1}\times\mathcal{N}(0,\sigma^2) + (1-I_{\beta = 1})\times \mathcal{N}(0,g*\sigma^2),}

where \eqn{\sigma^2} is determined by \code{prior.control$hypparams[3]}, \eqn{g} is determined by \code{ssvs.g}, and \eqn{I_{\beta = 1} = P(\beta = 1)} is an indicator function representing whether coefficient is included in the model. It is given a Bernoulli prior with probability of inclusion 0.5. After fitting, the posterior probability of \eqn{\beta} being included in the model is returned based on posterior mean of the indicator function \eqn{I_{\beta = 1}}. Note this is NOT the same as a \emph{p}-value seen in maximum likelihood estimation: a \emph{p}-value provides an indication of how much evidence there is against the null hypothesis of \eqn{\beta = 0}, while the posterior probability provides a measure of how likely it is for \eqn{\beta \ne 0} given the data.

SSVS can be applied at a grouped or individual coefficient level, and this is governed by \cr \code{prior.control$ssvs.index}: 
\itemize{
\item For elements of \code{ssvs.index} equal to -1, SSVS is not applied on the corresponding covariate of the \code{X}. 

\item For elements equal to 0, SSVS is applied to each individual coefficients of the corresponding covariate in \code{X}. That is, the fitted model will return posterior probabilities for this covariate, one for each column of \code{y}.

\item For elements taking positive integers {1,2,...}, SSVS is applied to each group of coefficients of the corresponding covariate in \code{X}. That is, the fitted model will return a single posterior probability for this covariate, indicating whether this covariate should be included for all columns of \code{y}; see O'Hara and Sillanpaa (2009) and Tenan et al. (2014) among many others for an discussion of Bayesian variable selection methods.
}


Note the last application of SSVS allows multiple covariates to be selected \emph{simultaneously}. For example, suppose \code{X} consists of five columns: the first two columns are environmental covariates, while the last three correspond to quadratic terms of the two covariates as well as their interaction. If we want to "test" whether any quadratic terms are required, then we can set \cr \code{prior.control$ssvs.index = c(-1,-1,1,1,1)}, so a single posterior probability of inclusion is returned for the last three columns of \code{X}. 

Finally, note that summaries such as posterior medians and HPD intervals of the coefficients, as well as performing residual analysis, from a boral model that has implemented SSVS may be problematic because the posterior distribution is by definition multi-modal. It may be advisable instead to separate out their application of SSVS and posterior inference.

\bold{SSVS on trait coefficients:}
If traits are included in boral, thereby leading to a fourth corner model (see \code{\link{about.traits}} for more details on this type of model), SSVS can also be performed on the associated trait coefficients. That is, in such model we have

\deqn{\beta_{0j} \sim N(\kappa_{01} + \bm{traits}^\top_j\bm{\kappa}_1, \sigma^2_1)}

for the column-specific intercepts, and 

\deqn{\beta_{jk} \sim N(\kappa_{0k} + \bm{traits}^\top_j\bm{\kappa}_k, \sigma^2_k)}

for \eqn{k = 1,\ldots,d} where \code{d = ncol(X)}. Then if the a particular index in the argument \cr \code{prior.control$ssvs.traitsindex} is set to 0, SSVS is performed on the corresponding element in \eqn{\bm{\kappa}_1} or \eqn{\bm{\kappa}_k}. For example, suppose \code{which.traits[[2]] = c(2,3)}, meaning that the \eqn{\beta_{j1}}'s are drawn from a normal distribution with mean depending only on the second and third columns of \code{traits}. Then \cr \code{prior.control$ssvs.traitsindex[[2]] = c(0,1)}, then a spike-and-slab prior is placed on the first coefficent in \eqn{\bm{\kappa}_2}, while the second coefficient is assigned the ``standard" prior governed by the \code{prior.control$hypparams}. That is, SSVS is performed on the first but not the second coefficient in \eqn{\bm{\kappa}_2}. 

Please keep in mind that because boral allows the user to manually decide which traits drive which covariates in \code{X}, then care must be taken when setting up both \code{which.traits} and \cr \code{prior.control$ssvs.traitsindex}. That is, when supplied then both objects should be lists of have the same length, and the length of the corresponding vectors comprising each element in the two lists should match as well e.g., \code{which.traits[[2]]} and \cr \code{prior.control$ssvs.traitsindex[[2]]} should be of the same length.

}


\references{
\itemize{
\item George, E. I. and McCulloch, R. E. (1993). Variable selection via Gibbs sampling. Journal of the American Statistical Association, 85, 398-409.

\item O'Hara, B., and Sillianpaa, M.J. (2009). A Review of Bayesian Variable Selection Methods: What, How and Which. Bayesian Analysis, 4, 85-118.

\item Tenan, S., et al. (2014). Bayesian model selection: The steepest mountain to climb. Ecological Modelling, 283, 62-69.
	}
}

\author{
\packageAuthor{boral}

Maintainer: \packageMaintainer{boral}
}

\section{Warnings}{
\itemize{
\item Summaries of the coefficients such as posterior medians and HPD intervals may also be problematic when SSVS is being used, since the posterior distribution will be multi-modal. 

\item If \code{save.model = TRUE}, the raw jags model is also returned. This can be quite very memory-consuming, since it indirectly saves all the MCMC samples.
}
}

\seealso{
\code{\link{boral}} for the main boral fitting function which implementing SSVS, and \code{\link{about.traits}} for how fourth corner models work before applying SSVS to them.
}

\examples{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
X <- scale(spider$x)
n <- nrow(y)
p <- ncol(y)

## NOTE: The two examples below and taken directly from the boral help file

example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
     n.thin = 1)

\dontrun{
## Example 3a - Extend example 2 to demonstrate grouped covariate selection
## on the last three covariates. 
example_prior_control <- list(type = c("normal","normal","normal","uniform"), 
     ssvs.index = c(-1,-1,-1,1,2,3))
spiderfit_nb2 <- boral(y, X = X, family = "negative.binomial", 
    mcmc.control = example_mcmc_control, prior.control = example_prior_control)
     
summary(spiderfit_nb2) 


## Example 3b - Extend example 2 to demonstrate individual covariate selection
## on the last three covariates. 
example_prior_control <- list(type = c("normal","normal","normal","uniform"), 
     ssvs.index = c(-1,-1,-1,0,0,0))
spiderfit_nb3 <- boral(y, X = X, family = "negative.binomial", 
    mcmc.control = example_mcmc_control, prior.control = example_prior_control)
summary(spiderfit_nb3) 


## Example 5a - model fitted to count data, no site effects, and
## two latent variables, plus traits included to explain environmental responses
data(antTraits)
y <- antTraits$abun
X <- as.matrix(scale(antTraits$env))
## Include only traits 1, 2, and 5
traits <- as.matrix(antTraits$traits[,c(1,2,5)])
example_which_traits <- vector("list",ncol(X)+1)
for(i in 1:length(example_which_traits)) 
     example_which_traits[[i]] <- 1:ncol(traits)
## Just for fun, the regression coefficients for the second column of X,
## corresponding to the third element in the list example_which_traits,
## will be estimated separately and not regressed against traits.
example_which_traits[[3]] <- 0

fit_traits <- boral(y, X = X, traits = traits, 
    which.traits = example_which_traits, family = "negative.binomial", 
    mcmc.control = example_mcmc_control, save.model = TRUE)

summary(fit_traits)


## Example 5b - perform selection on trait coefficients
ssvs_traitsindex <- vector("list",ncol(X)+1)
for(i in 1:length(ssvs_traitsindex)) 
     ssvs_traitsindex[[i]] <- rep(0,ncol(traits))
ssvs_traitsindex[[3]] <- -1
fit_traits <- boral(y, X = X, traits = traits, which.traits = example_which_traits, 
    family = "negative.binomial", mcmc.control = example_mcmc_control, 
    save.model = TRUE, prior.control = list(ssvs.traitsindex = ssvs_traitsindex))

summary(fit_traits)
}

}
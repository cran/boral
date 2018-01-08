\name{predict.boral}

\alias{predict.boral}

\title{Predict using a boral model}

\description{Obtain predictions and associated intervals (lower and upper limits) on the linear scale from a fitted boral object. Predictions can be made either conditionally on the predicted latent variables and any random row effects included in the model, or marginally (averaged) on the latent variables and any random effects included in the model.}

\usage{
\method{predict}{boral}(object, newX = NULL, newrow.ids = NULL, 
     predict.type = "conditional", est = "median", prob = 0.95, 
     lv.mc = 1000, ...)
     }

\arguments{
  \item{object}{An object of class "boral".}
  
  \item{newX}{An optional model matrix of covariates for extrapolation to the same sites (under different environmental conditions) or extrapolation to new sites. No intercept column should be included in \code{newX}. Defaults to \code{NULL}, in which case the model matrix of covariates is taken from the fitted boral object if found.}  

  \item{newrow.ids}{An optional matrix with the number of columns equal to the number of row effects to be included in the model. Element \eqn{(i,j)} indicates to the cluster ID of row \eqn{i} in \code{y} for random effect eqn{j}. Defaults to \code{NULL}, in which case row IDs are taken from the fitted boral object itself (if required) i.e., from \code{object$row.ids}.}  

  \item{predict.type}{The type of prediction to be made. Either takes value \code{"conditional"} in which case the prediction is made conditionally on the predicted latent variables and any random row effects in the model, or \code{"marginal"} in which case the prediction marginalizes (averages) over the latent variables and random row effects in the model. Defaults to \code{"conditional"}.}  

  \item{est}{A choice of either whether to print the posterior median (\code{est == "median"}) or posterior mean (\code{est == "mean"}) of the parameters.}

  \item{prob}{A numeric scalar in the interval (0,1) giving the target probability coverage of the intervals. Defaults to 0.95.}   

  \item{lv.mc}{If the predictions are made marginalizing over the latent variables, then number of Monte-Carlo samples to take when performing the relevant integration.}

     \item{...}{Not used.}
}


\details{
Due to the Bayesian MCMC framework, then predictive inference for boral models is based around the posterior predictive distribution, which is the integral of the quantity one wants to predict on, integrated or averaged over the posterior distribution of the parameters and latent variables. Currently, all predictions are made on the \emph{linear predictor scale} (although this might change in future updates) i.e., 

\deqn{\eta_{ij} = \alpha_i + \beta_{0j} + \bm{x}^\top_i\bm{\beta}_j + \bm{z}^\top_i\bm{\theta}_j; \quad i = 1,\ldots,n; j = 1,\ldots,p,}

where \eqn{\bm{z}_i} are a vector of latent variables included in the model, \eqn{\bm{\theta}_j} are the column-specific coefficients relating to these latent variables, \eqn{\bm{x}_i} are covariates included in the model, and \eqn{\bm{beta}_j} being the column-specific coefficients related to these covariates. The quantity \eqn{beta_{0j}} denotes the column-specific intercepts while \code{alpha_i} represents one or more optional row effects that may be treated as a fixed or random effect. 

Note that for the above to work, one must have saved the MCMC samples in the fitted boral object, that is, set \code{save.model = TRUE} when fitting. 

Two types of predictions are possible using this function:
\itemize{
     \item The first type is \code{predict.type = "conditional"}, meaning predictions are made conditionally on the predicted latent variables and any (random) row effects in the model. This is mainly used when predictions are made onto the \emph{same} set of sites that the boral model was fitted to, although a \code{newX} can be supplied in this case if we want to extrapolate on to the same set of sites but under different environmental conditions. 

     \item The second type of prediction is \code{predict.type = "marginal"}, meaning predictions are made marginally or averaging over the latent variables and any (random) row effects in the model. This is mainly used when predictions are made onto a \emph{new} set of sites where the latent variables and/or row effects are unknown. A \code{newX} and/or \code{newrow.ids} is often supplied since we are extrapolating to new sites. The integration over the latent variables and random row effects is done via Monte-Carlo integration. Please note however that, as mentioned before, the integration will be done on the linear predictor scale. 
     }

More information on conditional versus marginal predictions in latent variable models can be found in Warton et al., (2015). In both cases, the function returns a point prediction (either the posterior mean or median depending on \code{est}) and the lower and upper bounds of a 100\eqn{\alpha\%} interval of the posterior prediction. All of these quantities are calculated empirically based the MCMC samples e.g., the posterior mean is the average of the predictions across the MCMC samples, and the lower and upper bounds are based on quantiles. 
}

\value{
A list containing the following components:
     \item{linpred}{A matrix containing posterior point predictions (either posterior mean or median depending on \code{est}), on the linear predictor scale.}
     
     \item{lower}{A matrix containing the lower bound of the 100\code{alpha}\% interval of the posterior predictions, on the linear predictor scale.}

     \item{upper}{A matrix containing the upper bound of the 100\code{alpha}\% interval of the posterior predictions, on the linear predictor scale.}
}

\references{
\itemize{
\item Gelman et al. (2013) Bayesian Data Analysis. CRC Press.

\item Warton et al. (2015). So Many Variables: Joint Modeling in Community Ecology. Trends in Ecology and Evolution, 30, 766-779.
}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\section{Warnings}{
\itemize{
\item Marginal predictions can take quite a while to construct due to the need to perform Monte-Carlo integration to marginalize over the latent variables and any random row effects in the model.
}
}

\examples{
\dontrun{
library(mvabund) ## Load a dataset from the mvabund package
library(mvtnorm) 
data(spider)
y <- spider$abun
X <- scale(spider$x)
n <- nrow(y)
p <- ncol(y)

## NOTE: The values below MUST NOT be used in a real application;
## they are only used here to make the examples run quick!!!
example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
     n.thin = 1)

## Example 1 - model with two latent variables, random site effects, 
## 	and environmental covariates
spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
     row.eff = "random", num.lv = 2, 
     mcmc.control = example_mcmc_control, save.model = TRUE)


## Predictions conditional on predicted latent variables
getcondpreds <- predict(spiderfit_nb)

## Predictions marginal on latent variables, random row effects
## The intervals for these will generally be wider than the
##   conditional intervals.
getmargpreds <- predict(spiderfit_nb, predict.type = "marginal")


## Now suppose you extrpolate to new sites
newX <- rmvnorm(100, mean = rep(0,ncol(X)))

## Below won't work since conditional predictions are made to the same sites
getcondpreds <- predict(spiderfit_nb, newX = newX)

## Marginal predictions will work though provided newrow.ids is set up 
## properly. For example,
new_row_ids <- matrix(sample(1:28,100,replace=TRUE), 100, 1)
getmargpreds <- predict(spiderfit_nb, newX = newX, predict.type = "marginal", 
     newrow.ids = new_row_ids)

     
## Example 2 - simulate count data, based on a model with two latent variables, 
## no site variables, with two traits and one environmental covariates 
library(mvtnorm)

n <- 100; s <- 50
X <- as.matrix(scale(1:n))
colnames(X) <- c("elevation")

traits <- cbind(rbinom(s,1,0.5), rnorm(s)) 
## one categorical and one continuous variable
colnames(traits) <- c("thorns-dummy","SLA")

simfit <- list(true.lv = rmvnorm(n, mean = rep(0,2)), 
	lv.coefs = cbind(rnorm(s), rmvnorm(s, mean = rep(0,2)), 1), 
	traits.coefs = matrix(c(0.1,1,-0.5,0.1,0.5,0,-1,0.1), 2, byrow = TRUE))
rownames(simfit$traits.coefs) <- c("beta0","elevation")
colnames(simfit$traits.coefs) <- c("kappa0","thorns-dummy","SLA","sigma")

simy = create.life(true.lv = simfit$true.lv, lv.coefs = simfit$lv.coefs, X = X, 
	traits = traits, traits.coefs = simfit$traits.coefs, family = "normal") 


example_which_traits <- vector("list",ncol(X)+1)
for(i in 1:length(example_which_traits)) 
     example_which_traits[[i]] <- 1:ncol(traits)
fit_traits <- boral(y = simy, X = X, traits = traits, 
     which.traits = example_which_traits, family = "normal", 
     num.lv = 2, save.model = TRUE, mcmc.control = example_mcmc_control)	

     
## Predictions conditional on predicted latent variables   
getcondpreds <- predict(fit_traits)     
     
## Predictions marginal on latent variables
## The intervals for these will generally be wider than the
##   conditional intervals.
getmargpreds <- predict(fit_traits, predict.type = "marginal")
}

}
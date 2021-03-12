\name{predict.boral}

\alias{predict.boral}

\docType{package}

\title{Predict using a model}

\description{
\Sexpr[results=rd, stage=render]{lifecycle::badge("stable")}

Construct predictions and associated intervals (lower and upper limits) from a fitted boral object. Predictions can be made either conditionally on the predicted latent variables and any random row effects/response-specific random intercepts included in the model, or marginally (averaged) on the latent variables and these other effects (note integration is done on the linear predictor scale).}

\usage{
\method{predict}{boral}(object, newX = NULL, newrow.ids = NULL, newranef.ids = NULL,
     distmat =  NULL, predict.type = "conditional", scale = "link", 
     est = "median", prob = 0.95, lv.mc = 1000, return.alllinpred = FALSE, 
     ...)
     }

\arguments{
     \item{object}{An object of class "boral".}
  
     \item{newX}{An optional model matrix of covariates for extrapolation to the same sites (under different environmental conditions) or extrapolation to new sites. No intercept column should be included in \code{newX}. Defaults to \code{NULL}, in which case the model matrix of covariates is taken from the fitted boral object if found.}  

     \item{newrow.ids}{An optional matrix with the number of columns equal to the number of row effects included in the model. Element \eqn{(i,j)} indicates the cluster ID of row \eqn{i} in the response matrix for random effect \eqn{j}. Defaults to \code{NULL}, in which case row IDs are taken from the fitted boral object itself, if appropriate, i.e., from \code{object$row.ids}.}  

     \item{newranef.ids}{An optional matrix with the number of columns equal to the number of response-specific random intercepts included in the model. Element \eqn{(i,j)} indicates the cluster ID of row \eqn{i} in the response matrix for random intercept \eqn{j}. Defaults to \code{NULL}, in which case random intercept IDs are taken from the fitted boral object itself, if appropriate, i.e., from \code{object$ranef.ids}.}  

     \item{distmat}{A distance matrix required to calculate correlations across sites when a non-independence correlation structure on the latent variables is imposed.} 

     \item{predict.type}{The type of prediction to be made. Either takes value \code{"conditional"} in which case the prediction is made conditionally on the predicted latent variables and any random row effects in the model, or \code{"marginal"} in which case the prediction marginalizes (averages) over the latent variables and random row effects in the model. Defaults to \code{"conditional"}.}  

     \item{scale}{The type of prediction required. The default is on the scale of the linear predictors; the alternative \code{scale = "response"} is on the scale of the response variable. For example, if the binomial family is used, then the default predictions (and associated uncertainty intervals) provide probabilities on linear predictor scale, while \code{scale = "response"} gives the predicted probabilities between 0 and 1. 
  
     Note things are slightly more complicated for zero truncated distributions because the log-link connects the mean of the \emph{untruncated} distribution to the linear predictor. Therefore if \code{scale = "link"}, then the linear predictor is returned. But if \code{scale = "response"}, then actual predicted mean is returned.}

     \item{est}{A choice of either whether to print the posterior median (\code{est = "median"}) or posterior mean (\code{est = "mean"}) of the parameters.}

     \item{prob}{A numeric scalar in the interval (0,1) giving the target probability coverage of the intervals. Defaults to 0.95.}   

     \item{lv.mc}{If the predictions are made marginalizing over the latent variables, then number of Monte-Carlo samples to take when performing the relevant integration.}

     \item{return.alllinpred}{If \code{TRUE}, then the full array of predictions, on the linear predictor scale, across all MCMC samples is predicted. This is useful if the user wants to transform the predictions onto a different scale or for further manipulation, say. Defaults to \code{FALSE}.}

    \item{...}{Not used.}
     }


\details{
In the Bayesian MCMC framework, predictions are based around the posterior predictive distribution, which is the integral of the quantity one wants to predict on, integrated or averaged over the posterior distribution of the parameters and latent variables. For example, on the linear predictor scale, predictions are made as,

\deqn{\eta_{ij} = \alpha_i + \beta_{0j} + \bm{x}^\top_i\bm{\beta}_j + \bm{z}^\top_i\bm{b}_j + \bm{u}^\top_i\bm{\theta}_j; \quad i = 1,\ldots,n; j = 1,\ldots,p,}

where \eqn{\beta_{0j} + \bm{x}^\top_i\bm{\beta}_j} is the component of the linear predictor due to the covariates \eqn{\bm{X}} plus an intercept, \eqn{\bm{z}^\top_i\bm{b}_j} is the component due to response-specific random intercept, \eqn{\bm{u}^\top_i\bm{\theta}_j} is the component due to the latent variables, and \eqn{\alpha_i} is the component due to one or more fixed or random row effects. Not all of these components may be included in the model, and the above is just representing the general case. 

Note that for the above to work, one must have saved the MCMC samples in the fitted boral object, that is, set \code{save.model = TRUE} when fitting. 


Two types of predictions are possible using this function:
\itemize{
     \item The first type is \code{predict.type = "conditional"}, meaning predictions are made conditionally on the predicted latent variables and any (random) row effects and response-specific random intercepts included in the model. This is mainly used when predictions are made onto the \emph{same} set of sites that the model was fitted to, although a \code{newX} can be supplied in this case if we want to extrapolate on to the same set of sites but under different environmental conditions. 

     \item The second type of prediction is \code{predict.type = "marginal"}, meaning predictions are made marginally or averaging over the latent variables and any (random) row effects and response-specific random intercepts included in the model. This is mainly used when predictions are made onto a \emph{new} set of sites where the latent variables/row effects/response-specific random intercepts are unknown. Consequently, arguments \code{newX}, \code{newrow.ids} and \code{newranef.ids} are often supplied in such a setting since we are extrapolating to new observational units. The integration over the latent variables and random row effects is done via Monte-Carlo integration. Please note however that the integration is done on the linear predictor scale. 
     }

More information on conditional versus marginal predictions in latent variable models can be found in Warton et al., (2015). In both cases, and if \code{return.alllinpred = FALSE}, the function returns a point prediction (either the posterior mean or median depending on \code{est}) and the lower and upper bounds of a (100\eqn{\times}\code{prob}) \% interval of the posterior prediction. All of these quantities are calculated empirically based across the MCMC samples.
}

\value{
A list containing the following components:
    \item{linpred}{A matrix containing posterior point predictions (either posterior mean or median depending on \code{est}), on the linear predictor scale.}
     
    \item{lower}{A matrix containing the lower bound of the (100\eqn{\times}\code{prob}) \% interval of the posterior predictions, on the linear predictor scale.}

    \item{upper}{A matrix containing the upper bound of the (100\eqn{\times}\code{prob}) \% interval of the posterior predictions, on the linear predictor scale.}

    \item{all.linpred}{If \code{return.alllinpred = TRUE}, then only an array of predicted linear predictions across all MCMC samples.}

     }

\references{
\itemize{
\item Gelman et al. (2013) Bayesian Data Analysis. CRC Press.

\item Warton et al. (2015). So Many Variables: Joint Modeling in Community Ecology. Trends in Ecology and Evolution, 30, 766-779.
}
}

\author{
\packageAuthor{boral}

Maintainer: \packageMaintainer{boral}
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

testpath <- file.path(tempdir(), "jagsboralmodel.txt")


## Example 1 - model with two latent variables, random site effects, 
## 	and environmental covariates
spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
    row.eff = "random", lv.control = list(num.lv = 2), 
    mcmc.control = example_mcmc_control, save.model = TRUE,
    model.name = testpath)

## Predictions conditional on predicted latent variables
getcondpreds <- predict(spiderfit_nb)

## Predictions marginal on latent variables, random row effects
## The intervals for these will generally be wider than the
##   conditional intervals.
getmargpreds <- predict(spiderfit_nb, predict.type = "marginal")


## Now suppose you extrpolate to new sites
newX <- rmvnorm(100, mean = rep(0,ncol(X)))

## Below won't work since conditional predictions are made to the same sites
#getcondpreds <- predict(spiderfit_nb, newX = newX)

## Marginal predictions will work though, provided newrow.ids is set up 
## properly. For example,
new_row_ids <- matrix(sample(1:28,100,replace=TRUE), 100, 1)
while(length(table(new_row_ids)) != 28) { 
    new_row_ids <- matrix(sample(1:28,100,replace=TRUE), 100, 1)
    }
getmargpreds <- predict(spiderfit_nb, newX = newX, predict.type = "marginal", 
     newrow.ids = new_row_ids)

     
## Example 1b - Similar to 1 except with no random site effects, 
## 	and a non-independence correlation structure for the latent variables
##      based on a fake distance matrix
fakedistmat <- as.matrix(dist(1:n))
spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
     lv.control = list(type = "squared.exponential", num.lv = 2, 
          distmat = fakedistmat), model.name = testpath, 
     mcmc.control = example_mcmc_control, save.model = TRUE)

getmargpreds <- predict(spiderfit_nb, predict.type = "marginal", 
     distmat = fakedistmat)

## Now suppose you extrpolate to new sites
newfakedistmat <- as.matrix(dist(1:100))

getmargpreds <- predict(spiderfit_nb, newX = newX, 
     predict.type = "marginal", distmat = newfakedistmat)


## Example 1c - similar to 1 except instead of random site effects,
##   there are species-specific random intercepts at a so-called
##   "region" level
y <- spider$abun
X <- scale(spider$x)
spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
    lv.control = list(num.lv = 2), 
    ranef.ids = data.frame(region = rep(1:7,each=4)), 
    mcmc.control = example_mcmc_control, model.name = testpath, 
    save.model = TRUE) 

## Predictions conditional on predicted latent variables and 
##   random intercepts
getcondpreds <- predict(spiderfit_nb)

## Predictions marginal on latent variables, random intercepts
## The intervals for these will generally be wider than the
##   conditional intervals.
getmargpreds <- predict(spiderfit_nb, predict.type = "marginal")

## Now suppose you extrpolate to new sites
newX <- rmvnorm(100, mean = rep(0,ncol(X)))

## Marginal predictions will work though, provided newranef.ids is set up 
## properly. For example,
new_ranef_ids <- matrix(sample(1:7,100,replace=TRUE), 100, 1)
getmargpreds <- predict(spiderfit_nb, newX = newX, predict.type = "marginal",
     newranef.ids = new_ranef_ids)

     
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
     lv.control = list(num.lv = 2), save.model = TRUE, 
     mcmc.control = example_mcmc_control,
     model.name = testpath)

     
## Predictions conditional on predicted latent variables   
getcondpreds <- predict(fit_traits)     
     
## Predictions marginal on latent variables
## The intervals for these will generally be wider than the
##   conditional intervals.
getmargpreds <- predict(fit_traits, predict.type = "marginal")
}

}

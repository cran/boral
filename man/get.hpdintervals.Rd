\name{get.hpdintervals}

\alias{get.hpdintervals}

\docType{package}

\title{Highest posterior density intervals for a fitted model}

\description{
\Sexpr[results=rd, stage=render]{lifecycle::badge("stable")}

Calculates the lower and upper bounds of the highest posterior density intervals for parameters and latent variables in a fitted model.}

\usage{
get.hpdintervals(y, X = NULL, traits = NULL, row.ids = NULL, ranef.ids = NULL, 
	fit.mcmc, lv.control, prob = 0.95, num.lv = NULL)
}

\arguments{
  \item{y}{The response matrix that the model was fitted to.}
  
  \item{X}{The covariate matrix used in the model. Defaults to \code{NULL}, in which case it is assumed no model matrix was used.}  
  
  \item{traits}{The trait matrix used in the model. Defaults to \code{NULL}, in which case it is assumed no traits were included.}  

 \item{row.ids}{A matrix with the number of rows equal to the number of rows in the response matrix, and the number of columns equal to the number of row effects to be included in the model. Element \eqn{(i,j)} indicates the cluster ID of row \eqn{i} in the response matrix for random effect \eqn{j}; please see \code{\link{boral}} for details. Defaults to \code{NULL}, in which case iti assumed no random effects were included in the model.}

 \item{ranef.ids}{A matrix with the number of rows equal to the number of rows in the response matrix, and the number of columns equal to the number of random intercepts to be included in the model. Element \eqn{(i,j)} indicates the cluster ID of row \eqn{i} in the response matrix for random intercept \eqn{j}; please see \code{\link{about.ranefs}} for details. Defaults to \code{NULL}, in which case it is assumed no random intercepts are to be included in the model. If supplied, then response-specific random intercepts are assumed to come from a normal distribution with mean zero and unknown (response-specific) standard deviation.}

 \item{fit.mcmc}{All MCMC samples for the fitted model. These can be extracted by fitting a model using \code{boral} with \code{save.model = TRUE}, and then applying \code{get.mcmcsamples(fit)}.} 
  
 \item{lv.control}{A list (currently) with the following arguments:
    \itemize{
        \item \emph{num.lv:} which specifies the number of true latent variables to generate. Defaults to 0.
        
        \item \emph{type:} which specifies the type the correlation structure of the latent variables (across sites). Defaults to independence correlation structure.
        
        \item \emph{distmat:} which a distance matrix required to calculate correlations across sites when a non-independence correlation structure on the latent variables is imposed. 
        }
    Please see \code{\link{about.lvs}} for more information.
    }

    \item{prob}{A numeric scalar in the interval (0,1) giving the target probability coverage of the intervals. Defaults to 0.95.}   

    \item{num.lv}{Old argument superceded by \code{lv.control}. Defaults to \code{NULL} and ignored.} 
}

\details{
The function uses the \code{HPDinterval} function from the \code{coda} package to obtain the HPD intervals. See \code{HPDinterval} for details regarding the definition of the HPD interval. For interpreting the results, please check the dimension names of each of the components below to better ascertain what is being printed.
}

\value{
A list containing the following components, where applicable:

\item{lv.coefs}{An array giving the lower and upper bounds of the HPD intervals for the response-specific intercepts, latent variable coefficients, and dispersion parameters if appropriate.}

\item{lv}{An array giving the and upper bounds of the HPD intervals for the latent variables.}

\item{lv.covparams}{A matrix giving the lower and upper bounds of the HPD intervals for the parameters characterizing the correlation structure of the latent variables when they are assumed to be non-independent across rows.}

\item{row.coefs}{A list with each element being a matrix giving the lower and upper bounds of the HPD intervals for row effects. The number of elements in the list should equal the number of row effects included in the model i.e., \code{ncol(row.ids)}.}

\item{row.sigma}{A list with each element being a vector giving the lower and upper bounds of the HPD interval for the standard deviation of the normal distribution for the row effects. The number of elements in the list should equal the number of row effects included in the model i.e., \code{ncol(row.ids)}.}

\item{ranef.coefs}{A list with each element being a array giving the lower and upper bounds of the HPD intervals for response-specific random intercepts. The number of elements in the list should equal the number of row effects included in the model i.e., \code{ncol(ranef.ids)}.}

\item{ranef.sigma}{An array giving the lower and upper bounds of the HPD interval for the standard deviation of the normal distribution for the response-specific random intercepts. The number of elements in the list should equal the number of row effects included in the model i.e., \code{ncol(row.ids)}.}

\item{X.coefs}{An array giving the lower and upper bounds of the HPD intervals for coefficients relating to the covariate matrix.}

\item{traits.coefs}{An array giving the lower and upper of the HPD intervals for coefficients and standard deviation relating to the traits matrix.}

\item{cutoffs}{A matrix giving the lower and upper bounds of the HPD intervals for common cutoffs in proportional odds regression.}

\item{powerparam}{A vector giving the lower and upper bounds of the HPD interval for common power parameter in tweedie regression.}
}

\author{
\packageAuthor{boral}

Maintainer: \packageMaintainer{boral}
}

\note{
\code{\link{boral}} fits the model and returns the HPD intervals by default. 
}

\section{Warnings}{
\itemize{
\item HPD intervals tend to be quite wide, and inference is somewhat tricky with them. This is made more difficult by the multiple comparison problem due to the construction one interval for each parameter!

\item Be careful with interpretation of coefficients and HPD intervals if different columns of the response matrix have different distributions!

\item HPD intervals for the cutoffs in proportional odds regression may be poorly estimated for levels with few data.
}
}

\examples{
\dontrun{
library(mvabund) ## Load a dataset from the mvabund package
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
 
 
## Example 1 - model with two latent variables, site effects, 
## 	and no environmental covariates
spiderfit_nb <- boral(y, family = "negative.binomial", 
    lv.control = list(num.lv = 2), row.eff = "fixed", 
    mcmc.control = example_mcmc_control, model.name = testpath)

spiderfit_nb$hpdintervals

## Example 2a - model with no latent variables, no site effects, 
##      and environmental covariates
spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
    mcmc.control = example_mcmc_control, model.name = testpath)

spiderfit_nb$hpdintervals

## Example 2b - suppose now, for some reason, the 28 rows were
## 	sampled such into four replications of seven sites
## Let us account for this as a fixed effect
spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
    row.eff = "fixed", row.ids = matrix(rep(1:7,each=4),ncol=1),
    mcmc.control = example_mcmc_control, model.name = testpath)

spiderfit_nb$hpdintervals

## Example 2c - suppose now, for some reason, the 28 rows reflected
## 	a nested design with seven regions, each with four sub-regions
## We can account for this nesting as a random effect
spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
    row.eff = "random", 
    row.ids = cbind(1:n, rep(1:7,each=4)), 
    mcmc.control = example_mcmc_control, model.name = testpath)

spiderfit_nb$hpdintervals

## Example 2d - model with environmental covariates and 
##  two structured latent variables using fake distance matrix
fakedistmat <- as.matrix(dist(1:n))
spiderfit_lvstruc <- boral(y, X = X, family = "negative.binomial", 
    lv.control = list(num.lv = 2, type = "exponential", distmat = fakedistmat), 
     mcmc.control = example_mcmc_control, model.name = testpath)

spiderfit_nb$hpdintervals

## Example 2e - Similar to 2d, but we will species-specific random intercepts
##   for the seven regions (with row effects in the model)
spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
    ranef.ids = data.frame(region = rep(1:7,each=4)), 
    mcmc.control = example_mcmc_control, model.name = testpath) 

spiderfit_nb$hpdintervals
}
}

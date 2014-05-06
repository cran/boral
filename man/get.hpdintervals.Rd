\name{get.hpdintervals}
\alias{get.hpdintervals}
\title{Highest posterior density intervals for an boral model}
\description{Calculates the lower and upper bounds of the highest posterior density intervals for parameters and latent variables in a fitted boral model.}

\usage{
get.hpdintervals(y, X = NULL, fit.mcmc, num.lv)
}

\arguments{
  \item{y}{The response matrix that the boral model was fitted to.}
  
  \item{X}{The model matrix used in the boral model. Defaults to \code{NULL}, in which case it is assumed no model matrix was used.}  
  
  \item{fit.mcmc}{All MCMC samples for the fitted boral model, as obtained from JAGS. These can be extracted by fitting an boral model using \code{\link{boral}} with \code{save.model = TRUE}, and the applying \code{as.mcmc} on the \code{jags.model} component of the output.} 
  
  \item{num.lv}{The number of latent variables used in the boral model. If zero, then HPD intervals will not be produced for latent variables.}   
}

\details{
The function uses the \code{HPDinterval} function from the \code{coda} package to obtain the HPD intervals. See \code{HPDinterval} for details regarding the defintion of the HPD interval. 
}

\value{
\item{lv.coefs.hpd.lower/upper}{Two matrices corresponding to the lower and upper bounds of the HPD intervals for the column-specific intercepts, regression coefficients relating to the latent variables, and dispersion parameters.}
\item{lv.hpd.lower/upper}{Two matrices corresponding to the lower and upper bounds of the HPD intervals for the latent variables.}
\item{site.coefs.hpd.lower/upper}{Two matrices corresponding to the lower and upper bounds of the HPD intervals for row effects.}
\item{X.coefs.hpd.lower/upper}{Two matrices corresponding to the lower and upper bounds of the HPD intervals for regression coefficients relating to the model matrix \code{X}.}
\item{cutoffs.hpd.lower/upper}{Two vectors corresponding to the lower and upper bounds of the HPD intervals for common cutoffs in proportional odds regression.}
\item{powerparam.hpd.lower/upper}{Two scalars corresponding to the lower and upper bounds of the HPD interval for common power parameter in tweedie regression.}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\note{
\code{\link{boral}} fits the boral model and returns the HPD intervals by default. 
}

\section{Warnings}{
\itemize{
\item HPD intervals tend to be quite wide, and inference is somewhat tricky with them. This is made more difficult by the multiple comparison problem due to the construction one interval for each parameter!
\item Be very careful with interpretation of coefficients and HPD intervals if different columns of \eqn{y} have different distributions!
\item HPD intervals for the cutoffs in proportional odds regression may be poorly estimated for levels with few data.
}
}

\examples{
\dontrun{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
n <- nrow(y); p <- ncol(y); 
    
## Example 1 - model with two latent variables, site effects, 
## 	and no environmental covariates
spider.fit.nb <- boral(y, family = "negative.binomial", num.lv = 2, 
     site.eff = TRUE, save.model = TRUE)

## Returns a list with components corresponding to values described above.
spider.fit.nb$hpdintervals 

## Example 2 - model with one latent variable, site effects, 
## 	and environmental covariates
spider.fit.nb2 <- boral(y, X = spider$x, family = "negative.binomial", num.lv = 1, 
     site.eff = TRUE, save.model = TRUE)

## Returns a list with components corresponding to values described above.
spider.fit.nb2$hpdintervals 

}
}
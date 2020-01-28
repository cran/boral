\name{get.mcmcsamples}

\alias{get.mcmcsamples}

\docType{package}

\title{Extract MCMC samples from models}

\description{Extract the MCMC samples from fitted models, taking into account the burnin period and thinning.}

\usage{
get.mcmcsamples(object)
}

\arguments{
  \item{object}{An object for class "boral".}
}

\details{
For the function to work, the JAGS model file (containing the MCMC samples from the call to JAGS) has to have been saved when fitting the model, that is, \code{save.model = TRUE}. The function will throw an error if it cannot find the the JAGs model file.
}

\value{
A matrix containing the MCMC samples, with the number of rows equal to the number of MCMC samples after accounting the burnin period and thinning (i.e., number of rows = (n.iteration - n.burnin)/n.thin), and the number of columns equal to the number of parameters in the fitted model. 
}


\author{
\packageAuthor{boral}

Maintainer: \packageMaintainer{boral}
}

%\seealso{
%\code{\link{get.residual.cor}}, which calculates the residual correlation matrix for models involving latent variables.
%}

\examples{
\dontrun{
## NOTE: The values below MUST NOT be used in a real application;
## they are only used here to make the examples run quick!!!
example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
     n.thin = 1)
     
testpath <- file.path(tempdir(), "jagsboralmodel.txt")


library(mvabund) ## Load a dataset from the mvabund package
library(corrplot) ## For plotting correlations
data(spider)
y <- spider$abun
X <- scale(spider$x)
n <- nrow(y)
p <- ncol(y)
    
spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
     mcmc.control = example_mcmc_control, model.name = testpath, 
     save.model = TRUE)

mcmcsamps <- get.mcmcsamples(spiderfit_nb)
}

}

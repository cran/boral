\name{get.mcmcsamples}

\alias{get.mcmcsamples}

\title{Extract MCMC samples from boral models}

\description{Extract the MCMC samples from boral models, taking into account the burnin period and thinning.}

\usage{
get.mcmcsamples(object)
}

\arguments{
  \item{object}{An object for class "boral".}
}

\details{
For the function to work, the JAGS model file (containing the MCMC samples from the call to JAGS) has to have been saved when fitting the boral model, that is, \code{save.model = TRUE}. The function will throw an error if it cannot find the the JAGs model file.
}

\value{
A matrix containing the MCMC samples, with the number of rows equal to the number of MCMC samples after accounting the burnin period and thinning (i.e., number of rows = (n.iteration - n.burnin)/n.thin), and the number of columns equal to the number of parameters in the fitted boral model. 
}


\author{
Francis K.C. Hui \email{fhui28@gmail.com}; Wade Blanchard \email{wade.blanchard@anu.edu.au}
}

%\seealso{
%\code{\link{get.residual.cor}}, which calculates the residual correlation matrix for boral models involving latent variables.
%}

\examples{
\dontrun{
## NOTE: The values below MUST NOT be used in a real application;
## they are only used here to make the examples run quick!!!
example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
     n.thin = 1)
     
library(mvabund) ## Load a dataset from the mvabund package
library(corrplot) ## For plotting correlations
data(spider)
y <- spider$abun
X <- scale(spider$x)
n <- nrow(y)
p <- ncol(y)
    
spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
     save.model = TRUE, mcmc.control = example_mcmc_control, 
     save.model = TRUE)

mcmcsamps <- get.mcmcsamples(spiderfit_nb)
}

}
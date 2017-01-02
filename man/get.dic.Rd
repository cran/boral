\name{get.dic}

\alias{get.dic}

\title{Extract Deviance Information Criterion for boral model}

\description{Calculates the Deviance Information Criterion (DIC) for a boral model fitted using JAGS.}

\usage{
get.dic(jagsfit)
}

\arguments{
  \item{jagsfit}{The \code{jags.model} component of the output, from a model fitted using \code{boral} with \code{save.model = TRUE}.}
}

\details{
Details regarding the Deviance Information Criterion may be found in (Spiegelhalter et al., 2002; Ntzoufras, 2011; Gelman et al., 2013). The DIC here is based on the conditional log-likelihood i.e., the latent variables (and row effects if applicable) are treated as "fixed effects". A DIC based on the marginal likelihood is obtainable from \code{\link{get.more.measures}}, although this requires a much longer time to compute. For models with overdispered count data, conditional DIC may not perform as well as marginal DIC (Millar, 2009)
}

\value{
DIC value for the jags model.
}

\references{
\itemize{
\item Gelman et al. (2013). Bayesian data analysis. CRC press.

\item Millar, R. B. (2009). Comparison of hierarchical Bayesian models for overdispersed count data using DIC and Bayes' factors. Biometrics, 65, 962-969.

\item Ntzoufras, I. (2011). Bayesian modeling using WinBUGS (Vol. 698). John Wiley & Sons.

\item Spiegelhalter, et al. (2002). Bayesian measures of model complexity and fit. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 64, 583-639.
}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\note{
This function and consequently the DIC value is automatically returned when a boral model is fitted using \code{\link{boral}} with \code{calc.ics = TRUE}. 
}

\seealso{
\code{\link{get.measures}} and \code{\link{get.more.measures}} for other information criteria which could potentially be used for variable selection.
}

\examples{
\dontrun{
## NOTE: The values below MUST NOT be used in a real application;
## they are only used here to make the examples run quick!!!
example.mcmc.control <- list(n.burnin = 10, n.iteration = 100, 
     n.thin = 1)
     
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
n <- nrow(y); p <- ncol(y); 
    
spider.fit.nb <- boral(y, family = "negative.binomial", num.lv = 2, 
     save.model = TRUE, calc.ics = TRUE, 
     mcmc.control = example.mcmc.control)

spider.fit.nb$ics ## DIC returned as one of several information criteria.
}
}
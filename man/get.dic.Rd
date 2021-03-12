\name{get.dic}

\alias{get.dic}

\docType{package}

\title{Extract Deviance Information Criterion for a fitted model}

\description{
\Sexpr[results=rd, stage=render]{lifecycle::badge("defunct")}

Calculates the Deviance Information Criterion (DIC) for a model fitted using JAGS. WARNING: As of version 1.6, this function is no longer maintained (and probably doesn't work properly, if at all)!}

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

\item Spiegelhalter et al. (2002). Bayesian measures of model complexity and fit. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 64, 583-639.
}
}

\author{
\packageAuthor{boral}

Maintainer: \packageMaintainer{boral}
}

\note{
This function and consequently the DIC value is automatically returned when a model is fitted using \code{\link{boral}} with \code{calc.ics = TRUE}. 
}

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
     save.model = TRUE, calc.ics = TRUE, mcmc.control = example_mcmc_control,
     model.name = testpath)

spiderfit_nb$ics ## DIC returned as one of several information criteria.
}
}

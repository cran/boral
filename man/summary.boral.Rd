\name{summary.boral}

\alias{summary.boral}
\alias{print.summary.boral}

\docType{package}

\title{Summary of fitted boral object}

\description{
\Sexpr[results=rd, stage=render]{lifecycle::badge("stable")}

A summary of the fitted boral objects including the type of model fitted e.g., error distribution, number of latent variables parameter estimates, and so on.}

\usage{
\method{summary}{boral}(object, est = "median", ...)

\method{print}{summary.boral}(x,...)
}

\arguments{
  \item{object}{An object of class "boral".}

  \item{x}{An object of class "boral".}

  \item{est}{A choice of either whether to print the posterior median (\code{est = "median"}) or posterior mean (\code{est = "mean"}) of the parameters.}

  \item{...}{Not used.}
}

\value{
Attributes of the model fitted, parameter estimates, and posterior probabilities of including individual and/or grouped coefficients in the model based on SSVS if appropriate.
}

\author{
\packageAuthor{boral}

Maintainer: \packageMaintainer{boral}
}

\seealso{
\code{\link{boral}} for the fitting function on which \code{summary} is applied.
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

spiderfit_nb <- boral(y, family = "negative.binomial", 
     lv.control = list(num.lv = 2), row.eff = "fixed", 
     mcmc.control = example_mcmc_control, model.name = testpath)

summary(spiderfit_nb)
}
}

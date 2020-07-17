\name{fitted.boral}

\alias{fitted.boral}

\docType{package}

\title{Extract Model Fitted Values for an boral object}

\description{
\Sexpr[results=rd, stage=render]{lifecycle::badge("stable")}

Calculated the predicted mean responses based on the fitted model, by using the posterior medians or means of the parameters.}

\usage{ 
\method{fitted}{boral}(object, est = "median", linear.predictor = FALSE, ...) 
}

\arguments{
  \item{object}{An object of class "boral".}

  \item{est}{A choice of either the posterior median (\code{est == "median"}) or posterior mean (\code{est == "mean"}), which are then treated as estimates and the fitted values are calculated from. Default is posterior median.}

  \item{linear.predictor}{Determines the scale on which to return the fitted values. When set to \code{TRUE}, it returns the fitted mean values on the linear predictor scale. When set to \code{FALSE}, which is the defaul behavior, the fitted mean values are on the response scale. Note things are slightly more complicated for zero truncated distributions because, the log-link connects the mean of the \emph{untruncated} distribution to the linear predictor. Therefore if \code{linear.predictor = TRUE}, then the linear predictor is returned. But if \code{linear.predictor = FALSE}, then actual fitted mean value is returned.}

  \item{...}{Not used.}
}

\details{
This fitted values here are calculated based on a point estimate of the parameters, as determined by the argument \code{est}. A fully Bayesian approach would calculate the fitted values by averaging over the posterior distribution of the parameters i.e., ergodically average over the MCMC samples. For simplicity and speed though (to avoid generation of a large number of predicted values), this is not implemented.
}

\value{
A list containing \code{ordinal.probs} which is an array with dimensions (no. of rows of \code{y}) x (no. of rows of \code{y}) x (no. of levels) containing the predicted probabilities for ordinal columns, and \code{out} which is a matrix of the same dimension as the original response matrix \code{y} containing the fitted values. For ordinal columns, the "fitted values" are defined as the level/class that had the highest fitted probability. 
}

\author{
\packageAuthor{boral}

Maintainer: \packageMaintainer{boral}
}

\seealso{
\code{\link{plot.boral}} which uses the fitted values calculated from this function to construct plots for residual analysis; \code{\link{ds.residuals}} for calculating the Dunn-Smyth residuals for a fitted model. 
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

spiderfit_nb <- boral(y, family = "negative.binomial", lv.control = list(num.lv = 2),
     row.eff = "fixed", mcmc.control = example_mcmc_control, model.name = testpath)

fitted(spiderfit_nb)
}
}

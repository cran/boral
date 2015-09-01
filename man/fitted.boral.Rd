\name{fitted.boral}
\alias{fitted.boral}
\title{Extract Model Fitted Values for an boral object}
\description{Calculated the predicted mean responses based on the fitted boral model, by using the posterior medians or means of the parameters.}

\usage{ 
\method{fitted}{boral}(object, est = "median",...) 
}

\arguments{
  \item{object}{An object of class "boral".}

  \item{est}{A choice of either the posterior median (\code{est == "median"}) or posterior mean (\code{est == "mean"}), which are then treated as estimates and the fitted values are calculated from. Default is posterior median.}

  \item{...}{Not used.}
}

\details{
This fitted values here are calculated based on a point estimate of the parameters, as determined by the argument \code{est}. A fully Bayesian approach would calculate the fitted values by averaging over the posterior distribution of the parameters i.e., ergodically average over the MCMC samples. For simplicity and speed though (to avoid generation of a large number of predicted values), this is not implemented.
}

\value{
A list with potential \code{NULL} elements in it, containing \code{ordinal.probs} which is an array with dimensions (no. of rows of \code{y}) x (no. of rows of \code{y}) x (no. of levels) containing the predicted probabilities for ordinal columns, and \code{out} which is a matrix of the same dimension as the original response matrix \code{y} containing the fitted values.
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\seealso{
\code{\link{plot.boral}} which uses the fitted values calculated from this function to construct plots for residual analysis; \code{\link{ds.residuals}} for calculating the Dunn-Smyth residuals for a fitted boral model. 
}

\examples{
\dontrun{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun

spider.fit.nb <- boral(y, family = "negative.binomial", num.lv = 2, 
     row.eff = "fixed")

fitted(spider.fit.nb)
}
}
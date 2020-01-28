\name{ds.residuals}

\alias{ds.residuals}

\docType{package}

\title{Dunn-Smyth Residuals for a fitted model}

\description{Calculates the Dunn-Smyth residuals for a fitted model and, if some of the responses are ordinal, a confusion matrix between predicted and true levels.}

\usage{ds.residuals(object, est = "median")}

\arguments{
  \item{object}{An object for class "boral".}

  \item{est}{A choice of either the posterior median (\code{est == "median"}) or posterior mean (\code{est == "mean"}), which are then treated as parameter estimates and the residuals are calculated from. Default is posterior median.}
}

\details{
Details regarding Dunn-Smyth residuals, based on the randomized quantile residuals of Dunn and Smyth (1996), can be found in \code{plot.manyglm} function in the \code{mvabund} package (Wang et al., 2012) where they are implemented in all their glory. Due their inherent stochasticity, Dunn-Smyth residuals will be slightly different each time this function is run. As with other types of residuals, Dunn-Smyth residuals can be used in the context of residual analysis. 

For ordinal responses, a single confusion matrix between the predicted levels (as based on the class with the highest probability) and true levels is aso returned. The table pools the results over all columns assumed to be ordinal.

The Dunn-Smyth residuals are calculated based on a point estimate of the parameters, as determined by the argument \code{est}. A fully Bayesian approach would calculate the residuals by averaging over the posterior distribution of the parameters i.e., ergodically average over the MCMC samples. In general however, the results (as in the trends seen in residual analysis) from either approach should be very similar.
}

\value{
A list containing \code{agree.ordinal} which is a single confusion matrix for ordinal columns, and \code{residuals} which contains Dunn-Smyth residuals.
}

\references{
\itemize{
\item Dunn, P. K., and Smyth, G. K. (1996). Randomized quantile residuals. Journal of Computational and Graphical Statistics, 5, 236-244.

\item Wang, Y. et al. (2012). mvabund-an R package for model-based analysis of multivariate abundance data. Methods in Ecology and Evolution, 3, 471-474.
}
}

\author{
\packageAuthor{boral}

Maintainer: \packageMaintainer{boral}
}

\seealso{
\code{\link{plot.boral}} for constructing residual analysis plots directly; \code{\link{fitted.boral}} which calculated fitted values from a model.}

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

ds.residuals(spiderfit_nb) 
}
}

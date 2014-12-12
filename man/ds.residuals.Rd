\name{ds.residuals}
\alias{ds.residuals}
\title{Dunn-Smyth Residuals for a boral model}
\description{Calculates the Dunn-Smyth residuals for a fitted boral model or, if some of the responses are ordinal, a table of agreement between predicted and true levels.}

\usage{ds.residuals(object, est = "median")}

\arguments{
  \item{object}{An object for class "boral".}

  \item{est}{A choice of either the posterior median (\code{est == "median"}) or posterior mean (\code{est == "mean"}), which are then treated as parameter estimates and the residuals are calculated from. Default is posterior median.}
}

\details{
Details regarding Dunn-Smyth residuals, based on the randomized quantile residuals of Dunn and Smyth (1996), can be found in \code{plot.manyglm} function in the \code{mvabund} package (Wang et al., 2012) where they are implemented in all their glory. Due their inherent stochasticity, Dunn-Smyth residuals will be slightly different each time this function is run. As with other types of residuals, Dunn-Smyth residuals can be used in the context of residual analysis. 

For ordinal responses, a single table of agreement between the predicted levels (as based on the class with the highest probability) and true levels is returned. The table pools the results over all columns assumed to be ordinal.

The Dunn-Smyth residuals are calculated based on a point estimate of the parameters, as determined by the argument \code{est}. A fully Bayesian approach would calculate the residuals by averaging over the posterior distribution of the parameters i.e., ergodically average over the MCMC samples. In general however, the results (as in the trends seen in residual analysis) from either approach should be very similar.
}

\value{
A list with potentially \code{NULL} elements, containing \code{agree.ordinal} which is a single table of agreement for ordinal columns, and \code{residuals} which contains Dunn-Smyth residuals.
}

\references{
\itemize{
\item Dunn, P. K., and Smyth, G. K. (1996). Randomized quantile residuals. Journal of Computational and Graphical Statistics, 5, 236-244.
\item Wang, Y. et al. (2012). mvabund-an R package for model-based analysis of multivariate abundance data. Methods in Ecology and Evolution, 3, 471-474.
}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\seealso{
\code{\link{plot.boral}} for constructing residual analysis plots directly; \code{\link{fitted.boral}} which calculated fitted values from a boral model.}

\examples{
\dontrun{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun

spider.fit.nb <- boral(y, family = "negative.binomial", num.lv = 2, 
     row.eff = "fixed", n.burnin = 10, n.iteration = 100, 
     n.thin = 1, calc.ics = FALSE)

ds.residuals(spider.fit.nb) 
}
}
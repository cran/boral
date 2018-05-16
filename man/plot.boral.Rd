\name{plot.boral}

\alias{plot.boral}

\title{Plots of a fitted boral object}

\description{Produces four plots relating to the fitted boral object, which can be used for residual analysis. If some of the columns are ordinal, then a single confusion matrix is also produced.}

\usage{
\method{plot}{boral}(x, est = "median", jitter = FALSE, ...)
}

\arguments{
  \item{x}{An object of class "boral".}

  \item{est}{A choice of either the posterior median (\code{est == "median"}) or posterior mean (\code{est == "mean"}) of the parameters, which are then treated as parameter estimates and the fitted values/residuals used in the plots are calculated from. Default is posterior median.}

  \item{jitter}{If \code{jitter = TRUE}, then some jittering is applied so that points on the plots do not overlap exactly (which can often occur with discrete data). Please see \code{\link{jitter}} for its implementation.}

   \item{...}{Additional graphical options to be included in. These include values for \cr \code{cex, cex.lab, cex.axis, cex.main, lwd}, and so on.}
}

\details{
Four plots are provided:
\enumerate{
     \item Plot of Dunn-Smyth residuals against the linear predictors. This can be useful to assess whether the assumed mean-variance relationship is adequately satisfied, as well as to look for particular outliers. For ordinal responses things are more ambiguous due to the lack of single definition for "linear predictor". Therefore, instead of linear predictors the Dunn-Smyth residuals are plotted against the fitted values (defined as the level with the highest fitted probability). It is fully acknowledged that this makes things VERY hard to interpret if only some of your columns are ordinal. 
     \item Plot of Dunn-Smyth residuals against the row index/row names.  
     \item Plot of Dunn-Smyth residuals against the column index/column names. Both this and the previous plot are useful for assessing how well each row/column of the response matrix is being modeled. 
     \item A normal quantile plot of the Dunn-Smyth residuals, which can be used to assess the normality assumption and overall goodness of fit.  
     }
For ordinal responses, a single confusion matrix between the predicted levels (as based on the class with the highest probability) and true levels is aso returned. The table pools the results over all columns assumed to be ordinal.
}     

\note{
Due the inherent stochasticity, Dunn-Smyth residuals and consequently the plots will be slightly different time this function is run. Note also the fitted values and residuals are calculated from point estimates of the parameters, as opposed to a fully Bayesian approach (please see details in \code{\link{fitted.boral}} and \code{\link{ds.residuals}}). Consequently, it is recommended that this function is run several times to ensure that any trends observed in the plots are consistent throughout the runs.

As mentioned above, for ordinal responses things are much more challenging as there is no single definition for "linear predictor". Instead of linear predictors then, for the first plot the Dunn-Smyth residuals are plotted against the fitted values, defined as the level with the highest fitted probability. It is fully acknowledged that this makes things VERY hard to interpret if only some of your columns are ordinal though. Suggestions to improve this are welcome!!!
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}; Wade Blanchard \email{wade.blanchard@anu.edu.au}
}

\seealso{\code{\link{fitted.boral}} to obtain the fitted values, \code{\link{ds.residuals}} to obtain Dunn-Smyth residuals and details as to what they are.}

\examples{
\dontrun{
## NOTE: The values below MUST NOT be used in a real application;
## they are only used here to make the examples run quick!!!
example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
     n.thin = 1)

library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun

spider.fit.p <- boral(y, family = "poisson", num.lv = 2, 
	row.eff = "fixed", mcmc.control = example_mcmc_control)

par(mfrow = c(2,2))
plot(spider.fit.p) 
## A distinct fan pattern is observed in the plot of residuals 
## versus linear predictors plot. 


spiderfit_nb <- boral(y, family = "negative.binomial", num.lv = 2, 
     row.eff = "fixed", mcmc.control = example_mcmc_control)

par(mfrow = c(2,2))
plot(spiderfit_nb) 
## The fan shape is not as clear now, 
## and the normal quantile plot also suggests a better fit to the data 
}
}
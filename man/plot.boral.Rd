\name{plot.boral}
\alias{plot.boral}
\title{Plots of a fitted boral object}
\description{Produces four plots relating to the fitted boral object, which can be used for residual analysis.}

\usage{
\method{plot}{boral}(x, est = "median", jitter = FALSE, a = 1,...)
}

\arguments{
  \item{x}{An object of class "boral".}

  \item{est}{A choice of either the posterior median (\code{est == "median"}) or posterior mean (\code{est == "mean"}) of the parameters, which are then treated as parameter estimates and the fitted values/residuals used in the plots are calculated from. Default is posterior median.}

  \item{jitter}{If \code{jitter = TRUE}, then some jittering is applied so that points on the plots do not overlap exactly (which can often occur with discrete data). Please see \code{\link{jitter}} for its implementation.}

  \item{a}{Default parameter used in cex. Graphical options are then adjusted as \cr \code{par(ask = T, cex = a, cex.main = a, ...)}. Defaults to 1.}

   \item{...}{Additional graphical options to be included in \code{par}.}
}

\details{
Four types of plots are provided:
\enumerate{
     \item Plot of Dunn-Smyth residuals against the linear predictors. This can be useful to assess whether the assumed mean-variance relationship is adequately satisfied, as well as to look for particular outliers.  
     \item Plot of Dunn-Smyth residuals against the row index/row names.  
     \item Plot of Dunn-Smyth residuals against the column index/column names. Both this and the previous plot are useful for assessing how well each row/column of the response matrix is being modeled. 
     \item A normal quantile plot of the Dunn-Smyth residuals, which can be used to assess the normality assumption and overall goodness of fit.  
     }
}     

\note{
If all the columns of \code{y} were assumed to be ordinal, then this function is immediately stopped, as not residuals can be plotted in this case.

Due the inherent stochasticity, Dunn-Smyth residuals and consequently the plots will be slightly different time this function is run. Note also the fitted values and residuals are calculated from point estimates of the parameters, as opposed to a fully Bayesian approach (please see details in \code{\link{fitted.boral}} and \code{\link{ds.residuals}}). Consequently, it is recommended that this function is run several times to ensure that any trends observed in the plots are consistent throughout the runs.
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\seealso{\code{\link{fitted.boral}} to obtain the fitted values, \code{\link{ds.residuals}} to obtain Dunn-Smyth residuals and details as to what they are.}

\examples{
\dontrun{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun

spider.fit.p <- boral(y, family = "poisson", num.lv = 2, 
	row.eff = "fixed")

plot(spider.fit.p, ask = FALSE, mfrow = c(2,2)) 
## A distinct fan pattern is observed in the plot of residuals 
## versus linear predictors plot. 


spider.fit.nb <- boral(y, family = "negative.binomial", num.lv = 2, 
     row.eff = "fixed")

plot(spider.fit.nb, ask = FALSE, mfrow = c(2,2)) 
## The fan shape is not as clear now, 
## and the normal quantile plot also suggests a better fit to the data 
}
}
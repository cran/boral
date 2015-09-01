\name{get.residual.cor}
\alias{get.residual.cor}
\title{Extract residual correlations from boral models}
\description{Calculates the residual correlation from models that include latent variables.}

\usage{
get.residual.cor(object, est = "median", prob = 0.95)
}

\arguments{
  \item{object}{An object for class "boral".}
  
  \item{est}{A choice of either the posterior median (\code{est = "median"}) or posterior mean (\code{est = "mean"}), which are then treated as estimates and the fitted values are calculated from. Default is posterior median.}

  \item{prob}{A numeric scalar in the interval (0,1) giving the target probability coverage of the intervals, by which to determine whether the correlations are "significant". Defaults to 0.95.}   
}

\details{
In models with latent variables, the residual covariance matrix is calculated based on the matrix of latent variables regression coefficients formed by stacking the rows of \eqn{\bm{\theta}_j}. That is, if we denote \eqn{\bm{\Theta} = (\bm{\theta}_1 \ldots \bm{\theta}_p)'}, then the residual covariance and hence residual correlation matrix is calculated based on \eqn{\bm{\Theta}\bm{\Theta}'} (see also the help file for \code{\link{boral}}).

For multivariate abundance data, the inclusion of latent variables provides a parsimonious method of accounting for correlation between species. Specifically, the linear predictor,

\deqn{\beta_{0j} + \bm{x}^T_i\bm{\beta}_j + \bm{z}^T_i\bm{\theta}_j}

is normally distributed with a residual covariance matrix given by \eqn{\bm{\Theta}\bm{\Theta}'}. A strong residual covariance/correlation matrix between two species can then be interpreted as evidence of species interaction (e.g., facilitation or competition), missing covariates, as well as any additional species correlation not accounted for by shared environmental responses (see also Pollock et al., 2014, for residual correlation matrices in the context of Joint Species Distribution Models).

In addition to the residual correlation matrix, the median or mean point estimator of trace of the residual covariance matrix is returned, \eqn{\sum\limits_{j=1}^p [\bm{\Theta}\bm{\Theta}']_{jj}}. Often used in other areas of multivariate statistics, the trace may be interpreted as the amount of covariation explained by the latent variables. One situation where the trace may be useful is when comparing a pure LVM versus a model with latent variables and some predictors (correlated response models) -- the proportional difference in trace between these two models may be interpreted as the proportion of covariation between species explained by the predictors. Of course, the trace itself is random due to the MCMC sampling, and so it is not always guranteed to produce sensible answers =P
}

\value{
A list with the following components:
\item{cor}{A \eqn{p \times p} residual correlation matrix based on posteriori median or mean estimators of the latent variables and coefficients.}
\item{sig.cor}{A \eqn{p \times p} correlation matrix containing only the ``significant" correlations whose 95\% highest posterior interval does not contain zero. All non-significant correlations are zero to zero.}
\item{cov}{A \eqn{p \times p} covariance correlation matrix based on posteriori median or mean estimators of the latent variables and coefficients.}
\item{trace}{The median/mean point estimator of the trace (sum of the diagonal elements) of the residual covariance matrix.}
}

\references{
\itemize{
\item Pollock et al. (2014). Understanding co-occurrence by modelling species simultaneously with a Joint Species Distribution Model (JSDM). Methods in Ecology and Evolution, 5, 397-406.
}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\note{
Residual correlation matrices are reliably modeled only with two or more latent variables i.e., \code{num.lv > 1} when fitting the model using \code{boral}.  
}

\seealso{
\code{\link{get.enviro.cor}}, which calculates the correlation matrix due to similarities in the response to the explanatory variables (i.e., similarities due to a shared environmental response). 
}

\examples{
\dontrun{
library(mvabund) ## Load a dataset from the mvabund package
library(corrplot) ## For plotting correlations
data(spider)
y <- spider$abun
n <- nrow(y); p <- ncol(y); 
    
spider.fit.nb <- boral(y, X = spider$x, family = "negative.binomial", 
	num.lv = 2, save.model = TRUE)

res.cors <- get.residual.cor(spider.fit.nb)

corrplot(res.cors$sig.cor, title = "Residual correlations", 
	type = "lower", diag = FALSE, mar = c(3,0.5,2,1), tl.srt = 45)
}
}
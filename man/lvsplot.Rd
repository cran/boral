\name{lvsplot}
\alias{lvsplot}
\title{Plot the latent variables from an boral model} 
\description{Construct a 1-D index plot or 2-D scatterplot of the latent variables, and their corresponding coefficients (i.e. a biplot), from a fitted boral model.}

\usage{lvsplot(x, jitter = FALSE, a = 1, newplot = TRUE, biplot = TRUE, ind.spp = NULL, ...) }

\arguments{
  \item{x}{An object for class "boral".}

  \item{jitter}{If \code{jitter = TRUE}, then some jittering is applied so that points on the plots do not overlap exactly (which can often occur with discrete data, small sample sizes, and if some sites are identical in terms species co-occurence). Please see \code{\link{jitter}} for its implementation.}
     
  \item{a}{Default parameter used in cex. Graphical options are adjusted as \cr \code{par(cex = a, cex.axis = a, cex.lab = a+0.5, cex.main = a+0.5, ...)}. Defaults to 1.}
      
  \item{newplot}{Construct plots on a new graphics device? Defaults to TRUE.}
   
  \item{biplot}{If \code{biplot = TRUE}, then a biplot is construct such that both the latent variables \emph{and} their corresponding coefficients are plotted. Otherwise, only the latent variable scores are plotted. Defaults to \code{TRUE}.}

  \item{ind.spp}{Controls the number of latent variable coefficients to plot if \code{biplot = TRUE}. If \code{ind.spp} is an integer, then only the first \code{ind.spp} "most important" latent variable coefficients are included in the biplot, where "most important" means the latent variable coefficients with the largests L2-norms. Defaults to \code{NULL}, in which case all latent variable coefficients are included in the biplot.}
  
  \item{...}{Additional graphical options to be included in \code{par}.}
}

\details{
Two plots are produced, based on the posterior medians and posterior means of the latent variables respectively. The points in each plot are labeled using the row index of the response matrix \eqn{y}. 

If the fitted model does not include covariates (\code{X = NULL}), the plots can be interpreted in the exactly same manner as unconstrained ordination plots constructed from methods such as Nonmetric Multi-dimensional Scaling (NMDS, Kruskal, 1964) and Correspondence Analysis (CA, Hill, 1974). with multivariate abundance data for instance, where the response matrix \eqn{y} consists of \eqn{n} sites and \eqn{p} species, the ordination plots can be studied to look for possible clustering of sites, location and/or dispersion effects, an arch pattern indicative of some sort species succession over an environmental gradient, and so on.

If \code{biplot = TRUE}, then a biplot is constructed so that both the latent variables and their corresponding coefficients are included in their plot (Gabriel, 1971). The latent variable coefficients are shown in red, and are indexed by the column names of \eqn{y}. The number of latent variable coefficients to plot is controlled by \code{ind.spp}. In ecology for example, often we are only be interested in the "indicator" species, e.g. the species with most represent a particular set of sites or species with the strongest covariation (see Chapter 9, Legendre and Legendre, 2012, for additional discussion). In such case, we can then biplot only the \code{ind.spp} "most important" species, as indicated by the the L2-norm of their latent variable coefficients. Please note that in order to construct the biplot, a scaling is performed on the latent variable coefficients. We cannot ensure that this scaling "works", and it may sometimes produce some funny looking biplots =P In this or other such circumstances, the user may want to draw their own biplots by accessing the latent variables and their coefficients from the fit (please see the help file for the main \code{boral} fitting as to how to access these).
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\references{
\itemize{
\item Gabriel, K. R. (1971). The biplot graphic display of matrices with application to principal component analysis. Biometrika, 58, 453-467.
\item Hill, M. O. (1974). Correspondence analysis: a neglected multivariate method. Applied statistics, 23, 340-354.
\item Kruskal, J. B. (1964). Nonmetric multidimensional scaling: a numerical method. Psychometrika, 29, 115-129.
\item Legendre, P. and Legendre, L. (2012). Numerical ecology, Volume 20. Elsevier.
}
}

\examples{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
n <- nrow(y); p <- ncol(y); 

spider.fit.nb <- boral(y, family = "negative.binomial", num.lv = 2, 
     n.burnin = 10, n.iteration = 100, n.thin = 1,
	row.eff = "fixed", save.model = FALSE, calc.ics = FALSE)

## Three clusters are observable with site 25 roughly in the middle.
## The different species representative of each of the three clusers are 
## also shown in red.
lvsplot(spider.fit.nb) 
}
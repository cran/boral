\name{lvsplot}
\alias{lvsplot}
\title{Plot the latent variables from an boral model} 
\description{Construct a 1-D index plot or 2-D scatterplot of the latent variables, and their corresponding coefficients (i.e. a biplot), from a fitted boral model.}

\usage{lvsplot(x, jitter = FALSE, a = 1, biplot = TRUE, ind.spp = NULL, 
	alpha = 0.5, main = NULL, est = "median",...) 
	}
	
	
\arguments{
  \item{x}{An object for class "boral".}

  \item{jitter}{If \code{jitter = TRUE}, then some jittering is applied so that points on the plots do not overlap exactly (which can often occur with discrete data, small sample sizes, and if some sites are identical in terms species co-occurence). Please see \code{\link{jitter}} for its implementation.}
     
  \item{a}{Default parameter used in cex. Graphical options are adjusted as \cr \code{par(cex = a, cex.axis = a, cex.lab = a+0.5, cex.main = a+0.5, ...)}. Defaults to 1.}
        
  \item{biplot}{If \code{biplot = TRUE}, then a biplot is construct such that both the latent variables \emph{and} their corresponding coefficients are plotted. Otherwise, only the latent variable scores are plotted. Defaults to \code{TRUE}.}

  \item{ind.spp}{Controls the number of latent variable coefficients to plot if \code{biplot = TRUE}. If \code{ind.spp} is an integer, then only the first \code{ind.spp} "most important" latent variable coefficients are included in the biplot, where "most important" means the latent variable coefficients with the largests L2-norms. Defaults to \code{NULL}, in which case all latent variable coefficients are included in the biplot.}
  
  \item{alpha}{A numeric scalar between 0 and 1 that is used to control the relative scaling of the latent variables and their coefficients, when constructing a biplot. Defaults to 0.5, and we typically recommend between 0.45 to 0.55 so that the latent variables and their coefficients are on roughly the same scale.}

  \item{main}{Title for resulting ordination plot. Defaults to \code{NULL}, in which case a "standard" title is used.}

  \item{est}{A choice of either the posterior median (\code{est = "median"}) or posterior mean (\code{est = "mean"}), which are then treated as estimates and the ordinations based off. Default is posterior median.}

  \item{...}{Additional graphical options to be included in \code{par}.}
}

\details{
This function allows an ordination plot to be constructed, based on either the posterior medians and posterior means of the latent variables respectively depending on the choice of \code{est}. The latent variables are labeled using the row index of the response matrix \code{y}. 

If the fitted model did not contain any covariates, the ordination plot can be interpreted in the exactly same manner as unconstrained ordination plots constructed from methods such as Nonmetric Multi-dimensional Scaling (NMDS, Kruskal, 1964) and Correspondence Analysis (CA, Hill, 1974). With multivariate abundance data for instance, where the response matrix \code{y} consists of \eqn{n} sites and \eqn{p} species, the ordination plots can be studied to look for possible clustering of sites, location and/or dispersion effects, an arch pattern indicative of some sort species succession over an environmental gradient, and so on.

If the fitted model did include covariates, then a ``residual ordination" plot is produced, which can be interpreted can offering a graphical representation of the (main patterns of) residual covarations, i.e. covariations after accounting for the covariates. With multivariate abundance data for instance, these residual ordination plots represent could represent residual species co-occurrence due to phylogency, species competition and facilitation, missing covariates, and so on (Warton et al., 2015)

If \code{biplot = TRUE}, then a biplot is constructed so that both the latent variables and their corresponding coefficients are included in their plot (Gabriel, 1971). The latent variable coefficients are shown in red, and are indexed by the column names of \code{y}. The number of latent variable coefficients to plot is controlled by \code{ind.spp}. In ecology for example, often we are only be interested in the "indicator" species, e.g. the species with most represent a particular set of sites or species with the strongest covariation (see Chapter 9, Legendre and Legendre, 2012, for additional discussion). In such case, we can then biplot only the \code{ind.spp} "most important" species, as indicated by the the L2-norm of their latent variable coefficients. 

As with correspondence analysis, the relative scaling of the latent variables and the coefficients in a biplot is essentially arbitrary, and could be adjusted to focus on the sites, species, or put even weight on both (see Section 9.4, Legendre and Legendre, 2012). In \code{lvsplot}, this relative scaling is controlled by the \code{alpha} argument, which basically works by taking the latent variables to a power \code{alpha} and the latent variable coefficients to a power \code{1-alpha}. 

For latent variable models, we are generally interested in "symmetric plots" that place the latent variables and their coefficients on the same scale. In principle, this is achieved by setting \code{alpha = 0.5}, the default value, although sometimes this needs to be tweaked slighlty to a value between 0.45 and 0.55 (see also the \code{corresp} function in the \code{MASS} package that also produces symmetric plots, as well as Section 5.4, Borcard et al., 2011 for more details on scaling).
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\references{
\itemize{
\item Borcard et al. (2011). Numerical Ecology with R. Springer.
\item Gabriel, K. R. (1971). The biplot graphic display of matrices with application to principal component analysis. Biometrika, 58, 453-467.
\item Hill, M. O. (1974). Correspondence analysis: a neglected multivariate method. Applied statistics, 23, 340-354.
\item Kruskal, J. B. (1964). Nonmetric multidimensional scaling: a numerical method. Psychometrika, 29, 115-129.
\item Legendre, P. and Legendre, L. (2012). Numerical ecology, Volume 20. Elsevier.
\item Warton et al. (2015). So Many Variables: Joint Modeling in Community Ecology. Trends in Ecology and Evolution, in review.
}
}

\examples{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
n <- nrow(y); p <- ncol(y); 

spider.fit.nb <- boral(y, family = "negative.binomial", num.lv = 2, 
     n.burnin = 10, n.iteration = 100, n.thin = 1,
	row.eff = "fixed", calc.ics = FALSE)

lvsplot(spider.fit.nb) 
}
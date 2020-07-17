\name{about.lvs}

\alias{about.lvs}

\docType{package}

\title{Correlation structure for latent variables}

\description{
\Sexpr[results=rd, stage=render]{lifecycle::badge("stable")}

This help file provides more information how (non-independence) correlation structures can be assumed for latent variables.}

\details{
In the main boral function, when latent varaibles are included, the default option is to assume that the latent variables are independent across the rows (sites) of the response matrix i.e., \code{lv.type = "independent"}. That is, \eqn{\bm{z}_i \sim N(\bm{0},\bm{I}_d)} where \code{d = num.lv}. This is useful when we want to model between species correlations (is a parsimonious manner), but it does make an assumption that sites are independent. 

If one \emph{a-priori} believes that the sites are, in fact, correlated e.g., due to spatial correlation, and that it cannot be sufficiently well accounted for by row effects (see comment below), then we can account for this by assuming a non-independence correlation structure for the the latent variables across sites. Note however we continue to assume that the \eqn{d} latent variables are still independent of one another. That is, if we let \eqn{\bm{z}_i = (z_{i1}, \ldots, z_{id})}, then we assume that for \eqn{l = 1,\ldots,d},

\deqn{(z_{1l}, z_{2l}, \ldots, z_{nl}) \sim N(\bm{0}, \bm{\Sigma}),}

where \eqn{\bm{\Sigma}} is some correlation matrix. When \eqn{\bm{\Sigma} = \bm{I}_n} then we are back in the independence case. However, if we allow for the off-diagonals to be non-zero, then we the latent variables to be correlated, \eqn{\Sigma_{ij} = Cov(z_{il}, z_{jl})}. This in turn induces correlation across sites and species i.e., two species at two different sites are now correlated because of the correlation across sites. 

While there are fancier structures and attempts at accounting for correlations between sites (Cressie and Wikle, 2015), in boral we assume relatively simple structures. Specifically, we can assume that sites further away are less correlated, and so \eqn{\Sigma} can be characterized based on a distance matrix \code{distmat} and associated spatial covariance parameters which require estimation. Indeed, such simple spatial latent variable models have become rather popular in community ecology of late, at least as a first attempt at accounting for spatial (and also temporal) correlation e.g., Thorson et al., (2015, 2016); Ovaskainen et al., (2017). 

At the moment, several correlation structures are permitted. Let \eqn{D_{ij}} denote the distance between site \eqn{i} and {j} i.e., entry \eqn{(i,j)} in \code{distmat}. Also, let \eqn{(\vartheta_1,\vartheta_2)} denote the two spatial covariance parameters (noting that the second parameter is not required for some of structures). Then we have: 1) \code{lv.type = "exponential"} such that \eqn{\Sigma_{ij} = \exp(-D_{ij}/\vartheta_1)}; 2) \code{lv.type = "squared.exponential"}, such that \eqn{\Sigma_{ij} = \exp(-D_{ij}/\vartheta_1^2)}; 3) \code{lv.type = "power.exponential"}, such that \eqn{\Sigma_{ij} = \exp(-(D_{ij}/\vartheta_1)^{\vartheta_2})} where \eqn{\vartheta_1 \in (0,2]} ; 4) \code{lv.type = "spherical"}, such that \eqn{(D_{ij} < \vartheta_1)*(1 - 1.5*D_{ij}/\vartheta_1 + 0.5*(D_{ij}/\vartheta_1)^3)}. We refer the reader to the \code{geoR} and the function \code{cov.spatial} for more, simple information on spatial covariance functions (Ribeiro Jr and Diggle, 2016). 

It is important to keep in mind that moving away from an independence correlation structure for the latent variables \emph{massively} increases computation time for MCMC sampling (and indeed any estimation method for latent variable models). Given JAGS is not the fastest of methods when it comes to MCMC sampling, then one should be cautious about moving away from indepndence. For example, if you \emph{a-priori} have a nested experimental design which is inducing spatial correlation, then it is much faster and more effective to include (multiple) row effects in the model to account for this spatial correlation instead. 
}


\references{
\itemize{
\item Cressie, N. and Wikle, C. K. (2015) Statistics for spatio-temporal data. John Wiley & Sons.

\item Ovaskainen, O., Tikhonov, G. Norberg, A. Guillaume Blanchet, F. Duan, L. Dunson, D. Roslin, T., and Abrego, N. (2017). How to make more out of community data? A conceptual framework and its implementation as models and software. Ecology Letters, 20, 561-576.

\item Ribeiro Jr, P. J., and Diggle P. J., (2016). geoR: Analysis of Geostatistical Data. R package version 1.7-5.2. \url{https://CRAN.R-project.org/package=geoR}.

\item Thorson, J. T., Ianelli, J. N., Larsen, E. A., Ries, L., Scheuerell, M. D., Szuwalski, C., and Zipkin, E. F. (2016). Joint dynamic species distribution models: a tool for community ordination and spatio-temporal monitoring. Global Ecology and Biogeography, 25, 1144-1158

\item Thorson, J. T., Scheuerell, M. D., Shelton, A. O., See, K. E., Skaug, H. J., and Kristensen, K. (2015). Spatial factor analysis: a new tool for estimating joint species distributions and correlations in species range. Methods in Ecology and Evolution, 6, 627-63
}
}

\author{
\packageAuthor{boral}

Maintainer: \packageMaintainer{boral}
}

\seealso{
\code{\link{boral}} for the main boral fitting function.
}

\examples{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
X <- scale(spider$x)
n <- nrow(y)
p <- ncol(y)

## NOTE: The examples below and taken directly from the boral help file

example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
     n.thin = 1)

testpath <- file.path(tempdir(), "jagsboralmodel.txt")

\dontrun{
## Example 2d - model with environmental covariates and 
##  two structured latent variables using fake distance matrix
fakedistmat <- as.matrix(dist(1:n))
spiderfit_lvstruc <- boral(y, X = X, family = "negative.binomial", 
    lv.control = list(num.lv = 2, type = "exponential", distmat = fakedistmat), 
     mcmc.control = example_mcmc_control, model.name = testpath)

summary(spiderfit_lvstruc)
    }

}

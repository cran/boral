\name{get.residual.cor}
\alias{get.residual.cor}
\title{Extract residual correlations from boral models}
\description{Calculates the residual correlation from models that include latent variables.}

\usage{
get.residual.cor(y, fit.mcmc, est = "median")
}

\arguments{
  \item{y}{The response matrix that the boral model was fitted to.}
  
     \item{fit.mcmc}{All MCMC samples for the fitted boral model, as obtained from JAGS. These can be extracted by fitting an boral model using \code{\link{boral}} with \code{save.model = TRUE}, and then applying \code{as.mcmc} on the \code{jags.model} component of the output.} 

     \item{est}{A choice of either the posterior median (\code{est = "median"}) or posterior mean (\code{est = "mean"}), which are then treated as estimates and the fitted values are calculated from. Default is posterior median.}

}

\details{
In models with latent variables, the residual correlation between two columns \eqn{j} and \eqn{j'} of \eqn{\bm{y}} is calculated based on the linear predictors \eqn{\bm{z}'_i\bm{\theta}_j} and \eqn{\bm{z}'_i\bm{\theta}_{j'})}, where \eqn{\bm{z}_i} are the latent variables, and \eqn{\bm{\theta}_j} are column-specific coefficients relating to the latent variables (see also the help file for \code{\link{boral}}).

For multivariate abundance data, the inclusion of latent variables provides a parsimonious method of accounting for correlation between species. Specfically, the linear predictor,

\deqn{\theta_{0j} + \bm{x}'_i\bm{\beta}_j + \bm{z}'_i\bm{\theta}_j}

is normally distributed with a residual covariance/correlation matrix induced by \eqn{\bm{z}'_i\bm{\theta}_j}. A strong residual covariance/correlation matrix between two species can then be interpreted as evidence of species interaction (e.g., facilitation or competition), missing covariates, as well as any additional species correlation not accounted for by shared environmental responses (see also Pollock et al., 2014, for residual correlation matrices in the context of Joint Species Distribution Models).
}

\value{
A \eqn{p \times p} residual correlation matrix matrix based on posterian median or mean estimators of the latent variables and coefficients.
}

\references{
\itemize{
\item Pollock, L. J., Tingley, R., Morris, W. K., Golding, N., O'Hara, R. B., Parris, K. M., Vesk, P. A., and McCarthy, M. A. (2014). Understanding co-occurrence by modelling species simultaneously with a Joint Species Distribution Model (JSDM). Methods in Ecology and Evolution, 5, 397-406.
}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\note{
Residual correlatiom matrices are reliably modeled only with two or more latent variables i.e., \code{num.lv > 1} when fitting the model using \code{boral}.  
}

\examples{
library(mvabund) ## Load a dataset from the mvabund package
library(corrplot) ## For plotting correlations
data(spider)
y <- spider$abun
n <- nrow(y); p <- ncol(y); 
    
spider.fit.nb <- boral(y, X = spider$X, family = "negative.binomial", num.lv = 2, 
     site.eff = FALSE, n.burnin = 10, n.iteration = 100, 
     n.thin = 1, save.model = TRUE, calc.ics = FALSE)

res.cors <- get.residual.cor(y, fit.mcmc = as.mcmc(spider.fit.nb$jags.model)[[1]])

corrplot(res.cors, title = "Residual correlations", type = "lower")
}
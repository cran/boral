\name{get.residual.cor}
\alias{get.residual.cor}
\title{Extract residual correlations from boral models}
\description{Calculates the residual correlation from models that include latent variables.}

\usage{
get.residual.cor(object, est = "median")
}

\arguments{
  \item{object}{An object for class "boral".}
  
  \item{est}{A choice of either the posterior median (\code{est = "median"}) or posterior mean (\code{est = "mean"}), which are then treated as estimates and the fitted values are calculated from. Default is posterior median.}

}

\details{
In models with latent variables, the residual correlation matrix is calculated based on the matrix of latent variables regression coefficients formed by stacking the rows of \eqn{\bm{\theta}_j}. That is, if we denote \eqn{\bm{\Theta} = (\bm{\theta}_1 \ldots \bm{\theta}_p)'}, then the residual correlation matrix is calculated based on \eqn{\bm{\Theta}\bm{\Theta}'} (see also the help file for \code{\link{boral}}).

For multivariate abundance data, the inclusion of latent variables provides a parsimonious method of accounting for correlation between species. Specifically, the linear predictor,

\deqn{\beta_{0j} + \bm{x}'_i\bm{\beta}_j + \bm{z}'_i\bm{\theta}_j}

is normally distributed with a residual covariance given by \eqn{\bm{\Theta}\bm{\Theta}'}. A strong residual covariance/correlation matrix between two species can then be interpreted as evidence of species interaction (e.g., facilitation or competition), missing covariates, as well as any additional species correlation not accounted for by shared environmental responses (see also Pollock et al., 2014, for residual correlation matrices in the context of Joint Species Distribution Models).
}

\value{
A \eqn{p \times p} residual correlation matrix matrix based on posteriori median or mean estimators of the latent variables and coefficients.
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
    
spider.fit.nb <- boral(y, X = spider$x, family = "negative.binomial", num.lv = 2, 
     row.eff = FALSE, n.burnin = 10, n.iteration = 100, 
     n.thin = 1, save.model = TRUE, calc.ics = FALSE)

res.cors <- get.residual.cor(spider.fit.nb)

corrplot(res.cors, title = "Residual correlations", type = "lower", diag = FALSE)
}
}
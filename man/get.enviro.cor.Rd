\name{get.enviro.cor}
\alias{get.enviro.cor}
\title{Extract covariances and correlations due to shared environmental responses from boral models}
\description{Calculates the correlation between columns of the response matrix, due to similarities in the response to explanatory variables (i.e., shared environmental response)}

\usage{
get.enviro.cor(object, est = "median", prob = 0.95)
}

\arguments{
  \item{object}{An object for class "boral".}

  \item{est}{A choice of either the posterior median (\code{est = "median"}) or posterior mean (\code{est = "mean"}), which are then treated as estimates and the fitted values are calculated from. Default is posterior median.}

  \item{prob}{A numeric scalar in the interval (0,1) giving the target probability coverage of the intervals, by which to determine whether the correlations are "significant". Defaults to 0.95.}   
}

\details{
In both independent response and correlated response models, where the each of the columns of the response matrix \code{y} are fitted to a set of explanatory variables given by \code{X}, the covariance and thus between two columns \eqn{j} and \eqn{j'} due to similarities in their response to the model matrix is calculated based on the linear predictors \eqn{\bm{x}^T_i\bm{\beta}_j} and \eqn{\bm{x}^T_i\bm{\beta}_{j'})}, where \eqn{\bm{\beta}_j} are column-specific coefficients relating to the explanatory variables (see also the help file for \code{\link{boral}}).

For multivariate abundance data, the correlation calculated by this function can be interpreted as the correlation attributable to similarities in the environmental response between species. Such correlation matrices are discussed and found in Ovaskainen et al., (2010), Pollock et al., 2014.
}

\value{
A list with the following components:
\item{cor}{A \eqn{p \times p} correlation matrix based on model matrix and the posterior or mean estimators of the associated regression coefficients.}
\item{sig.cor}{A \eqn{p \times p} correlation matrix containing only the ``significant" correlations whose 95\% highest posterior interval does not contain zero. All non-significant correlations are zero to zero.}
\item{cov}{A \eqn{p \times p} covariance matrix based on model matrix and the posterior or mean estimators of the associated regression coefficients.}
}

\references{
\itemize{
\item Ovaskainen et al. (2010). Modeling species co-occurrence by multivariate logistic regression generates new hypotheses on fungal interactions. Ecology, 91, 2514-2521.
\item Pollock et al. (2014). Understanding co-occurrence by modelling species simultaneously with a Joint Species Distribution Model (JSDM). Methods in Ecology and Evolution, 5, 397-406.
}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\seealso{
\code{\link{get.residual.cor}}, which calculates the residual correlation matrix for boral models involving latent variables.
}

\examples{
\dontrun{
library(mvabund) ## Load a dataset from the mvabund package
library(corrplot) ## For plotting correlations
data(spider)
y <- spider$abun
X <- scale(spider$x)
n <- nrow(y); p <- ncol(y); 
    
spider.fit.nb <- boral(y, X = X, family = "negative.binomial", 
     save.model = TRUE)

enviro.cors <- get.enviro.cor(spider.fit.nb)

corrplot(enviro.cors$sig.cor, title = "Shared response correlations", 
	type = "lower", diag = FALSE, mar = c(3,0.5,2,1), tl.srt = 45)
}
}
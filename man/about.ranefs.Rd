\name{about.ranefs}

\alias{about.ranefs}

\docType{package}

\title{Including response-specific random intercepts in boral}

\description{
\Sexpr[results=rd, stage=render]{lifecycle::badge("stable")}

This help file provides more information regarding the how response-specific random intercepts can be included to account for sampling design, induce correlation betweens observational units that are different to each response etc...}

\details{
As of version 2.0, it is now possible to include response-specific random intercepts in a model. There may be a number of reasons why a user may be to do this, but the most common reasons (at least in community ecology) would be to account sampling design on a response-specific basis, and more generally if there \emph{a-priori} knowledge of clustering between observational units and so random intercepts are to be included to account for such potential within-cluster correlation. 

Alternatively, if you are familiar with including random row effects in a model, then one may think of response-specific random intercepts as similar to this, except the row effects are now response-specific specific rather than a common across all responses. In doing so, response-specific random intercepts are more flexible and allow the induced correlations (and thus the standard deviations) to be on a per-response basis, although this comes at the expense of some random intercepts being poorly estimates for responses with relatively little information in their observations e.g., rare species. 

In for example, the formulation for the correlated response model 

\deqn{g(\mu_{ij}) = \alpha_i + \beta_{0j} + \bm{x}^\top_i\bm{\beta}_j + \bm{z}^\top_i\bm{b}_j; \quad i = 1,\ldots,n; j = 1,\ldots,p,}

the \eqn{\bm{z}^\top_i\bm{b}_j} denote response-specific random intercepts included, with \eqn{bm{b}_j} denoting the response-specific random intercepts, and \eqn{\bm{z}\top_i} denoting the corresponding design vector for observational unit \eqn{i}, and are "constructed" based on the input \code{ranef.ids}. 

Akin to other packages such as \code{lme4} or \code{glmmTMB}, all random intercepts are assumed to be normally distributed with the corresponding variance components are assumed to be response-specific. In fact, if we consider the simpler independent response model with no row effects 

\deqn{g(\mu_{ij}) = \beta_{0j} + \bm{x}^\top_i\bm{\beta}_j + \bm{z}^\top_i\bm{b}_j; \quad i = 1,\ldots,n; j = 1,\ldots,p,}

then the above would be equivalent to fitting a generalized linear mixed model (GLMM) to each response separately, for which packages such as \code{lme4} and \code{glmmTMB} are well suited for. In that sense, \code{boral} can fit some models, but can also incorporate row effects and/or latent variables to induce correlation between responses.

Perhaps not surprisingly, the way response-specific random intercepts are included is very similar to how row effects are included in the model. Specifically, the argument \code{ranef.ids} identifies the number of random intercepts to be included and how each observational unit maps to a cluster. \code{ranefs.ids} is a matrix with the number of rows equal to the number of rows in the response matrix, and the number of columns equal to the number of random intercepts to be included in the model. Element \eqn{(i,j)} indicates the cluster ID of row \eqn{i} in the response matrix for random intercept eqn{j}. Examples of its use are provided in the help file below.

After the model is fitted, for each set of random intercepts included (i.e., each column of \code{ranef.ids}), estimates along with HPD intervals of the response-specific standard deviations of the corresponding random effects distributions, as well of the random intercept predictions are returned. These can then be visualized for example using the \code{\link{ranefsplot}} function, or analyzed as appropriately by the user.
}


\section{Warnings}{
\itemize{
\item It is usually not recommended to have both random row effects and response-specific random intercepts simultaneously in the same model, unless they are were at different levels of of the data)
}
}

\author{
\packageAuthor{boral}

Maintainer: \packageMaintainer{boral}
}


\seealso{
\code{\link{boral}} for the main boral fitting function, 
\code{\link{ranefsplot}} for horizontal line or "caterpillar plot" of the response-specific random effects predictons (if applicable).
}

\examples{
\dontrun{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
X <- scale(spider$x)
n <- nrow(y)
p <- ncol(y)

## NOTE: The example below is taken directly from the boral help file

example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
     n.thin = 1)

testpath <- file.path(tempdir(), "jagsboralmodel.txt")

## Example 2e - Similar to 2c, but we will species-specific random intercepts
##   for the seven regions (with row effects in the model)
spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
    ranef.ids = data.frame(region = rep(1:7,each=4)), 
    mcmc.control = example_mcmc_control, model.name = testpath) 

spiderfit_nb$ranef.coefs.median
spiderfit_nb$ranef.sigma.median
}
}


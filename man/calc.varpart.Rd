\name{calc.varpart}

\alias{calc.varpart}

\docType{package}

\title{Variance partitioning for a latent variable model}

\description{For each response (species), partition the variance of the linear predictor into components associated with (groups of) the covariates, the latent variables, any row effects. If traits are also included in the model, then it also calculates an R-squared value for the proportion of the variance in the environmental response (due to the covariates) which can be explained by traits.}

\usage{
calc.varpart(object, groupX = NULL)
}

\arguments{
  \item{object}{An object of class "boral".}

  \item{groupX}{A vector of group indicator variables, which allows the variance partitioning to be done for groups of covariates (including the intercept) i.e., how much of the total variation does a certain subset of the covariates explain. Defaults to \code{NULL}, in whih case all the covariates are treated as single group.}
}


\details{
As an alternative to looking at differences in trace of the residual covariance matrix (Hui et al., 2014; Warton et al., 2015), an alternative way to quantify the amount of variance explained by covariates, traits, row effects, is to perform a variance decomposition of the linear predictor of a latent variable model. In particular, for a general model the linear predictor for response \eqn{j = 1,\ldots,p} at row \eqn{i = 1,\ldots,n} is given by

\deqn{\eta_{ij} = \alpha_i + \beta_{0j} + \bm{x}^\top_i\bm{\beta}_j, + \bm{z}^\top_i\bm{\theta}_j,}

where \eqn{\beta_{0j} + \bm{x}^\top_i\bm{\beta}_j} is the component of the linear predictor due to the covariates \code{X} plus an intercept, \eqn{\bm{z}^\top_i\bm{\theta}_j} is the component due to the latent variables, and \eqn{\alpha_i} is the component due to one or more fixed or random row effects. The regression coefficients \eqn{\bm{\beta}_j} may be further as random effects and regressed against traits; please see \code{\link{about.traits}} for further information on this. 

For the response, a variation partitioning of the linear is performed by calculating the variance due to component in \eqn{\eta_{ij}} and then rescaling them to ensure that they sum to one. The general details of this type of variation partitioning is given in Ovaskainen et al., (2017); see also Nakagawa and Schielzeth (2013) for R-squared and proportion of variance explained in the case of generalized linear mixed model. In brief, for response \eqn{j = 1,\ldots,p}: 1) the variance due to the \code{X} covariates and intercept is given by the variance of \eqn{\beta_{0j} + \bm{x}^\top_i\bm{\beta}_j} calculated across the \eqn{n} rows; 2) the variance due the latent variables is given by the diagonal elements of \eqn{\bm{\theta}^\top_j\bm{\theta}_j}; 3) the variance due to (all) the random effects is given by variance of \eqn{\alpha_i} calculated across the \eqn{n} rows for fixed row effects (\code{row.eff = "fixed"}), and given by the (sum of the) variance \eqn{\sigma^2_{\alpha}} for random row effects (\code{row.eff = "random"}). After scaling, we can then obtain the proportion of variance for each response which is explained by the variance components. These proportions are calculated for each MCMC sample and then acrossed them to calculate a posterior mean variance partitioning. 

If \code{groupX} is supplied, the variance due to the \code{X} covariates is done based on subsets of \code{X} variables (including the intercept) as identified by code{groupX}, and then rescaled correspondingly. This is useful if one was to, for example, quantify the proportion of variation in each species which is explained by each \code{X} covariate.

If a fitted model also containing traits, which are included to help explain/mediate differences in species environmental responses, then the function calculates \eqn{R^2} value for the proportion of variance in the \code{X} variables which is explained by the traits. In brief, this is calculated based the correlation between \eqn{\beta_{0j} + \bm{x}^\top_i\bm{\beta}_j} and \eqn{\tau_{0j} + \bm{x}^\top_i\bm{\tau}_j}, where \eqn{\tau_{0j}} and \eqn{\bm{\tau}_j} are the ``predicted" values of the species coefficients based on values i.e., \eqn{\tau_{0j} = \kappa_{01} + \bm{traits}^\top_j\bm{\kappa}_1} and \eqn{\tau_{jk} = \kappa_{0k} + \bm{traits}^\top_j\bm{\kappa}_k} for element \eqn{k} in \eqn{\bm{\tau}_j}. 
}

\value{
A list containing the following components if applicable:
     \item{varpart.X}{Vector containing the proportion of variance (in the linear predictor) for each response which is explained by \code{X}.}
     
     \item{varpart.lv}{Vector containing the proportion of variance (in the linear predictor) for each response which is explained by the latent variables.}
     
     \item{varpart.row}{Vector containing the proportion of variance (in the linear predictor) for each response which is explained by the row effects.}
     
     \item{R2.traits}{Vector containing the proportion of variance due to the covariates for each response, which can be explained by traits for each response.}
}

\references{
\itemize{
\item Nakagawa, S., and Schielzeth, H. (2013). A general and simple method for obtaining R2 from generalized linear mixed-effects models. Methods in Ecology and Evolution 4, 133-142.

\item Ovaskainen, et al. (2017). How to make more out of community data? A conceptual framework and its implementation as models and software. Ecology Letters 20, 561-576.

\item Hui et al. (2014). Model-based approaches to unconstrained ordination. Methods in Ecology and Evolution, 6, 399-411.

\item Warton et al. (2015). So Many Variables: Joint Modeling in Community Ecology. Trends in Ecology and Evolution, 30, 766-779.
}
}

\author{
\packageAuthor{boral}

Maintainer: \packageMaintainer{boral}
}

\section{Warnings}{
There is considerable controversy over exactly what quantities such as R-squared and proportion of variance explained are in the case mixed models and latent variable models, and how they can interpreted e.g., what is considered a high value for the proportion of variance by \code{X} variables, is it consistent with whether the coefficients are significantly different from zero or not; see for instance \href{https://stats.stackexchange.com/questions/111150/calculating-r2-in-mixed-models-using-nakagawa-schielzeths-2013-r2glmm-me}{R2 controversy}.

When reporting these values, researchers should be at least aware of this and that there are multiple ways of manufacturing such quantities with no single best approach e.g., using relative changes in trace of the residual covariance matrix, relative changes in marginal and conditional log-likelihoods are other possible approaches. 
}

\examples{
\dontrun{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
X <- scale(spider$x)
n <- nrow(y)
p <- ncol(y)

## NOTE: The values below MUST NOT be used in a real application;
## they are only used here to make the examples run quick!!!
example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
     n.thin = 1)

testpath <- file.path(tempdir(), "jagsboralmodel.txt")


## Example 1 - model with X variables, two latent variables, and no row effects
spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
     lv.control = list(num.lv = 2), 
     save.model = TRUE, mcmc.control = example_mcmc_control,
     model.name = testpath)

## Partition variance for each species into that explained by covariates 
## and by the latent variables
dovar <- calc.varpart(spiderfit_nb)

## Consider the intercept and first two covariates in X as one group, 
## and remaining four covariates in X as another group, 
## then partition variance for each species based on these groups.
dovar <- calc.varpart(spiderfit_nb, groupX = c(1,1,1,2,2,2,2))


## Example 2 - model fitted to count data, no site effects, and
## two latent variables, plus traits included to explain environmental responses
data(antTraits)
y <- antTraits$abun
X <- as.matrix(scale(antTraits$env))
## Include only traits 1, 2, and 5
traits <- as.matrix(antTraits$traits[,c(1,2,5)])
example_which_traits <- vector("list",ncol(X)+1)
for(i in 1:length(example_which_traits)) 
     example_which_traits[[i]] <- 1:ncol(traits)
## Just for fun, the regression coefficients for the second column of X,
## corresponding to the third element in the list example_which_traits,
## will be estimated separately and not regressed against traits.
example_which_traits[[3]] <- 0

fit_traits <- boral(y, X = X, traits = traits, which.traits = example_which_traits, 
    family = "negative.binomial", mcmc.control = example_mcmc_control, 
    save.model = TRUE, model.name = testpath)

## Partition variance for each species due to covariates in X 
## and latent variables. Also calculate proportion of variance 
## due to the covariates which can be explained by traits 
dovar <- calc.varpart(fit_traits)
}

}

\name{boral}
\alias{boral}
\alias{boral.default}
\alias{print.boral}
\title{Fitting boral (Bayesian Ordination and Regression AnaLysis) models}
\description{Bayesian ordination and regression models for analyzing multivariate data in ecology. Three ``types" of models may be fitted: 1) With covariates and no latent variables, \code{boral} fits independent column GLMs such that the columns of \eqn{y} are assumed to be independent; 2) With no covariates, \code{boral} fits a purely latent variable model to perform model-based unconstrained ordination (Hui et al., 2014); 3) With covariates and latent variables, \code{boral} fits correlated column GLMs, with latent variables accounting for correlations between the columns of \eqn{y}.}

\usage{
boral(y, ...)

\method{boral}{default}(y, X = NULL, family, trial.size = NULL, num.lv = 0, site.eff = FALSE, 
n.burnin = 4000, n.iteration = 24000, n.thin = 5, save.model = FALSE, 
	seed = 123, calc.ics = TRUE, hypparams = c(100,100), ...)

\method{print}{boral}(x, ...)
}

\arguments{
  \item{y}{A \eqn{n} times \eqn{p} response matrix of multivariate data e.g., counts, binomial or Bernoulli responses, continuous response, and so on. In ecology for instance, \eqn{n} is the number of sites and \eqn{p} is the number of species. Any categorical (multinomial) responses \bold{must} be converted to integer values. For ordinal data, the minumum level of \eqn{y} must be 1 instead of 0.}
  
  \item{X}{An \eqn{n} times \eqn{d} model matrix of covariates, which can be included as part of the boral model. Defaults to \code{NULL}, in which case no model matrix was used. No intercept column should be included in \code{X}.}  

  \item{x}{An object for class "boral".}

\item{family}{Either a single element, or a vector of length equal to the number of columns in \eqn{y}. The former assumes all columns of \eqn{y} come from this distribution. The latter option allows for different distributions for each column of \eqn{y}. Elements can be one of "binomial" (with logit link), "poisson" (with log link), "negative.binomial" (with log link), "normal" (with identity link), "lnormal" for lognormal (with log link), "tweedie" (with log link), "ordinal" which performs cumulative logit regression. For the negative binomial distribution, the variance is parameterized as \eqn{Var(y) = \mu + \phi \mu^2}, where \eqn{\phi} is the column-specific overdispersion parameter. For the normal distribution, the variance is parameterized as \eqn{Var(y) = \phi}, where \eqn{\phi} is the column-specific variance. For the tweedie distribution, the variance is parameterized as \eqn{Var(y) = \phi \mu^p} where \eqn{\phi} is the column-specific dispersion parameter and \eqn{p} is a power parameter common to all columns assumed to be tweedie, with \eqn{1 < p < 2}. 

All columns assumed to have ordinal responses are constrained to have the same cutoffs points, with a column-specific intercept to account for differences between the columns (please see \emph{Details} for formulation). 
}

 \item{trial.size}{Either equal to NULL, a single element, or a vector of length equal to the number of columns in \eqn{y}. If a single element, then all columns assumed to be binomially distributed will have trial size set to this. If a vector, different trial sizes are allowed in each column of \eqn{y}.}

 \item{num.lv}{Number of latent variables to fit. Can take any non-negative integer value. Defaults to 0.} 

 \item{site.eff}{A logical value indicating whether to include row effects or not. Defaults to FALSE.}    

 \item{n.burnin}{Length of burnin i.e., the number of iterations to discard at the beginning of the MCMC sampler. Default is 4000.}

 \item{n.iteration}{Number of iterations including burnin. Default is 24000}

 \item{n.thin}{Thinning rate. Must be a positive integer. Default is 4. With the default values of \code{n.burin, n.iteration} and \code{n.thin}, this leads to a total of 5000 MCMC samples.}

 \item{save.model}{A logical value indicating whether to save the JAGS model file as a text file (called "jagsboralmodel.txt") in the current working directory, as well as the MCMC samples from the call to JAGS. If saved, various functions available in the \code{coda} package can be applied to the MCMC samples. Note MCMC samples can take up a lot of memory. Defaults to \code{FALSE}.}

 \item{seed}{seed for JAGS sampler. Within \code{boral}, a \code{set.seed(seed)} command is run immediately before starting the MCMC sampler. Defaults to the value 123.}

\item{calc.ics}{A logical values indicating whether to return various information criteria values, which could be used to perform model selection (see \code{\link{get.measures}}). Defaults to \code{TRUE}.}

 \item{hypparams}{Vector of 2 hyperparameters in the Bayesian estimation. The first hyperparameter is the variance for the normal priors of all coefficients i.e., row effects, column-specific intercepts, column-specific regression coefficients corresponding to all latent variables and covariates, cutoff points for ordinal data. The second hyperparameter is the maximum of the uniform prior used for the dispersion parameters i.e, for the negative binomial, normal, and tweedie distributions. Defaults to 100 for both elements. The common power parameter in the tweedie distributon is assumed to have uniform prior from 1 to 2.}

 \item{...}{Not used.}
}

\details{
The boral package is designed to fit three types models which may be useful in ecology (and probably outside of ecology as well). 

\bold{Independent column models:} \code{boral} also allows explanatory variables to be entered into the model via a model matrix \code{X}. This model matrix can contain anything the user wants, provided factors have been parameterized as dummy variables. It should \emph{not} include an intercept column. 

Without latent variables i.e., \code{num.lv = 0}, \code{boral} fits separate GLMs to each column of \eqn{y}, where the columns are assumed to be independent. That is,

\deqn{g(\mu_{ij}) = \alpha_i + \theta_{0j} + \bm{x}'_i\bm{\beta}_j; \quad i = 1,\ldots,n; j = 1,\ldots,p,}

where the mean response for element \eqn{(i,j)}, denoted as \eqn{mu_{ij}}, is regressed against the covariates \eqn{\bm{x}_i} via a link function \eqn{g(\cdot)}. \eqn{\theta_{0j}} and \eqn{\bm{\beta}_j} are the column-specific intercepts and regression coefficients respectively, while \eqn{\alpha_i} is an optional row effect (if \code{site.eff = TRUE}). Fitting the above type of model is thus more or less a Bayesian analog of the \code{manyglm} function in the \code{mvabund} package (Wang et al., 2013). Unlike \code{manyglm} though, different distributional assumptions are allowed however for each column of \eqn{y}, and row effects can be added as a type of "row-standardization". 

\bold{A not-so-brief tangent on distributions:} In the event different responses are collect for different columns, e.g., some columns of \eqn{y} are counts, while other columns are presence-absence, one can specify different distributions for each column. Aspects such as variable selection, residual analysis, and plotting of the latent variables are not affected by having different distributions. Naturally though, one has to be more careful with interpretation of the row effects \eqn{\alpha_i} and latent variables \eqn{\bm{z}_i}, with different scalings due to different link functions used. A situation where different distributions may prove useful is when \eqn{y} is a species-traits matrix, where each row is a species and each column is a trait such as specific leaf area, presence of thorns and so on. In this case, traits could be of different response types, and the goal perhaps is to perform unconstrained ordination to look for patterns between species on an underlying trait surface e.g., a defense index for a species (Moles et al., 2013; see the discussion below on how to perform model-based unconstrained ordination).

For multivariate abundance data in ecology (also known as community composition data, Legendre and Gallagher, 2001), species counts are often overdispersed. Using a negative binomial distribution (\code{family = "negative.binomial"}) to model the counts usually helps to account for this overdispersion. For non-negative continuous data such as biomass, the lognormal and tweedie distributions may be used (Foster and Bravington, 2013). Note though a a common power paramter is used for tweedie columns -- there almost always insufficient information to model column-specific power parameters. Normal responses are also implemented, just in case you encounter some normal stuff in ecology (pun intended)!  

For ordinal response columns, cumulative logit regression is used (also known as proportional odds regression, Agresti, 2010). \code{boral} assumes all ordinal columns are measured using the same scale i.e., all columns have the same number of theoretical levels. The number of levels is then assumed to be given by the maximum value from all the ordinal columns of \eqn{y}. Thus, the parameterization of proportional odds regression is such that all ordinal columns have the \emph{same} cutoff points, \eqn{\bm{\tau}}, while the column-specific intercept effect, \eqn{\theta_{0j}}, allows for deviations away from these common cutoffs. That is, 

\deqn{logit(P(y_{ij} \le k)) = \tau_k + \theta_{0j} + \ldots,}

where \eqn{logit(\cdot)} is the logit function, \eqn{P(y_{ij} \le k)} is the cumulative probability of element \eqn{y_{ij}} being less than or equal to level \eqn{k}, \eqn{\tau_k} is the cutoff linking levels \eqn{k} and \eqn{k+1} (and which are increasing in \eqn{k}), \eqn{\theta_{0j}} are the column effects, and "\eqn{\ldots}" generically represents what else is included in the model e.g., latent variables and related coefficients. A sum-to-zero constraint is imposed on the \eqn{\theta_{0j}}'s of all ordinal columns, to ensure model identifiability. 

The parameterization above is useful for modeling multivariate abundance data in ecology. When ordinal responses are recorded, usually the same ordinal scale is applied to all species e.g., level 1 = not there, level 2 = a bit there, level 3 = lots there, level 4 = everywhere! The quantity \eqn{\tau_k} can be interpreted as this common ordinal scale, while \eqn{\theta_{0j}} allows for deviations away from these to account for differences in species prevalence. Admittedly, the current implementation of \code{boral} for ordinal data can be rather slow. 

\bold{Purely latent variable models:} If no explantory variables are included and \code{num.lv} > 0, \code{boral} fits a purely latent variable model to perform model-based unconstrained ordination (Hui et al., 2014),

\deqn{g(\mu_{ij}) = \alpha_i + \theta_{0j} + \bm{z}'_i\bm{\theta}_j,}

where instead of measured covariates, we now have a vector of latent variables \eqn{\bm{z}_i} with \eqn{\bm{\theta}_j} being the column-specific coefficients relating to these latent variables. Unconstrained ordination is used for visualizing multivariate data in a low-dimensional space, without reference to explanatory variables (see Chapter 9, Legendre and Legendre, 2012). Typically, \code{num.lv} = 1 or 2 latent variables is used, allowing the latent variables to plotted (using \code{\link{lvsplot}}, for instance). The resulting plot can be interpreted in the same manner as plots from Nonmetric Multi-dimensional Scaling (NMDS, Kruskal, 1964) and Correspondence Analysis (CA, Hill, 1974), for example. For instance, with multivariate abundance data, unconstrained ordination is used to visualize the relationships between sites in terms of species abundance or composition. The column-specific intercept, \eqn{\theta_{0j}}, accounts for differences between species prevalence, while the row effect, \eqn{\alpha_i}, is included to account for differences in site total abundance. The inclusion of the latter ensures the ordination is in terms of species composition. If \eqn{\alpha_i} is omitted from the model i.e., \code{site.eff = FALSE}, the ordination will be in terms of species abundance.

\bold{Correlated column models:} When one or more latent variables are included in conjunction with covariates i.e., \code{X} is non-empty and \code{num.lv} > 1, \code{boral} fits separate GLMs to each column of \eqn{y} while allowing for any potential correlation between columns via the latent variables. This is useful for multivariate abundance data in ecology, where a separate GLM is fitted to species (modeling its response against environmental covariates, say), while accounting for the fact species at a site are correlated. The model fitted takes the following form,

\deqn{g(\mu_{ij}) = \theta_{0j} + \bm{x}'_i\bm{\beta}_j, + \bm{z}'_i\bm{\theta}_j,}

which we can see is a mash of the first two types of models. The linear predictor \eqn{\bm{z}'_i\bm{\theta}_j} induces a ``residual covariance/correlation"" between the columns of \eqn{y}, which is of rank \code{num.lv}. For multivariate abundance data, this leads to a parsimonious method of accounting for correlation between species not due to the shared environmental responses. After fitting the model, the residual correlation matrix can be obtained using the \code{\link{get.residual.cor}} function. Note \code{num.lv > 1} is necessarily in order to flexibly model the residual correlations (see also Pollock et al., 2014, for residual correlation matrices in the context of Joint Species Distribution Models).

\bold{Using information criteria at your own risk:} Using information criteria from \code{calc.ics = TRUE} for model selection should be done with extreme caution, for two reasons: 1) The implementation of some of these criteria is heuristic and experimental. 2) Deciding what model to fit should also be driven by the science. For example, it may be the case that criteria suggests a model with 3 or 4 latent variables is more appropriate. However, if we are interested in visualizing the data for ordination purposes, then models with 1 or 2 latent variables are more appropriate. As another example, whether or not we include row effects when ordinating multivariate abundance data depends on if we are interested in differences between sites in terms of species abundance (\code{site.eff = FALSE}) or species composition (\code{site.eff = TRUE}).  

Estimation for \code{boral} is performed using Bayesian Markov Chain Monte Carlo (MCMC) methods via JAGS (Plummer, 2003). Uninformative priors are used for all parameters i.e., normal priors with mean zero and variance given by \code{hypparams[1]} for all coefficients, row and column effects, and cutoffs for ordinal responses. For the negative binomial, normal, and tweedie distributions, uniform priors from 0 to \eqn{1/hypparams[2]} are used on \eqn{\phi}.  
}

\value{
An object of class "boral" is returned, being a list containing the following components where applicable:
     \item{call}{The matched call.}
     \item{lv.coefs.mean/median/sd/iqr}{Matrices containing the mean/median/standard deviation/interquartile range of the posterior distributions of the column-specific coefficients relating to the latent variables. This also includes the column-specific intercepts and dispersion parameters.}
     \item{lv.mean/median/sd/iqr}{A matrix containing the mean/median/standard deviation/interquartile range of the posterior distributions of the latent variables.}
     \item{X.coefs.mean/median/sd/iqr}{Matrices containing the mean/median/standard deviation/interquartile range of the posterior distributions of the column-specific coefficients relating to the model matrix \code{X}.}
     \item{cutoffs.mean/median/sd/iqr}{Vectors containing the mean/median/standard deviation/interquartile range of the posterior distributions of the common cutoffs for ordinal responses.}
     \item{powerparam.mean/median/sd/iqr}{Scalars for the mean/median/standard deviation/interquartile range of the posterior distributions of the common power parameter for tweedie responses.}
     \item{site.coefs.mean/median/sd/iqr}{Matrices containing the mean/median/standard deviation/interquartile range of the posterior distributions of the row-specific effects.}
     \item{hpdintervals}{A list containing components which correspond to the lower and upper bounds of highest posterior density (HPD) intervals for all the parameters indicated above. Please see \code{\link{get.hpdintervals}} for more details.}
     \item{ics}{If \code{calc.ics = TRUE}, then a vector of different information criteria values for the model. Please see \code{\link{get.measures}} for details on the criteria.}
     \item{jags.model}{If \code{save.model = TRUE}, the raw jags model fitted is returned. This can be quite large! }
     \item{n, p, family, trial.size, num.lv, ...}{Various attributes of the model fitted, including the dimension of \code{y}, the response and model matrix used, distributional assumptions and trial sizes, number of latent variables, the number of covariates, whether information criteria values were calculated, hyper-parameters used in the Bayesian estimation, and the number of levels for ordinal responses.}
}

\references{
\itemize{
\item Agresti, A. (2010). Analysis of Ordinal Categorical Data. Wiley.
\item Foster, S. D. & Bravington, M. V. (2013). A Poisson-Gamma model for analysis of ecological non-negative continuous data. Journal of Environmental and Ecological Statistics, 20, 533-552.
\item Hill, M. O. (1974). Correspondence analysis: a neglected multivariate method. Applied statistics, 23, 340-354.
\item Kruskal, J. B. (1964). Nonmetric multidimensional scaling: a numerical method. Psychometrika, 29, 115-129.
\item Legendre, P. & Gallagher, E. D. (2001). Ecologically meaningful transformations for ordination of species data. Oecologia, 129, 271-280.
Numerical ecology, Volume 20. Elsevier.
\item Legendre, P. & Legendre, L. (2012). Numerical ecology, Volume 20. Elsevier.
\item  Moles et al. (2013). Correlations between physical and chemical defences in plants: Trade-offs,
syndromes, or just many different ways to skin a herbivorous cat? New Phytologist 198, 252-263.
\item Plummer, M. (2003). JAGS: A program for analysis of Bayesian graphical models using Gibbs sampling. In Proceedings of the 3rd International Workshop on Distributed Statistical Computing. March (pp. 20-22).
\item Pollock, L. J., Tingley, R., Morris, W. K., Golding, N., O'Hara, R. B., Parris, K. M., Vesk, P. A., and McCarthy, M. A. (2014). Understanding co-occurrence by modelling species simultaneously with a Joint Species Distribution Model (JSDM). Methods in Ecology and Evolution, 5, 397-406.
\item Skrondal, A., & Rabe-Hesketh, S. (2004). Generalized latent variable modeling: Multilevel, longitudinal, and structural equation models. CRC Press.
\item Warton, D. I., Wright, S. T., & Wang, Y. (2012). Distance-based multivariate analyses confound location and dispersion effects. Methods in Ecology and Evolution, 3, 89-101.
\item Yi Wang, Ulrike Naumann, Stephen Wright and David Warton (2013). mvabund: statistical methods for analysing
  multivariate abundance data. R package version 3.8.4.}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\section{Warnings}{
\itemize{
\item \emph{No} intercept column is required in \code{X}. Column-specific intercepts are estimated automatically and given by the first column of \code{lv.coefs}.
\item If num.lv > 5, a warning is printed asking whether you really want to fit an boral with more than five latent variables.
\item For models including both explanatory covariates and latent variables, one requires \code{num.lv > 1} to allow flexible modeling of the residual correlation matrix. 
\item MCMC can take a long time to run, especially with if the response matrix \code{y} is large! The calculation of information criteria (\code{calc.ics = TRUE}) can also take a while.
\item MCMC with lots of ordinal columns take an especially long time to run! Moreover, estimates for the cutoffs in proportional odds regression may be poor for levels with little data. 
\item As discussed in the details, the use of information criteria should be done so with caution. What model to select should be first and foremost driven by the question of interest.
\item If \code{save.model = TRUE}, the raw jags model is also returned. This can be quite very memory-consuming, since it indirectly saves all the MCMC samples.
}
}

\seealso{
\code{\link{lvsplot}} for a scatter plot of the latent variables when \code{num.lv = 1} or \code{2}, \code{\link{summary.boral}} for a summary of the fitted boral model, \code{\link{get.measures}} and \code{\link{get.more.measures}} for information criteria from the fitted boral model, \code{\link{get.residual.cor}} for calculating the residual correlation matrix.
}

\examples{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
n <- nrow(y); p <- ncol(y); 

## Example 1 - model with two latent variables, site effects, 
## 	and no environmental covariates
spider.fit.nb <- boral(y, family = "negative.binomial", num.lv = 2, 
     site.eff = TRUE, n.burnin = 10, n.iteration = 100, 
     n.thin = 1, calc.ics = FALSE)

summary(spider.fit.nb)

plot(spider.fit.nb) ## Plots used in residual analysis, 
## Used to check if assumptions such an mean-variance relationship 
## are adequately satisfied.

lvsplot(spider.fit.nb) ## Plots a 2-D scatterplot of the latent variables, 
## which can be interpreted in the same manner as an ordination plot.

## Example 2 - model where the first six columns are assumed to be Poisson, 
## and the second six columns are assumed to be negative binomial...why not?
spider.fit.mix <- boral(y, family = rep(c("poisson","negative.binomial"),each=p/2), 
     num.lv = 2, site.eff = TRUE, n.burnin = 10, n.iteration = 100, 
     n.thin = 1, calc.ics = FALSE)

summary(spider.fit.mix)

## Example 3 - model with no latent variables, no site effects, 
## 	and environmental covariates
spider.fit.nb <- boral(y, X = spider$x, family = "negative.binomial", num.lv = 0, 
     site.eff = FALSE, n.burnin = 10, n.iteration = 100, n.thin = 1,
     calc.ics = FALSE)

summary(spider.fit.nb) 
## The results can be compared with the default example from 
## the manyglm() function in mvabund. Hopefully they are similar =D

}

\name{boral}
\alias{boral}
\alias{boral.default}
\alias{print.boral}
\title{Fitting boral (Bayesian Ordination and Regression AnaLysis) models}
\description{Bayesian ordination and regression models for analyzing multivariate data in ecology. Three ``types" of models may be fitted: 1) With covariates and no latent variables, \code{boral} fits independent column GLMs such that the columns of \eqn{y} are assumed to be independent; 2) With no covariates, \code{boral} fits a purely latent variable model to perform model-based unconstrained ordination (Hui et al., 2014); 3) With covariates and latent variables, \code{boral} fits correlated column GLMs, with latent variables accounting for correlations between the columns of \eqn{y}. Model selection is performed using either information criterion, or via Stochastic Search Variable Selection (SSVS, George and McCulloch, 1993) }

\usage{
boral(y, ...)

\method{boral}{default}(y, X = NULL, family, trial.size = NULL, num.lv = 0, 
	row.eff = FALSE, n.burnin = 5000, n.iteration = 35000, n.thin = 5, 
	save.model = FALSE, seed = 123, calc.ics = TRUE, 
	hypparams = c(100,100,100,100), ssvs.index = -1, ...)

\method{print}{boral}(x, ...)
}

\arguments{
  \item{y}{A \eqn{n} times \eqn{p} response matrix of multivariate data e.g., counts, binomial or Bernoulli responses, continuous response, and so on. In ecology for instance, \eqn{n} is the number of sites and \eqn{p} is the number of species. Any categorical (multinomial) responses \bold{must} be converted to integer values. For ordinal data, the minimum level of \eqn{y} must be 1 instead of 0.}
  
  \item{X}{An \eqn{n} times \eqn{d} model matrix of covariates, which can be included as part of the boral model. Defaults to \code{NULL}, in which case no model matrix was used. No intercept column should be included in \code{X}.}  

  \item{x}{An object for class "boral".}

\item{family}{Either a single element, or a vector of length equal to the number of columns in \eqn{y}. The former assumes all columns of \eqn{y} come from this distribution. The latter option allows for different distributions for each column of \eqn{y}. Elements can be one of "binomial" (with logit link), "poisson" (with log link), "negative.binomial" (with log link), "normal" (with identity link), "lnormal" for lognormal (with log link), "tweedie" (with log link), "exponential" (with log link), "gamma" (with log link), "beta" (with logit link), "ordinal" which performs cumulative logit regression. For the negative binomial distribution, the variance is parameterized as \eqn{Var(y) = \mu + \phi \mu^2}, where \eqn{\phi} is the column-specific overdispersion parameter. For the normal distribution, the variance is parameterized as \eqn{Var(y) = \phi}, where \eqn{\phi} is the column-specific variance. For the tweedie distribution, the variance is parameterized as \eqn{Var(y) = \phi \mu^p} where \eqn{\phi} is the column-specific dispersion parameter and \eqn{p} is a power parameter common to all columns assumed to be tweedie, with \eqn{1 < p < 2}. For the gamma distribution, the variance is parameterized as \eqn{Var(y) = \mu/\phi} where \eqn{\phi} is the column-specific rate (henceforth referred to also as dispersion parameter). For the beta distribution, the parameterization is in terms of the mean \eqn{\mu} and sample size \eqn{\phi} (henceforth referred to also as dispersion parameter), so that the two shape parameters are given by \eqn{a = \mu\phi} and \eqn{b = (1-\mu)\phi}.

All columns assumed to have ordinal responses are constrained to have the same cutoffs points, with a column-specific intercept to account for differences between the columns (please see \emph{Details} for formulation). 
}

 \item{trial.size}{Either equal to NULL, a single element, or a vector of length equal to the number of columns in \eqn{y}. If a single element, then all columns assumed to be binomially distributed will have trial size set to this. If a vector, different trial sizes are allowed in each column of \eqn{y}.}

 \item{num.lv}{Number of latent variables to fit. Can take any non-negative integer value. Defaults to 0.} 

 \item{row.eff}{A logical value indicating whether to include row effects or not. Defaults to FALSE.}    

 \item{n.burnin}{Length of burnin i.e., the number of iterations to discard at the beginning of the MCMC sampler.}

 \item{n.iteration}{Number of iterations including burnin.}

 \item{n.thin}{Thinning rate. Must be a positive integer. With the default values of \code{n.burnin} \cr, \code{n.iteration} and \code{n.thin}, this leads to a total of 6000 MCMC samples.}

 \item{save.model}{A logical value indicating whether to save the JAGS model file as a text file (called "jagsboralmodel.txt") in the current working directory, as well as the MCMC samples from the call to JAGS. If saved, various functions available in the \code{coda} package can be applied to the MCMC samples. Note MCMC samples can take up a lot of memory. Defaults to \code{FALSE}.}

 \item{seed}{seed for JAGS sampler. Within \code{boral}, a \code{set.seed(seed)} command is run immediately before starting the MCMC sampler. Defaults to the value 123.}

\item{calc.ics}{A logical values indicating whether to return various information criteria values, which could be used to perform model selection (see \code{\link{get.measures}}). Defaults to \code{TRUE}.}

 \item{hypparams}{Vector of 4 hyperparameters used in the set up of prior distributions. The first element is the variance for the normal priors of all (column-specific) intercepts, row effects, cutoff points for ordinal data. The second element is the variance for the normal priors of all column-specific coefficients relating to the latent variables (ignored if \code{num.lv = 0}). The third element is the variance for the normal priors of all column-specific coefficients relating to the model matrix \code{X} (ignored if \code{X = NULL}). The third element is the maximum of the uniform prior used for the dispersion parameters i.e, for the negative binomial, normal, lognormal, tweedie, gamma, beta (the dispersion parameter here refers to the sample size) distributions. Defaults to 100 for all three elements. Note the common power parameter in the tweedie distribution is assumed to have uniform prior from 1 to 2.}

 \item{ssvs.index}{Indices to be used for Stochastic Search Variable Selection (SSVS, George and McCulloch, 1993). Either a single element or a vector with length equal to the number of columns in \code{X}. Each element can take values of -1 (no SSVS is performed on this covariate), 0 (SSVS is performed on individual coefficients for this covariate, or any integer exceeding 1 (SSVS is performed on collectively all coefficients on this covariate/s.) Defaults to -1, in which case no model selection is performed on the fitted model at all. Please see the \code{details} for more information regarding the implementation of SSVS.}
 
 \item{...}{Not used.}
}

\details{
The boral package is designed to fit three types models which may be useful in ecology (and probably outside of ecology as well). 

\bold{Independent column models:} \code{boral} allows explanatory variables to be entered into the model via a model matrix \code{X}. This model matrix can contain anything the user wants, provided factors have been parameterized as dummy variables. It should \emph{not} include an intercept column. 

Without latent variables i.e., \code{num.lv = 0}, \code{boral} fits separate GLMs to each column of \eqn{y}, where the columns are assumed to be independent. That is,

\deqn{g(\mu_{ij}) = \alpha_i + \theta_{0j} + \bm{x}'_i\bm{\beta}_j; \quad i = 1,\ldots,n; j = 1,\ldots,p,}

where the mean response for element \eqn{(i,j)}, denoted as \eqn{mu_{ij}}, is regressed against the covariates \eqn{\bm{x}_i} via a link function \eqn{g(\cdot)}. \eqn{\theta_{0j}} and \eqn{\bm{\beta}_j} are the column-specific intercepts and coefficients respectively, while \eqn{\alpha_i} is an optional row effect (if \code{row.eff = TRUE}). Fitting the above type of model is thus more or less a Bayesian analog of the \code{manyglm} function in the \code{mvabund} package (Wang et al., 2013). Unlike \code{manyglm} though, different distributional assumptions are allowed however for each column of \eqn{y}, and row effects can be added as a type of "row-standardization". 

\bold{A not-so-brief tangent on distributions:} In the event different responses are collect for different columns, e.g., some columns of \eqn{y} are counts, while other columns are presence-absence, one can specify different distributions for each column. Aspects such as variable selection, residual analysis, and plotting of the latent variables are not affected by having different distributions. Naturally though, one has to be more careful with interpretation of the row effects \eqn{\alpha_i} and latent variables \eqn{\bm{z}_i}, with different scalings due to different link functions used. A situation where different distributions may prove useful is when \eqn{y} is a species-traits matrix, where each row is a species and each column is a trait such as specific leaf area, presence of thorns and so on. In this case, traits could be of different response types, and the goal perhaps is to perform unconstrained ordination to look for patterns between species on an underlying trait surface e.g., a defense index for a species (Moles et al., 2013; see the discussion below on how to perform model-based unconstrained ordination).

For multivariate abundance data in ecology (also known as community composition data, Legendre and Gallagher, 2001), species counts are often overdispersed. Using a negative binomial distribution (\code{family = "negative.binomial"}) to model the counts usually helps to account for this overdispersion. For non-negative continuous data such as biomass, the lognormal and tweedie distributions may be used (Foster and Bravington, 2013). Note however that a common power parameter is used for tweedie columns -- there almost always insufficient information to model column-specific power parameters. Normal responses are also implemented, just in case you encounter some normal stuff in ecology (pun intended)! 

The beta distribution can be used to model data between values between but not including 0 and 1. In principle, this would make it quite useful for percent cover data in ecology, if it not were for the fact that percent cover data is commonly characterized by having lots of zeros (which are not permitted for beta regression). An \emph{ad-hoc} fix to this would be to add a very small value to pull the data away from exact zeros and/ones. This is however quite heuristic, and pulls the model towards producing conservative results (see Smithson and Verkuilen, 2006, for a detailed discussion on beta regression, and Korhonen et al., 2007, for an example of an application to forest canopy cover data). Note the parameterization of the beta distribution used here is directly in terms of the mean \eqn{\mu} and sample size (referred to here as the dispersion parameter) \eqn{\phi}. In terms of the two shape parameters, this is equivalent to \eqn{shape1 = a = \mu\phi} and \eqn{shape2 = b = (1-\mu)\phi}.

For ordinal response columns, cumulative logit regression is used (also known as proportional odds regression, Agresti, 2010). \code{boral} assumes all ordinal columns are measured using the same scale i.e., all columns have the same number of theoretical levels. The number of levels is then assumed to be given by the maximum value from all the ordinal columns of \eqn{y}. Thus, the parameterization of proportional odds regression is such that all ordinal columns have the \emph{same} cutoff points, \eqn{\bm{\tau}}, while the column-specific intercept effect, \eqn{\theta_{0j}}, allows for deviations away from these common cutoffs. That is, 

\deqn{logit(P(y_{ij} \le k)) = \tau_k + \theta_{0j} + \ldots,}

where \eqn{logit(\cdot)} is the logit function, \eqn{P(y_{ij} \le k)} is the cumulative probability of element \eqn{y_{ij}} being less than or equal to level \eqn{k}, \eqn{\tau_k} is the cutoff linking levels \eqn{k} and \eqn{k+1} (and which are increasing in \eqn{k}), \eqn{\theta_{0j}} are the column effects, and "\eqn{\ldots}" generically represents what else is included in the model e.g., latent variables and related coefficients. A sum-to-zero constraint is imposed on the \eqn{\theta_{0j}}'s of all ordinal columns, to ensure model identifiability. 

The parameterization above is useful for modeling multivariate abundance data in ecology. When ordinal responses are recorded, usually the same ordinal scale is applied to all species e.g., level 1 = not there, level 2 = a bit there, level 3 = lots there, level 4 = everywhere! The quantity \eqn{\tau_k} can be interpreted as this common ordinal scale, while \eqn{\theta_{0j}} allows for deviations away from these to account for differences in species prevalence. Admittedly, the current implementation of \code{boral} for ordinal data can be quite slow. 

\bold{Purely latent variable models:} If no explantory variables are included and \code{num.lv} > 0, \code{boral} fits a purely latent variable model to perform model-based unconstrained ordination (Hui et al., 2014),

\deqn{g(\mu_{ij}) = \alpha_i + \theta_{0j} + \bm{z}'_i\bm{\theta}_j,}

where instead of measured covariates, we now have a vector of latent variables \eqn{\bm{z}_i} with \eqn{\bm{\theta}_j} being the column-specific coefficients relating to these latent variables. Unconstrained ordination is used for visualizing multivariate data in a low-dimensional space, without reference to explanatory variables (see Chapter 9, Legendre and Legendre, 2012). Typically, \code{num.lv} = 1 or 2 latent variables is used, allowing the latent variables to plotted (using \code{\link{lvsplot}}, for instance). The resulting plot can be interpreted in the same manner as plots from Nonmetric Multi-dimensional Scaling (NMDS, Kruskal, 1964) and Correspondence Analysis (CA, Hill, 1974), for example. For instance, with multivariate abundance data, unconstrained ordination is used to visualize the relationships between sites in terms of species abundance or composition. The column-specific intercept, \eqn{\theta_{0j}}, accounts for differences between species prevalence, while the row effect, \eqn{\alpha_i}, is included to account for differences in site total abundance. The inclusion of the latter ensures the ordination is in terms of species composition. If \eqn{\alpha_i} is omitted from the model i.e., \code{row.eff = FALSE}, the ordination will be in terms of relative species abundance.

\bold{Correlated column models:} When one or more latent variables are included in conjunction with covariates i.e., \code{X} is non-empty and \code{num.lv} > 1, \code{boral} fits separate GLMs to each column of \eqn{y} while allowing for any potential correlation between columns via the latent variables. This is useful for multivariate abundance data in ecology, where a separate GLM is fitted to species (modeling its response against environmental covariates, say), while accounting for the fact species at a site are correlated. The model fitted takes the following form,

\deqn{g(\mu_{ij}) = \theta_{0j} + \bm{x}'_i\bm{\beta}_j, + \bm{z}'_i\bm{\theta}_j,}

which we can see is a mash of the first two types of models. The linear predictor \eqn{\bm{z}'_i\bm{\theta}_j} induces a ``residual covariance/correlation"" between the columns of \eqn{y}, which is of rank \code{num.lv}. For multivariate abundance data, this leads to a parsimonious method of accounting for correlation between species not due to the shared environmental responses. After fitting the model, the residual correlation matrix can be obtained using the \code{\link{get.residual.cor}} function. Note \code{num.lv > 1} is necessarily in order to flexibly model the residual correlations (see also Pollock et al., 2014, for residual correlation matrices in the context of Joint Species Distribution Models).

\bold{Estimation:} For \code{boral} models, estimation is performed using Bayesian Markov Chain Monte Carlo (MCMC) methods via JAGS (Plummer, 2003). By default, uninformative priors are used for all parameters. That is, normal priors with mean zero and variance given by \code{hypparams[1]} are assigned to all intercepts, coefficients relating to latent variables, cutoffs for ordinal responses. Normal priors with mean zero and variance given by \code{hypparams[2]} are assigned to all coefficients relating to covariates in \code{X}. For the negative binomial, normal, lognormal, and tweedie distributions, uniform priors from 0 to \eqn{hypparams[3]} are used on the dispersion parameters.  

\bold{Using information criteria at your own risk:} Using information criterion from \code{calc.ics = TRUE} for model selection should be done with extreme caution, for two reasons: 1) The implementation of some of these criteria is heuristic and experimental. 2) Deciding what model to fit should also be driven by the science. For example, it may be the case that a criterion suggests a model with 3 or 4 latent variables is more appropriate. However, if we are interested in visualizing the data for ordination purposes, then models with 1 or 2 latent variables are more appropriate. As another example, whether or not we include row effects when ordinating multivariate abundance data depends on if we are interested in differences between sites in terms of relative species abundance (\code{site.eff = FALSE}) or species composition (\code{site.eff = TRUE}).  

\bold{SSVS:} As an alternative to using information criterion for model selection, Stochastic Search Variable Selection (SSVS, George and McCulloch, 1993) is also implemented for the column-specific coefficients relating to the model matrix \code{X}, \eqn{\bm{\beta}_j}. Basically, SSVS works by placing a spike-and-slab priors on these coefficients, such that the spike is a narrow normal distribution concentrated around zero and the spike is a normal distribution with a large variance.

\deqn{\rho(\beta) = I_{\beta = 1}\times\mathcal{N}(0,\sigma^2) + (1-I_{\beta = 1})\times \mathcal{N}(0,0.0001*\sigma^2),}

where \eqn{\sigma^2} is determined by \code{hypparams[2]} (see section on \bold{Estimation} above) and \eqn{I_{\beta = 1} = P(\beta = 1)} is an indicator function representing whether coefficient is included in the model, and is given a Bernoulli prior with probability of inclusion 0.5. After fitting, the posterior probability of \eqn{\beta} being included in the model is returned based on posterior mean of the indicator function \eqn{I_{\beta = 1}}. Note this is NOT the same as a \emph{p}-value seen in maximum likelihood estimation -- a \emph{p}-value provides an indication of how much evidence there is against the null hypothesis of \eqn{\beta = 0}, while the posterior probability provides a (direct) measure of how likely it is for \eqn{\beta \ne 0} given the data.

In \code{boral}, SSVS can be applied at a grouped or individual coefficient level, and this is governed by \code{ssvs.index}. For elements of \code{ssvs.index} equal to -1, SSVS is not applied on the corresponding covariate of the model matrix \code{X}. For elements equal to 0, SSVS is applied to each individual coefficient of the corresponding covariate in \code{X}. That is, the fitted model will return \eqn{p} posterior probabilities for this covariate, one for each column of \eqn{\bm{y}}. For elements taking positive integers {1,2,...}, SSVS is applied to each group of coefficients of the corresponding covariate in \code{X}. That is, the fitted model will return a single posterior probability for this covariate, indicating whether this covariate should be included for all columns of \eqn{\bm{y}}. Please see O'Hara and Sillanpaa (2009) for an discussion of Bayesian variable selection methods.

Note the last application of SSVS allows multiple covariates to be tested \emph{simultaneously}. For example, suppose \code{X} consists of five columns -- the first two columns are environmental covariates, while the last three correspond to quadratic terms of the two covariates as well as their interaction. If we want to ``test" whether any quadratic terms are required, then we can set \code{ssvs.index = c(-1,-1,1,1,1)}, so a single posterior probability of inclusion is returned for the last three columns of \code{X}. 

Finally, note the use of information criterion (and possibly residual analysis) in the presence of variable selection being performed using SSVS is questionable, and so the it maybe advised to separate out their use e.g., choose the explanatory variables first using SSVS, and then use information criterion to select the number of latent variables?
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
     \item{ssvs.indcoefs.mean/ssvs.indcoefs.sd}{Matrices containing the SSVS posterior probabilities and associated standard deviation of including individual coefficients in the model.}
     \item{ssvs.gpcoefs.mean/ssvs.gpcoefs.sd}{Matrices containing the SSVS posterior probabilities and associated standard deviation of including grouped coefficients in the model.}
     \item{hpdintervals}{A list containing components which correspond to the lower and upper bounds of highest posterior density (HPD) intervals for all the parameters indicated above. Please see \code{\link{get.hpdintervals}} for more details.}
     \item{ics}{If \code{calc.ics = TRUE}, then a list of different information criteria values for the model calculated using \code{\link{get.measures}} is run. Please see help file for \code{\link{get.measures}} regarding details on the criteria. Also, please note the ics returned are based on \code{\link{get.measures}} with \code{more.measures = FALSE}.}
     \item{jags.model}{If \code{save.model = TRUE}, the raw jags model fitted is returned. This can be quite large! }
     \item{n, p, family, trial.size, num.lv, ...}{Various attributes of the model fitted, including the dimension of \code{y}, the response and model matrix used, distributional assumptions and trial sizes, number of latent variables, the number of covariates, whether information criteria values were calculated, hyperparameters used in the Bayesian estimation, indices for SSVS, and the number of levels for ordinal responses.}
}

\references{
\itemize{
\item Agresti, A. (2010). Analysis of Ordinal Categorical Data. Wiley.
\item Foster, S. D. and Bravington, M. V. (2013). A Poisson-Gamma model for analysis of ecological non-negative continuous data. Journal of Environmental and Ecological Statistics, 20, 533-552.
\item George, E. I. and McCulloch, R. E. (1993). Variable selection via Gibbs sampling. Journal of the American Statistical Association, 85, 398-409.
\item Hill, M. O. (1974). Correspondence analysis: a neglected multivariate method. Applied statistics, 23, 340-354.
\item Korhonen, L., et al. (2007). Local models for forest canopy cover with beta regression. Silva Fennica, 41, 671-685.
\item Kruskal, J. B. (1964). Nonmetric multidimensional scaling: a numerical method. Psychometrika, 29, 115-129.
\item Legendre, P. and Gallagher, E. D. (2001). Ecologically meaningful transformations for ordination of species data. Oecologia, 129, 271-280.
Numerical ecology, Volume 20. Elsevier.
\item Legendre, P. and Legendre, L. (2012). Numerical ecology, Volume 20. Elsevier.
\item  Moles et al. (2013). Correlations between physical and chemical defences in plants: Trade-offs,
syndromes, or just many different ways to skin a herbivorous cat? New Phytologist, 198, 252-263.
\item O'Hara, B. and Sillianpaa, M.J. (2009). A Review of Bayesian Variable Selection Methods: What, How and Which. Bayesian Analysis, 4, 85-118.
\item Plummer, M. (2003). JAGS: A program for analysis of Bayesian graphical models using Gibbs sampling. In Proceedings of the 3rd International Workshop on Distributed Statistical Computing. March (pp. 20-22).
\item Pollock, L. J., et al. (2014). Understanding co-occurrence by modelling species simultaneously with a Joint Species Distribution Model (JSDM). Methods in Ecology and Evolution, 5, 397-406.
\item Skrondal, A., and Rabe-Hesketh, S. (2004). Generalized latent variable modeling: Multilevel, longitudinal, and structural equation models. CRC Press.
\item Smithson, M., and Verkuilen, J. (2006). A better lemon squeezer? Maximum-likelihood regression with beta-distributed dependent variables. Psychological methods, 11, 54-71.
\item Warton et al. (2012). Distance-based multivariate analyses confound location and dispersion effects. Methods in Ecology and Evolution, 3, 89-101.
\item Wang et al. (2013). mvabund: statistical methods for analysing multivariate abundance data. R package version 3.8.4.}
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
\item As discussed in the details, the use of information criterion should be done so with caution. What model to select should be first and foremost driven by the question of interest. Also, the use of information criterion in the presence of model seelction using SSVS is questionable.
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
     row.eff = TRUE, n.burnin = 10, n.iteration = 100, 
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
     num.lv = 2, row.eff = TRUE, n.burnin = 10, n.iteration = 100, 
     n.thin = 1, calc.ics = FALSE)

summary(spider.fit.mix)

## Example 3 - model with no latent variables, no site effects, 
## 	and environmental covariates
spider.fit.nb <- boral(y, X = spider$x, family = "negative.binomial", num.lv = 0, 
     row.eff = FALSE, n.burnin = 10, n.iteration = 100, n.thin = 1,
     calc.ics = FALSE)

summary(spider.fit.nb) 
## The results can be compared with the default example from 
## the manyglm() function in mvabund. Hopefully they are similar =D

## Example 4 - extend example 3 with grouped covariate selection
spider.fit.nb2 <- boral(y, X = spider$x, family = "negative.binomial", num.lv = 0, 
     row.eff = FALSE, n.burnin = 10, n.iteration = 100, n.thin = 1,
     calc.ics = FALSE, ssvs.index = c(-1,-1,-1,1,2,3))

summary(spider.fit.nb2) 

## Example 5 - extend example 3 with individual and grouped covariate selection
spider.fit.nb2 <- boral(y, X = spider$x, family = "negative.binomial", num.lv = 0, 
     row.eff = FALSE, n.burnin = 10, n.iteration = 100, n.thin = 1,
     calc.ics = FALSE, ssvs.index = c(-1,-1,-1,1,0,0))
}

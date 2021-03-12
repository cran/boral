\name{about.distributions}

\alias{about.distributions}

\docType{package}

\title{
Distributions available in boral
}

\description{
\Sexpr[results=rd, stage=render]{lifecycle::badge("stable")}

This help file provides more information regarding the distributions i.e., the \code{family} argument, available in the boral package, to handle various responses types.
}


\details{
A variety of families are available in boral, designed to accommodate multivariate abundance data of varying response types. Please see the \code{family} argument in the \code{\link{boral}} which lists all distributions that are currently available. 

For multivariate abundance data in ecology, species counts are often overdispersed. Using a negative binomial distribution (\code{family = "negative.binomial"}) to model the counts usually helps to account for this overdispersion. Please note the variance for the negative binomial distribution is parameterized as \eqn{Var(y) = \mu + \phi\mu^2}, where \eqn{\phi} is the dispersion parameter. 

For non-negative continuous data such as biomass, the lognormal, Gamma, and tweedie distributions may be used (Foster and Bravington, 2013). For the gamma distribution, the variance is parameterized as \eqn{Var(y) = \mu/\phi} where \eqn{\phi} is the response-specific rate (henceforth referred to also as dispersion parameter). 

For the tweedie distribution, a common power parameter is across all columns with this family, because there is almost always insufficient information to model response-specific power parameters. Specifically, the variance is parameterized as \eqn{Var(y) = \phi \mu^p} where \eqn{\phi} is the response-specific dispersion parameter and \eqn{p} is a power parameter common to all columns assumed to be tweedie, with \eqn{1 < p < 2}. 

Normal responses are also implemented, just in case you encounter normal stuff in ecology (pun intended)! For the normal distribution, the variance is parameterized as \eqn{Var(y) = \phi^2}, where \eqn{\phi} is the response-specific standard deviation. 

The beta distribution can be used to model data between values between but \emph{not} including 0 and 1. In principle, this would make it useful for percent cover data in ecology, if it not were for the fact that percent cover is commonly characterized by having lots of zeros (which are not permitted for beta regression). An \emph{ad-hoc} fix to this would be to add a very small value to shift the data away from exact zeros and/or ones. This is however heuristic, and pulls the model towards producing conservative results (see Smithson and Verkuilen, 2006, for a detailed discussion on beta regression, and Korhonen et al., 2007, for an example of an application to forest canopy cover data). Note the parameterization of the beta distribution used here is directly in terms of the mean \eqn{\mu} and the dispersion parameter \eqn{\phi} (more commonly know as the "sample size"). In terms of the two shape parameters, if we denote the two shape parameters as the vector \eqn{(a,b)}, his is equivalent to \eqn{a = \mu\phi} and \eqn{b = (1-\mu)\phi}.

For ordinal response columns, cumulative probit regression is used (Agresti, 2010). boral assumes all ordinal columns are measured using the same scale i.e., all columns have the same number of theoretical levels, even though some levels for some species may not be observed. The number of levels is then assumed to be given by the maximum value from all the ordinal columns of the response matrix. Because of this, all ordinal columns then assumed to have the \emph{same} cutoffs, \eqn{\bm{\tau}}, while the response-specific intercept, \eqn{\beta_{0j}}, allows for deviations away from these common cutoffs. That is, 

\deqn{\Phi(P(y_{ij} \le k)) = \tau_k + \beta_{0j} + \ldots,}

where \eqn{\Phi(\cdot)} is the probit function, \eqn{P(y_{ij} \le k)} is the cumulative probability of element \eqn{y_{ij}} being less than or equal to level \eqn{k}, \eqn{\tau_k} is the cutoff linking levels \eqn{k} and \eqn{k+1} (and which are increasing in \eqn{k}), \eqn{\beta_{0j}} are the column effects, and \eqn{\ldots} denotes what else is included in the model, e.g. latent variables and related coefficients. To ensure model identifiability, and also because they are interpreted as response-specific deviations from the common cutoffs, the \eqn{\beta_{0j}}'s are treated as random effects and drawn from a normal distribution with mean zero and unknown standard deviation.

The parameterization above is useful for modeling ordinal in ecology. When ordinal responses are recorded, usually the same scale is applied to all species e.g., level 1 = not there, level 2 = a bit there, level 3 = lots there, level 4 = everywhere! The quantity \eqn{\tau_k} can thus be interpreted as this common scale, while \eqn{\beta_{0j}} allows for deviations away from these to account for differences in species prevalence. Admittedly, the current implementation of boral for ordinal data can be quite slow. 

For count distributions where zeros are not permitted, then the zero truncated Poisson (\code{family = "ztpoisson"}) and zero truncated negative binomial distributions (\code{family = "ztnegative.binomial"}) are possible. Note for these two distributions, and as is commonly implemented in other regression models e.g., the \code{countreg} package (Zeileis and Kleiber, 2018), the models are set up such that a log-link connects the mean of the \emph{untruncated} distribution to the linear predictor. While not necessarily useful on its own, the zero truncated distributions may often be used in conjunction with an model for modeling presence-absence data, and together they can be used to construct the hurdle class of models (noting direct implementation of hurdle models is currently not available). 

Finally, in the event different responses are collected for different columns, e.g., some columns of the response matrix are counts, while other columns are presence-absence, one can specify different distributions for each column. Aspects such as variable selection, residual analysis, and plotting of the latent variables are, in principle, not affected by having different distributions. Naturally though, one has to be more careful with interpretation of the row effects \eqn{\alpha_i} and latent variables \eqn{\bm{u}_i}, as different link functions will be applied to each column of the response matrix. A situation where different distributions may prove useful is when the response matrix is a species--traits matrix, where each row is a species and each column a trait such as specific leaf area. In this case, traits could be of different response types, and the goal perhaps is to perform unconstrained ordination to look for patterns between species on an underlying trait surface e.g., a defense index for a species (Moles et al., 2013).

}


\references{
\itemize{
\item Agresti, A. (2010). Analysis of Ordinal Categorical Data. Wiley.

\item Foster, S. D. and Bravington, M. V. (2013). A Poisson-Gamma model for analysis of ecological non-negative continuous data. Journal of Environmental and Ecological Statistics, 20, 533-552.

\item Korhonen, L., et al. (2007). Local models for forest canopy cover with beta regression. Silva Fennica, 41, 671-685.

\item  Moles et al. (2013). Correlations between physical and chemical defences in plants: Trade-offs, syndromes, or just many different ways to skin a herbivorous cat? New Phytologist, 198, 252-263.

\item Smithson, M., and Verkuilen, J. (2006). A better lemon squeezer? Maximum-likelihood regression with beta-distributed dependent variables. Psychological methods, 11, 54-71.

\item Zeileis, A., and Kleiber C. (2018). countreg: Count Data Regression. R package version 0.2-1. URL http://R-Forge.R-project.org/projects/countreg/
}
}

\author{
\packageAuthor{boral}

Maintainer: \packageMaintainer{boral}
}

\section{Warnings}{
\itemize{
\item MCMC with lots of ordinal columns take an especially long time to run! Moreover, estimates for the cutoffs in cumulative probit regression may be poor for levels with little data. Major apologies for this advance =(
}
}

\seealso{
\code{\link{boral}} for the main boral fitting function. 
}

\examples{
## Please see main boral function for examples. 
}

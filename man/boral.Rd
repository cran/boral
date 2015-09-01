\name{boral}
\alias{boral}
\alias{boral.default}
\alias{print.boral}
\title{Fitting boral (Bayesian Ordination and Regression AnaLysis) models}
\description{Bayesian ordination and regression models for analyzing multivariate data in ecology. Three "types" of models may be fitted: 1) With covariates and no latent variables, boral fits independent response GLMs; 2) With no covariates, boral fits a purely latent variable model; 3) With covariates and latent variables, boral fits correlated response GLMs.}

\usage{
boral(y, ...)

\method{boral}{default}(y, X = NULL, traits = NULL, which.traits = NULL, family, 
	trial.size = 1, num.lv = 0, row.eff = "none", n.burnin = 10000, 
	n.iteration = 40000, n.thin = 30, save.model = FALSE, seed = 123, 
	calc.ics = TRUE, hypparams = c(100,20,100,50), ssvs.index = -1, 
	do.fit = TRUE, model.name = NULL, ...)

\method{print}{boral}(x, ...)
}

\arguments{
  \item{y}{A response matrix of multivariate data e.g., counts, binomial or Bernoulli responses, continuous response, and so on. With multivariate abundance data ecology for instance, rows correspond to sites and columns correspond to species. Any categorical (multinomial) responses \bold{must} be converted to integer values. For ordinal data, the minimum level of \code{y} must be 1 instead of 0.}
  
  \item{X}{A model matrix of covariates, which can be included as part of the boral model. Defaults to \code{NULL}, in which case no model matrix was used. No intercept column should be included in \code{X}.}  

  \item{x}{An object for class "boral".}

  \item{traits}{A model matrix of species covariates, which can be included as part of the boral model. Defaults to \code{NULL}, in which case no matrix was used. An intercept column should be included in \code{traits} if appropriate (usually is).}  

  \item{which.traits}{A list of length equal to (number of columns in \code{X} + 1), informing which columns of \code{traits} the column-specific intercepts and each of the column-specific regression coefficients should be regressed against. The first element in the list applies to the column-specific intercept, while the remaining elements apply to the regression coefficients. Each element of \code{which.traits} is a vector indicating which traits are to be used. For example, if \code{which.traits[[2]] = c(2,3)}, then the regression coefficients corresponding to the first column in \code{X} are regressed against the second and third columns of \code{traits}. If \code{which.traits[[2]] = 0}, then the regression coefficients are treated as independent. Please see help file below for more details.
  
  Defaults to \code{NULL}, in conjunction with \code{traits = NULL}).}

  \item{family}{Either a single element, or a vector of length equal to the number of columns in \code{y}. The former assumes all columns of \code{y} come from this distribution. The latter option allows for different distributions for each column of \code{y}. Elements can be one of "binomial" (with probit link), "poisson" (with log link), "negative.binomial" (with log link), "normal" (with identity link), "lnormal" for lognormal (with log link), "tweedie" (with log link), "exponential" (with log link), "gamma" (with log link), "beta" (with logit link), "ordinal" (cumulative probit regression). 

  For the negative binomial distribution, the variance is parameterized as \eqn{Var(y) = \mu + \phi\mu^2}, where \eqn{\phi} is the column-specific dispersion parameter. For the normal distribution, the variance is parameterized as \eqn{Var(y) = \phi^2}, where \eqn{\phi} is the column-specific standard deviation. For the tweedie distribution, the variance is parameterized as \eqn{Var(y) = \phi \mu^p} where \eqn{\phi} is the column-specific dispersion parameter and \eqn{p} is a power parameter common to all columns assumed to be tweedie, with \eqn{1 < p < 2}. For the gamma distribution, the variance is parameterized as \eqn{Var(y) = \mu/\phi} where \eqn{\phi} is the column-specific rate (henceforth referred to also as dispersion parameter). For the beta distribution, the parameterization is in terms of the mean \eqn{\mu} and sample size \eqn{\phi} (henceforth referred to also as dispersion parameter), so that the two shape parameters are given by \eqn{a = \mu\phi} and \eqn{b = (1-\mu)\phi}.

 All columns assumed to have ordinal responses are constrained to have the same cutoffs points, with a column-specific intercept to account for differences between the columns (please see \emph{Details} for formulation). 
}

 \item{trial.size}{Either equal to a single element, or a vector of length equal to the number of columns in y. If a single element, then all columns assumed to be binomially distributed will have trial size set to this. If a vector, different trial sizes are allowed in each column of y. The argument is ignored for all columns not assumed to be binomially distributed. Defaults to 1, i.e. Bernoulli distribution.}

 \item{num.lv}{Number of latent variables to fit. Can take any non-negative integer value. Defaults to 0.} 

 \item{row.eff}{Single element indicating whether row effects are included as fixed effects ("fixed"), random effects ("random") or not included ("none") in the boral model. If random effects, they are drawn from a normal distribution with mean zero and unknown variance, analogous to a random intercept in mixed models. Defaults to "none".}    

 \item{n.burnin}{Length of burnin i.e., the number of iterations to discard at the beginning of the MCMC sampler.}

 \item{n.iteration}{Number of iterations including burnin.}

 \item{n.thin}{Thinning rate. Must be a positive integer. With the default values of \code{n.burnin} \cr, \code{n.iteration} and \code{n.thin}, this leads to a final of 1000 MCMC samples.}

 \item{save.model}{A logical value indicating whether to save the JAGS model file as a text file (with name based on \code{model.name}) in the current working directory, as well as the MCMC samples from the call to JAGS. If saved, various functions available in the \code{coda} package can be applied to the MCMC samples. Note MCMC samples can take up a lot of memory. Defaults to FALSE.}

 \item{seed}{Seed for JAGS sampler. A \code{set.seed(seed)} command is run immediately before starting the MCMC sampler. Defaults to the value 123.}

\item{calc.ics}{A logical values indicating whether to return various information criteria values, which could be used to perform model selection (see \code{\link{get.measures}}). Defaults to TRUE.}

 \item{hypparams}{Vector of four hyperparameters used in the set up of prior distributions. The first element is the variance for the normal priors of all column-specific intercepts, row effects, and cutoff points for ordinal data. It also controls the maximum of the uniform prior for the standard deviation of the random effects normal distribution, if \code{row.eff = "random"}. The second element is the variance for the normal priors of all latent variable coefficients (ignored if \code{num.lv = 0}). The third element is the variance for the normal priors of all column-specific coefficients relating to the model matrix \code{X} (ignored if \code{X = NULL}). When traits are included in the model, it also controls the maximum of the uniform prior for the standard deviation of the normally distributed random effects (please see section on \emph{Including species traits} below). The fourth element controls the maximum of the uniform prior used for dispersion parameters, \eqn{\phi}. Note the common power parameter in the tweedie distribution is assumed to have uniform prior from 1 to 2.}

 \item{ssvs.index}{Indices to be used for Stochastic Search Variable Selection (SSVS, George and McCulloch, 1993). Either a single element or a vector with length equal to the number of columns in \code{X}. Each element can take values of -1 (no SSVS is performed on this covariate), 0 (SSVS is performed on individual coefficients for this covariate), or any integer exceeding 1 (SSVS is performed on collectively all coefficients on this covariate/s.) Defaults to -1, in which case no model selection is performed on the fitted model at all. Please see the \code{details} for more information regarding the implementation of SSVS.}
 
 \item{do.fit}{A logical value indicating whether to actually fit the boral model. If set to FALSE, then only the JAGS model file is written to the current working directly (as text file with name based on \code{model.name}), no MCMC sampling is performed, and \emph{nothing else} is returned. Defaults to TRUE.}    

 \item{model.name}{Name of the text file that the JAGS model is written to. Defaults to NULL, in which case the default of "jagsboralmodel.txt" is used.}

 \item{...}{Not used.}
}


\details{
The boral package is designed to fit three types models which may be useful in ecology (and probably outside of ecology as well =D). 

\bold{Independent response models:} boral allows explanatory variables to be entered into the model via a model matrix \code{X}. This model matrix can contain anything the user wants, provided factors have been parameterized as dummy variables. It should NOT include an intercept column. 

Without latent variables, i.e. \code{num.lv = 0}, boral fits separate GLMs to each column of the \eqn{n \times p} matrix \code{y}, where the columns are assumed to be independent. 

\deqn{g(\mu_{ij}) = \alpha_i + \beta_{0j} + \bm{x}^T_i\bm{\beta}_j; \quad i = 1,\ldots,n; j = 1,\ldots,p,}

where the mean response for element (i,j), denoted as \eqn{mu_{ij}}, is regressed against the covariates \eqn{\bm{x}_i} via a link function \eqn{g(\cdot)}. The quantities \eqn{beta_{0j}} and \eqn{\bm{beta}_j} denote the column-specific intercepts and coefficients respectively, while \code{alpha_i} is an optional row effect that may be treated as a fixed or random effect. The latter assumes the row effects are drawn from a normal distribution with unknown variance \eqn{\phi^2}.

Fitting the above type of model is sort of like a Bayesian analog of the \code{manyglm} function in the \code{mvabund} package (Wang et al., 2013). Unlike \code{manyglm} though, row effects can be added easily as a type of "row-standardization". Also, a wider range of assumed distributions (families) are possible, as discussed below.

\bold{A not-so-brief tangent on distributions:} In the event different responses are collect for different columns, e.g., some columns of \code{y} are counts, while other columns are presence-absence, one can specify different distributions for each column. Aspects such as variable selection, residual analysis, and plotting of the latent variables are, in principle, not affected by having different distributions. Naturally though, one has to be more careful with interpretation of the row effects \eqn{\alpha_i} and latent variables \eqn{\bm{z}_i}, as different link functions will be applied to each column of \code{y}. A situation where different distributions may prove useful is when \code{y} is a species-traits matrix, where each row is a species and each column a trait such as specific leaf area. In this case, traits could be of different response types, and the goal perhaps is to perform unconstrained ordination to look for patterns between species on an underlying trait surface e.g., a defense index for a species (Moles et al., 2013; see also the discussion below on how to perform model-based unconstrained ordination).

For multivariate abundance data in ecology (also known as community composition data, Legendre and Gallagher, 2001), species counts are often overdispersed. Using a negative binomial distribution (\code{family = "negative.binomial"}) to model the counts usually helps to account for this overdispersion. Please note the variance for the negative binomial distribution is parameterized as \eqn{Var(y) = \mu + \phi\mu^2}, where \eqn{\phi} is the dispersion parameter. 

For non-negative continuous data such as biomass, the lognormal and tweedie distributions may be used (Foster and Bravington, 2013). Note however that a common power parameter is used for tweedie columns -- there is almost always insufficient information to model column-specific power parameters. Normal responses are also implemented, just in case you encounter normal stuff in ecology (pun intended)! 

The beta distribution can be used to model data between values between but \emph{not} including 0 and 1. In principle, this would make it useful for percent cover data in ecology, if it not were for the fact that percent cover is commonly characterized by having lots of zeros (which are not permitted for beta regression). An \emph{ad-hoc} fix to this would be to add a very small value to shift the data away from exact zeros and/or ones. This is however heuristic, and pulls the model towards producing conservative results (see Smithson and Verkuilen, 2006, for a detailed discussion on beta regression, and Korhonen et al., 2007, for an example of an application to forest canopy cover data). Note the parameterization of the beta distribution used here is directly in terms of the mean \eqn{\mu} and the dispersion parameter \eqn{\phi} (more commonly know as the "sample size"). In terms of the two shape parameters, this is equivalent to \eqn{shape1 = a = \mu\phi} and \eqn{shape2 = b = (1-\mu)\phi}.

For ordinal response columns, cumulative probit regression is used (Agresti, 2010). boral assumes all ordinal columns are measured using the same scale i.e., all columns have the same number of theoretical levels. The number of levels is then assumed to be given by the maximum value from all the ordinal columns of \code{y}. Because of this, all ordinal columns then assumed to have the \emph{same} cutoffs, \eqn{\bm{\tau}}, while the column-specific intercept effect, \eqn{\beta_{0j}}, allows for deviations away from these common cutoffs. That is, 

\deqn{probit(P(y_{ij} \le k)) = \tau_k + \beta_{0j} + \ldots,}

where \eqn{probit(\cdot)} is the probit function, \eqn{P(y_{ij} \le k)} is the cumulative probability of element \eqn{y_{ij}} being less than or equal to level \eqn{k}, \eqn{\tau_k} is the cutoff linking levels \eqn{k} and \eqn{k+1} (and which are increasing in \eqn{k}), \eqn{\beta_{0j}} are the column effects, and \eqn{\ldots} denotes what else is included in the model, e.g. latent variables and related coefficients. A sum-to-zero constraint is imposed on the \eqn{\beta_{0j}}'s of all ordinal columns, to ensure model identifiability. 

The parameterization above is useful for modeling ordinal in ecology. When ordinal responses are recorded, usually the same scale is applied to all species e.g., level 1 = not there, level 2 = a bit there, level 3 = lots there, level 4 = everywhere! The quantity \eqn{\tau_k} can thus be interpreted as this common scale, while \eqn{\beta_{0j}} allows for deviations away from these to account for differences in species prevalence. Admittedly, the current implementation of boral for ordinal data can be quite slow. 


\bold{Purely latent variable models:} If no explantory variables are included and \code{num.lv} > 0, boral fits a purely latent variable model to perform model-based unconstrained ordination (Hui et al., 2014),

\deqn{g(\mu_{ij}) = \alpha_i + \beta_{0j} + \bm{z}^T_i\bm{\theta}_j,}

where instead of measured covariates, we now have a vector of latent variables \eqn{\bm{z}_i} with \eqn{\bm{\theta}_j} being the column-specific coefficients relating to these latent variables. The column-specific intercept, beta_{0j}, accounts for differences between species prevalence, while the row effect, \eqn{alpha_i}, is included to account for differences in site total abundance (typically assuming a fixed effect, \code{row.eff = "fixed"}, although see Jamil and ter Braak, 2013, for a motivation for using random site effects), so that the ordination is then in terms of species composition. If \eqn{\alpha_i} is omitted from the model i.e., \code{row.eff = FALSE}, then the ordination will be in terms of relative species abundance.

Unconstrained ordination is used for visualizing multivariate data in a low-dimensional space, without reference to covariates (Chapter 9, Legendre and Legendre, 2012). Typically, \code{num.lv} = 1 to 3 latent variables is used, allowing the latent variables to plotted (using \code{\link{lvsplot}}, for instance). The resulting plot can be interpreted in the same manner as plots from Nonmetric Multi-dimensional Scaling (NMDS, Kruskal, 1964) and Correspondence Analysis (CA, Hill, 1974), for example. A biplot can also be constructed by setting \code{biplot = TRUE} when using \code{\link{lvsplot}}, so that both the latent variables and their corresponding coefficients are plotted. For instance, with multivariate abundance data, biplots are used to visualize the relationships between sites in terms of species abundance or composition, as well as the indicator species for the sites. 


\bold{Correlated response models:} When one or more latent variables are included in conjunction with covariates i.e., \code{X} is given and \code{num.lv} > 1, boral fits separate GLMs to each column of \code{y} while allowing for residual correlation between columns via the latent variables. This is quite useful for multivariate abundance data in ecology, where a separate GLM is fitted to species (modeling its response against environmental covariates), while accounting for the fact species at a site are likely to be correlated for reason other than similiarites in environmental responses, e.g. biotic interaction, phylogeny, and so on. Correlated response model take the following form,

\deqn{g(\mu_{ij}) = \beta_{0j} + \bm{x}^T_i\bm{\beta}_j, + \bm{z}^T_i\bm{\theta}_j.}

This model is thus a mash of the first two types of models. The linear predictor \eqn{\bm{z}^T_i\bm{\theta}_j} induces a residual covariance between the columns of \code{y} (which is of rank \code{num.lv}). For multivariate abundance data, this leads to a parsimonious method of accounting for correlation between species not due to the shared environmental responses. After fitting the model, the residual correlation matrix then can be obtained via the \code{\link{get.residual.cor}} function. Note \code{num.lv > 1} is necessarily in order to flexibly model the residual correlations; see Pollock et al. (2014) for residual correlation matrices in the context of Joint Species Distribution Models, and Warton et al. (2015) for an overview of latent variable models in multivariate ecology.


\bold{Including species traits:} When covariates \code{X} are included (i.e. both the independent and correlated response models), one has the option of also including traits to help explain differences in species environmental responses to these covariates. Specifically, when \code{traits} and \code{which.traits} are supplied, then the \eqn{\beta_{0j}}'s and \eqn{\bm{\beta}_j}'s are then regarded as random effects drawn from a normal distribution. For the species-specific intercepts, we have

\deqn{\beta_{0j} \sim N(\kappa_{01} + \bm{traits}^T_j\bm{\kappa}_1, \sigma^2_1),}

where \eqn{(\kappa_{01},\bm{\kappa}_1)} are the regression coefficients relating to the traits to the intercepts and \eqn{\sigma_1} is the error standard deviation. These are now the "parameters" in the model, in the sense that priors are assigned to them and MCMC sampling is used to estimate them (see the next section on estimation). Please note that in order of \eqn{\kappa_{01}} to be included in the model, an intercept column MUST be included in \code{traits}.

In an analogous manner, each of the elements in \eqn{\bm{\beta}_j = (\beta_{j1},\ldots,\beta_{jd})} are now drawn as random effects from a normal distribution. That is, for \eqn{k = 1,\ldots,d} where \code{d = ncol(X)}, we have,

\deqn{\beta_{jk} \sim N(\kappa_{0k} + \bm{traits}^T_j\bm{\kappa}_k, \sigma^2_k),}

Which traits are to included (regressed) in the mean of the normal distributions is determined by the list \code{which.traits}. The first element in the list applies to \eqn{beta_{0j}}, while the remaining elements apply to the the \eqn{\bm{\beta}_j}. Each element of \code{which.traits} is a vector indicating which traits are to be used. For example, if \code{which.traits[[2]] = c(2,3)}, then the \eqn{\beta_{j1}}'s are drawn from a normal distribution with mean depending only on the second and third columns of \code{traits}. If \code{which.traits[[2]] = 0}, then the regression coefficients are treated as independent, i.e. the values of \eqn{\beta_{j1}} are given their own priors and estimated separately from each other. 

Including species traits in the model can be regarded as a method of simplifying the model -- rather than each to estimates \eqn{p} sets of species-specific coefficients, we instead say that these coefficients are linearly related to the corresponing values of their traits (Warton et al., 2015). That is, we are using trait data to help explain similarities/differences in species responses to the environment. This idea has close relations to the fourth corner problem in ecology (Brown et al., 2014). Unlike the models of Brown et al. (2014) however, which treat the \eqn{\beta_{0j}}'s and \eqn{\beta_{jk}}'s are fixed effects and fully explained by the traits, boral adopts a random effects approach similar to Jamil et al. (2013) to "soak up" any additional between species differences in environmental responses not explained by traits.


\bold{Estimation:} For boral models, estimation is performed using Bayesian Markov Chain Monte Carlo (MCMC) methods via JAGS (Plummer, 2003). Please note that only \emph{one} MCMC chain in run -- this point is discussed further in this help file. Regarding prior distributions, the default settings are as follows: 
\itemize{
\item Normal priors with mean zero and variance given by \code{hypparams[1]} are assigned to all intercepts, cutoffs for ordinal responses, and row effects. If the row effects are assumed to random, then the standard deviation of the normal random effect is assigned a uniform prior with maximum \code{hypparams[1]},
\item Normal priors with mean zero and variance given by \code{hypparams[2]} are assigned coefficients relating to latent variables, \eqn{\bm{\theta}_j}, 
\item Normal priors with mean zero and variance given by \code{hypparams[3]} are assigned to all coefficients relating to covariates in \eqn{\bm{\beta}_j}. If traits are included, the same normal priors are assigned to the \eqn{\kappa}'s, and the standard deviation \eqn{\sigma_k} are assigned uniform priors with maximum equal to \code{hypparams[4]}.
\item For the negative binomial, normal, lognormal, and tweedie distributions, uniform priors with maximum equal to \code{hypparams[4]} are used on the dispersion parameters. Please note that for the normal and lognormal distributions, these uniform priors are assigned to the standard deviations \eqn{\phi} (see Gelman, 2006).
}

With the default values of \code{hypparams}, all parameters are given uninformative prior distributions except for the priors of the latent variable coefficients \eqn{\bm{\theta}_j}. We recommend such a ``semi-informative" prior for the latent variable coefficients, as this tends to be produce more stable MCMC sampling particularly if the response matrix is large and sparse. 


\bold{Using information criteria at your own risk:} Using information criterion from \code{calc.ics = TRUE} for model selection should be done with extreme caution, for two reasons: 1) The implementation of some of these criteria is heuristic and experimental, 2) Deciding what model to fit should also be driven by the science. For example, it may be the case that a criterion suggests a model with 3 or 4 latent variables is more appropriate. However, if we are interested in visualizing the data for ordination purposes, then models with 1 or 2 latent variables are more appropriate. As another example, whether or not we include row effects when ordinating multivariate abundance data depends on if we are interested in differences between sites in terms of relative species abundance (\code{row.eff = "none"}) or species composition (\code{row.eff = "fixed"}).  

We also make the important point that if traits are included in the model, then the regression coefficients \eqn{\beta_{0j}, \bm{\beta}_j} are now random effects. However, currently the calculation of all information criteria do not take this into account! 


\bold{SSVS:} As an alternative to using information criterion for model selection, Stochastic Search Variable Selection (SSVS, George and McCulloch, 1993) is also implemented for the column-specific coefficients \eqn{\bm{\beta}_j}. Basically, SSVS works by placing a spike-and-slab priors on these coefficients, such that the spike is a narrow normal distribution concentrated around zero and the spike is a normal distribution with a large variance.

\deqn{\rho(\beta) = I_{\beta = 1}\times\mathcal{N}(0,\sigma^2) + (1-I_{\beta = 1})\times \mathcal{N}(0,0.0001*\sigma^2),}

where \eqn{\sigma^2} is determined by \code{hypparams[3]} (see section on estimation above) and \eqn{I_{\beta = 1} = P(\beta = 1)} is an indicator function representing whether coefficient is included in the model. It is given a Bernoulli prior with probability of inclusion 0.5. After fitting, the posterior probability of \eqn{\beta} being included in the model is returned based on posterior mean of the indicator function \eqn{I_{\beta = 1}}. Note this is NOT the same as a \emph{p}-value seen in maximum likelihood estimation -- a \emph{p}-value provides an indication of how much evidence there is against the null hypothesis of \eqn{\beta = 0}, while the posterior probability provides a measure of how likely it is for \eqn{\beta \ne 0} given the data.

In boral, SSVS can be applied at a grouped or individual coefficient level, and this is governed by \code{ssvs.index}. For elements of \code{ssvs.index} equal to -1, SSVS is not applied on the corresponding covariate of the model matrix \code{X}. For elements equal to 0, SSVS is applied to each individual coefficient of the corresponding covariate in \code{X}. That is, the fitted model will return \eqn{p} posterior probabilities for this covariate, one for each column of \code{y}. For elements taking positive integers {1,2,...}, SSVS is applied to each group of coefficients of the corresponding covariate in \code{X}. That is, the fitted model will return a single posterior probability for this covariate, indicating whether this covariate should be included for all columns of \code{y}; see O'Hara and Sillanpaa (2009) for an discussion of Bayesian variable selection methods.

Note the last application of SSVS allows multiple covariates to be tested \emph{simultaneously}. For example, suppose \code{X} consists of five columns -- the first two columns are environmental covariates, while the last three correspond to quadratic terms of the two covariates as well as their interaction. If we want to "test" whether any quadratic terms are required, then we can set \cr \code{ssvs.index = c(-1,-1,1,1,1)}, so a single posterior probability of inclusion is returned for the last three columns of \code{X}. 

Finally, note using information criterion (and possibly residual analysis) should probably not be done at the same as when SSVS is used, and it is advised to separate out their applications e.g., choose the explanatory variables first using SSVS, and then use information criterion to select the number of latent variables???
}


\section{Why is only one MCMC chain run?}{
Much like the \code{MCMCfactanal} function in the \code{MCMCpack} package (Martin et al., 2011) for conducting factor analysis, which is a special case of the purely latent variable model with Gaussian responses, boral deliberately runs only one MCMC chain. This runs contrary to the recommendation of most Bayesian analyses, where the advice is to run multiple MCMC chains and check convergence using (most commonly) the Gelman-Rubin statistic or ``Rhat'' (Gelman et al., 2013). The main reason for this is that, in the context of MCMC sampling, the latent variable model is invariant to a switch of the sign, i.e. \eqn{\bm{z}^T_i\bm{\theta}_j = (-\bm{z})^T_i(-\bm{\theta}_j)}, and so is actually unidentifiable. This is similar to well-known problem of label switching that occurs during the course of MCMC sampling for mixture models (see for instance, Section 4.9, McLachlan and Peel, 2004), and is due to the fact that the sign of the latent variables (ordination coordinates) is inherently arbitrary.

As a result of this sign-switching problem, it means that different MCMC chains can produce latent variables and corresponding coefficients values that, while having similar magnitudes, will be different in sign. Consequently, combining MCMC chains and checking Rhats, computing posterior means and medians etc...becomes inappropriate (in principle, one way to resolve this problem would be to post-process the MCMC chains and deal with sign switching, but this is really hard!). Therefore, to alleviate this issue together, boral chooses to only run one MCMC chain.

What does this mean for the user? 
\itemize{
\item For checking convergence, we recommend you look at trace plots of the MCMC chains. Using the \code{coda} package, which is automatically loaded when the \code{boral} package is loaded, try something like \code{traceplot(fit$jags.model, ask = T)}. You could also try \code{geweke.diag} for Geweke's convergence diagnostic, although no promises this necessarily does what is meant it!
\item If you have a lot of data, e.g. lots of sites compared to species, sign-switching tends to be less of problem and pops up less often.
\item IMPORTANTLY, if the goal of your analysis is to inference while account for residual correlations between the columns of \code{y}, and not for model-based ordination, then the sign-switching problem is not a problem at all! This is because while the signs of the latent variables and associated coefficients may switch, the correlation and their signs are unaffected. In other words, looking the point estimates and credible intervals of regression coefficients \eqn{\bm{\beta}_j}, and functions like \code{get.residual.cor} are unaffected by sign-switching. 
}

}

\value{
An object of class "boral" is returned, being a list containing the following components where applicable:
     \item{call}{The matched call.}
     
     \item{lv.coefs.mean/median/sd/iqr}{Matrices containing the mean/median/standard deviation/interquartile range of the posterior distributions of the latent variable coefficients. This also includes the column-specific intercepts, and dispersion parameters if appropriate.}
     
     \item{lv.mean/median/sd/iqr}{A matrix containing the mean/median/standard deviation/interquartile range of the posterior distributions of the latent variables.}
     
     \item{X.coefs.mean/median/sd/iqr}{Matrices containing the mean/median/standard deviation/interquartile range of the posterior distributions of the column-specific coefficients relating to the model matrix \code{X}.}
     
     \item{traits.coefs.mean/median/sd/iqr}{Matrices containing the mean/median/standard deviation/interquartile range of the posterior distributions of the coefficients and standard deviation relating to the species traits (please see the section on including traits above).}

     \item{cutoffs.mean/median/sd/iqr}{Vectors containing the mean/median/standard deviation/interquartile range of the posterior distributions of the common cutoffs for ordinal responses (please see the not-so-brief tangent on distributions above).}
     
     \item{powerparam.mean/median/sd/iqr}{Scalars for the mean/median/standard deviation/interquartile range of the posterior distributions of the common power parameter for tweedie responses (please see the not-so-brief tangent on distributions above).}
     
     \item{row.coefs.mean/median/sd/iq}{Vectors containing the mean/median/standard deviation/interquartile range of the posterior distributions of the row effects.}
     
	\item{row.sigma.mean/median/sd/iqr}{Scalars containing the mean/median/standard deviation/interquartile range of the posterior distributions of the standard deviation for the row random effects normal distribution.}
	
     \item{ssvs.indcoefs.mean/ssvs.indcoefs.sd}{Matrices containing the SSVS posterior probabilities and associated standard deviation of including individual coefficients in the model (please see the section on SSVS above). }
     
     \item{ssvs.gpcoefs.mean/ssvs.gpcoefs.sd}{Matrices containing the SSVS posterior probabilities and associated standard deviation of including grouped coefficients in the model (please see the section on SSVS above). }
     
     \item{hpdintervals}{A list containing components which correspond to the lower and upper bounds of highest posterior density (HPD) intervals for all the parameters indicated above. Please see \code{\link{get.hpdintervals}} for more details.}
     
     \item{ics}{If \code{calc.ics = TRUE}, then a list of different information criteria values for the model calculated using \code{\link{get.measures}} is run. Please see help file for \code{\link{get.measures}} regarding details on the criteria. Also, please note the ics returned are based on \code{\link{get.measures}} with \code{more.measures = FALSE}.}
     
     \item{jags.model}{If \code{save.model = TRUE}, the raw jags model fitted is returned. This can be quite large!}
     
     \item{n, p, family, trial.size, num.lv, ...}{Various attributes of the model fitted, including the dimension of \code{y}, the response and model matrix used, distributional assumptions and trial sizes, number of latent variables, the number of covariates and traits, whether information criteria values were calculated, hyperparameters used in the Bayesian estimation, indices for SSVS, the number of levels for ordinal responses, and n.burin, n.iteration and n.thin.}
}

\references{
\itemize{
\item Agresti, A. (2010). Analysis of Ordinal Categorical Data. Wiley.
\item Brown, et al. (2014). The fourth-corner solution - using predictive models to understand how species traits interact with the environment. Methods in Ecology and Evolution 5, 344-352.
\item Foster, S. D. and Bravington, M. V. (2013). A Poisson-Gamma model for analysis of ecological non-negative continuous data. Journal of Environmental and Ecological Statistics, 20, 533-552.
\item Gelman et al. (2013) Bayesian Data Analysis. CRC Press.
\item Gelman A. (2006) Prior distributions for variance parameters in hierarchical models. Bayesian Analysis 1, 515-533.
\item George, E. I. and McCulloch, R. E. (1993). Variable selection via Gibbs sampling. Journal of the American Statistical Association, 85, 398-409.
\item Hui et al. (2014). Model-based approaches to unconstrained ordination. Methods in Ecology and Evolution, 6, 399-411.
\item Hill, M. O. (1974). Correspondence analysis: a neglected multivariate method. Applied statistics, 23, 340-354.
\item Jamil, T., and ter Braak, C.J.F. (2013). Generalized linear mixed models can detect unimodal species-environment relationships. PeerJ 1: e95.
\item Jamil, T. et al. (2013). Selecting traits that explain species-environment relationships: a generalized linear mixed model approach. Journal of Vegetation Science 24, 988-1000
\item Korhonen, L., et al. (2007). Local models for forest canopy cover with beta regression. Silva Fennica, 41, 671-685.
\item Kruskal, J. B. (1964). Nonmetric multidimensional scaling: a numerical method. Psychometrika, 29, 115-129.
\item Legendre, P. and Gallagher, E. D. (2001). Ecologically meaningful transformations for ordination of species data. Oecologia, 129, 271-280.
Numerical ecology, Volume 20. Elsevier.
\item Legendre, P. and Legendre, L. (2012). Numerical ecology, Volume 20. Elsevier.
\item Martin et al. (2011). MCMCpack: Markov Chain Monte Carlo in R. Journal of Statistical Software, 42, 1-21. URL: http://www.jstatsoft.org/v42/i09/.
\item McLachlan, G., and Peel, D. (2004). Finite Mixture Models. Wiley.
\item  Moles et al. (2013). Correlations between physical and chemical defences in plants: Trade-offs,
syndromes, or just many different ways to skin a herbivorous cat? New Phytologist, 198, 252-263.
\item O'Hara, B., and Sillianpaa, M.J. (2009). A Review of Bayesian Variable Selection Methods: What, How and Which. Bayesian Analysis, 4, 85-118.
\item Plummer, M. (2003). JAGS: A program for analysis of Bayesian graphical models using Gibbs sampling. In Proceedings of the 3rd International Workshop on Distributed Statistical Computing. March (pp. 20-22).
\item Pollock, L. J. et al. (2014). Understanding co-occurrence by modelling species simultaneously with a Joint Species Distribution Model (JSDM). Methods in Ecology and Evolution, 5, 397-406.
\item Skrondal, A., and Rabe-Hesketh, S. (2004). Generalized latent variable modeling: Multilevel, longitudinal, and structural equation models. CRC Press.
\item Smithson, M., and Verkuilen, J. (2006). A better lemon squeezer? Maximum-likelihood regression with beta-distributed dependent variables. Psychological methods, 11, 54-71.
\item Warton et al. (2015). So Many Variables: Joint Modeling in Community Ecology. Trends in Ecology and Evolution, in review.
\item Warton et al. (2012). Distance-based multivariate analyses confound location and dispersion effects. Methods in Ecology and Evolution, 3, 89-101.
\item Wang et al. (2013). \code{mvabund}: statistical methods for analysing multivariate abundance data. R package version 3.8.4.}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\section{Warnings}{
\itemize{
\item \emph{No} intercept column is required in \code{X}. Column-specific intercepts are estimated automatically and given by the first column of \code{lv.coefs}.
\item If num.lv > 5, a warning is printed asking whether you really want to fit an boral with more than five latent variables. A warning is also printed if num.lv == 1, as this is not going to be successful in modeling between the correlation between columns.
\item For models including both explanatory covariates and latent variables, one requires \code{num.lv > 1} to allow flexible modeling of the residual correlation matrix. 
\item MCMC can take a long time to run, especially with if the response matrix \code{y} is large! The calculation of information criteria (\code{calc.ics = TRUE}) can also take a while. Apologies for this advance =(
\item MCMC with lots of ordinal columns take an especially long time to run! Moreover, estimates for the cutoffs in cumulative probit regression may be poor for levels with little data. Major apologies for this advance =(
\item As discussed in the details, the use of information criterion should be done so with caution. What model to select should be first and foremost driven by the question of interest. Also, the use of information criterion in the presence of model seelction using SSVS is questionable.
\item If \code{save.model = TRUE}, the raw jags model is also returned. This can be quite very memory-consuming, since it indirectly saves all the MCMC samples.
}
}

\seealso{
\code{\link{lvsplot}} for a scatter plot of the latent variables (and their coefficients if applicable) when \code{num.lv = 1} or \code{2}, \code{\link{summary.boral}} for a summary of the fitted boral model, \code{\link{get.measures}} and \cr \code{\link{get.more.measures}} for information criteria from the fitted boral model, \code{\link{get.residual.cor}} for calculating the residual correlation matrix. 
}

\examples{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
n <- nrow(y); p <- ncol(y); 

## Example 1 - model with two latent variables, site effects, 
## 	and no environmental covariates
spider.fit.nb <- boral(y, family = "negative.binomial", num.lv = 2, 
     row.eff = "fixed", n.burnin = 10, n.iteration = 100, 
     n.thin = 1, calc.ics = FALSE)

summary(spider.fit.nb)

plot(spider.fit.nb, ask = FALSE, mfrow = c(2,2)) ## Plots used in residual analysis, 
## Used to check if assumptions such an mean-variance relationship 
## are adequately satisfied.

lvsplot(spider.fit.nb) ## Biplot of the latent variables, 
## which can be interpreted in the same manner as an ordination plot.


## Example 2 - model with no latent variables, no site effects, 
## 	and environmental covariates
X <- scale(spider$x)
spider.fit.nb <- boral(y, X = X, family = "negative.binomial", 
	num.lv = 0, n.burnin = 10, n.iteration = 100, n.thin = 1)

summary(spider.fit.nb) 
## The results can be compared with the default example from 
## the manyglm() function in mvabund. Hopefully they are similar =D


## Example 3 - Extend example 2 to demonstrate grouped and individual
## covariate selection
spider.fit.nb2 <- boral(y, X = X, family = "negative.binomial", 
	num.lv = 0, n.burnin = 10, n.iteration = 100, n.thin = 1,
	calc.ics = FALSE, ssvs.index = c(-1,-1,-1,0,1,2))
     
summary(spider.fit.nb2) 


## Example 3 - model fitted to presence-absence data, no site effects, and
## two latent variables
data(tikus)
y <- tikus$abun
y[y > 0] <- 1
y <- y[1:20,] ## Consider only years 1981 and 1983
y <- y[,apply(y > 0,2,sum) > 2] ## Consider only spp with more than 2 presences
     
tikus.fit <- boral(y, family = "binomial", num.lv = 2, 
	n.burnin = 10, n.iteration = 100, n.thin = 1, calc.ics = FALSE)
     
lvsplot(tikus.fit, biplot = FALSE) 
## A strong location between the two sampling years 


## Example 4 - model fitted to count data, no site effects, and
## two latent variables, plus traits included to explain environmental responses
data(antTraits)
y <- antTraits$abun
X <- as.matrix(scale(antTraits$env))
## Include only traits 1, 2, and 5
traits <- as.matrix(cbind(1,antTraits$traits[,c(1,2,5)]))
which.traits <- vector("list",ncol(X)+1)
for(i in 1:length(which.traits)) which.traits[[i]] <- 1:ncol(traits)
## Just for fun, the regression coefficients for the second column of X
## will be estimated separately and not regressed against traits.
which.traits[[3]] <- 0

fit.traits <- boral(y, X = X, traits = traits, which.traits = which.traits, 
	family = "negative.binomial", num.lv = 2, n.burnin = 10, n.iteration = 100, 
	n.thin = 1, calc.ics = FALSE)

summary(fit.traits)
}
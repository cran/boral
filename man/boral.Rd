\name{boral}

\alias{boral}
\alias{boral.default}
\alias{print.boral}

\title{Fitting boral (Bayesian Ordination and Regression AnaLysis) models}

\description{Bayesian ordination and regression models for analyzing multivariate data in ecology. Three "types" of models may be fitted: 1) With covariates and no latent variables, boral fits independent response GLMs; 2) With no covariates, boral fits a pure latent variable model; 3) With covariates and latent variables, boral fits correlated response GLMs.}

\usage{
boral(y, ...)

\method{boral}{default}(y, X = NULL, traits = NULL, which.traits = NULL, 
	family, trial.size = 1, num.lv = 0, row.eff = "none", row.ids = NULL,
	offset = NULL, save.model = FALSE, calc.ics = FALSE, 
     mcmc.control = list(n.burnin = 10000, n.iteration = 40000, 
     n.thin = 30, seed = 123), 
	prior.control = list(type = c("normal","normal","normal","uniform"), 
	hypparams = c(10, 10, 10, 30), ssvs.index = -1, ssvs.g = 1e-6, 
	ssvs.traitsindex = -1), do.fit = TRUE, model.name = NULL, ...)

\method{print}{boral}(x, ...)
}

\arguments{
  \item{y}{A response matrix of multivariate data e.g., counts, binomial or Bernoulli responses, continuous response, and so on. With multivariate abundance data ecology for instance, rows correspond to sites and columns correspond to species. Any categorical (multinomial) responses \bold{must} be converted to integer values. For ordinal data, the minimum level of \code{y} must be 1 instead of 0.}
  
  \item{X}{A model matrix of covariates, which can be included as part of the boral model. Defaults to \code{NULL}, in which case no model matrix was used. No intercept column should be included in \code{X}.}  

  \item{x}{An object for class "boral".}

  \item{traits}{A model matrix of species covariates, which can be included as part of the boral model. Defaults to \code{NULL}, in which case no matrix was used. No intercept column should be included in \code{traits}, as it is included automatically.}  

  \item{which.traits}{A list of length equal to (number of columns in \code{X} + 1), informing which columns of \code{traits} the column-specific intercepts and each of the column-specific regression coefficients should be regressed against. The first element in the list applies to the column-specific intercept, while the remaining elements apply to the regression coefficients. Each element of \code{which.traits} is a vector indicating which traits are to be used.
  
  For example, if \code{which.traits[[2]] = c(2,3)}, then the regression coefficients corresponding to the first column in \code{X} are regressed against the second and third columns of \code{traits}. If \code{which.traits[[2]][1] = 0}, then the regression coefficients for each column are treated as independent. Please see \code{\link{about.traits}} for more details.
  
  Defaults to \code{NULL}, and used in conjunction with \code{traits} and \cr \code{prior.control$ssvs.traitsindex}.}

  \item{family}{Either a single element, or a vector of length equal to the number of columns in \code{y}. The former assumes all columns of \code{y} come from this distribution. The latter option allows for different distributions for each column of \code{y}. Elements can be one of "binomial" (with probit link), "poisson" (with log link), "negative.binomial" (with log link), "normal" (with identity link), "lnormal" for lognormal (with log link), "tweedie" (with log link), "exponential" (with log link), "gamma" (with log link), "beta" (with logit link), "ordinal" (cumulative probit regression). 

  Please see \code{\link{about.distributions}} for information on distributions available in boral overall.
  }

 \item{trial.size}{Either equal to a single element, or a vector of length equal to the number of columns in y. If a single element, then all columns assumed to be binomially distributed will have trial size set to this. If a vector, different trial sizes are allowed in each column of y. The argument is ignored for all columns not assumed to be binomially distributed. Defaults to 1, i.e. Bernoulli distribution.}

 \item{num.lv}{Number of latent variables to fit. Can take any non-negative integer value. Defaults to 0.} 

 \item{row.eff}{Single element indicating whether (multiple) row effects are included as fixed effects ("fixed"), random effects ("random") or not included ("none") in the boral model. If random effects, they are drawn from a normal distribution with mean zero and unknown variance, analogous to a random intercept in mixed models. Defaults to "none".}    

 \item{row.ids}{A matrix with the number of rows equal to the number of rows in \code{y}, and the number of columns equal to the number of row effects to be included in the model. Element \eqn{(i,j)} indicates to the cluster ID of row \eqn{i} in \code{y} for random effect eqn{j}. This matrix is useful if one wants to specify more complicated row effect structures beyond a single, row effect unique to each row; please see details below as well as examples below. Whether these row effects are included as fixed or random effects is governed by \code{row.eff}. Defaults to \code{NULL}, so that if \code{row.eff = "none"} then the argument is ignored, otherwise if \code{row.eff = "fixed"} or \code{"random"}, then \cr \code{row.ids = matrix(1:nrow(y), ncol = 1)} i.e., a single, row effect unique to each row.}    

 \item{offset}{A matrix with the same dimensions as the response matrix \code{y}, specifying an a-priori known component to be included in the linear predictor during fitting. Defaults to \code{NULL}.}

 \item{save.model}{If \code{save.model = TRUE}, then the JAGS model file is saved as a text file (with name given by \code{model.name}) in the current working directory as well as the MCMC samples, which themselves can be extracted using the \code{get.mcmcsamples} function. Various functions available in the \code{coda} package can be applied to the MCMC samples for diagnosing convergence. Note MCMC samples can take up a lot of memory. Defaults to \code{FALSE}.}

\item{calc.ics}{If \code{calc.ics = TRUE}, then various information criteria values are also returned, which could be used to perform model selection (see \code{\link{get.measures}}). Defaults to \code{FALSE}. WARNING: As of version 1.5, functions to calculate information criteria will no longer be updated...use at your own peril!!!}

\item{mcmc.control}{A list of parameters for controlling the MCMC sampling. Values not set will assume default values. These include:
	\itemize{
		\item \emph{n.burnin:} Length of burnin i.e., the number of iterations to discard at the beginning of the MCMC sampler. Defaults to 10000.

		\item \emph{n.iteration:} Number of iterations including burnin. Defaults to 40000.
	
		\item \emph{n.thin:} Thinning rate. Must be a positive integer. Defaults to 30.

		\item \emph{seed:} Seed for JAGS sampler. A \code{set.seed(seed)} command is run immediately before starting the MCMC sampler. Defaults to the value 123.
		}
	}
	
\item{prior.control}{A list of parameters for controlling the prior distributions. Values not set will assume default values. These include:
	\itemize{
		\item \emph{type:} Vector of four strings indicating the type of prior distributions to use. In order, these are: 1) priors for all column-specific intercepts, row effects, and cutoff points for ordinal data; 2) priors for the latent variable coefficients. This is ignored if \code{num.lv = 0}; 3) priors for all column-specific coefficients relating to \code{X} (ignored if \code{X = NULL}). When traits are included in the model, this is also the prior for the trait regression coefficients (please see \code{\link{about.traits}} for more information); 4) priors for any dispersion parameters and variance (standard deviation, to be precise) parameters in the model.
				
		For elements 1-3, the prior distributions currently available include: I) ``normal", which is a normal prior with the variance controlled by elements 1-3 in \code{hypparams}; II) ``cauchy", which is a Cauchy prior with variance controlled by elements 1-3 in \code{hypparams}. Gelman, et al. (2008) considers using Cauchy priors with variance \eqn{2.5^2} as weakly informative priors for coefficients in logistic and potentially other generalized linear models; III) ``uniform", which is a symmetric uniform prior with minimum and maximum values controlled by element 1-3 in \code{hypparams}. 
		
		For element 4, the prior distributions currently available include: I) ``uniform", which is uniform prior with minimum zero and maximum controlled by element 4 in \code{hypparmas}; II) ``halfnormal", which is half-normal prior with variance controlled by \code{hypparams}; III) ``halfcauchy", which is a half-Cauchy prior with variance controlled by element 4 in \code{hypparams}.
		
		Defaults to the vector \code{c("normal","normal","normal","uniform")}. 

		\item \emph{hypparams:} Vector of four hyperparameters used in the set up of prior distributions. In order, these are: 1) affects the prior distribution for all column-specific intercepts, row effects, and cutoff points for ordinal data; 2) affects the prior distribution for all latent variable coefficients. This is ignored if \code{num.lv = 0}; 3) affects the prior distribution for column-specific coefficients relating to \code{X} (ignored if \code{X = NULL}). When traits are included in the model, it also affects the prior distribution for the trait regression coefficients; 4) affects the prior distribution for any dispersion parameters, as well as the prior distributions for the standard deviation of the random effects normal distribution if \code{row.eff = "random"}, the standard deviation of the column-specific random intercepts for these columns if more than two of the columns are ordinal, and the standard deviation of the random effects normal distribution for trait regression coefficients when traits are included in the model.
		
		Defaults to the vector \code{c(10, 10, 10, 30)}. The use of normal distributions with mean zero and variance 10 as priors is seen as one type of (very) weakly informative prior, according to \href{https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations}{Prior choice recommendations}.
		
		\item \emph{ssvs.index:} Indices to be used for stochastic search variable selection (SSVS, George and McCulloch, 1993). Either a single element or a vector with length equal to the number of columns in \code{X}. Each element can take values of -1 (no SSVS is performed on this covariate), 0 (SSVS is performed on individual coefficients for this covariate), or any integer greater than 0 (SSVS is performed on collectively all coefficients on this covariate/s.) 
		
		Please see \code{\link{about.ssvs}} for more information regarding the implementation of SSVS. Defaults to -1, in which case SSVS is not performed on \code{X} variables. 

		\item \emph{ssvs.g:} Multiplicative, shrinkage factor for SSVS, which controls the strength of the "spike" in the SSVS mixture prior. In summary, if the coefficient is included in the model, the "slab" prior is a normal distribution with mean zero and variance given by element 3 in \code{hypparams}, while if the coefficient is not included in the model, the "spike" prior is normal distribution with mean zero and variance given by element 3 in \code{hypparams} multiplied by \code{ssvs.g}. Please see \code{\link{about.ssvs}} for more information regarding the implementation of SSVS. Defaults to 1e-6.  		
		
		\item \emph{ssvs.traitsindex:} Used in conjunction with \code{traits} and \code{which.traits}, this is a list of indices to be used 
		for performing SSVS on the trait coefficients. Should be a list with the same length as \code{which.traits}, and with each element a vector of indices with the same length as the corresponding element in \code{which.traits}. Each index either can take values of -1 (no SSVS on this trait coefficient) or 0 (no SSVS on this trait coefficient). 

		Please see \code{\link{about.ssvs}} for more information regarding the implementation of SSVS. Defaults to -1, in which case SSVS is not performed on any of the trait coefficients, if they are included in the model.
		}
	}
 
 \item{do.fit}{If \code{do.fit = FALSE}, then only the JAGS model file is written to the current working directly (as text file with name based on \code{model.name}). No MCMC sampling is performed, and \emph{nothing else} is returned. Defaults to \code{TRUE}.}    

 \item{model.name}{Name of the text file that the JAGS model is written to. Defaults to \code{NULL}, in which case the default of "jagsboralmodel.txt" is used.}

% \item{new.format}{If \code{new.format = TRUE}, then the new (and hopefully simpler) format is used for the output, as implemented in versions \eqn{>= 1.2}. If \code{FALSE}, then the old outpt format is used. Defaults to \code{TRUE}.}

 \item{...}{Not used.}
}


\details{
The boral package is designed to fit three types models which may be useful in ecology (and probably outside of ecology as well =D). 

\bold{Independent response models:} boral allows explanatory variables to be entered into the model via \code{X}. This model matrix can contain anything the user wants, provided factors have been parameterized as dummy variables. It should NOT include an intercept column. 

Without latent variables, i.e. \code{num.lv = 0}, boral fits separate GLMs to each column of the \eqn{n \times p} matrix \code{y}, where the columns are assumed to be independent. 

\deqn{g(\mu_{ij}) = \alpha_i + \beta_{0j} + \bm{x}\top_i\bm{\beta}_j; \quad i = 1,\ldots,n; j = 1,\ldots,p,}

where the mean response for element (i,j), denoted as \eqn{mu_{ij}}, is regressed against the covariates \eqn{\bm{x}_i} via a link function \eqn{g(\cdot)}. The quantities \eqn{beta_{0j}} and \eqn{\bm{beta}_j} denote the column-specific intercepts and coefficients respectively, while \code{alpha_i} is an optional row effect that may be treated as a fixed or random effect. The latter assumes the row effects are drawn from a normal distribution with unknown variance \eqn{\phi^2}. One reason we might want to include row effects is to account differences in sampling intensity between sites: these can lead to differences in site total abundance, and so by including fixed effects they play the same role as an offset to account for these differences. 

Note boral can also handle multiple, hierarchical row effects, which may be useful to account for sampling design. This is controlled using the \code{row.ids} argument. For example, if the first five rows of \eqn{y} correspond to replications from site 1, the next five rows correspond to replications from site 2, and so on, then one can set \code{row.ids = matrix(c(1,1,1,1,1,2,2,2,2,2,...), ncol = 1)} to take this in account. While this way of coding row effects via the \code{row.ids} argument takes some getting used to, it has been done this way partly to force the user to think more carefully about exactly the structure of the data i.e., with great power comes great responsibility...

If \code{offset} is supplied, then it will be included in the linear predictor below (and all linear predictors below where appropriate).

Without row effects, the above independent response model is basically a Bayesian analog of the \code{manyglm} function in the \code{mvabund} package (Wang et al., 2013). 

\bold{Pure latent variable models:} If no explanatory variables are included and \code{num.lv} > 0, boral fits a pure latent variable model to perform model-based unconstrained ordination (Hui et al., 2014),

\deqn{g(\mu_{ij}) = \alpha_i + \beta_{0j} + \bm{z}\top_i\bm{\theta}_j,}

where instead of measured covariates, we now have a vector of latent variables \eqn{\bm{z}_i} with \eqn{\bm{\theta}_j} being the column-specific coefficients relating to these latent variables. The column-specific intercept, beta_{0j}, accounts for differences between species prevalence, while the row effect, \eqn{alpha_i}, is included to account for differences in site total abundance (typically assuming a fixed effect, \code{row.eff = "fixed"}, although see Jamil and ter Braak, 2013, for a motivation for using random site effects), so that the ordination is then in terms of species composition. If \eqn{\alpha_i} is omitted from the model i.e., \code{row.eff = FALSE}, then the ordination will be in terms of relative species abundance. As mentioned previously, one reason for including fixed row effects is to account of any potential differences in sampling intensity between sites.  

Unconstrained ordination is used for visualizing multivariate data in a low-dimensional space, without reference to covariates (Chapter 9, Legendre and Legendre, 2012). Typically, \code{num.lv} = 1 to 3 latent variables is used, allowing the latent variables to plotted (using \code{\link{lvsplot}}, for instance). The resulting plot can be interpreted in the same manner as plots from Nonmetric Multi-dimensional Scaling (NMDS, Kruskal, 1964) and Correspondence Analysis (CA, Hill, 1974), for example. A biplot can also be constructed by setting \code{biplot = TRUE} when using \code{\link{lvsplot}}, so that both the latent variables and their corresponding coefficients are plotted. For instance, with multivariate abundance data, biplots are used to visualize the relationships between sites in terms of species abundance or composition, as well as the indicator species for the sites. 


\bold{Correlated response models:} When one or more latent variables are included in conjunction with covariates i.e., \code{X} is given and \code{num.lv} > 1, boral fits separate GLMs to each column of \code{y} while allowing for residual correlation between columns via the latent variables. This is quite useful for multivariate abundance data in ecology, where a separate GLM is fitted to species (modeling its response against environmental covariates), while accounting for the fact species at a site are likely to be correlated for reason other than similarities in environmental responses, e.g. biotic interaction, phylogeny, and so on. Correlated response model take the following form,

\deqn{g(\mu_{ij}) = \alpha_i + \beta_{0j} + \bm{x}^\top_i\bm{\beta}_j, + \bm{z}^\top_i\bm{\theta}_j.}

This model is thus a mash of the first two types of models. The linear predictor \eqn{\bm{z}^\top_i\bm{\theta}_j} induces a residual covariance between the columns of \code{y} (which is of rank \code{num.lv}). For multivariate abundance data, this leads to a parsimonious method of accounting for correlation between species not due to the shared environmental responses. After fitting the model, the residual correlation matrix then can be obtained via the \code{\link{get.residual.cor}} function. Note \code{num.lv > 1} is necessarily in order to flexibly model the residual correlations; see Pollock et al. (2014) for residual correlation matrices in the context of Joint Species Distribution Models, and Warton et al. (2015) for an overview of latent variable models in multivariate ecology.


\bold{Distributions:} A variety of families are available in boral, designed to handle multivariate abundance data of varying response types. Please see \code{\link{about.distributions}} for more information on this.


\bold{Including species traits:} When covariates \code{X} are included i.e. both the independent and correlated response models, one has the option of also including traits to help explain differences in species environmental responses to these covariates. Please see \code{\link{about.traits}} for more information on this.


\bold{Estimation:} For boral models, estimation is performed using Bayesian Markov Chain Monte Carlo (MCMC) methods via JAGS (Plummer, 2003). Please note that only \emph{one} MCMC chain in run: this point is discussed later in this help file. Regarding prior distributions, the default settings, based on the \code{prior.control} argument, are as follows: 

\itemize{
\item Normal priors with mean zero and variance given by element 1 in \code{hypparams} are assigned to all intercepts, cutoffs for ordinal responses, and row effects. 

\item Normal priors with mean zero and variance given by element 2 in \code{hypparams} are assigned coefficients relating to latent variables, \eqn{\bm{\theta}_j}.

\item Normal priors with mean zero and variance given by element 3 \code{hypparams} are assigned to all coefficients relating to covariates in \eqn{\bm{\beta}_j}. If traits are included, the same normal priors are assigned to the \eqn{\kappa}'s, and the standard deviation \eqn{\sigma_k} are assigned uniform priors with maximum equal to element 4 in \code{hypparams}.

\item For the negative binomial, normal, lognormal, and tweedie distributions, uniform priors with maximum equal to element 4 in \code{hypparams} are used on the dispersion parameters. Please note that for the normal and lognormal distributions, these uniform priors are assigned to the standard deviations \eqn{\phi} (see Gelman, 2006). If there are any variance (standard deviation, to be precise) parameters in the model, then these are also assigned uniform priors with maximum equal to element 4 in \code{hypparams} e.g., standard deviation of the normal random effect if the row effects are assumed to random, the standard deviation of the normal random column-specific intercepts if more than two columns are ordinal responses etc...
}


\bold{Using information criteria at your own risk:} Using information criterion from \code{calc.ics = TRUE} for model selection should be done with extreme caution, for two reasons: 1) The implementation of some of these criteria is heuristic and experimental, 2) Deciding what model to fit should also be driven by the science. For example, it may be the case that a criterion suggests a model with 3 or 4 latent variables is more appropriate. However, if we are interested in visualizing the data for ordination purposes, then models with 1 or 2 latent variables are more appropriate. As another example, whether or not we include row effects when ordinating multivariate abundance data depends on if we are interested in differences between sites in terms of relative species abundance (\code{row.eff = "none"}) or species composition (\code{row.eff = "fixed"}). We also make the important point that if traits are included in the model, then the regression coefficients \eqn{\beta_{0j}, \bm{\beta}_j} are now random effects. However, currently the calculation of all information criteria do not take this into account! 


\bold{SSVS:} Stochastic search variable selection (SSVS, George and McCulloch, 1993) is also implemented for the column-specific coefficients \eqn{\bm{\beta}_j}. Please see \code{\link{about.ssvs}} for more information on this approach.
}


\section{Why is only one MCMC chain run?}{
Much like the \code{MCMCfactanal} function in the \code{MCMCpack} package (Martin et al., 2011) for conducting factor analysis, which is a special case of the pure latent variable model with Gaussian responses, boral deliberately runs only one MCMC chain. This runs contrary to the recommendation of most Bayesian analyses, where the advice is to run multiple MCMC chains and check convergence using (most commonly) the Gelman-Rubin statistic (Gelman et al., 2013). The main reason for this is that, in the context of MCMC sampling, the latent variable model is invariant to a switch of the sign, i.e. \eqn{\bm{z}^\top_i\bm{\theta}_j = (-\bm{z}^\top_i(-\bm{\theta}_j)}, and so is actually unidentifiable. 

As a result of sign-switching, different MCMC chains can produce latent variables and corresponding coefficients values that, while having similar magnitudes, will be different in sign. Consequently, combining MCMC chains and checking Rhats, computing posterior means and medians etc...becomes complicated (in principle, one way to resolve this problem would be to post-process the MCMC chains and deal with sign switching, but this is really hard!). Therefore, to alleviate this issue together, boral chooses to only run one MCMC chain.

What does this mean for the user? 
\itemize{
\item boral automatically calculates the Geweke convergence diagnostic (Geweke, 1992), which is a diagnostic applicable with only one MCMC chain; please see the help file \code{geweke.diag} in the \code{coda} package for more information. The output is a list containing Z-scores for the appropriate parameters in the model, and each score can be interpreted in the same manner as the test statistic from conducting a Z-test i.e., if the score exeeds roughly 1.96 then the p-value is less than 0.05, and there is evidence the MCMC chain (for this particular parameter) has not converged. 

The output from boral also provides the proportion of Z-scores whose corresponding p-values are less than 0.05. Of course, because there are a large number of parameters in the model, then there are large number of Z-scores, and boral does not make any multiple comparison adjustment for this when calculating the number of ``significant" Z-scores. If you do indeed want to use this diagnostic to formally check for convergence, then we recommend you conduct some adjustment e.g., using Holm's method, by doing something such as \cr \code{gew.pvals <- 2*pnorm(abs(unlist(fit$geweke.diag[[1]])), lower.tail = FALSE)} and then \code{p.adjust(gew.pvals, method = "holm")}.
 
\item For checking convergence, we recommend you look at trace plots of the MCMC chains. Using the \code{coda} package, which is automatically loaded when the \code{boral} package is loaded, try something like \code{plot(get.mcmcsamples(fit))}. 

\item If you have a lot of data, e.g. lots of sites compared to species, sign-switching tends to be less of problem and pops up less often.

%\item IMPORTANTLY, if the goal of your analysis is to inference while account for residual correlations between the columns of \code{y}, and not for model-based ordination, then the sign-switching problem is not a problem at all! This is because while the signs of the latent variables and associated coefficients may switch, the correlation and their signs are unaffected. In other words, looking the point estimates and credible intervals of regression coefficients \eqn{\bm{\beta}_j}, and functions like \code{get.residual.cor} are unaffected by sign-switching. 
}

}

\value{
An object of class "boral" is returned, being a list containing the following components where applicable:
     \item{call}{The matched call.}
     
     \item{lv.coefs.mean/median/sd/iqr}{Matrices containing the mean/median/standard deviation/interquartile range of the posterior distributions of the latent variable coefficients. This also includes the column-specific intercepts, and dispersion parameters if appropriate.}
     
     \item{lv.mean/median/sd/iqr}{A matrix containing the mean/median/standard deviation/interquartile range of the posterior distributions of the latent variables.}
     
     \item{X.coefs.mean/median/sd/iqr}{Matrices containing the mean/median/standard deviation/interquartile range of the posterior distributions of the column-specific coefficients relating to \code{X}.}
     
     \item{traits.coefs.mean/median/sd/iqr}{Matrices containing the mean/median/standard deviation/interquartile range of the posterior distributions of the coefficients and standard deviation relating to the species traits; please see \code{\link{about.traits}}.}

     \item{cutoffs.mean/median/sd/iqr}{Vectors containing the mean/median/standard deviation/interquartile range of the posterior distributions of the common cutoffs for ordinal responses (please see the not-so-brief tangent on distributions above).}
     
	\item{ordinal.sigma.mean/median/sd/iqr}{Scalars containing the mean/median/standard deviation/interquartile range of the posterior distributions of the standard deviation for the random intercept normal distribution corresponding to the ordinal response columns.}

     \item{powerparam.mean/median/sd/iqr}{Scalars for the mean/median/standard deviation/interquartile range of the posterior distributions of the common power parameter for tweedie responses (please see the not-so-brief tangent on distributions above).}
     
     \item{row.coefs.mean/median/sd/iq}{A list with each element containing the vectors of mean/median/standard deviation/interquartile range of the posterior distributions of the row effects. The length of the list is equal to the number of row effects included i.e., \code{ncol(row.ids)}.}
     
	\item{row.sigma.mean/median/sd/iqr}{A list with each element containing the mean/median/standard deviation/interquartile range of the posterior distributions of the standard deviation for the row random effects normal distribution. The length of the list is equal to the number of row effects included i.e., \code{ncol(row.ids)}.}
	
     \item{ssvs.indcoefs.mean/ssvs.indcoefs.sd}{Matrices containing posterior probabilities and associated standard deviation for individual SSVS of coefficients in \code{X}.}
     
     \item{ssvs.gpcoefs.mean/ssvs.gpcoefs.sd}{Matrices containing posterior probabilities and associated standard deviation for group SSVS of coefficients in \code{X}.}
     
     \item{ssvs.traitscoefs.mean/ssvs.traitscoefs.sd}{Matrices containing posterior probabilities and associated standard deviation for individual SSVS of coefficients relating to species traits.}

     \item{hpdintervals}{A list containing components which correspond to the lower and upper bounds of highest posterior density (HPD) intervals for all the parameters indicated above. Please see \code{\link{get.hpdintervals}} for more details.}
     
     \item{ics}{If \code{calc.ics = TRUE}, then a list of different information criteria values for the model calculated using \code{\link{get.measures}} is run. Please see \code{\link{get.measures}} for details regarding the criteria. Also, please note the ics returned are based on \code{\link{get.measures}} with \code{more.measures = FALSE}.}
     
     \item{jags.model}{If \code{save.model = TRUE}, the raw jags model fitted is returned. This can be quite large!}
     
     \item{geweke.diag}{A list with two elements. The first element is itself a list containing the Geweke convergence diagnostic (Z-scores) for all appropriate parameters in the model. The second element contains the proportion of Z-scores that whose corresponding p-value is less than 0.05. No adjustment is made for multiple comparison on the p-values. Please see the section \emph{Why is only one MCMC chain run?} for more information on this diagnostic.}

     \item{n, p, family, trial.size, num.lv, ...}{Various attributes of the model fitted, including the dimension of \code{y}, the response and model matrix used, distributional assumptions and trial sizes, number of latent variables, the number of covariates and traits, hyperparameters used in the Bayesian estimation, indices for SSVS, the number of levels for ordinal responses, and n.burnin, n.iteration and n.thin.}
}

\references{
\itemize{
\item Gelman A. (2006) Prior distributions for variance parameters in hierarchical models. Bayesian Analysis 1, 515-533.

\item Gelman, et al. (2008). A weakly informative default prior distribution for logistic and other regression models. The Annals of Applied Statistics, 2, 1360-1383.

\item Gelman et al. (2013) Bayesian Data Analysis. CRC Press.

\item George, E. I. and McCulloch, R. E. (1993). Variable selection via Gibbs sampling. Journal of the American Statistical Association, 85, 398-409.

\item Geweke, J. (1992) Evaluating the accuracy of sampling-based approaches to calculating posterior moments. In Bayesian Statistics 4 (editors JM Bernado, JO Berger, AP Dawid and AFM Smith). Clarendon Press.

\item Hui et al. (2014). Model-based approaches to unconstrained ordination. Methods in Ecology and Evolution, 6, 399-411.

\item Hill, M. O. (1974). Correspondence analysis: a neglected multivariate method. Applied statistics, 23, 340-354.

\item Jamil, T., and ter Braak, C.J.F. (2013). Generalized linear mixed models can detect unimodal species-environment relationships. PeerJ 1: e95.

\item Kruskal, J. B. (1964). Nonmetric multidimensional scaling: a numerical method. Psychometrika, 29, 115-129.

\item Legendre, P. and Legendre, L. (2012). Numerical ecology, Volume 20. Elsevier.

\item Martin et al. (2011). MCMCpack: Markov Chain Monte Carlo in R. Journal of Statistical Software, 42, 1-21. URL: http://www.jstatsoft.org/v42/i09/.

\item McLachlan, G., and Peel, D. (2004). Finite Mixture Models. Wiley.

\item Plummer, M. (2003). JAGS: A program for analysis of Bayesian graphical models using Gibbs sampling. In Proceedings of the 3rd International Workshop on Distributed Statistical Computing. March (pp. 20-22).

\item Pollock, L. J. et al. (2014). Understanding co-occurrence by modelling species simultaneously with a Joint Species Distribution Model (JSDM). Methods in Ecology and Evolution, 5, 397-406.

\item Skrondal, A., and Rabe-Hesketh, S. (2004). Generalized latent variable modeling: Multilevel, longitudinal, and structural equation models. CRC Press.

\item Warton et al. (2015). So Many Variables: Joint Modeling in Community Ecology. Trends in Ecology and Evolution, 30, 766-779.

\item Warton et al. (2012). Distance-based multivariate analyses confound location and dispersion effects. Methods in Ecology and Evolution, 3, 89-101.

\item Wang et al. (2013). \code{mvabund}: statistical methods for analysing multivariate abundance data. R package version 3.8.4.}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}

\section{Warnings}{
\itemize{
\item \emph{No} intercept column is required in \code{X}. Column-specific intercepts are estimated automatically and given by the first column of \code{lv.coefs}. Similarly, \emph{no} intercept column is required in \code{traits}, as it is included automatically.

\item If num.lv > 5, a warning is printed asking whether you really want to fit an boral with more than five latent variables. A warning is also printed if num.lv == 1, as this is not going to be successful in modeling between the correlation between columns.

\item For models including both explanatory covariates and latent variables, one requires \code{num.lv > 1} to allow flexible modeling of the residual correlation matrix. 

\item As of version 1.5, functions to calculate information criteria will no longer be updated...use \code{calc.ics = TRUE} at your own peril!!!

\item MCMC with lots of ordinal columns take an especially long time to run! Moreover, estimates for the cutoffs in cumulative probit regression may be poor for levels with little data. Major apologies for this advance =(

\item Summaries of the coefficients such as posterior medians and HPD intervals may also be problematic when SSVS is being used, since the posterior distribution will be multi-modal. 

\item If \code{save.model = TRUE}, the raw jags model is also returned. This can be quite very memory-consuming, since it indirectly saves all the MCMC samples.
}
}

\seealso{
\code{\link{lvsplot}} for a scatter plot of the latent variables (and their coefficients if applicable) when \code{num.lv = 1} or \code{2}, \code{\link{coefsplot}} for horizontal line or "caterpillar plot" of the regression coefficients corresponding to \code{X} (if applicable), \code{\link{summary.boral}} for a summary of the fitted boral model, \code{\link{get.enviro.cor}} and \code{\link{get.residual.cor}} for calculating the correlation matrix between the explanatory variables in \code{X} and the residual correlation matrix respectively, and \code{\link{calc.varpart}} to calculate variance partitioning of the explanatory variables in \code{X}. 
}

\examples{
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

     
## Example 1 - model with two latent variables, site effects, 
## 	and no environmental covariates
spiderfit_nb <- boral(y, family = "negative.binomial", num.lv = 2, 
     row.eff = "fixed", mcmc.control = example_mcmc_control)

summary(spiderfit_nb)

par(mfrow = c(2,2))
plot(spiderfit_nb) ## Plots used in residual analysis, 
## Used to check if assumptions such an mean-variance relationship 
## are adequately satisfied.

lvsplot(spiderfit_nb) ## Biplot of the latent variables, 
## which can be interpreted in the same manner as an ordination plot.

\dontrun{
## Example 2a - model with no latent variables, no site effects, 
## 	and environmental covariates
spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
	num.lv = 0, mcmc.control = example_mcmc_control)

summary(spiderfit_nb) 
## The results can be compared with the default example from 
## the manyglm() function in mvabund. 

## Caterpillar plots for the coefficients
par(mfrow=c(2,3), mar = c(5,6,1,1))
sapply(colnames(spiderfit_nb$X), coefsplot, x = spiderfit_nb)


## Example 2b - suppose now, for some reason, the 28 rows were
## 	sampled such into four replications of seven sites
## Let us account for this as a fixed effect
spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
	num.lv = 0, row.eff = "fixed", row.ids = matrix(rep(1:7,each=4),ncol=1),
	mcmc.control = example_mcmc_control)

spiderfit_nb$row.coefs

## Example 2c - suppose now, for some reason, the 28 rows reflected
## 	a nested design with seven regions, each with four sub-regions
## We can account for this nesting as a random effect
spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
	num.lv = 0, row.eff = "random", 
	row.ids = cbind(1:n, rep(1:7,each=4)), 
	mcmc.control = example_mcmc_control)

spiderfit_nb$row.coefs


## Example 3a - Extend example 2 to demonstrate grouped covariate selection
## on the last three covariates. 
example_prior_control <- list(type = c("normal","normal","normal","uniform"), 
     ssvs.index = c(-1,-1,-1,1,2,3))
spiderfit_nb2 <- boral(y, X = X, family = "negative.binomial", 
	num.lv = 0, mcmc.control = example_mcmc_control,
	prior.control = example_prior_control)
     
summary(spiderfit_nb2) 


## Example 3b - Extend example 2 to demonstrate individual covariate selection
## on the last three covariates. 
example_prior_control <- list(type = c("normal","normal","normal","uniform"), 
     ssvs.index = c(-1,-1,-1,0,0,0))
spiderfit_nb3 <- boral(y, X = X, family = "negative.binomial", 
	num.lv = 0, mcmc.control = example_mcmc_control,
	prior.control = example_prior_control)
summary(spiderfit_nb3) 

	
## Example 4 - model fitted to presence-absence data, no site effects, and
## two latent variables
data(tikus)
y <- tikus$abun
y[y > 0] <- 1
y <- y[1:20,] ## Consider only years 1981 and 1983
y <- y[,apply(y > 0,2,sum) > 2] ## Consider only spp with more than 2 presences
     
tikus.fit <- boral(y, family = "binomial", num.lv = 2, 
	mcmc.control = example_mcmc_control)
     
lvsplot(tikus.fit, biplot = FALSE) 
## A strong location between the two sampling years 


## Example 5a - model fitted to count data, no site effects, and
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

fit_traits <- boral(y, X = X, traits = traits, 
     which.traits = example_which_traits, family = "negative.binomial",
     mcmc.control = example_mcmc_control, save.model = TRUE)

summary(fit_traits)


## Example 5b - perform selection on trait coefficients
ssvs_traitsindex <- vector("list",ncol(X)+1)
for(i in 1:length(ssvs_traitsindex)) 
     ssvs_traitsindex[[i]] <- rep(0,ncol(traits))
ssvs_traitsindex[[3]] <- -1
fit_traits <- boral(y, X = X, traits = traits, which.traits = example_which_traits, 
	family = "negative.binomial", mcmc.control = example_mcmc_control, 
	save.model = TRUE, prior.control = list(ssvs.traitsindex = ssvs_traitsindex))

summary(fit_traits)


## Example 6 - simulate Bernoulli data, based on a model with two latent variables, 
## no site variables, with two traits and one environmental covariates 
## This example is a proof of concept that traits can used to 
## explain environmental responses 
library(mvtnorm)

n <- 100; s <- 50
X <- as.matrix(scale(1:n))
colnames(X) <- c("elevation")

traits <- cbind(rbinom(s,1,0.5), rnorm(s)) 
## one categorical and one continuous variable
colnames(traits) <- c("thorns-dummy","SLA")

simfit <- list(true.lv = rmvnorm(n, mean = rep(0,2)), 
	lv.coefs = cbind(rnorm(s), rmvnorm(s, mean = rep(0,2))), 
	traits.coefs = matrix(c(0.1,1,-0.5,1,0.5,0,-1,1), 2, byrow = TRUE))
rownames(simfit$traits.coefs) <- c("beta0","elevation")
colnames(simfit$traits.coefs) <- c("kappa0","thorns-dummy","SLA","sigma")

simy = create.life(true.lv = simfit$true.lv, lv.coefs = simfit$lv.coefs, X = X, 
	traits = traits, traits.coefs = simfit$traits.coefs, family = "binomial") 


example_which_traits <- vector("list",ncol(X)+1)
for(i in 1:length(example_which_traits)) 
     example_which_traits[[i]] <- 1:ncol(traits)
fit_traits <- boral(y = simy, X = X, traits = traits, 
     which.traits = example_which_traits, family = "binomial", 
     num.lv = 2, save.model = TRUE, mcmc.control = example_mcmc_control)

}

}
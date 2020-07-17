\name{make.jagsboralnullmodel}

\alias{make.jagsboralnullmodel}

\docType{package}

\title{Write a text file containing a model for use into JAGS}

\description{
\Sexpr[results=rd, stage=render]{lifecycle::badge("stable")}

This function is designed to write models with no latent variables i.e., so-called "null" models.}

\usage{make.jagsboralnullmodel(family, num.X = 0, X.ind = NULL, num.traits = 0, 
   which.traits = NULL, row.eff = "none", row.ids = NULL,
   offset = NULL, trial.size = 1, n, p, model.name = NULL, 
   prior.control = list(type = c("normal","normal","normal","uniform"), 
	hypparams = c(10, 10, 10, 30), ssvs.index = -1, ssvs.g = 1e-6,
	ssvs.traitsindex = -1))
	}

\arguments{

  \item{family}{Either a single element, or a vector of length equal to the number of columns in \code{y}. The former assumes all columns of \code{y} come from this distribution. The latter option allows for different distributions for each column of \code{y}. Elements can be one of "binomial" (with probit link), "poisson" (with log link), "negative.binomial" (with log link), "normal" (with identity link), "lnormal" for lognormal (with log link), "tweedie" (with log link), "exponential" (with log link), "gamma" (with log link), "beta" (with logit link), "ordinal" (cumulative probit regression), "ztpoisson" (zero truncated Poisson with log link), "ztnegative.binomial" (zero truncated negative binomial with log link). 

  Please see \code{\link{about.distributions}} for information on distributions available in boral overall.
  }

\item{num.X}{Number of columns in \code{X}. Defaults to 0, in which case it is assumed that no covariates are included in the model. Recall that no intercept should be included in \code{X}.}

  \item{X.ind}{An matrix of 1s and 0s, indicating whether a particular covariate should be included (1) or excluded (0) in the mean structure of a particular response. The matrix should the number of rows equal to the number of columns in \code{y}, and the number of columns equal to the number of columns in \code{X}. Defaults to \code{NULL}, in which case it is assumed that all covariates are included in the mean structure of all responses i.e., all 1s.}  

  \item{num.traits}{Number of columns in the model matrix \code{traits}. Defaults to 0, in which case it is assumed no traits are included in model. Recall that no intercept should be included in \code{traits}.}

  \item{which.traits}{A list of length equal to (number of columns in \code{X} + 1), informing which columns of \code{traits} the column-specific intercepts and each of the column-specific regression coefficients should be regressed against. The first element in the list applies to the column-specific intercept, while the remaining elements apply to the regression coefficients. Each element of \code{which.traits} is a vector indicating which traits are to be used. 
  
  For example, if \code{which.traits[[2]] = c(2,3)}, then the regression coefficients corresponding to the first column in \code{X} are regressed against the second and third columns of \code{traits}. If \code{which.traits[[2]][1] = 0}, then the regression coefficients for each column are treated as independent. Please see \code{\link{about.traits}} for more details.
  
  Defaults to \code{NULL}, and used in conjunction with \code{traits} and \cr \code{prior.control$ssvs.traitsindex}.}
  
  \item{row.eff}{Single element indicating whether row effects are included as fixed effects ("fixed"), random effects ("random") or not included ("none") in the fitted model. If fixed effects, then for parameter identifiability the first row effect is set to zero, which analogous to acting as a reference level when dummy variables are used. If random effects, they are drawn from a normal distribution with mean zero and unknown standard deviation. Defaults to "none". } 

 \item{row.ids}{A matrix with the number of rows equal to the number of rows in \code{y}, and the number of columns equal to the number of row effects to be included in the model. Element \eqn{(i,j)} indicates to the cluster ID of row \eqn{i} in \code{y} for random effect eqn{j}; please see \code{\link{boral}} for more details. for details. Defaults to \code{NULL}, so that if \code{row.eff = "none"} then the argument is ignored, otherwise if \cr \code{row.eff = "fixed"} or \code{"random"}, \cr then \code{row.ids = matrix(1:nrow(y), ncol = 1)} i.e., a single, row effect unique to each row.}    

 \item{offset}{A matrix with the same dimensions as the response matrix \code{y}, specifying an a-priori known component to be included in the linear predictor during fitting. Defaults to \code{NULL}.}

 \item{trial.size}{Either equal to a single element, or a vector of length equal to the number of columns in y. If a single element, then all columns assumed to be binomially distributed will have trial size set to this. If a vector, different trial sizes are allowed in each column of y. The argument is ignored for all columns not assumed to be binomially distributed. Defaults to 1, i.e. Bernoulli distribution.}

 \item{n}{The number of rows in the response matrix \code{y}.}

  \item{p}{The number of columns in the response matrix \code{y}.}

 \item{model.name}{Name of the text file that the JAGS model is written to. Defaults to \code{NULL}, in which case the default of "jagsboralmodel.txt" is used.}

 \item{prior.control}{A list of parameters for controlling the prior distributions. These include:
	\itemize{
		\item \emph{type:} Vector of four strings indicating the type of prior distributions to use. In order, these are: 1) priors for all column-specific intercepts, row effects, and cutoff points for ordinal data; 2) priors for the latent variable coefficients and correlation parameters. This is ignored for this function; 3) priors for all column-specific coefficients relating to \code{X} (ignored if \code{X = NULL}). When traits are included in the model, this is also the prior for the trait regression coefficients (please see \code{\link{about.traits}} for more information); 4) priors for any dispersion parameters and variance (standard deviation, to be precise) parameters in the model.
				
		For elements 1-3, the prior distributions currently available include: I) ``normal", which is a normal prior with the variance controlled by elements 1-3 in \code{hypparams}; II) ``cauchy", which is a Cauchy prior with variance controlled by elements 1-3 in \code{hypparams}. Gelman, et al. (2008) considers using Cauchy priors with variance \eqn{2.5^2} as weakly informative priors for coefficients in logistic and potentially other generalized linear models; III) ``uniform", which is a symmetric uniform prior with minimum and maximum values controlled by element 1-3 in \code{hypparams}. 
		
		For element 4, the prior distributions currently available include: I) ``uniform", which is uniform prior with minimum zero and maximum controlled by element 4 in \code{hypparmas}; II) ``halfnormal", which is half-normal prior with variance controlled by \code{hypparams}; III) ``halfcauchy", which is a half-Cauchy prior with variance controlled by element 4 in \code{hypparams}.
		
		Defaults to the vector \code{c("normal","normal","normal","uniform")}. 

		\item \emph{hypparams} Vector of four hyperparameters used in the set up of prior distributions. In order, these are: 1) affects the prior distribution for all column-specific intercepts, row effects, and cutoff points for ordinal data; 2) affects the prior distribution for all latent variable coefficients and correlation parameters. This is ignored for this function; 3) affects the prior distribution for column-specific coefficients relating to \code{X} (ignored if \code{X = NULL}). When traits are included in the model, it also affects the prior distribution for the trait regression coefficients; 4) affects the prior distribution for any dispersion parameters, as well as the prior distributions for the standard deviation of the random effects normal distribution if \code{row.eff = "random"}, the standard deviation of the column-specific random intercepts for these columns if more than two of the columns are ordinal, and the standard deviation of the random effects normal distribution for trait regression coefficients when traits are included in the model.

		Defaults to the vector \code{c(10, 10, 10, 30)}. The use of normal distributions with mean zero and variance 10 as priors is seen as one type of (very) weakly informative prior, according to \href{https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations}{Prior choice recommendations}.
		
		\item \emph{ssvs.index:} Indices to be used for stochastic search variable selection (SSVS, George and McCulloch, 1993). Either a single element or a vector with length equal to the number of columns in the implied model matrix \code{X}. Each element can take values of -1 (no SSVS is performed on this covariate), 0 (SSVS is performed on individual coefficients for this covariate), or any integer greater than 0 (SSVS is performed on collectively all coefficients on this covariate/s.) 
		
		Please see \code{\link{about.ssvs}} for more information regarding the implementation of SSVS. Defaults to -1, in which case SSVS is not performed on \code{X} variables. 

		\item \emph{ssvs.g:} Multiplicative, shrinkage factor for SSVS, which controls the strength of the "spike" in the SSVS mixture prior. In summary, if the coefficient is included in the model, the "slab" prior is a normal distribution with mean zero and variance given by element 3 in \code{hypparams}, while if the coefficient is not included in the model, the "spike" prior is normal distribution with mean zero and variance given by element 3 in \code{hypparams} multiplied by \code{ssvs.g}. Please see \code{\link{about.ssvs}} for more information regarding the implementation of SSVS. Defaults to 1e-6.  		
		
		\item \emph{ssvs.traitsindex:} Used in conjunction with \code{traits} and \code{which.traits}, this is a list of indices to be used 
		for performing SSVS on the trait coefficients. Should be a list with the same length as \code{which.traits}, and with each element a vector of indices with the same length as the corresponding element in \code{which.traits}. Each index either can take values of -1 (no SSVS on this trait coefficient) or 0 (no SSVS on this trait coefficient). 

		Please see \code{\link{about.ssvs}} for more information regarding the implementation of SSVS. Defaults to -1, in which case SSVS is not performed on any of the trait coefficients, if they are included in the model.
		}
	}

 }

\details{
This function is automatically executed inside \code{\link{boral}}, and therefore does not need to be run separately before fitting the model. It can however be run independently if one is: 1) interested in what the actual JAGS file for a particular model looks like, 2) wanting to modify a basic JAGS model file to construct more complex model e.g., include environmental variables. 

Please note that \code{\link{boral}} currently does not allow the user to manually enter a script to be run. 

When running the main function \code{\link{boral}}, setting \code{save.model = TRUE} which automatically save the JAGS model file as a text file (with name based on the \code{model.name}) in the current working directory.
}

\value{
A text file is created, containing the JAGS model to be called by the boral function for entering into jags. This file is automatically deleted once boral has finished running unless \code{save.model = TRUE}.}

\references{
\itemize{
\item Gelman, et al. (2008). A weakly informative default prior distribution for logistic and other regression models. The Annals of Applied Statistics, 2, 1360-1383.
}
}

\author{
\packageAuthor{boral}

Maintainer: \packageMaintainer{boral}
}

\seealso{
\code{\link{make.jagsboralmodel}} for writing JAGS scripts for models with one or more latent variables.}

\examples{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
n <- nrow(y)
p <- ncol(y)

testpath <- file.path(tempdir(), "jagsboralmodel.txt")


## Create a "null" model JAGS script, where distributions alternative 
## between Poisson and negative distributions 
##   across the rows of y.
make.jagsboralnullmodel(family = rep(c("poisson","negative.binomial"),length=p), 
    num.X = ncol(spider$x), row.eff = "fixed", n = n, p = p,
    model.name = testpath)

     
## Create a "null" model JAGS script, where distributions are all negative 
## 	binomial distributions and covariates will be included!
make.jagsboralnullmodel(family = "negative.binomial", 
    num.X = ncol(spider$x), n = n, p = p, 
    model.name = testpath)

     
## Have a look at the JAGS model file for a model involving traits,
## based on the ants data from mvabund.
library(mvabund)
data(antTraits)

y <- antTraits$abun
X <- as.matrix(antTraits$env)
## Include only traits 1, 2, and 5, plus an intercept
traits <- as.matrix(antTraits$traits[,c(1,2,5)])
## Please see help file for boral regarding the use of which.traits
example_which_traits <- vector("list",ncol(X)+1)
for(i in 1:length(example_which_traits)) 
     example_which_traits[[i]] <- 1:ncol(traits)

\dontrun{
## NOTE: The values below MUST NOT be used in a real application;
## they are only used here to make the examples run quick!!!
example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
     n.thin = 1)

fit_traits <- boral(y, X = X, traits = traits, which.traits = example_which_traits, 
    family = "negative.binomial", model.name = testpath, 
     mcmc.control = example_mcmc_control, do.fit = FALSE)
}
}

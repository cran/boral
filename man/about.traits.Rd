\name{about.traits}

\alias{about.traits}

\title{Including species traits in boral}

\description{This help file provides more information regarding the how species can be included to help mediate environmental responses, analogous to the so-called fourth corner problem.}

\details{
In the main boral function, when covariates \code{X} are included i.e. both the independent and correlated response models, one has the option of also including traits to help explain differences in species environmental responses to these covariates. Specifically, when \code{traits} and \code{which.traits} are supplied, then the \eqn{\beta_{0j}}'s and \eqn{\bm{\beta}_j}'s are then regarded as random effects drawn from a normal distribution. For the species-specific intercepts, we have

\deqn{\beta_{0j} \sim N(\kappa_{01} + \bm{traits}^T_j\bm{\kappa}_1, \sigma^2_1),}

where \eqn{(\kappa_{01},\bm{\kappa}_1)} are the regression coefficients relating to the traits to the intercepts and \eqn{\sigma_1} is the error standard deviation. These are now the "parameters" in the model, in the sense that priors are assigned to them and MCMC sampling is used to estimate them (see the next section on estimation). 

In an analogous manner, each of the elements in \eqn{\bm{\beta}_j = (\beta_{j1},\ldots,\beta_{jd})} are now drawn as random effects from a normal distribution. That is, for \eqn{k = 1,\ldots,d} where \code{d = ncol(X)}, we have,

\deqn{\beta_{jk} \sim N(\kappa_{0k} + \bm{traits}^T_j\bm{\kappa}_k, \sigma^2_k),}

Which traits are to included (regressed) in the mean of the normal distributions is determined by the list argument \code{which.traits} in the main boral function. The first element in the list applies to \eqn{beta_{0j}}, while the remaining elements apply to the the \eqn{\bm{\beta}_j}. Each element of \code{which.traits} is a vector indicating which traits are to be used. For example, if \code{which.traits[[2]] = c(2,3)}, then the \eqn{\beta_{j1}}'s are drawn from a normal distribution with mean depending only on the second and third columns of \code{traits}. If \code{which.traits[[2]] = 0}, then the regression coefficients are treated as independent, i.e. the values of \eqn{\beta_{j1}} are given their own priors and estimated separately from each other. 

Including species traits in the model can be regarded as a method of simplifying the model: rather than each to estimates \eqn{p} sets of species-specific coefficients, we instead say that these coefficients are linearly related to the corresponding values of their traits (Warton et al., 2015). That is, we are using trait data to help explain similarities/differences in species responses to the environment. This idea has close relations to the fourth corner problem in ecology (Brown et al., 2014). Unlike the models of Brown et al. (2014) however, which treat the \eqn{\beta_{0j}}'s and \eqn{\beta_{jk}}'s are fixed effects and fully explained by the traits, boral adopts a random effects approach similar to Jamil et al. (2013) to "soak up" any additional between species differences in environmental responses not explained by traits.
}


\references{
\itemize{
\item Brown, et al. (2014). The fourth-corner solution - using predictive models to understand how species traits interact with the environment. Methods in Ecology and Evolution 5, 344-352.

\item Jamil, T. et al. (2013). Selecting traits that explain species-environment relationships: a generalized linear mixed model approach. Journal of Vegetation Science 24, 988-1000

\item Warton et al. (2015). So Many Variables: Joint Modeling in Community Ecology. Trends in Ecology and Evolution, to appear.
}
}

\section{Warnings}{
\itemize{
\item \emph{No} intercept column should be required in \code{traits}, as it is included automatically.
}
}

\author{
Francis K.C. Hui \email{fhui28@gmail.com}
}


\seealso{
\code{\link{boral}} for the main boral fitting function. 
}

\examples{
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
X <- scale(spider$x)
n <- nrow(y); p <- ncol(y); 

## NOTE: The two examples below and taken directly from the boral help file

example.mcmc.control <- list(n.burnin = 10, n.iteration = 100, 
     n.thin = 1)

\dontrun{
## Example 5 - model fitted to count data, no site effects, and
## two latent variables, plus traits included to explain environmental responses
data(antTraits)
y <- antTraits$abun
X <- as.matrix(scale(antTraits$env))
## Include only traits 1, 2, and 5
traits <- as.matrix(antTraits$traits[,c(1,2,5)])
which.traits <- vector("list",ncol(X)+1)
for(i in 1:length(which.traits)) which.traits[[i]] <- 1:ncol(traits)
## Just for fun, the regression coefficients for the second column of X,
## corresponding to the third element in the list which.traits,
## will be estimated separately and not regressed against traits.
which.traits[[3]] <- 0

fit.traits <- boral(y, X = X, traits = traits, which.traits = which.traits, 
	family = "negative.binomial", calc.ics = FALSE, 
	mcmc.control = example.mcmc.control)

summary(fit.traits)


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


which.traits <- vector("list",ncol(X)+1)
for(i in 1:length(which.traits)) which.traits[[i]] <- 1:ncol(traits)
fit.traits <- boral(y = simy, X = X, traits = traits, which.traits = which.traits, 
	family = "binomial", num.lv = 2, save.model = TRUE, calc.ics = FALSE, 
	mcmc.control = example.mcmc.control)
}

}
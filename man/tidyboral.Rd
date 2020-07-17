\name{tidyboral}

\alias{tidyboral}

\docType{package}

\title{Reformats output from a boral fit}

\description{
\Sexpr[results=rd, stage=render]{lifecycle::badge("maturing")}

Reformats estimates from a fitted boral model to be in a slightly tidier format, meaning everything is presented as a long data frame.}

\usage{
tidyboral(object)
}

\arguments{
  \item{object}{An object of class "boral".}
}

\details{
Formatting the output into a long data frame maybe useful if one wishes to take the estimated parameters (namely the posterior mean/median/interquartile range/standard deviation, and the lower and upper limits
of the HPD intervals), and subsequently wrangle them for presentation purposes using packages from the tidyverse e.g., Wickham and Henry (2018), construct plots from them using the ggplot2 package (Wickham, 2016), and so on.

It is important to note that this function is solely designed to format output from a fitted boral model relating to the estimated parameters. It does not an additional information like model call and MCMC samples. 
Please do NOT erase the original fitted model in place of this!
}


\value{
A a list containing the following components where applicable:
     \item{lv.coefs}{A long format data frame containing the mean/median/standard deviation/interquartile range of the posterior distributions of the latent variable coefficients. This also includes the column-specific intercepts, and dispersion parameters if appropriate.}
     
     \item{lv}{A long format data frame containing the mean/median/standard deviation/interquartile range of the posterior distributions of the latent variables.}
     
     \item{lv.covparams}{A long format data frame containing the mean/median/standard deviation/interquartile range of the posterior distributions for the parameters characterizing the correlation structure of the latent variables when they are assumed to be non-independent across rows.}

     \item{X.coefs}{A long format data frame containing the mean/median/standard deviation/interquartile range of the posterior distributions of the column-specific coefficients relating to \code{X}.}
     
     \item{traits.coefs}{A long format data frame containing the mean/median/standard deviation/interquartile range of the posterior distributions of the coefficients and standard deviation relating to the species traits; please see \code{\link{about.traits}}.}

     \item{cutoffs}{A long format data frame containing the mean/median/standard deviation/interquartile range of the posterior distributions of the common cutoffs for ordinal responses (please see the not-so-brief tangent on distributions above).}
     
	\item{ordinal.sigma}{A long format data frame containing the mean/median/standard deviation/interquartile range of the posterior distributions of the standard deviation for the random intercept normal distribution corresponding to the ordinal response columns.}

     \item{powerparam}{A long format data frame containing the mean/median/standard deviation/interquartile range of the posterior distributions of the common power parameter for tweedie responses (please see the not-so-brief tangent on distributions above).}
     
     \item{row.coefs}{A list with each element being a long format data frame containing the mean/median/standard deviation/interquartile range of the posterior distributions of the row effects. The length of the list is equal to the number of row effects included i.e., \code{ncol(row.ids)}.}
     
	\item{row.sigma}{A list with each element being a long format data frame containing the mean/median/standard deviation/interquartile range of the posterior distributions of the standard deviation for the row random effects normal distribution. The length of the list is equal to the number of row effects included i.e., \code{ncol(row.ids)}.}
	
     \item{ssvs.indcoefs}{A long format data frame containing posterior probabilities and associated standard deviation for individual SSVS of coefficients in \code{X}.}
     
     \item{ssvs.gpcoefs}{A long format data frame containing posterior probabilities and associated standard deviation for group SSVS of coefficients in \code{X}.}
     
     \item{ssvs.traitscoefs}{A long format data frame containing posterior probabilities and associated standard deviation for individual SSVS of coefficients relating to species traits.}

     \item{hpdintervals}{A list containing long format data frames corresponding to the lower and upper bounds of highest posterior density (HPD) intervals for all the parameters indicated above. Please see \code{\link{get.hpdintervals}} for more details.}
     }

\references{
\itemize{
\item Wickham, H. (2016). ggplot2: elegant graphics for data analysis. Springer.

\item Wickham, H., & Henry, L. (2017). Tidyr: Easily tidy data with `spread ()' and `gather ()' functions.}
}

\author{
\packageAuthor{boral}

Maintainer: \packageMaintainer{boral}
}

\section{Warnings}{
\itemize{
\item This function is solely designed to format output from a fitted boral model relating to the estimated parameters. It does not an additional information like model call and MCMC samples. Please do NOT erase the original fitted model in place of this!
}
}

\seealso{
\code{\link{boral}} for the fitting function on which \code{tidyboral} can be applied.
}

\examples{
\dontrun{
## NOTE: The values below MUST NOT be used in a real application;
## they are only used here to make the examples run quick!!!
example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
    n.thin = 1)

testpath <- file.path(tempdir(), "jagsboralmodel.txt")


library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun

spiderfit_nb <- boral(y, family = "negative.binomial", 
     lv.control = list(num.lv = 2), row.eff = "fixed", 
     mcmc.control = example_mcmc_control, model.name = testpath)

spiderfit_nb_tidy <- tidyboral(spiderfit_nb)
}
}

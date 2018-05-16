\name{coefsplot}

\alias{coefsplot}

\title{Caterpillar plots of the regression coefficients from a boral model} 

\description{Constructs horizontal line plot (point estimate and HPD intervals), otherwise known as "caterpillar plots", for the column-specific regression coefficients corresponding to a covariate in \code{X} fitted in the boral model.}

\usage{
coefsplot(covname, x, labely = NULL, est = "median", ...) 
	}
	
\arguments{
  \item{covname}{The name of one of the covariates fitted in the boral model. That is, it must be a character vector corresponding to one of the elements in \code{colnames(x)$X.coefs.median}.}

  \item{x}{An object for class "boral".}
             
  \item{labely}{Controls the labels on the y-axis for the line plot. If it is not \code{NULL}, then it must be a vector either of length 1 or the same length as the number of columns in the \code{y} in the fitted boral object. In the former, it is treated as the y-axis label. In the latter, it is used in place of the column names of \code{y} to label each line. Defaults to \code{NULL}, in which the each line in the plot is labeled according to the columns of \code{y}, or equivalently \code{rownames(x$X.coefs.median)}.}

  \item{est}{A choice of either the posterior median (\code{est = "median"}) or posterior mean (\code{est = "mean"}), which are then used as the point estimates in the lines. Default is posterior median.}

   \item{...}{Additional graphical options to be included in. These include values for \cr \code{cex, cex.lab, cex.axis, cex.main, lwd}, and so on.}
}

\details{
For each species (column of \code{y}), the horizontal line or "caterpillar" is constructed by first marking the point estimate (posterior mean or median) with an "x" symbol. Then the line is construed based on the lower and upper limits of the highest posterior density (HPD) intervals as found in \code{x$hpdintervals}. By default these intervals of 95\% HPD intervals. To complete the plot, a vertical dotted line is drawn to denote the zero value. All HPD intervals that include zero are colored gray, while HPD intervals that exclude zero are colored black. 

The graph is probably better explained by, well, plotting it using the toy example below =P

Thanks to Robert O'Hara for suggesting and providing the original code for this function.
}


\value{
	If SSVS was applied individually to each coefficient of \code{X} when fitting the boral model, then the posterior probabilities of including the specified covariate are printed out i.e., \cr those from \code{x$ssvs.indcoefs.mean}.
	}

\seealso{
\code{caterplot} from the \code{mcmcplots} package for other, sexier caterpillar plots. 
}


\author{
Francis K.C. Hui \email{fhui28@gmail.com}; Wade Blanchard \email{wade.blanchard@anu.edu.au}
}


\examples{
\dontrun{
## NOTE: The values below MUST NOT be used in a real application;
## they are only used here to make the examples run quick!!!
example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
     n.thin = 1)

library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
n <- nrow(y)
p <- ncol(y)

X <- scale(spider$x)
spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
	num.lv = 2, mcmc.control = example_mcmc_control)


## Do separate line plots for all the coefficients of X
par(mfrow=c(2,3), mar = c(5,6,1,1))
sapply(colnames(spiderfit_nb$X), coefsplot, 
	spiderfit_nb)
}

}
\name{ranefsplot}

\alias{ranefsplot}

\docType{package}

\title{Caterpillar plots of response-specific random effects from a fitted model} 

\description{
\Sexpr[results=rd, stage=render]{lifecycle::badge("stable")}

Constructs horizontal line plot (point estimate and HPD intervals), otherwise known as "caterpillar plots", for the response-specific random intercept predictions in the fitted model.}

\usage{
ranefsplot(sel.spp, object, ordered = FALSE, est = "median", ...) 
	}
	
\arguments{
     \item{sel.spp}{A vector selecting which response' random intercept predictions are to be plotted. It can either be a numeric vector indexing the columns of \code{object$y} to be plotted, or a character vector with each being an element in \code{colnames(object$y)} indexing the names of the responses of the plotted. Default is \code{NULL}, in which case plots are done for all responses, one response at a time.}

     \item{object}{An object for class "boral".}
               
     \item{ordered}{If set to \code{TRUE}, then the random intercept predictions in each caterpillar plot are plotted from smallest to largest. Defaults to \code{FALSE}, in which case the caterpillar plot is simply ordered as per the rows of object{object$ranef.ids}.}

     \item{est}{A choice of either the posterior median (\code{est = "median"}) or posterior mean (\code{est = "mean"}), which are then used as the point estimates in the lines. Default is posterior median.}

     \item{...}{Additional graphical options to be included in. These include values for \cr \code{cex, cex.lab, cex.axis, cex.main, lwd}, and so on.}
}

\details{
For each response (column of the response matrix) and random intercept, the horizontal line or "caterpillar" is constructed by first marking the point prediction (posterior mean or median) with an "x" symbol. Then the line is construed based on the lower and upper limits of the highest posterior density (HPD) intervals as found in \code{object$hpdintervals}. By default, these are 95\% HPD intervals. To complete the plot, a vertical dotted line is drawn to denote the zero value. All HPD intervals that include zero are colored gray, while HPD intervals that exclude zero are colored black. 

By defaults, the plots are constructed one response at a time. That is, for a particular response, caterpillar plots of all the random intercepts (noting the number of plots is detemined by the number of random intercepts included in the model i.e., \code{ncol(object$ranef.ids}) are constructed on a single page. Then it moves onto the next response, and so on. If the user is only interested in plots from a subset of responses, then please make use of the \code{sel.spp} argument. This may be recommended especially since, in community ecology for example, the number of responses may be very large and so plotting all graphs may take a lot time.

The graph is probably better explained by, well, plotting it using the toy example below! 
}


\seealso{
\code{\link{coefsplot}} for horizontal line or "caterpillar plot" of the regression coefficients corresponding to the covariate matrix (if applicable),
the help file for \code{ranef} function in the \code{lme4} package, for other examples of caterpillar plots of random effect predictions,
\code{caterplot} from the \code{mcmcplots} package, as well as the \code{ggpubr} package, for other, sexier caterpillar plots. 
}


\author{
\packageAuthor{boral}

Maintainer: \packageMaintainer{boral}
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

spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
    ranef.ids = data.frame(region = rep(1:7,each=4)), 
    mcmc.control = example_mcmc_control, model.name = testpath) 

ranefsplot(sel.spp = 1:5, object = spiderfit_nb) 

ranefsplot(sel.spp = 1:5, object = spiderfit_nb, ordered = TRUE) 

ranefsplot(sel.spp = c("Alopacce","Zoraspin"), object = spiderfit_nb, ordered = TRUE) 
}
}

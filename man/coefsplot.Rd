\name{coefsplot}

\alias{coefsplot}

\docType{package}

\title{Caterpillar plots of the regression coefficients from a fitted model} 

\description{
\Sexpr[results=rd, stage=render]{lifecycle::badge("stable")}

Constructs horizontal line plot (point estimate and HPD intervals), otherwise known as "caterpillar plots", for the response-specific regression coefficients corresponding to the covariate in the fitted model. If a fourth-corner model is fitted, then "caterpillar plots" can optionally be produced for all the fourth-corner regression coefficients.}

\usage{
coefsplot(covname, object, fourthcorner = FALSE, labely = NULL, est = "median", ...) 
	}
	
\arguments{
  \item{covname}{The name of one of the covariates in the fitted model. That is, it must be a character vector corresponding to one of the elements in \code{colnames(object$X.coefs.median)}. 
  
  If \code{fourthcorner = TRUE}, then this argument is ignored.}

  \item{object}{An object for class "boral".}
             
  \item{fourthcorner}{If set to \code{TRUE}, then a caterpillar plot of the fourth-corner regression coefficients, as given by \code{object$traits.coefs.median} say, are plotted instead, assuming a model involving traits is fitted. Defaults to \code{FALSE}, in which case the plot is of response-specific regression coefficients. If \code{fourthcorner = TRUE}, then both the \code{covname} and \code{labely} arguments are ignored.}

  \item{labely}{Controls the labels on the y-axis for the line plot. If it is not \code{NULL}, then it must be a vector either of length 1 or the same length as the number of columns in the response matrix in the fitted boral object. In the former, it is treated as the y-axis label. In the latter, it is used in place of the column names of the response matrix to label each line. Defaults to \code{NULL}, in which the each line in the plot is labeled according to the columns of the response matrix, or equivalently \code{rownames(object$X.coefs.median)}.
  
  If \code{fourthcorner = TRUE}, then this argument is ignored.}

  \item{est}{A choice of either the posterior median (\code{est = "median"}) or posterior mean (\code{est = "mean"}), which are then used as the point estimates in the lines. Default is posterior median.}

     \item{...}{Additional graphical options to be included in. These include values for \cr \code{cex, cex.lab, cex.axis, cex.main, lwd}, and so on.}
}

\details{
For each response (column of the response matrix, the horizontal line or "caterpillar" is constructed by first marking the point estimate (posterior mean or median) with an "x" symbol. Then the line is construed based on the lower and upper limits of the highest posterior density (HPD) intervals as found in \code{object$hpdintervals}. By default, these are 95\% HPD intervals. To complete the plot, a vertical dotted line is drawn to denote the zero value. All HPD intervals that include zero are colored gray, while HPD intervals that exclude zero are colored black. 

For plots of fourth-corner regression coefficients, the coefficients are labelled such that on the left vertical axis the names of the covariates includes in the model are given, while on the right vertical axis the names of traits included in the model are given. Please see the \code{\link{about.traits}} for more about fourth-corner models where traits are included to help explain differences in species environmental responses to covariates. 

The graph is probably better explained by, well, plotting it using the toy example below! Thanks to Robert O'Hara for suggesting and providing the original code for this function.
}


\value{
	If SSVS was applied individually to each coefficient of the covariate matrix when fitting the model, then the posterior probabilities of including the specified covariate are printed out i.e., \cr those from \code{object$ssvs.indcoefs.mean}.
	
	For fourth-corher models, if SSVS was applied individually to fourth-corner coefficients when fitting the model, then the posterior probabilities of including the specified coefficient are printed out i.e., \cr those from \code{object$ssvs.traitscoefs.mean}.
	}

\seealso{
\code{\link{ranefsplot}} for horizontal line or "caterpillar plot" of the response-specific random effects predictons (if applicable),
\code{caterplot} from the \code{mcmcplots} package, as well as the \code{ggpubr} package, for other, sexier caterpillar plots. 
}


\author{
\packageAuthor{boral}

Maintainer: \packageMaintainer{boral}
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
n <- nrow(y)
p <- ncol(y)

X <- scale(spider$x)
spiderfit_nb <- boral(y, X = X, family = "negative.binomial", 
    lv.control = list(num.lv = 2), mcmc.control = example_mcmc_control, 
    model.name = testpath)


## Do separate line plots for all the coefficients of X
par(mfrow=c(2,3), mar = c(5,6,1,1))
sapply(colnames(spiderfit_nb$X), coefsplot, 
    spiderfit_nb)
    
        
## Consider a model based on Example 5a in the main boral help file
## The model is fitted to count data, no site effects, two latent variables, 
## plus traits included to explain environmental responses
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
                    mcmc.control = example_mcmc_control, model.name = testpath,
                    save.model = TRUE)

summary(fit_traits)

par(mar = c(3,10,2,10))
coefsplot(object = fit_traits, fourthcorner = TRUE)
}

}

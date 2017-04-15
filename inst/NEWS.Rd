\name{boral}
\title{News for package 'boral'}
\encoding{UTF-8}

\section{CHANGES IN boral VERSION 1.3}{
This list below is written in terms of what is deemed to be most important to least important changes =P 

  \itemize{
     \item A \code{predict.boral} function is now available, which can be used for constructing point predictions and associated intervals, on a linear predictor scale, from a fitted boral model. The predictions can be made either conditionally or marginally on the latent variables and any random row effects in the model. 

     \item Prior distributions for all variance (standard deviation to be precise!) parameters in a boral model are now (also) controlled by the fourth element of prior.control\$type and prior.control\$hypparams. Previously, the priors for variance parameters were all set to uniform distributions with the maximum controlled by different values of prior.control\$hypparams. Now, they are controlled by the same parameters that govern the priors for all dispersion parameters (e.g., the variance in the normal responses, overdispersion parameters in negative binomial responses). 
     
     This change allows for more flexibility in constructing the prior distribuion for the variance parameters, in particular, it brings in weakly informative priors such as half-normal and half-cauchy distributions.

     \item The Geweke convergence diagnostic is now calculated and outputed automatically as part of every boral model fitted. In particular, a new output \code{fit\$geweke.diag} produces a list with the two elements. The first element is itself a list containing the Geweke convergence diagnostic (Z-scores) for all appropriate parameters in the model. The second element contains the proportion of Z-scores that whose corresponding p-value is less than 0.05. No adjustment is made for multiple comparison on the p-values. 
     
     As written in the help file for boral, if you do indeed want to use this diagnostic to formally check for convergence, then we recommend you conduct some adjustment e.g., using Holm's method, by doing something such as \code{gew.pvals <- 2*pnorm(abs(unlist(fit\$geweke.diag[[1]])), lower.tail = FALSE); p.adjust(gew.pvals, method = "holm")}.

     \item Both \code{calc.logLik.lv0} and \code{calc.marglogL}, and hence information criteria calculated using these functions downstream in \code{get.measures} and \code{get.more.measures}, now permit multiple random row effects. 

     \item An \code{offset} argument is now available in the main boral fitting function as well as all relevant functions upstream and and downstream. This permits users to enter a matrix of a-priori effects to be included in the linear predictor. The offset must of the same dimensions as the response matrix \code{y}. 

     \item The names of some of the parameters, when constructing the JAGS scripts using \code{makejagsboralmodel} and \code{makejagsboralnullmodel}, have been changed. The purpose of this is to offer more consistency and transparency with the output from boral. In particular, the following name changes have been implemented: 
     \itemize{
          \item \code{all.params} => \code{lv.coefs}, but keep in mind that in the boral output, \code{lv.coefs} also contains the column-specific intercepts and dispersion parameters if appropriate. 
          \item \code{lvs} becomes \code{lv}
          \item \code{X.params} becomes \code{X.coefs}
          \item \code{row.params} becomes \code{row.coefs}. Note though the functions for calculating conditional and marginal log-likelihoods, \code{create.life}, and \code{simulate}, all still have \code{row.params} as an input argument (to allow the flexibility of fixed and random row effects)
          \item \code{row.ranef.sigma} becomes \code{row.sigma} 
          \item \code{probindX} becomes \code{ssvs.indcoefs}
          \item \code{probGpX} becomes \code{ssvs.gpcoefs}
          \item \code{alpha} becomes \code{cutoffs}
          }
     
     \item When possible, the marginal log-likelihood evaluated at the posterior median is now also outputed as part of the \code{get.measures} function and hence for the main \code{boral} fitting function (when \code{calc.ics = TRUE}). 

     }

	
}

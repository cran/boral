\name{boral}
\title{News for package 'boral'}
\encoding{UTF-8}

\section{CHANGES IN boral VERSION 1.4}{
This list below is written in terms of what is deemed to be most important to least important changes =P 

  \itemize{
     \item A \code{get.mcmcsamples} function has been created to allow users to simply extract MCMC samples from fitted boral objects. Apart for simplifying access of MCMC samples for assessing diagnostics, it also benefits downstream functions such as \code{get.measures} and \code{get.more.measures}. Thanks to Patrick Ewing for this idea!

     \item Inverse correlation, otherwise known as precision matrices, and significant precision matrices, are now available as part of the output for \code{get.residual.cor}.

     \item A \code{calc.varpart} function has been created, which attempts a variance partitioning on the linear predictor scale, averaged with respect to MCMC samples from the posterior distribution. Please see the associated help file for details.

     \item The priors on the elements of latent variable coefficient matrix (lv.coefs) which have sign contraints is now broadened to behave like with the choice of prior.control\$type[2]. Specifically, it is set to normal then a half normal distribution is used, if it is set to uniform[-c,c] then a uniform[0,c] is used, and if it is set to Cauchy then a half Cauchy prior is used. 
     
     \item Default \code{hypparams} have changed: For normal priors, the default variance is now set at 10, while for uniform priors, the default maximum value is now set to 30. The former has been done to echo the ideas from Gelman et al. (2008) in that that we do not expect regression coefficients on the linear predictor scale to go beyond -5 to 5; see also \url{https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations}.

     \item A citation file has been added. 

     \item Fixed a couple of bugs in the \code{predict.boral} example and code itself. Thanks to Chongliang Zhang for picking this up!
     }

}

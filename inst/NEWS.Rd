\name{boral}
\title{News for package 'boral'}
\encoding{UTF-8}

\section{CHANGES IN boral VERSION 1.8}{
This list below is written in terms of what is deemed to be most important to least important changes.

  \itemize{
     \item Fixed a change in 1.8.1 that was required as part of the latest developmental version of R
  
     \item Both \code{get.enviro.cor} and \code{get.residual.cor} functions now return highest posterior density intervals for each element in the relevant correlation and precision matrix. The probability defining the interval is given by a new prob argument.

     \item The \code{predict.boral} function now has a \code{scale} argument which allows prediction on the response scale (although note actual predicted responses are not simulated). Previously, and by default, predictions are done the linear predictor scale.

     \item The \code{predict.boral} function now allows marginal predictions extrapolate to new X and row IDs. This works as you would expect (hopefully!), by sampling new row effects from relevant normal random effects distribution.
     
     \item The get.residual.cor function now makes use of the \code{cor2pcor} function from \code{corpcor} package to calculate the inverse correlation or precision matrix. 

     \item A \code{tidyboral} function has been created to reformat output for boral so that instead of a separate element for mean/median/IQR/sd, a long data frame is outputed. This may be useful for researchers wanting to wrangle/ggplot some of the estimated parameters from a fitted boral model 

     \item For analyzing presence-absence responses, the parameterization of the LVM has been changed to use a direct probit link function rather than the step function with normal auxilary random variable parameterization. Hopefully this better handles the issue of calculating residual correlations correctly.

     \item Fixed an identifiablity constraint issue when row effects are fixed i.e., the first element in each fixed row effect is now constrained to be zero

     \item Fixed an issue found by co-maintainder Wade Blanchard to check that the dimensions of the trait matrix supplied to boral is consistant with other elements supplied.
     }

}

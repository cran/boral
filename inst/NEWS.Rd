\name{boral}
\title{News for package 'boral'}
\encoding{UTF-8}

\section{CHANGES IN boral VERSION 1.9}{
This list below is written in terms of what is deemed to be most important to least important changes.

  \itemize{
     \item Zero truncated Poisson and zero truncated negative binomial distributions are now available, and fingers crossed they work! They are accessed by setting \code{family = "ztpoisson"} and \code{family = "ztnegative.binomial"}. Please note that, as is commonly used in other regression models, the models are set up such that a log-link function connects the mean of the \emph{untruncated} distribution to the linear predictor.

     \item Creator, author, and contributor have been more formally defined now. Furthermore. lifecycle stickers have also been included for each function, courtesy of the \code{lifecycle} package although, it is all a bit ``opinionated" in the way it is applied. This are most clearly visible when using \code{boral} in RStudio.
       
     \item The function \code{fitted.boral} now has an additional argument \code{linear.predictor}. When set to \code{TRUE}, it allows one to return the fitted values on the linear predictor scale. When set to \code{FALSE}, which is the defaul behavior, the fitted mean values are on the response scale. Note things are slightly more complicated for the two newly implemented zero truncated distributions because, as mentioned above, the log-link function connects the mean of the \emph{untruncated} distribution to the linear predictor. Therefore if \code{linear.predictor = TRUE}, then the linear predictor is returned. But if \code{linear.predictor = FALSE}, then actual fitted mean value is returned. 

     \item The default for value for \code{mcmc.control$seed} has been changed from \code{123} to \code{NULL}. Thus the user must now specify a seed if they wish to "seed" the MCMC fitting algorithm. 

     \item An issue pertaining to the calculation of the Dunn-Smyth residuals for Tweedie responses has been fixed. Thaks to an ongoing collobration with Adrien Chevallier to leading to this issue being rectified.
     }

}

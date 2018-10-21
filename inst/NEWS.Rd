\name{boral}
\title{News for package 'boral'}
\encoding{UTF-8}

\section{CHANGES IN boral VERSION 1.7}{
This list below is written in terms of what is deemed to be most important to least important changes =P 

  \itemize{
	\item The syntax for how latent variables are included in the model has been revamped, with a new lv.control list argument controlling this. This was done as now structured latent variables are available, such that latent variables are allowed to be, say, spatially correlated based on some pre-defined distance matrix. Please see the about.lvs help file for more information.
  
     \item A X.ind argument is now implement in the main boral function, which allows users to manually force particular covariates to be excluded for particular species. The argument takes the form of a matrix with ones and zeros, where a zero indicates that a particular covariate is excluded for a particular species. Defaults to NULL, in which case all covariates are included for all species.

     \item The predict.boral function now has the option of returning all the full array of linear predictors, which the user can then manipulate to their heart's content!

     }

}

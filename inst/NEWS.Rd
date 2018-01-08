\name{boral}
\title{News for package 'boral'}
\encoding{UTF-8}

\section{CHANGES IN boral VERSION 1.5}{
This list below is written in terms of what is deemed to be most important to least important changes =P 

  \itemize{
     \item Stochastic Search Variable Selection (SSVS) is now available for the trait coefficients, by supplying a \code{prior.control\$ssvs.traitsindex} argument. This argument should have the same set up as \code{which.traits}, except with elements equal to either -1 (no SSVS) or 0 (SSVS). Please see the help files for the main \code{boral} function, \code{about.ssvs}, and \code{about.traits} for more information.

     \item The \code{create.life} and \code{simulate.boral} functions both now have the option of generating new latent variables for simulation, in addition to the previously available (default) option of supplying or conditioning the latent variables from the true model. This is controlled via the \code{lv.control} and \code{new.lvs argument} in the two functions respectively. Please see the help file for the \code{create.life} function for more information.

     \item \code{calc.ics = FALSE} is now the default argument in the main \code{boral} function.

     \item Both the \code{make.jagsboralmodel} and \code{make.jagsboralnullmodel} functions now print the \code{boral} version and the time of file creation as part of JAGS script.

     \item Some of the objects inside the function and in the examples in the help files have been altered to use ``_" instead ".". This has been done to reduce potential confusion with functions. Note however that arguments to and outputs from functions all still use ``." when required.

     \item Fixed a bug in the \code{plot.boral} code. Thanks to Matthew Helmus for picking this up!
     }

}

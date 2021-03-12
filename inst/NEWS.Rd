\name{boral}
\title{News for package 'boral'}
\encoding{UTF-8}

\section{CHANGES IN boral VERSION 2.0}{
This list below is written in terms of what is deemed to be most important to least important changes.

  \itemize{
     \item Response-specific random intercepts can now be included in a model. If you are familiar with including random row effects in a model, then one may think of response-specific random intercepts as similar to this, except the row effects are now response-specific specific rather than a common across all responses e.g., species-specific random intercepts to accounted for sampling design on a per-species basis. This is done via the \code{ranef.ids} argument in the main function \code{boral} with downstream uses elsewhere also. Please see the help file for \code{boral}, as well as the new \code{about.ranefs} help file for more information.

     \item Quite a bit of editing has been done to the help files, as well as some very minor names changes to returned objects. For help files, note that the mathematical notation for latent variables has been changed from \eqn{z_i} to \eqn{u_i}, while species-specifc random intercept "covariates" are now represented by \eqn{z_i}.

     \item The \code{coefsplot} function now implements plotting for fourth-corner coefficients, if a fourth corner model is fitted and \code{fourthcorner = TRUE} is set. Please see the help file for \code{coefsplot} for more information.

     \item A new argument \code{formula.X} is now available in the main \code{boral} function. Like other regression functions such as \code{lm} and \code{glm}, the argument \code{formula.X} should be an object of class "formula", which represents a symbolic description of the covariate matrix one wishes to include.  Please see the help file for \code{boral} for more details. 
       
     \item ALL functions to calculate information criteria and log-likelihood values of any sort are now classed as defuncted and no longer maintained...and probably do not work properly!
     }

}

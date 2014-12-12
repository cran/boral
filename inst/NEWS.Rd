\name{boral}
\title{News for package 'boral'}
\encoding{UTF-8}

\section{CHANGES IN boral VERSION 0.6 (2014-12-12)}{
This list below is written in terms of what is deemed to be most important to least important changes =D
  \itemize{
  	\item The \code{row.eff} argument in the \code{boral} fitting function now allows for random effect
  	e.g., in multivariate abundance data, if you want to assume the sites as a random sample. This is done
  	by setting \code{row.eff = "random"} in which case the row effects are assumed to be drawn from a normal
  	distribution with the mean and variance to be estimated. In particular, priors are placed on both 
  	the mean and variance, with the first element of \code{hypparams} partly controlling this. 
  	
  	Fixed row effects can still be assumed if \code{row.eff = "fixed"}. Defaults to \code{row.eff = "none"},
  	in which case no row effects are assumed. 
  	
  	This update should still be compatible with version 0.5. That is, \code{row.eff = FALSE} will set 
  	\code{row.eff = "none"}, and \code{row.eff = TRUE} will set \code{row.eff = "fixed"}.
  	
  	\item \code{get.residual.cor} now calculates the residual correlation matrix differently.
  	Whereas in version 0.5 it was done based on the linear predictor \eqn{\bm{z}'_i\bm{\theta}_j},
  	is is now corrected based on the \eqn{p \times p} matrix \eqn{\bm{\Theta}\bm{\Theta}'} where
  	\eqn{\bm{\Theta}} is the matrix formed by stacking the rows of the latent variable coefficients
  	\eqn{\bm{\theta}_j} together. 

	\item The negative binomial family has been slightly reparameterized, but the related arguments
	are still the same i.e., the mean-variance relationship is still \eqn{Var(y) = \mu + \phi \mu^2}
	and where the prior for \eqn{\phi} is controlled by the fourth element in \code{hypparams}. Feedback
	from users previously found the parameterization in version 0.5 unstable, and often were focued to
	select the poisson family for many of the columns. Hopefully this new parametrization is more reliable!
	
  	\item A \code{do.fit} argument is now available in the \code{boral} fitting function. 
  	If set to \code{FALSE}, then only the model file is created and saved in the current
  	working directly i.e., no MCMC sampling via JAGS is performed. This may be useful if 
  	users want to simply create the JAGS file and then alter it and run it themselves. 
  	Defaults to \code{TRUE}.
  	
  	\item A \code{model.name} argument has been added to the \code{boral} fitting function, as
  	well as \code{make.jagsboralmodel} and \code{make.jagsboralnullmodel}. This argument allows
  	the user to given a custom name to the text file where the JAGS model is to be written to.
  	Defaults to \code{NULL}, in which case the standard name of "jagsboralmodel.txt" is adopted.
  	Of course, the model file is only saved in \code{boral} if \code{save.model} is set to \code{TRUE}.
  	
  	\item \code{get.enviro.cor} and \code{get.residual.cor} have been rewritten so that they now 
  	take an fitted boral model directly, instead of requiring the fitted MCMC samples and 
  	several other arguments. 
  	
  	\item Estimation of the Compound Laplace-Metrpolis estimator in \code{get.measures}
  	has been corrected.
  	
  	\item \code{n.burnin, n.thin, n.iterations} now returned as part of the output list 
  	from fitting a model using \code{boral}.
  	
  	\item Some minor errata in the help files have now been corrected. This mainly concern
  	the extraction of MCMC samples from the boral model. In version 0.5, this was done using
  	\code{as.mcmc} function. However this was found to be unreliable in the newer version of
  	the \code{coda} package, and so alterative method of extracting the MCMC samples is used.
  	Please see the updated help file for \code{calc.condlogLik} for instance.
   	
  	\item A default of \code{trial.size = 1} has now been set to all applicable 
  	functions. 
  	
   	\item "theta0", which was previously used to denote the column-specific intercepts \eqn{\theta_{0j}},
   	has now been changed to "beta0" and \eqn{\beta_{0j}}, to make it consistent with standard GLM and
   	regression notation.
	}
}

%
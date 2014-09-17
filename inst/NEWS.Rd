\name{boral}
\title{News for package 'boral'}
\encoding{UTF-8}

\section{CHANGES IN boral VERSION 0.5 (2014-09-17)}{
  \itemize{
  	\item The \code{site.eff} argument in the \code{boral} function has now been 
  	altered to \code{row.eff}, with default \code{row.eff = FALSE}.
  	
  	\item The default arguments for \code{n.burnin} and \code{n.iteration} in the 
  	\code{bora} function have modified to \code{n.burnin = 5000} and 
  	\code{n.iteration = 35000} i.e., both higher than in version 0.4
  	
  	\item \code{get.hpd.intervals} now takes a \code{prob} argument, to allow
  	control of percent coverage of HPD intervals. Defaults to 95\%
  	
  	\item New options for \code{family} argument include ``exponential" (with
  	log link function), ``gamma" (with log link function) such as the the dispersion
  	parameter returned in \code{lv.coefs.mean/median/iqr/sd} corresponds to the rate
  	parameter (please see \code{?rgamma} for details on parameterization), 
  	and ``beta" (with logit link function) for regression of data inside the interval 
  	but not including 0 and 1
  	
  	\item A new function \code{get.enviro.cor} is available to calculate the correlation
  	between the columns of the response matrix \code{y} which can be attributed to 
	similarities in the response to the explanatory variables in \code{X}
	
	\item The \code{get.more.measures} function can now be run within the \code{get.measures}
	function by setting \code{more.measure = TRUE}. Help files still exists separately for
	both functions
	
	\item The \code{X.eff} argument for \code{make.jagsboralmodel} and 
	\code{make.jagsboralnullmodel} has been replaced by with \code{num.X} argument, which 
	requires the user to input the number of columns in the model matrix \code{X}. Defaults
	to \code{num.X = 0} in which case it assumes no covariates are included in the model
	i.e., corresponds to \code{X.eff = FALSE}
	
	\item The \code{hypparams} argument has been extended to take four elements: 1) variance 
	for the normal prior of all intercepts, row effects, and ordinal cutoff points, 
	2) variance for the normal prior of all coefficients corresponding to the latent variables
	(required by ignored if \code{num.lv = 0}), 3) variance for the normal prior of all 
	coefficients corresponding the model matrix \code{X}, 4) upper limit for uniform prior 
	of the dispersion parameters
	
	\item Stochastic search variable selection (SSVS) now implemented, as alternative to 
	using information criteria to perform variable selection, and also as a method of
	``hypothesis testing". Please see the help file for \code{boral} for details. In brief,
	SSVS can return the posterior probability of a coefficient/covariate/group of covariates 
	belonging in a model. 
	}
}

%
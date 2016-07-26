\name{boral}
\title{News for package 'boral'}
\encoding{UTF-8}

\section{CHANGES IN boral VERSION 1.0 and 1.0.1}{
This list below is written in terms of what is deemed to be most important to least important changes =D
  \itemize{
	\item In the main \code{boral} function, the \code{trait} argument should now be a matrix which does NOT include an intercept column. This makes it now consistent with the \code{X} argument. An if statement has been created within the main function to check this.

	\item In the main \code{boral} function, the arguments \code{hypparams} and \code{ssvs.index} are now nested under a new argument called \code{prior.control}, which is a list argument. This argument also includes an element called \code{type}, which allows for different prior distributions, as discussed in the next point.

	\item A new function \code{coefsplot} has been created, which can be used to construct horizontal line or "caterpillar plots" for column-specific regression coefficients corresponding to covariates \code{X} (if applicable). These plots consist of a symbol marking a point estimate (posterior median by default), and a line based on the limits of the HPD intervals. 

	\item Both the \code{lvsplot} and \code{plot.boral} functions no longer have an argument for \code{a}. Instead, all graphical options are now controlled externally and/or through the \code{...} argument.

	\item As requested by several users, the range of prior distributions has been expanded, slightly. Specifically, there is now a \code{type} argument nested, within the \code{prior.control} list argument, that requires a vector of four elements. These four elements are closely related to the four respective elements in \code{prior.control$hypparams}. By default, \code{prior.control$type = c("normal","normal","normal","uniform")}, which coincides with the priors used in all version of boral up until now. 	
	The first three elements of \code{prior.control$type} can also take string values \code{"uniform"}, for a uniform prior, and \code{"cauchy"} for a Cauchy prior. 
	The fourth element of \code{prior.control$type} can now also take string values \code{"halfnormal"}, for a half-normal prior, and \code{"halfnormal"} for a half-Cauchy prior. Please see the help file for the \code{boral} function for more detail on these priors.
	
	\item In addition to the changes above, there now a \code{ssvs.g} argument that is also available with \code{prior.control} list argument. This is used to control the shrinkage factor, i.e. strength of the spike, in Stochastic Search Variable Selection (SSVS). Please see the SSVS section in the help file for the \code{boral} function for more detail. Defaults to \code{prior.control$ssvs.g = 1e-6}.

	\item In the main \code{boral} function, the arguments \code{n.burnin, n.iteration, n.thin} and \code{seed} are now all nested under a new argument called \code{mcmc.control}, which is a list argument (not to be confused with the \code{prior.control} argument!). The default values for these four arguments remain the same as in the previous version of boral.
	
	\item Dunn-Smyth residuals are now available for ordinal responses as well. This is now produced in addition to the single table of agreement ("confusion matrix") in the \code{ds.residuals} function. Fitted values are also produced for ordinal responses (using the \code{fitted} function); these fitted values are simply defined as the level with the highest fitted probability. 
	
	\item Both the above modifications to Dunn-Smyth residuals and fitted values have a flow-on effect to the \code{plot} function, which now also performs residual analysis for ordinal responses. 

	\item The function \code{create.life} now has an \code{save.params} argument, which allows the user to (additionally) return the parameters of the "true model" used for simulating the multivariate response matrix. This argument is useful mainly when \code{traits} and \code{traits.coefs} are supplied, as \code{lv.coefs[,1]} and \code{X.coefs} are overwritten in such a scenario; please see the help file for \code{create.life} for more information. When \code{save.params = TRUE}, then a list is returned containing the simulated response matrix as well as the parameters of the true model. Defaults to FALSE.

	\item A \code{num.lv} argument is now required for \code{make.jagsboralmodel}. Note the function is typically not used directly by the user. Defaults to \code{num.lv = 2}. 
	
	\item Fixed bug in \code{get.enviro.cor} and \code{get.residual.cor} functions. Thanks to Marcio Leite for picking this up!
	
	}

	
}

%
\name{boral}
\title{News for package 'boral'}
\encoding{UTF-8}

\section{CHANGES IN boral VERSION 1.0}{
This list below is written in terms of what is deemed to be most important to least important changes =D
  \itemize{
	\item Unit testing added to the code, which hopefully will continue with future updates.

 	\item All normal random effects distributions, e.g. for the row effects if \code{row.eff = "random"} and
 	if traits are included in the model (please see the Section on \code{Including species traits} in
 	the \code{boral} help file)), are now parameterized in terms of their standard deviation instead of their
 	variance. Importantly, the standard deviation is now assigned a uniform prior distribution. This modification
 	is on the back of Gelman (2006) and other papers since, which have recommended placing a uniform prior on the
 	standard deviations rather than variances. The output for parameter estimates from boral in now also in terms 
 	of standard deviation.
 	
 	\item On a similar note to the above, for any columns of \code{y} that are assigned a normal or lognormal distributions,
 	they are now parameterized in terms of their standard deviations, with a uniform prior \code{dunif(0,hypparams[4])} now
 	assigned to the standard deviation for normal responses and standard deviation on the log scale for lognormal responses.
 		
 	\item Both \code{get.enviro.cor} and \code{get.residual.cor} now return covariance (\code{cov}) and 
 	correlation (\code{cor}) matrices as part of a list output. In addition, both functions will also print 
 	out a \code{sig.cor} matrix containing only the "signficant" correlations, that is, correlations whose 
 	highest posterior density interval does not include zero. All non-significant correlations are set 
 	to zero in this matrix. The target coverage of the intervals is determined by a \code{prob} argument 
 	in both functions, similar to how the \code{get.hpdintervals} function works, with a default of 0.95 or 95\%.
 
 	\item In \code{lvsplot}, the \code{new.plot} argument has been removed and replaced with a \code{est} argument,
 	so that the plots are constructed based on either the posterior mean or median. This is analogous to 
 	other functions like \code{get.residual.cor} that also accept an \code{est} argument. Please note this implies
 	that only ONE ordination plot is now returned instead of two.
 	
 	\item In \code{lvsplot}, there is now an \code{alpha} parameter that gives the user more ability to adjust and 
 	"tweak" the relative scaling of the latent variables and their coefficients, when constructing a biplot. 
 	The notion of an \code{alpha} is akin to the different type of scalings available for correspondence analysis, and 
 	the notion of symmetric and asymmetric plots. We generally prefer to have alpha around 0.5, the default 
 	value, which is designed to place the latent variables and their coefficients on the same scale...in principle. 
 
 	\item In \code{lvsplot}, there is now a \code{main} argument to allow a custom title for the ordination plot. 
 	The default to \code{NULL}, in which case a "standard" title in used. 
 	
 	\item The default MCMC sampling settings in \code{boral} have been increased to 
 	\code{n.burnin = 10000, n.iteration = 40000, n.thin = 30}.
 	
 	\item The normal quantile plot in \code{plot.boral} now includes a 1-to-1 line for better assessment
 	of the normality of residuals.
			
	}

	
}

%
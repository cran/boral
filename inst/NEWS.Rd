\name{boral}
\title{News for package 'boral'}
\encoding{UTF-8}

\section{CHANGES IN boral VERSION 0.7 (2015-07-08)}{
This list below is written in terms of what is deemed to be most important to least important changes =D
  \itemize{
	\item \code{lvsplot} now has a \code{biplot} argument, such if \code{biplot = T}, then both of the latent
	variables AND their coefficients are presented on a single plot. Such biplots are common in multivariate analysis
	in ecology, for instance, to visualize both patterns between sites and the indicator species associated 
	with these sites. In addition, there is \code{ind.spp} argument in \code{lvsplot}, that allows the user to 
	include only the \code{ind.spp} "most important" latent variable coefficients (as determined by the L2-norm).
	This can be used to avoid overly messy plots and just focus on the "indicator" columns (species). 
	Please note that a necessary scaling is performed in order to construct the biplot, which is somewhat 
	experimental in that it make not always get it "right" =P
	
	\item Tougher priors are now used on for the latent variable coefficients. In particular, the \code{hypparams} 
  	argument in the \code{boral, make.jagsboralmodel} and \code{make.jagsboralnullmodel} functions is now set to
  	\code{hypparams = c(100,20,100,100)}. This implies that a normal prior with mean zero and variance 20 
  	instead of 100 is used for the latent variable coefficients \eqn{\bm{\theta}_j}. We have found that 
  	a slightly more informative prior tends to stablize the MCMC sampling and produce more believable results
  	for the ordinations and residual correlations.

	\item As well as returning a residual correlation matrix, \code{get.residual.cor} also now returns the 
	mean/median point estimator of the trace of residual covariance matrix. This number may be interpreted as 
	the amount of covariation explained by the latent variables. Where this comes in handy is if you are comparing 
	a pure LVM versus a model with latent variables and predictor variables. The proportional difference in 
	the trace between these two models is *somewhat* analogous to the amount of covariation explained by the predictors. 
  	
  	\item If \code{row.eff = "random"}, row effects are now sampled from a normal distribution with MEAN ZERO. 
  	In boral 0.6, it was sampled from a normal distribution with non-zero mean. This correction ensures that any
  	differences in row total abundance are "transferred" over to the column-specific coefficients. It is also makes
  	it consistent with the notion of a row effect being considered a random intercept. 
  	
  	\item Several users have noted that the negative binomial family tends to have issues sampling the overdispersion 
  	parameter. While this is dataset specific, and it is hard to offer a universal solution to this, a warning to 
  	now printed out explicity informing the user that such a problem has arisen (when it occurs of course), and that
  	they should consider enforcing a stricter prior on the overdispersion parameter or switch to a Poisson family 
  	for some of the responses.

  	\item Corrected issues for running JAGS when the tweedie family is used in combination with other families.
  	This previously straight out failed in JAGS, and should now be rectified.

	}
}

%
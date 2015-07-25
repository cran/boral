\name{boral}
\title{News for package 'boral'}
\encoding{UTF-8}

\section{CHANGES IN boral VERSION 0.8 (2015-07-28)}{
This list below is written in terms of what is deemed to be most important to least important changes =D
  \itemize{
	\item There was a MAJOR error for in the specification of the normal distribution when 
	\code{ran.eff = "random"} is version 0.7. This has now been corrected!
	
	\item The probit link function is now used with the \code{binomial} and \code{ordinal} families. 
	Previously, the logit link function was applied. While the logit function is the default choice
	in regression analysis, given it is the canonical link for the binomial family and its convenient interpretation
	as log odds, things have been switched over to the probit link because: 1) it offers the advantage that 
	it can be readily intepreted in terms of having a auxilary ("latent") variable having a standard 
	normal distribution, with cutoffs dividing up the distribution into categories. For instance with Bernoulli
	responses, if the auxilary variable is greater than zero than this corresponds to a presence, otherwise it
	corresponds to an absence. Regression coefficients can thus be interpreted on this auxilary variable scale,
	2) Related to point 1 is the fact that by regarding the probit link as having a standard normal auxilary 
	variable, this can be exploit in MCMC sampling, and often tends to produce more stable sampling. 

	\item The \code{negative binomial} family has been reparameterized such that the variance is now written 
	as \eqn{Var(y) = \mu + \mu^2/\phi}, where \eqn{\phi} is now the "new" dispersion parameter. I recognize this
	will annoy a lot of people (ecologists?) who are used to the \eqn{\phi\mu^2} way of thinking, and 
	I do apologize for this. The main reason for the change is: 1) it is now in line with the \code{nbinom} family
	in base \code{R} as well as in JAGS, which parameterizes things as \eqn{Var(y) = \mu + \mu^2/(size)}, 2) the 
	sampling algorithm in JAGS seems to work more stably with this parameterization. 
	
	\item Species traits can now be included in \code{boral}. Details can be found in the help file, but in summary, 
	both the column-specific intercepts \eqn{\beta_{0j}} and regression coefficients \eqn{\bm{\beta}_j} may be
	regressed against species traits in a linear model type manner. This makes the intercepts and coefficients
	random effects. Consequently, this option may be considered as a form of dimension reduction -- rather than 
	separate species coefficients, these coefficients are now linearly related to the value of the species traits 
	themselves. This also has links to the so-called fourth corner problem in multivariate ecology. 
	
	\item In the help file for \code{boral}, there is now a note about why only one MCMC chain is sampled. In brief,
	this is due to the sign switching problem -- different chains can have latent variables and their
	coefficients that are similar in magnitude but opposite in sign. Consequently, checking convergence and combining chains
	to produce do inference etc...becomes inappropriate. I recognize this is a cop out, and that there are ways to
	deal with this through post-processing, but these are complicated to implement =P

	\item For reasons of computation and increased customizability, both the functions \code{make.jagsboralmodel} 
	and \code{make.jagsboralnullmodel} now incoporate an argument for \code{trial.size}, which is specified 
	in the same manner as in \code{boral}. The default value is \code{trial.size = 1}, as is the case with \code{boral}.

	\item Due to its instability (often produces NA values), the compound Laplace-Metropolis ICs has been 
	as a potential information criterion available from \code{get.measures}.
	
	\item Both the plot functions, \code{lvsplot} and \code{plot.boral}, now reset the \code{par} arguments 
	after constructing the plots
	
	\item If all columns of \code{y} are assumed to be Bernoulli responses, i.e. \code{family = "binomial", trial.size = 1}, then, given the use of the probit link, the model is now parameterized directly in terms of auxilary variable having standard normal distributions, with success \eqn{y=1} corresponding to if this auxilary variable is greater than zero. This parameterization tends rectifies a previous issue regarding the scale of the latent variables not being identifiable when all the responses are Bernoulli. 

	\item To ease computation slightly, dispersion parameters are now no longer sampled in JAGS when they 
	don't need to be, i.e. if all columns of \code{y} are binomial, exponential, multinomial, ordinal, and so on 
	}

	
}

%
\name{boral}
\title{News for package 'boral'}
\encoding{UTF-8}

\section{CHANGES IN boral VERSION 1.2}{
This list below is written in terms of what is deemed to be most important to least important changes =P Thanks to Wade Blanchard from the ANU for suggesting the first two changes.

  \itemize{
	\item{The output for \code{hpdintervals} has been reformatted to simplify the structure of what is returned. For example, consider the lower and upper limits of the credible intervals for the latent variables: previously the function returned was \code{lv.lower} and \code{lv.upper} as two separate elements. In the new version, the function returns a single element called \code{lv}, which is a three dimensional array with \code{lv[,,"lower"]} and \code{lv[,,"upper"]} for the limits. Overall, the new structure should be fairly straightforward to understand.}

	\item{The help file for the main \code{boral} function has been formatted by moving some of the sections into their own separate help files. In particular, 1) stochastic search variable selection (SSVS) now has its own help file as \code{about.ssvs}; 2) Information on the distributions (families) available now has its own help file as \code{about.distributions}; 3) Information on including species traits to help mediate environmental responses now has its own help file as \code{about.traits}.
	
	Please note that these new help files are NOT directly linked to any, newly created functions. 
	}

	\item{In \code{plot.boral}, the 0--1 line in the normal quantile plot of Dunn-Smyth residuals has now been removed, and replaced with an estimated straight line based on the \code{qqline}. For the purposes of residual analysis, if the fitted model is appropriate then the Dunn-Smyth residuals should lie approximately on a straight line. Any deviations away from indicate skew, heavy tails etc... in the residuals, which in turn suggest outliers and/or the model assumptions are not satisfied. 
	
	Of course, the normal quantile plot is just one plot to check. Arguably the more important plot to check is the residuals versus fitted values!}

	\item{Traits are now allowed for ordinal regression, used to mediate environmental responses in the analogous manner to what is done with other response types; see \code{about.traits} help file. The only thing to keep in mind is that, for parameter identifiability reasons, if there is only one ordinal column, then the intercept for this is forced to be zero and therefore not regressed against traits.}

	\item{All examples in help files, even those enclosed with in a \code{dontrun}, now use a \code{example.mcmc.control} to speed things up just in case users want to run them and get an idea of how they work!}

	\item{Fixed a bug associated with the use of uniform distribution as priors.}
	
	\item{Fixed a bug associated with calculating residual correlations for binary response matrices (Thanks to Mirko Di Febbraro for finding this typo!)}
	}

	
}

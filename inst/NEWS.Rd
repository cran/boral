\name{boral}
\title{News for package 'boral'}
\encoding{UTF-8}

\section{CHANGES IN boral VERSION 1.1}{
This list below is written in terms of what is deemed to be most important to least important changes =P
  \itemize{
	\item{The way row effects can be incorporated has now been expanded slightly. Whereas all previous versions of \code{boral} only permitted for a single, unique site effect for each row, one can use the \code{row.ids} argument to specify site effects that cover multiple rows. For instance, if the first five rows in the dataset are replications from a site 1, the next five rows are replications from site 2 and so on, then one can specify \code{row.ids = matrix(c(1,1,1,1,1,2,2,2,2,2,...), ncol = 1)} to include a ``site effect". 
	
	Also, multiple row effects can now be included to reflect a nested or hierarchical structures in the rows of the data. This nesting is reflected in how the user specifies \code{row.ids}, in particular, by having more than one column in the matrix. For instance, if \code{y} is structured such that the first five rows are creeks nested within region 1, the second five rows are creeks nested within region 2 and so on, and one wanted to included a effect effect of creek nested within region, then one can specify \code{row.ids = cbind(c(1,2,3,4,5,6,7,8,...), c(1,1,1,1,1,2,2,2,2,2,...))}, such that the first column is the effect for creek and the second is the effect of region.
	
	One still has to specify whether to include the effects as fixed and random via the \code{row.eff} argument. Arguably though, all nested designs should be using random effects. Please note that at the moment, if multiple random row effects are included, then one can only set them all to be fixed effects or all to be random.}
	
	\item{Even though it is permitted as outlined above, the \code{calc.logLik.lv0} and \code{calc.marglogLik} currently does not permit the calculation of log-likelihoods when the row effects are random and are more complicated than a single, unique random effect for each row. Put simply, it is just too complicated to marginalize over random effects structures in general. Sorry! Please note the \code{get.measures} and \code{get.more.measures} functions for calculating information criteria have now also been modified to reflect this.}

	\item{The function \code{lvsplot} now has two additional arguments. The first is called \code{return.vals}. If \code{return.vals = TRUE}, then the \emph{scaled} latent variables scores and corresponding scaled coefficients are returned (based on the value of \code{alpha} used). This is useful is the user wants to construct their own custom model-based ordinations.}

	\item{The second new argument added to \code{lvsplot} is called \code{which.lvs}. This now allows model-based ordinations in cases where the fitted boral model contains more than two latent variables. In which case, \code{which.lvs} can be used to specify which latent variables (ordination axes to plot). By default, the first two axes are plotted i.e., \code{which.lvs = c(1,2)}. Please note this argument is ignored if the fitted boral model only contains one latent variables.}

	\item{If the number of ordinal responses exceeds 2, the species intercepts \eqn{\beta_{0j}} are now drawn as a random intercept from a normal distribution with mean zero and standard deviation given by \code{ordinal.ranef.sigma}. The random intercepts represent species-specific deviations from the common cutoffs i.e., the \eqn{\tau_k}'s. Previously, the intercepts were estimated as fixed effect but with a sum to zero constraint in these case of ordinal response. The change to a random intercept is in preparation for future updates when traits will be allowed for ordinal responses.}
	
	\item{A loophole has been closed! Specifically, the functions \code{get.measures} and \code{get.more.measures} now formally prevent calculation of marginal log-likelihood based information criteria in general cases where the random effects are complicated e.g., traits are included such that the regression coefficients \eqn{\beta_{0j}, \bm{\beta}_j} are random effects, when more than two columns are ordinal responses, or when more than one row effect is included (such that \code{row.ids} is more than just a one column matrix). The conditional log-likelihood based information criteria from \code{get.measures} are still calculated in such case though.}
	
	\item{Help files for main boral function, \code{make.jagsboralmodel} and \code{make.jagsboralnullmodel} have been corrected/update to specifically state the types of priors now available.}

	\item{The error message in the case where there are issues with MCMC sampling for the negative binomial overdispersion parameter has been updated, and is hopefully more informative now.}

	\item{The \code{more.measures} argument has now been scrapped from the \code{get.measures} function for computing some information criteria.}

	\item{Fixed a bug in the \code{calc.logLik.lv0} function.}

	\item{The \code{calc.condlogLik} function now accepts \code{lv = NULL} i.e., fitted boral models that do not have any latent variables included.}

	\item{For future clarity, the \code{X.mc} argument in \code{calc.marglogLik} function has been renamed \code{lv.mc}.}
	}

	
}

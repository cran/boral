##############
## Latent Variable Ordination and Regression using MCMC 
## Site effects fitted as fixed effects
## NB parameterized as V = mu + phi*mu^2
## Ordinal data handled as propotional odds regression, with same cutoff points for all spp. 
## Multinomial regression is available, but CURRENTLY NOT ALLOWED DUE TO COMPUTATION TIME
## This should make sense if the same ordinal scale is applied to all species

## Changes from v0.8
## 1) qqplot in plot.boral now adds a 1:1 line for better diagnosis
## 3) get.enviro.cor and ger.residual.cor now return covariance and correlation matrices
## 4) get.residual.cor and ger.enviro.cor now also print out correlation matrices that contains only significant correlations (HPDintervals not containing zero). As a result, there is a prob argument now included in both, with a default of 0.95
## 5) Default settings for n.burnin, n.iteration, and n.thin have been ramped up
## 6) In lvsplot, the new.plot argument has been removed and replaced with a est argument to plot either posterior mean or median. Also, a main argument has been included to allow custom titles
## 7) For lvsplot, include an alpha parameter that allows the user to tweak the scaling themselves. This is akin to the different type of scalings in CA. For biplots, we generally prefer adjust alpha until the LVs and coefs are the same scale. This is usually around the value of alpha = 0.5 to 0.55
## 8) After some complaints from users not the swtich to size for the negative.binomial family, it has now been switched back to V = mu + phi*mu^2. However,...
## 9) All random effects distributions are now parametrized in terms of their standard deviation, with a uniform prior now placed on the standard deviation sigma2 instead of on the variance. This is following advice by Gelman, 2006, Prior distributions for variance parameters in hierarchical models
## 10) Similarly, all normal and lognormal responses are now parameterized in terms of standard deviation instead of variance.

## Changes from v0.9 and v0.9.1
## 1) arguments for hypparams and ssvs.index are now nested under a prior.control "list" argument. This argument also includes an element called type, which allows for different prior distributions. Please see next point. 
## 2) uniform and cauchy priors now available for elements 1-3 of prior.control$type, and half-Cauchy and half-normal prior now available for dispersion parameters
## 3) A ssvs.g argument is available within the prior.control list argument, to control the strength of the spike is SSVS. Defaults to 1e-6
## 4) argument for n.burnin, n.iteration, n.thin and seed are now nested under a mcmc.control argument.
## 5) The trait argument should now be a matrix that does NOT include an intercept column. An if loop has been created to stop this
## 6) A hidden function called setup.resp.families.lv has been created, which is used as part of make.jagsboralmodel and make.jagsboralnullmodel. Basically, this simplifies the writing the JAGS script such that if all columns of y assume the same distribution, then a loop for(j in 1:p) is used instead of writing each column out individually, as is one if columns may take different responses. Note that as a result, there is now an num.lv argument used for make.jagsboralmodel.
## 7) create.life now has an additional argument save.params, which allows the user to (additionally) return the parameters of the "true model" used for simulating the multivariate response matrix. This argument is used mainly when traits.coefs is supplied, as lv.coefs[,1] and X.coefs are overwritten in such case.
## 8) The argument a is now removed from the lvsplot and plot.boral arguments, with all graphical options now controlled externally and/or through the ... argumnet
## 9) A coefsplot function has been created which makes caterpillar plots of the separate spp coefficients corresponding to X.
## 10) Calculate Dunn-Symth residuals for ordinal data, which in turn are used in residual analysis


## TODO: 
## 1) correct problems with increasing number of LVs causing increase in covariance and variance i.e., scale problem?; 
## 2) think about how to parallelize JAGS? 
## 3) Edit lvsplot function so that it actually returns scaled loadings and lvs as well as plotting them. Perhaps have an argument for plot = T or F
## 4) draw coefficients as random effects, or reduce rank them? HARD!!!
## 5) add a row id to allow a random effect based on that!
## 6) Allows traits for ordinal regression
## 7) allow for multiple chains, but don't check convergence on the LVs and their loadings (the code is close to be complete there); 
## 9) Extend lvsplot to include a which argument for which axis if lv > 2 (see Suggestions from BobOHara - email.txt and R code).



##############
# library(R2jags); 
# library(mvtnorm); 
# library(mvabund); 
# library(coda); 
# library(MASS)
# library(fishMod)
#source("/media/fhui/LEXAR/Maths/PHD statistics/Miscellaneous/Rpackage-boral/boral10/R/auxilaryfunctions.R")


# n = 60; p <- 30
# X <- matrix(rnorm(n*2),n,2); beta <- cbind(matrix(rnorm(p*3),p,3),runif(p,0,5)); true.power <- 1.6
# mu <- exp(cbind(1,X)%*%t(beta[,1:3]))
# y <- matrix(NA,n,p)
# for(j in 1:ncol(y)) { y[,j] <- rTweedie(nrow(y), mu = mu[,j], phi = beta[j,4], p = true.power) }
# family = "tweedie"
# num.lv = 0; row.eff = "none"; n.burnin = 10000; n.iteration = 40000; n.thin = 30; save.model = TRUE; seed = 1; calc.ics = TRUE; trial.size = NULL; hypparams = c(50,20,50,50); 
# 
boral <- function(y, ...) UseMethod("boral")


# y = simy; X = as.matrix(X[,-1]); traits = traits[,-1]; family = "binomial"; num.lv = 2; save.model = TRUE; calc.ics = FALSE; mcmc.control = list(n.burnin = 200, n.thin = 10, n.iteration = 1500); row.eff = "none"; trial.size = 1; prior.control = list(type = c("normal","normal","normal","uniform"), hypparams = c(100, 20, 100, 50), ssvs.index = -1, ssvs.g = 1e-6); do.fit = TRUE; model.name = NULL

## Model is g(mu_{ij}) = row + beta0 + LV_i*theta_j + X_i*beta_j
boral.default <- function (y, X = NULL, traits = NULL, which.traits = NULL, family, trial.size = 1, num.lv = 0, row.eff = "none", save.model = FALSE, calc.ics = TRUE, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123),
prior.control = list(type = c("normal","normal","normal","uniform"), hypparams = c(100, 20, 100, 50), ssvs.index = -1, ssvs.g = 1e-6), do.fit = TRUE, model.name = NULL, ...) {

	if(is.null(dim(y))) { 
		message("Converting y into a one column matrix."); y <- matrix(y, ncol = 1) }
	if(!is.null(X) & is.null(dim(X))) { 
		message("Converting X into a one column matrix."); X <- matrix(X, ncol = 1) }
	if(!is.null(traits) & is.null(dim(traits))) { 
		message("Converting traits into a one column matrix."); traits <- matrix(traits, ncol = 1) }
	if(!is.null(X)) { 
		if(!is.matrix(X)) X <- as.matrix(X)
		if(any(apply(X,2,function(x) all(x == 1)))) { stop("No intercept column should be included in X") } 
		}
	
	
	if(!("type" %in% names(prior.control))) prior.control$type <- c("normal","normal","normal","uniform")
	if(!("hypparams" %in% names(prior.control))) prior.control$hypparams <- c(100, 20, 100, 50)
	if(!("ssvs.index" %in% names(prior.control))) prior.control$ssvs.index <- -1		
	if(!("ssvs.g" %in% names(prior.control))) prior.control$ssvs.g <- 1e-6
	
	if(length(prior.control$hypparams) != 4 || length(prior.control$type) != 4) 
		stop("prior.control$type and prior.control$hypparams must be a vector of four elements. Please see boral help file as to what the elements correspond to.\n")
	if(!all(prior.control$type[-4] %in% c("normal","uniform","cauchy")))	
		stop("At least one of the first three elements of prior.control$type is not supported in current version of boral...sorry!")
	if(!(prior.control$type[4] %in% c("uniform","halfcauchy","halfnormal")))	
		stop("The fourth element of prior.control$type is not supported in current version of boral...sorry!")
		
 	if(!is.null(which.traits) & any(prior.control$ssvs.index > -1)) 
		stop("Current version of boral only supports ssvs.index = -1 when traits are supplied...sorry!")
	if(!is.null(traits)) { 
		if(!is.matrix(traits)) traits <- as.matrix(traits) 
		if(any(apply(traits,2,function(x) all(x == 1)))) { stop("No intercept column should be included in traits. It will be included automatically") } 
		}
    
    
	if(num.lv == 1) warning("We won't stop you, but one latent variable is unlikely to be successful in capturing between column correlation!")
	if(num.lv > 5) warning("We won't stop you, but please consider if you really want more than five latent variables in the model!")
	
	
	if(length(family) != ncol(y) & length(family) != 1) { stop("Number of elements in family must either one or the # of columns in y") }
	if(length(family) == 1) family <- rep(family, ncol(y))
	if(!all(family %in% c("negative.binomial", "poisson", "binomial", "normal", "lnormal", "tweedie", "ordinal", "exponential", "gamma", "beta"))) stop("At least one of the elements in family is not supported in current version of boral...sorry!")
	if(any(family == "ordinal")) {
		if(sum(y[, family == "ordinal"] == 0) > 0) stop("For ordinal data, please shift minimum level to 1.")
		message("It is assumed all ordinal columns have the same number of levels -- please see help file as to the motivation behind this.")
		message("boral may take a ridiculously long time to fit ordinal data models. Apologies in advance!") 
		if(!is.null(traits)) stop("Current version of boral does not allow traits for ordinal responses. Sorry!")
		}
	
	
	if(row.eff == FALSE) row.eff <- "none"; 
	if(row.eff == TRUE) row.eff <- "fixed"
	if(!(row.eff %in% c("none", "fixed", "random"))) stop("row.eff must be one of none/fixed/random.")
	
	if(!is.null(X)) { num.X <- ncol(X) } else { num.X <- 0 }
	if(!is.null(traits)) { num.traits <- ncol(traits) } else { num.traits <- 0 }
	
	if(num.X == 0 & num.traits > 0) 
		stop("num.traits > 0 suggests traits are to be regressed against covariates X, so please supply X.") 
 	if(num.traits > 0 & is.null(which.traits)) 
		stop("If traits are supplied, then please also supply which.traits to inform what traits are regressed against which covariates.") 
 	if(!is.null(which.traits) & ((num.X+1) != length(which.traits))) 
		stop("which.traits should have equal to 1+length(ncol(X))") 
 	if(!is.null(which.traits)) { 
		if(any(sapply(which.traits,length) > num.traits)) stop("Each element in the list which.traits should have at most ncol(traits) elements.") 
		}
 	if(!is.null(which.traits)) { 
		if(any(sapply(which.traits, function(x) any(x > ncol(traits))))) stop("The values contained in the list which.traits can be takes from 1 to ncol(traits).") 
		}
 	#if(is.null(which.traits)) { which.traits <- vector("list",num.X+1); for(k in 1:length(num.X+1)) which.traits[[k]] <- 0 } 

	
	if(!(length(prior.control$ssvs.index) %in% c(1, ncol(X)))) 
		stop("Number of elements in prior.control$ssvs.index must either be one or the # of columns in X.")
	if(length(prior.control$ssvs.index) == 1 & num.X > 0) prior.control$ssvs.index <- rep(prior.control$ssvs.index, ncol(X))
	if(any(prior.control$ssvs.index < -1)) 
		stop("Elements of prior.control$ssvs.index can only take values in -1, 0, or any positive integer; please see help file for more information.")
	
	
	if(any(family == "binomial") & !(length(trial.size) %in% c(1, length(family)))) 
        stop("trial.size needs to be specified if any columns are binomially distributed; can either be a single element or a vector equal to the # of columns in y. The latter will assume the specified trial size for all rows labelled binomial in the family argument.")
	if(any(family == "binomial") & length(trial.size) == 1) {
		complete.trial.size <- rep(0, ncol(y))
		complete.trial.size[which(family == "binomial")] <- trial.size }
	if(any(family == "binomial") & length(trial.size) == ncol(y)) { complete.trial.size <- trial.size }
	if(all(family != "binomial")) { complete.trial.size <- rep(0, ncol(y)) }
	if(all(family == "binomial") & all(complete.trial.size == 1)) { family <- rep("bernoulli",ncol(y)) }
	
	
	if(all(family != "ordinal")) { num.ord.levels <- 0; }
	if(any(family == "ordinal")) { num.ord.levels <- max(y[, family == "ordinal"]); }
	if(all(family != "multinom")) { num.multinom.levels <- 0; index.multinom.cols <- NULL }
# 	if(any(family == "multinom")) { 
# 		num.multinom.levels <- apply(y[, family == "multinom"], 2, max)
# 		index.multinom.cols <- which(family == "multinom") 
# 		}


	if(!("n.burnin" %in% names(mcmc.control))) mcmc.control$n.burnin <- 10000
	if(!("n.iteration" %in% names(mcmc.control))) mcmc.control$n.iteration <- 40000
	if(!("n.thin" %in% names(mcmc.control))) mcmc.control$n.thin <- 30
	if(!("seed" %in% names(mcmc.control))) mcmc.control$seed <- 123		


	
	n <- nrow(y); p <- ncol(y)
 	n.chains <- 1; ## Run one chain only to avoid arbitrary rotation problems
 	
 	if(num.lv > 0) 
		make.jagsboralmodel(family = family, num.X = num.X, num.traits = num.traits, which.traits = which.traits, num.lv = num.lv, row.eff = row.eff, trial.size = complete.trial.size, n = n, p = p, model.name = model.name, prior.control = prior.control)
	if(num.lv == 0)  
		make.jagsboralnullmodel(family = family, num.X = num.X, num.traits = num.traits, which.traits = which.traits, row.eff = row.eff, trial.size = complete.trial.size, n = n, p = p, model.name = model.name, prior.control = prior.control)
 	if(!do.fit) { 
		message("JAGS model file created only. Thank you, come again!")
		return() 
		}
	
		
	jags.data <- list("y", "n", "p", "num.lv", "num.X", "num.traits", "num.ord.levels", "num.multinom.levels")
	if(num.X > 0) jags.data <- c("X", jags.data)
	if(num.traits > 0) jags.data <- c("traits", jags.data)
	if(any(family == "ordinal")) { ones <- matrix(1, n, p); jags.data <- c(jags.data, "ones") }
	
	
	jags.params <- c("all.params")
	if(num.lv > 0) jags.params <- c(jags.params, "lvs")
	if(row.eff != "none") jags.params <- c(jags.params, "row.params")
	if(row.eff == "random") jags.params <- c(jags.params, "row.ranef.sigma")
	if(num.X > 0 & any(family != "multinom")) jags.params <- c(jags.params, "X.params")
	#if(num.X > 0 & any(family == "multinom")) jags.params <- c(jags.params, "X.multinom.params")
	if(num.traits > 0) jags.params <- c(jags.params, "traits.int", "traits.params", "sigma.trait")
	if(any(family == "tweedie")) jags.params <- c(jags.params, "powerparam")
	if(any(family == "ordinal")) jags.params <- c(jags.params, "alpha")
	if(any(prior.control$ssvs.index == 0)) jags.params <- c(jags.params, paste("probindX", which(prior.control$ssvs.index == 0), sep = ""))
	if(any(prior.control$ssvs.index > 0)) jags.params <- c(jags.params, paste("probGpX", unique(prior.control$ssvs.index[prior.control$ssvs.index > 0]), sep = ""))
	
	
	jags.inits <- function() {
		initial.list <- list()
		#if(num.lv > 0) initial.list$lvs <- matrix(rnorm(n*num.lv,0,0.1),n,num.lv)
#  		if(num.lv > 0) { 
#  			initial.list$all.params <- matrix(0,p,num.lv+1) 
# 			if(!all(family %in% c("poisson", "binomial", "ordinal", "multinom", "exponential","bernoulli"))) { initial.list$all.params <- cbind(initial.list$all.params,0.01) }
# 			}
		if(any(family %in% "tweedie")) initial.list$numfish = matrix(1, n, sum(family=="tweedie"))
		if(any(family %in% "ordinal")) initial.list$alpha0 <- seq(-1, 1, length = num.ord.levels - 1)
		
		if(all(family %in% "bernoulli")) {
			Tau <- rWishart(1,p+1,diag(p))[,,1]
			Sigma <- solve(Tau)
			Z <- abs(t(rmvnorm(n,rep(0,p),Sigma)))
			Z <- ifelse(as.matrix(y), Z, -1 * Z)
			initial.list$Z <- Z }
			
		return(initial.list)
		}

	set.seed(mcmc.control$seed)
	actual.filename <- model.name
	if(is.null(actual.filename)) actual.filename <- "jagsboralmodel.txt"


	## The fit! ##
 	jagsfit <- try(suppressWarnings(jags(data = jags.data, inits = jags.inits, parameters.to.save = jags.params, model.file = actual.filename, n.iter = mcmc.control$n.iteration, n.burnin = mcmc.control$n.burnin, n.chains = 1, n.thin = mcmc.control$n.thin)),silent=TRUE)
#   	jagsfit <- jags(data = jags.data, inits = jags.inits, parameters.to.save = jags.params, model.file = actual.filename, n.iter = mcmc.control$n.iteration, n.burnin = mcmc.control$n.burnin, n.chains = 1, n.thin = mcmc.control$n.thin)
    
	#print(jagsfit)
	if(inherits(jagsfit,"try-error")) {
		lookfornegbinerror <- grep("Slicer stuck at value with infinite density", jagsfit[[1]])
		if(any(family == "negative.binomial") & length(lookfornegbinerror) == 1) { 
			message("MCMC fitting through JAGS failed. This is likely due to the prior on the dispersion (size) parameter of the negative binomial distribution been too uninformative (see below). Please consider a tougher prior or switch to a Poisson family for those response that don't appear to actually be overdispersed.")
			print(jagsfit) }

		else {
			message("MCMC fitting through JAGS failed:")
			print(jagsfit) }

		message("boral fit failed...Exiting. Sorry!") 
		return()
		}
    
    
	## Format into big matrix; 
	fit.mcmcBase <- jagsfit$BUGSoutput
	fit.mcmc <- mcmc(fit.mcmcBase$sims.matrix, start = 1, thin = mcmc.control$n.thin) ## Thanks to Guilliaume Blanchet for this!
	if(n.chains == 1) {
		combined.fit.mcmc <- fit.mcmc
		}		
		
# 	if(n.chains > 1) {
# 		combined.fit.mcmc <- as.mcmc(fit.mcmcBase$sims.matrix)
# 
#    		fit.rhats <- (rhats(jagsfit, asc = FALSE))
#    		make.rhatslist <- list(lv.coefs = matrix(fit.rhats[grep("all.params", rownames(fit.rhats))], nrow = p))
#    		if(num.lv > 0) { make.rhatslist$lv.coefs <- as.matrix(make.rhatslist$lv.coefs[,-c(2:(num.lv+1))]) } ## Drop check on LV coefs
#    		rownames(make.rhatslist$lv.coefs) <- colnames(y); 
#    		if(ncol(make.rhatslist$lv.coefs) > 1) { colnames(make.rhatslist$lv.coefs) <- c("beta0","Disperson") } else { colnames(make.rhatslist$lv.coefs) <- c("beta0") }
#    		#if(num.lv > 0) { fit.rhats <- fit.rhats[-grep("lvs",rownames(fit.rhats)),] } ## Drop check on LVs
# #   	
#    		if(row.eff != "none") {
#    			make.rhatslist$row.coefs <- fit.rhats[grep("row.params", rownames(fit.rhats))]		
#    			names(make.rhatslist$row.coefs) <- rownames(y) 
#    			}
#    		if(row.eff == "random") {
#    			make.rhatslist$row.sigma <- fit.rhats[grep("row.ranef.sigma", rownames(fit.rhats))]
#   			names(make.rhatslist$row.sigma) <- c("Row random effects sigma") 
#    			}
# #   		
#   		if(num.X > 0) {
#    			make.rhatslist$X.coefs <- matrix(fit.rhats[grep("X.params", rownames(fit.rhats))], nrow = p)
#    			rownames(make.rhatslist$X.coefs) <- colnames(y); colnames(make.rhatslist$X.coefs) <- colnames(X) 
#    			}
# #   		
#    		if(num.traits > 0) {
#    			make.rhatslist$traits.coefs <- cbind(matrix(fit.rhats[grep("traits.params", rownames(fit.rhats))], nrow = ncol(X)+1),fit.rhats[grep("sigma.trait", rownames(fit.rhats))])
#    			rownames(make.rhatslist$traits.coefs) <- c("beta0",colnames(X)); colnames(make.rhatslist$traits.coefs) <- c(colnames(traits),"sigma")
#    			}
# 
#   		if(any(family == "ordinal")) {
#   			make.rhatslist$cutoffs <- fit.rhats[grep("alpha", rownames(fit.rhats))]
#   			names(make.rhatslist$cutoffs) <- paste(1:(num.ord.levels - 1), "|", 2:num.ord.levels, sep = "") 
#   			}
# #  
#   		if(any(family == "tweedie")) {
#   			make.rhatslist$powerparam <- fit.rhats[grep("powerparam", rownames(fit.rhats))]
#   			names(make.rhatslist$powerparam) <- "Common power parameter" 
#   			}

# 		exceed.rhatcutoff <- sum(sapply(make.rhatslist, function(x) sum(x > rhat.cutoff)))
# 		message("There were", exceed.rhatcutoff, "(", 100*exceed.rhatcutoff/sum(sapply(make.rhatslist,length)), "%) parameters whose Rhat exceeded the prespecified cutoff of", rhat.cutoff, "\n")		
# 		rm(fit.rhats)
# 		}
	
	rm(fit.mcmc)
    
    
#   	## Flip dispersion parameters?
#  	sel.thetas <- grep("all.params", colnames(combined.fit.mcmc))
#  	sel.thetas2 <- as.numeric(sel.thetas[(length(sel.thetas) - p + 1):length(sel.thetas)])		
#  	if(any(family %in% c("negative.binomial"))) {
#  		combined.fit.mcmc[, sel.thetas2[family %in% c("negative.binomial")]] <- 1/combined.fit.mcmc[, sel.thetas2[family %in% c("negative.binomial")]] 
#  		}

		
#   	## For any multinomial columns, set the corresponding rows in X.coefs to zero
# 	if(any(family == "multinom") & num.X > 0) {
# 		for(k in index.multinom.cols) {
# 			sel.multinom.col <- grep(paste("X.params\\[", k, ",+", sep = ""), colnames(combined.fit.mcmc))
# 			combined.fit.mcmc[, sel.multinom.col] <- 0 }
# 		}

		
 	## Make output beautiful
	if(is.null(colnames(y))) colnames(y) <- 1:ncol(y); if(is.null(rownames(y))) rownames(y) <- 1:nrow(y)
	if(num.X > 0) { if(is.null(colnames(X))) colnames(X) <- 1:ncol(X); if(is.null(rownames(X))) rownames(X) <- 1:nrow(X) }
 	if(num.traits > 0) { 
 		if(is.null(colnames(traits))) colnames(traits) <- 1:ncol(traits); if(is.null(rownames(traits))) rownames(traits) <- 1:nrow(traits) }

		
	out.fit <- list(lv.coefs.median = matrix(apply(combined.fit.mcmc[, grep("all.params", colnames(combined.fit.mcmc))], 2, median), nrow = p), lv.coefs.iqr = matrix(apply(combined.fit.mcmc[, grep("all.params", colnames(combined.fit.mcmc))], 2, IQR), nrow = p), lv.coefs.mean = matrix(apply(combined.fit.mcmc[, grep("all.params", colnames(combined.fit.mcmc))], 2, mean), nrow = p), lv.coefs.sd = matrix(apply(combined.fit.mcmc[, grep("all.params", colnames(combined.fit.mcmc))], 2, sd), nrow = p))
	rownames(out.fit$lv.coefs.median) <- rownames(out.fit$lv.coefs.iqr) <- rownames(out.fit$lv.coefs.mean) <- rownames(out.fit$lv.coefs.sd) <- colnames(y)

	
	if(num.lv > 0) { 
		out.fit$lv.median = matrix(apply(combined.fit.mcmc[, grep("lvs", colnames(combined.fit.mcmc))], 2, median), nrow = n)
		out.fit$lv.iqr = matrix(apply(combined.fit.mcmc[, grep("lvs", colnames(combined.fit.mcmc))], 2, IQR), nrow = n)
		out.fit$lv.mean = matrix(apply(combined.fit.mcmc[, grep("lvs", colnames(combined.fit.mcmc))], 2, mean), nrow = n)
		out.fit$lv.sd = matrix(apply(combined.fit.mcmc[, grep("lvs", colnames(combined.fit.mcmc))], 2, sd), nrow = n)
		rownames(out.fit$lv.median) <- rownames(out.fit$lv.iqr) <- rownames(out.fit$lv.mean) <- rownames(out.fit$lv.sd) <- rownames(y)
		colnames(out.fit$lv.median) <- colnames(out.fit$lv.iqr) <- colnames(out.fit$lv.mean) <- colnames(out.fit$lv.sd) <- paste("LV", 1:num.lv, sep = "")
		
		if(ncol(out.fit$lv.coefs.median) == (num.lv+2)) colnames(out.fit$lv.coefs.median) <- colnames(out.fit$lv.coefs.iqr) <- colnames(out.fit$lv.coefs.mean) <- colnames(out.fit$lv.coefs.sd) <- c("beta0", paste("theta", 1:num.lv, sep = ""), "Dispersion") 
		if(ncol(out.fit$lv.coefs.median) == (num.lv+1)) colnames(out.fit$lv.coefs.median) <- colnames(out.fit$lv.coefs.iqr) <- colnames(out.fit$lv.coefs.mean) <- colnames(out.fit$lv.coefs.sd) <- c("beta0", paste("theta", 1:num.lv, sep = "")) 
		}
	if(num.lv == 0) {
		if(ncol(out.fit$lv.coefs.median) == 2) colnames(out.fit$lv.coefs.median) <- colnames(out.fit$lv.coefs.iqr) <- colnames(out.fit$lv.coefs.mean) <- colnames(out.fit$lv.coefs.sd) <- c("beta0", "Dispersion") 
		if(ncol(out.fit$lv.coefs.median) == 1) colnames(out.fit$lv.coefs.median) <- colnames(out.fit$lv.coefs.iqr) <- colnames(out.fit$lv.coefs.mean) <- colnames(out.fit$lv.coefs.sd) <- c("beta0") 
		}	
	
	
	if(row.eff != "none") {
		out.fit$row.coefs.median <- apply(combined.fit.mcmc[, grep("row.params", colnames(combined.fit.mcmc))], 2, median)
		out.fit$row.coefs.iqr <- apply(combined.fit.mcmc[, grep("row.params", colnames(combined.fit.mcmc))], 2, IQR)
		out.fit$row.coefs.mean <- apply(combined.fit.mcmc[, grep("row.params", colnames(combined.fit.mcmc))], 2, mean)
		out.fit$row.coefs.sd <- apply(combined.fit.mcmc[, grep("row.params", colnames(combined.fit.mcmc))], 2, sd)
		
		names(out.fit$row.coefs.median) <- names(out.fit$row.coefs.iqr) <- names(out.fit$row.coefs.mean) <- names(out.fit$row.coefs.sd) <- rownames(y)
	
		if(row.eff == "random") {
			out.fit$row.sigma.median <- median(combined.fit.mcmc[, grep("row.ranef.sigma", colnames(combined.fit.mcmc))])
			out.fit$row.sigma.iqr <- IQR(combined.fit.mcmc[, grep("row.ranef.sigma", colnames(combined.fit.mcmc))])
			out.fit$row.sigma.mean <- mean(combined.fit.mcmc[, grep("row.ranef.sigma", colnames(combined.fit.mcmc))])
			out.fit$row.sigma.sd <- sd(combined.fit.mcmc[, grep("row.ranef.sigma", colnames(combined.fit.mcmc))])
            
			names(out.fit$row.sigma.median) <- names(out.fit$row.sigma.iqr) <- names(out.fit$row.sigma.mean) <- names(out.fit$row.sigma.sd) <- c("Row random effects sigma") }
		}

		
	if(num.X > 0) {
		out.fit$X.coefs.median <- matrix(apply(combined.fit.mcmc[, grep("X.params", colnames(combined.fit.mcmc))], 2, median), nrow = p)
		out.fit$X.coefs.iqr <- matrix(apply(combined.fit.mcmc[, grep("X.params", colnames(combined.fit.mcmc))], 2, IQR), nrow = p)
		out.fit$X.coefs.mean <- matrix(apply(combined.fit.mcmc[, grep("X.params", colnames(combined.fit.mcmc))], 2, mean), nrow = p)
		out.fit$X.coefs.sd <- matrix(apply(combined.fit.mcmc[, grep("X.params", colnames(combined.fit.mcmc))], 2, sd), nrow = p)
		rownames(out.fit$X.coefs.median) <- rownames(out.fit$X.coefs.iqr) <- rownames(out.fit$X.coefs.mean) <- rownames(out.fit$X.coefs.sd) <- colnames(y)
		colnames(out.fit$X.coefs.median) <- colnames(out.fit$X.coefs.iqr) <- colnames(out.fit$X.coefs.mean) <- colnames(out.fit$X.coefs.sd) <- colnames(X)
		
		if(any(prior.control$ssvs.index == 0) & num.traits == 0) {
			out.fit$ssvs.indcoefs.mean <- matrix(apply(combined.fit.mcmc[, grep("probindX", colnames(combined.fit.mcmc))], 
                2, mean), nrow = p)
			rownames(out.fit$ssvs.indcoefs.mean) <- colnames(y)
			colnames(out.fit$ssvs.indcoefs.mean) <- colnames(X)[which(prior.control$ssvs.index == 0)]
			out.fit$ssvs.indcoefs.sd <- matrix(apply(combined.fit.mcmc[, grep("probindX", colnames(combined.fit.mcmc))], 2, sd), nrow = p)
			rownames(out.fit$ssvs.indcoefs.sd) <- colnames(y)
			colnames(out.fit$ssvs.indcoefs.sd) <- colnames(X)[which(prior.control$ssvs.index == 0)]
			}
		if(any(prior.control$ssvs.index == 0) & num.traits > 0) {
			out.fit$ssvs.indcoefs.mean <- matrix(apply(combined.fit.mcmc[, grep("probindX", colnames(combined.fit.mcmc))], 2, mean), ncol = ncol(traits), byrow = TRUE)
			out.fit$ssvs.indcoefs.mean[out.fit$ssvs.indcoefs.mean == 2] <- NA 
			rownames(out.fit$ssvs.indcoefs.mean) <- colnames(X)
			colnames(out.fit$ssvs.indcoefs.mean) <- colnames(traits)

			out.fit$ssvs.indcoefs.sd <- matrix(apply(combined.fit.mcmc[, grep("probindX", colnames(combined.fit.mcmc))], 2, sd), ncol = ncol(traits), byrow = TRUE)
			out.fit$ssvs.indcoefs.sd[out.fit$ssvs.indcoefs.mean == 0] <- NA 
			rownames(out.fit$ssvs.indcoefs.sd) <- colnames(X)
			colnames(out.fit$ssvs.indcoefs.sd) <- colnames(traits)
			}
		if(any(prior.control$ssvs.index > 0)) {
			out.fit$ssvs.gpcoefs.mean <- apply(as.matrix(combined.fit.mcmc[, grep("probGpX", colnames(combined.fit.mcmc))]), 2, mean)
			names(out.fit$ssvs.gpcoefs.mean) <- paste("Gp", unique(prior.control$ssvs.index[prior.control$ssvs.index > 0]), sep = "")
			out.fit$ssvs.gpcoefs.sd <- apply(as.matrix(combined.fit.mcmc[, grep("probGpX", colnames(combined.fit.mcmc))]), 2, sd)
			names(out.fit$ssvs.gpcoefs.sd) <- paste("Gp", unique(prior.control$ssvs.index[prior.control$ssvs.index > 0]), sep = "")
			}
		}
		
		
	if(num.traits > 0) {
		out.fit$traits.coefs.median <- cbind(apply(combined.fit.mcmc[, grep("traits.int", colnames(combined.fit.mcmc))], 2, median), matrix(apply(combined.fit.mcmc[, grep("traits.params", colnames(combined.fit.mcmc))], 2, median), nrow = num.X+1), apply(combined.fit.mcmc[, grep("sigma.trait", colnames(combined.fit.mcmc))], 2, median))
				
		out.fit$traits.coefs.iqr <- cbind(apply(combined.fit.mcmc[, grep("traits.int", colnames(combined.fit.mcmc))], 2, IQR), matrix(apply(combined.fit.mcmc[, grep("traits.params", colnames(combined.fit.mcmc))], 2, IQR), nrow = num.X+1), apply(combined.fit.mcmc[, grep("sigma.trait", colnames(combined.fit.mcmc))], 2, IQR))
		out.fit$traits.coefs.mean <- cbind(apply(combined.fit.mcmc[, grep("traits.int", colnames(combined.fit.mcmc))], 2, mean), matrix(apply(combined.fit.mcmc[, grep("traits.params", colnames(combined.fit.mcmc))], 2, mean), nrow = num.X+1), apply(combined.fit.mcmc[, grep("sigma.trait", colnames(combined.fit.mcmc))], 2, mean))
		out.fit$traits.coefs.sd <- cbind(apply(combined.fit.mcmc[, grep("traits.int", colnames(combined.fit.mcmc))], 2, sd), matrix(apply(combined.fit.mcmc[, grep("traits.params", colnames(combined.fit.mcmc))], 2, sd), nrow = num.X+1), apply(combined.fit.mcmc[, grep("sigma.trait", colnames(combined.fit.mcmc))], 2, sd))

		rownames(out.fit$traits.coefs.median) <- rownames(out.fit$traits.coefs.iqr) <- rownames(out.fit$traits.coefs.mean) <- rownames(out.fit$traits.coefs.sd) <- c("beta0",colnames(X))
		colnames(out.fit$traits.coefs.median) <- colnames(out.fit$traits.coefs.iqr) <- colnames(out.fit$traits.coefs.mean) <- colnames(out.fit$traits.coefs.sd) <- c("kappa0",colnames(traits),"sigma")
		}

#   	if(num.X > 0 & any(family == "multinom")) {
#   		out.fit$X.multinom.coefs.median <- array(apply(combined.fit.mcmc[,grep("X.multinom.params", colnames(combined.fit.mcmc))],2,median),dim=c(length(index.multinom.cols),num.X,num.multinom.levels))
#   		out.fit$X.multinom.coefs.iqr <- array(apply(combined.fit.mcmc[,grep("X.multinom.params", colnames(combined.fit.mcmc))],2,IQR),dim=c(length(index.multinom.cols),num.X,num.multinom.levels))
#   		out.fit$X.multinom.coefs.mean <- array(apply(combined.fit.mcmc[,grep("X.multinom.params", colnames(combined.fit.mcmc))],2,mean),dim=c(length(index.multinom.cols),num.X,num.multinom.levels))
#   		out.fit$X.multinom.coefs.sd <- array(apply(combined.fit.mcmc[,grep("X.multinom.params", colnames(combined.fit.mcmc))],2,sd),dim=c(length(index.multinom.cols),num.X,num.multinom.levels))
#   
#   		dimnames(out.fit$X.multinom.coefs.median) <- dimnames(out.fit$X.multinom.coefs.iqr) <- dimnames(out.fit$X.multinom.coefs.mean) <- dimnames(out.fit$X.multinom.coefs.sd) <- list("1" = index.multinom.cols, "2" = colnames(X), "level" = 1:num.multinom.levels)
#   		}

	if(any(family == "ordinal")) {
		out.fit$cutoffs.median <- apply(combined.fit.mcmc[, grep("alpha", colnames(combined.fit.mcmc))], 2, median)
		out.fit$cutoffs.iqr <- apply(combined.fit.mcmc[, grep("alpha", colnames(combined.fit.mcmc))], 2, IQR)
		out.fit$cutoffs.mean <- apply(combined.fit.mcmc[, grep("alpha", colnames(combined.fit.mcmc))], 2, mean)
		out.fit$cutoffs.sd <- apply(combined.fit.mcmc[, grep("alpha", colnames(combined.fit.mcmc))], 2, sd)
		names(out.fit$cutoffs.median) <- names(out.fit$cutoffs.iqr) <- names(out.fit$cutoffs.mean) <- names(out.fit$cutoffs.sd) <- paste(1:(num.ord.levels - 1), "|", 2:num.ord.levels, sep = "") }

	if(any(family == "tweedie")) {
		out.fit$powerparam.median <- median(combined.fit.mcmc[, grep("powerparam", colnames(combined.fit.mcmc))])
		out.fit$powerparam.iqr <- IQR(combined.fit.mcmc[, grep("powerparam", colnames(combined.fit.mcmc))])
		out.fit$powerparam.mean <- mean(combined.fit.mcmc[, grep("powerparam", colnames(combined.fit.mcmc))])
		out.fit$powerparam.sd <- sd(combined.fit.mcmc[, grep("powerparam", colnames(combined.fit.mcmc))])
		names(out.fit$powerparam.median) <- names(out.fit$powerparam.iqr) <- names(out.fit$powerparam.mean) <- names(out.fit$powerparam.sd) <- "Common power parameter"
		}

	#print(out.fit$lv.coefs.mean)
	get.hpds <- get.hpdintervals(y, X, traits, combined.fit.mcmc, num.lv)
	out.fit$hpdintervals <- get.hpds
	if(calc.ics) {
		message("Calculating Information criteria...")
		if(num.traits > 0) 
			print("Please note that AIC and BIC at post median are currently calculated without taking into account the traits, i.e. that the column-specific coefficients are now random effects!")
		get.ics <- get.measures(y, X, family, complete.trial.size, row.eff, num.lv, combined.fit.mcmc, more.measures = FALSE)
		ics <- c(get.dic(jagsfit), get.ics$waic, get.ics$eaic, get.ics$ebic, get.ics$aic.median, get.ics$bic.median)#, get.ics$comp.lm)
		names(ics) <- c("Conditional DIC", "WAIC", "EAIC", "EBIC", "AIC at post. median", "BIC at post. median")#, "Compound L-M at post. median")
		out.fit$ics <- ics
		}
			
	if(save.model) { out.fit$jags.model <- jagsfit }

	out.fit$call <- match.call()
	out.fit$n <- n; out.fit$p <- p
	out.fit$X <- X
	out.fit$traits <- traits
	out.fit$y <- y

	out.fit$family <- family; if(all(family == "bernoulli")) out.fit$family <- rep("binomial",p)
	out.fit$num.lv <- num.lv
	out.fit$num.X <- num.X; out.fit$num.traits <- num.traits
	out.fit$which.traits <- which.traits
	out.fit$row.eff <- row.eff
	out.fit$calc.ics <- calc.ics
	out.fit$trial.size <- complete.trial.size
	out.fit$prior.control <- prior.control	
	out.fit$num.ord.levels <- num.ord.levels
	out.fit$mcmc.control <- mcmc.control
	#out.fit$n.chain <- out.fit$n.chains; 
	#if(n.chains == 1) out.fit$geweke.diag <- gdiagnostic
	
	class(out.fit) <- "boral"
	if(!save.model) { if(file.exists(actual.filename)) file.remove(actual.filename) }

	return(out.fit) }
 	

 	
################	
coefsplot <- function(covname, x, labely = NULL, est = "median", ...) {
	if(!is.null(labely)) if(!(length(labely) == nrow(x$X.coefs.median) || length(labely)==1)) 
		stop("If labely is not NULL, then it must be either of length one or a vector as long as the number of rows in x$X.coefs.median (number of species). Thanks!")
	if(!(covname %in% colnames(x$X.coefs.mean))) 
		stop("covname not found among the covariates in the boral object x.")
	
	col.seq <- rep("black", length(x$hpdintervals$X.coefs.upper[,covname]))
	col.seq[x$hpdintervals$X.coefs.lower[,covname] < 0 & x$hpdintervals$X.coefs.upper[,covname] > 0] <- "grey"
	
	At.y <- rev(1:nrow(x$X.coefs.median)) ## So that the graphs plots in the same order as rownames of x$X.coefs.median

	#.pardefault <- par(no.readonly = TRUE)
	#par(cex = a, cex.axis = a, cex.lab = a+0.5, mar = c(5,6,1,1), las = 1, cex.main = a+0.8, ...) 

	if(est == "median")
		plot(x = x$X.coefs.median[,covname], y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = covname, xlim = c(min(x$hpdintervals$X.coefs.lower[,covname]), max(x$hpdintervals$X.coefs.upper[,covname])), pch = "x", ...)
	if(est == "mean")
		plot(x = x$X.coefs.mean[,covname], y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = covname, xlim = c(min(x$hpdintervals$X.coefs.lower[,covname]), max(x$hpdintervals$X.coefs.upper[,covname])), pch = "x", ...)

	segments(x0 = x$hpdintervals$X.coefs.lower[,covname], y0 = At.y, x1 = x$hpdintervals$X.coefs.upper[,covname], y1 = At.y, col = col.seq, ...)  
	abline(v=0, lty=3)
	
	if(is.null(labely)) { axis(2, at=At.y, labels = rownames(x$X.coefs.mean), las=1, ...) } 
	if(!is.null(labely)) {
		if(length(labely) == nrow(x$X.coefs.mean)) axis(2, at=At.y, labels=labely, las=1, ...) 
		if(length(labely) == 1) mtext(text = labely, side = 2, line = 3, las = 3, ...)
		} 

 	#par(.pardefault) 	
	
	if(exists("ssvs.indcoefs.mean", x)) {
		message("Posterior probabilities of inclusion for ", covname, ":")
		print(round(x$ssvs.indcoefs.mean[,covname],3))
		message()
		}
	}

		
lvsplot <- function(x, jitter = FALSE, biplot = TRUE, ind.spp = NULL, alpha = 0.5, main = NULL, est = "median", ...) {
 	if(x$num.lv > 2) stop("Manual plotting required for plotting beyond 2 latent variables.")
 	if(x$num.lv == 0) stop("No latent variables to plot.")
 
 	n <- nrow(x$lv.median); p <- nrow(x$lv.coefs.median)
 	if(!is.null(ind.spp)) { if(ind.spp > p) { ind.spp <- p } }
	if(biplot == TRUE & !is.null(ind.spp)) { 
		message("Only the first ", ind.spp, " ``most important'' latent variable coefficients included in biplot.") }
	if(biplot == TRUE & is.null(ind.spp)) { 
		ind.spp <- p; message("All latent variable coefficients included in biplot.") }

	#.pardefault <- par(no.readonly = TRUE)
	#par(cex = a, cex.axis = a, cex.lab = a+0.5, mar = c(5,5,3,1), las = 1, cex.main = a+0.8, ...) 
 
 	if(x$num.lv == 1) {
		choose.x <- x$lv.median; 
		main <- "Plot of latent variable posterior medians"
		if(est == "mean") { choose.x <- x$lv.mean; main <- "Plot of latent variable posterior means" }
		plot(1:n, choose.x, xlab = "Row index", ylab = "Latent variable 1", main = main, type = "n", ...)
		text(x = 1:n, y = x$lv.median, label = 1:n, ...)
		}


 	if(x$num.lv == 2) {
# 		## Scale by L2norms
# 		x$lv.median2 <- scale(x$lv.median,center=TRUE,scale=FALSE) #*matrix(1/sqrt(colSums(x$lv.median^2)),n,2,byrow=TRUE)
# 		x$lv.coefs.median2 <- scale(x$lv.coefs.median[,2:3]*matrix(sqrt(colSums(x$lv.median2^2))/sqrt(colSums(x$lv.coefs.median[,2:3]^2)),p,2,byrow=TRUE),center=TRUE,scale=FALSE) 

#  		get.lv.norms <- sqrt(colSums(x$lv.median^2))
#  		get.lv.coefs.norms <- sqrt(colSums(x$lv.coefs.median[,2:3]^2))
# 		D <- diag(x=get.lv.norms,2,2)*diag(x=get.lv.coefs.norms,2,2)
#  		x$lv.median2 <- (x$lv.median*matrix(1/get.lv.norms,x$n,2,byrow=TRUE))%*%D^alpha
#  		x$lv.coefs.median2 <- (x$lv.coefs.median[,2:3]*matrix(1/get.lv.coefs.norms,x$p,2,byrow=TRUE))%*%D^(1-alpha)
 		
   		testcov <- x$lv.median%*%t(x$lv.coefs.median[,2:3])
		if(est == "mean") { testcov <- x$lv.mean%*%t(x$lv.coefs.mean[,2:3]) }

		do.svd <- svd(testcov,x$num.lv,x$num.lv)   		
   		choose.lvs <- scale(do.svd$u*matrix(do.svd$d[1:x$num.lv]^alpha,nrow=x$n,ncol=2,byrow=T),center=T, scale = F)
   		choose.lv.coefs <- scale(do.svd$v*matrix(do.svd$d[1:x$num.lv]^(1-alpha),nrow=x$p,ncol=2,byrow=T),center=T, scale = F)
   		
#  		## Scale by L2norms
# 		x$lv.mean2 <- scale(x$lv.mean,center=TRUE,scale=FALSE) #*matrix(sqrt(colSums(x$lv.coefs.mean[,2:3]^2))/sqrt(colSums(x$lv.mean^2)),n,2,byrow=TRUE) 
# 		x$lv.coefs.mean2 <- scale(x$lv.coefs.mean[,2:3]*matrix(sqrt(colSums(x$lv.mean2^2))/sqrt(colSums(x$lv.coefs.mean[,2:3]^2)),p,2,byrow=TRUE),center=FALSE,scale=FALSE) 

		if(!biplot) {
			if(is.null(main) & est == "median") { main = "Plot of latent variable posterior medians" }
			if(is.null(main) & est == "mean") { main = "Plot of latent variable posterior means" }
			plot(choose.lvs, xlab = "Latent variable 1", ylab = "Latent variable 2", main = main, type = "n", ...)
			if(!jitter) text(choose.lvs, label = 1:n, ...)
			if(jitter) text(jitter(choose.lvs[,1]), jitter(choose.lvs[,2]), label = 1:n, ...)
			}

		if(biplot) {
			if(is.null(main) & est == "median") { main = "Biplot of latent variable posterior medians" }
			if(is.null(main) & est == "mean") { main = "Biplot of latent variable posterior means" }
			largest.lnorms <- order(rowSums(choose.lv.coefs^2),decreasing=TRUE)[1:ind.spp]	
			
			plot(rbind(choose.lvs,choose.lv.coefs), xlab = "Latent variable 1", ylab = "Latent variable 2", main = main, type = "n", xlim = 1.1*range(rbind(choose.lvs,choose.lv.coefs)[,1]), ylim = 1.1*range(rbind(choose.lvs,choose.lv.coefs)[,2]), ...)
			if(!jitter) text(choose.lvs, label = 1:n, ...)
			if(jitter) text(jitter(choose.lvs[,1]), jitter(choose.lvs[,2]), label = 1:n, ...)
			text(choose.lv.coefs[largest.lnorms,], label = rownames(x$lv.coefs.mean[largest.lnorms,]), col = "red", ...)	
			}
 		}	

 	#par(.pardefault) 
 	}


print.boral <- function(x, ...) {
 	message("Call:")
 	print(x$call)
 	message()
 	message("Response matrix attributes\n \t# of rows:", x$n, "\n\t# of columns:", x$p) 
 	message("Model attributes\n \tColumn families used:", unique(x$family), "\n\t# of latent variables:", x$num.lv, "\n\tRow effect included (none/fixed/random)?", x$row.eff, "\n") 
 	if(any(x$family == "binomial")) message("Trial sizes used (columns with binomial families):", x$trial.size)
 	if(any(x$family == "ordinal")) message("Number of levels for ordinal data:", x$num.ord.levels)
 	if(x$num.X > 0) message("Model matrix with", x$num.X, "covariates also fitted\n")
 	if(x$num.traits > 0) message("Trait matrix with", x$num.traits, "traits also included\n")
 	if(any(x$prior.control$ssvs.index > -1)) message("SSVS performed on covariates with indices", x$prior.control$ssvs.index, "\n")
# 	message("Output attributes\n")
# 	print(attributes(x))
 	}

 	
print.summary.boral <- function(x, ...) {
 	message("Call:\n")
 	print(x$call)
 	message()
 	
 	if(x$est == "median") { message("Median point estimates\n\n LV coefficients (thetas) and dispersion parameter (if applicable)"); print(x$coefficients); message() }
 	if(x$est == "mean") { message("Mean point estimates\n\n LV coefficients (thetas) and dispersion parameter (if applicable)"); print(x$coefficients); message() }	
 	
 	if(!is.null(x$row.coefficients)) { message("Row coefficients\n"); print(x$row.coefficients); message() }
 	if(!is.null(x$X.coefficients)) { message("X coefficients (betas)\n"); print(x$X.coefficients); message() }
 	if(!is.null(x$X.multinom.coefficients)) { message("There are also coefficients corresponding to multinomial columns which have not been printed"); }
 	if(!is.null(x$traits.coefficients)) { message("Trait coefficients"); print(x$traits.coefficients); message() }
 	
 	if(any(x$family == "ordinal")) { message("Proportional odds (Cumulative probit) regression intercepts"); print(x$cutoffs); message() }
 	if(any(x$family == "tweedie")) { message("Tweedie power parameter"); print(x$powerparam); message() }
 
 	if(x$calc.ics) {
 		message("DIC (pD = var(deviance)/2):", as.vector(unlist(x$ics[1])))
 		message("WAIC:", as.vector(unlist(x$ics[2])))
 		message("EAIC:", as.vector(unlist(x$ics[3])))
 		message("EBIC:", as.vector(unlist(x$ics[4])))
 		message("AIC at posterior median:", as.vector(unlist(x$ics[5])))
 		message("BIC at posterior median:", as.vector(unlist(x$ics[6])))
		}

	if(!is.null(x$ssvs.indcoefs.prob)) { message("SSVS probabilities on individual coefficients"); print(x$ssvs.indcoefs.prob); message() }
	if(!is.null(x$ssvs.gpcoefs.prob)) { message("SSVS probabilities on groups of coefficients"); print(x$ssvs.gpcoefs.prob); message() }			
	}	
		
		
summary.boral <- function(object, est = "median", ...) {
	if(est == "median") {
 		gather.output <- list(call = object$call, coefficients = round(object$lv.coefs.median,3))
 		if(object$row.eff != "none") gather.output$row.coefficients = round(object$row.coefs.median,3)
 		if(object$num.X > 0) gather.output$X.coefficients = round(object$X.coefs.median,3)
 		if(object$num.traits > 0) gather.output$traits.coefficients = round(object$traits.coefs.median,3)
 		if(any(object$family == "ordinal")) gather.output$cutoffs = round(object$cutoffs.median,3)
 		if(any(object$family == "tweedie")) gather.output$powerparam = round(object$powerparam.median,3)
 		if(!is.null(object$X.multinom.coefs.median)) gather.output$X.multinom.coefficients = round(object$X.multinom.coefs.median,3) }
 
 	if(est == "mean") {
 		gather.output <- list(call = object$call, coefficients = round(object$lv.coefs.mean,3))
 		if(object$row.eff != "none") gather.output$row.coefficients = round(object$row.coefs.mean,3)
 		if(object$num.X > 0) gather.output$X.coefficients = round(object$X.coefs.mean,3)
 		if(object$num.traits > 0) gather.output$traits.coefficients = round(object$traits.coefs.mean,3)
 		if(any(object$family == "ordinal")) gather.output$cutoffs = round(object$cutoffs.mean,3)
 		if(any(object$family == "tweedie")) gather.output$powerparam = round(object$powerparam.mean,3)
 		if(!is.null(object$X.multinom.coefs.mean)) gather.output$X.multinom.coefficients = round(object$X.multinom.coefs.mean,3) }
 
 
 	gather.output$est <- est
 	gather.output$calc.ics <- object$calc.ics
	gather.output$trial.size <- object$trial.size
 	gather.output$num.ord.levels <- object$num.ord.levels
 	gather.output$prior.control$ssvs.index <- object$prior.control$ssvs.index 
 
 
	if(any(object$prior.control$ssvs.index == 0)) gather.output$ssvs.indcoefs.prob <- round(object$ssvs.indcoefs.mean,3)
	if(any(object$prior.control$ssvs.index > 0)) gather.output$ssvs.gpcoefs.prob <- round(object$ssvs.gpcoefs.mean,3) 

	if(object$calc.ics) gather.output$ics <- object$ics
 	class(gather.output) <- "summary.boral"
 	gather.output 
 	}
 			
 			
plot.boral <- function(x, est = "median", jitter = FALSE, ...) {
 	#if(all(x$family %in% c("ordinal","multinom"))) stop("Residuals are not defined, and therefore residual analysis cannot be performed, if all columns of y are ordinal.")
 	if(any(x$family %in% c("ordinal","multinom"))) 
		warning("For all columns of y that are ordinal, the first plot constructed is of Dunn-Smyth residuals against fitted values (i.e., the level with the highest predicted probability). Note this can make things very confusing to interpret if only SOME of the columns in the response matrix are ordinal.")
 	get.mus <- fitted.boral(x, est = est)$out
 	get.etas <- get.mus
 	get.ds.res <- ds.residuals(object = x, est = est)
 	print(get.ds.res$agree.ordinal)
 	get.ds.res <- get.ds.res$residuals
 	
	for(j in 1:ncol(x$y)) {
 		if(x$family[j] %in% c("binomial","beta")) get.etas[,j] <- log((get.mus[,j]+1e-5)/(1-get.mus[,j]+1e-5))
 		if(x$family[j] %in% c("poisson","lnormal","negative.binomial","tweedie","exponential")) get.etas[,j] <- log(get.mus[,j]+1e-5)
 		if(x$family[j] == "normal") get.etas[,j] <- (get.mus[,j]) 
 		if(x$family[j] == "ordinal") { } ## Fitted values are the class with highest probability, which is already attained from fitted.boral
 		}
 
 
	#.pardefault <- par(no.readonly = TRUE)	
	#par(ask = TRUE, cex = a, mar = c(5,5,2,1), cex.lab = 0.8*a, cex.main = a, las = 1, ...) 
 	palette(rainbow(ncol(get.etas)))
 
 	matplot(get.etas, get.ds.res, ylab = "Dunn-Smyth Residuals", xlab = "Linear Predictors", type="n", ...)
 	for(i in 1:ncol(get.etas)) { points(get.etas[,i], get.ds.res[,i], col=palette()[i], ...) }
 	abline(h=0, lty = 2, lwd = 2)
	
# 	matplot(get.mus, get.ds.res, ylab = "Dunn-Smyth Residuals", xlab = "Fitted Values", type="n")
# 	for(i in 1:ncol(get.mus)) { points(get.mus[,i], get.ds.res[,i], col=palette()[i]) }
# 	abline(h=0, lty = 2, lwd = 2)father

	matplot(get.ds.res, ylab = "Dunn-Smyth Residuals", xlab = "Row index",type="n", xaxt = "n", ...)
 	axis(side = 1, at = 1:nrow(x$y), labels = rownames(x$lv.mean), ...)
 	for (i in 1:ncol(get.mus)) { points(seq(1,nrow(x$y)),get.ds.res[,i], col=palette()[i], ...) }
 	abline(0,0,lty=2)
 
 	matplot(t(get.ds.res), ylab = "Dunn-Smyth Residuals", xlab = "Column index", type="n", xaxt = "n", ...)
 	axis(side = 1, at = 1:ncol(x$y), labels = rownames(x$coefs.mean), ...)
 	for(i in 1:ncol(get.mus)) { points(rep(i,nrow(get.etas)), get.ds.res[,i], col=palette()[i], ...) }
 	abline(h=0, lty = 2, lwd = 2)
 
	get.ds.res2 <- as.vector(unlist(get.ds.res))
 	qqnorm(get.ds.res2[is.finite(get.ds.res2)], main = "Normal Quantile Plot", ...)
 	abline(0,1, ...)
 	
 	palette("default")
 	#par(.pardefault) 	
 	}

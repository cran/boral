##############
## Latent Variable Ordination and Regression using MCMC 
## NB parameterized as V = mu + phi*mu^2
## Ordinal data handled as propotional odds regression, with common cutoff points for all spp, but spp intercepts included as random effects; This should make sense if the same ordinal scale is applied to all species
## Multinomial regression is available, but CURRENTLY NOT ALLOWED DUE TO COMPUTATION TIME

## Changes from v1.1 (news files to be updated!)
## 1) 0-1 line on qqplots for Dunn Smyth residuals is now replaced with a general qqline
## 2) Help files for boral reformatted to move some to a separate sub help file, analogous to how mgcv does things =D
## 3) hpdintervals reformatted to simplify output
## 4) traits are now allows for ordinal regression. The exception is if you have a single ordinal column, in which case the intercept for this is forced to be zero by parameter identifiability reasons.


## TODO: 
## 0) Reformat output for boral so that instead of a separate element for mean,median,IQR,sd; present it as a 4 column array or vector. Without forcing all users to adopt this, the only simple way of doing this is to create a function that reformats boral output?
## 1) correct problems with increasing number of LVs causing increase in covariance and variance i.e., scale proborla outblem?; 
## 2) Construct predict functions; see predictfunctions.R 
## 3) Consider switching to weakly informative priors by default?
## 4) Reduce rank species coefficients to get constrained ordination? HARD!!!
## 5) Draw species coefficients as random effects. In principle, the current function could allow this if you require users to manually include an intercept into traits. Then if only an intercept column is included for traits, the species coefs for each covariate are draw from a common mean. 
## 7) allow for multiple chains, but don't check convergence on the LVs and their loadings. Also cannot combined chains for LV and loadings unless you post process them, which is annoying. Thanks about how to parallelize JAGS if you can get this going though
## 8) Allow SSVS for selecting traits?

## Make sure you test all dontrun examples before building package!!!
## example(topic = "about.ssvs", package = "boral", run.dontrun=TRUE)
## example(topic = "about.traits", package = "boral", run.dontrun=TRUE)
## example(topic = "boral", package = "boral", run.dontrun=TRUE)
## example(topic = "calc.condlogLik", package = "boral", run.dontrun=TRUE)
## etc...


##############
# library(R2jags); 
# library(mvtnorm); 
# library(mvabund); 
# library(coda); 
# library(MASS)
# library(fishMod)
# library(abind)
# source("auxilaryfunctions.R")
# source("makejagsboralmodel.R")
# source("makejagsboralnullmodel.R")
# source("unseenfunctions.R")


# n = 60; p <- 30
# X <- matrix(rnorm(n*2),n,2); beta <- cbind(matrix(rnorm(p*3),p,3),runif(p,0,5)); true.power <- 1.6
# mu <- exp(cbind(1,X)%*%t(beta[,1:3]))
# y <- matrix(NA,n,p)
# for(j in 1:ncol(y)) { y[,j] <- rTweedie(nrow(y), mu = mu[,j], phi = beta[j,4], p = true.power) }
# family = "tweedie"
# num.lv = 0; row.eff = "none"; n.burnin = 10000; n.iteration = 40000; n.thin = 30; save.model = TRUE; seed = 1; calc.ics = TRUE; trial.size = NULL; hypparams = c(50,20,50,50); 
# 
# y = sim.y$resp; family = "ordinal"; num.lv = 2; mcmc.control = example.mcmc.control; traits = NULL; which.traits = NULL; save.model = TRUE; calc.ics = TRUE; row.eff = "none"; trial.size = 1; prior.control = list(type = c("normal","normal","normal","uniform"), hypparams = c(100, 20, 100, 50), ssvs.index = -1, ssvs.g = 1e-6); do.fit = TRUE; model.name = NULL; row.ids = NULL


boral <- function(y, ...) UseMethod("boral")


## Model is g(mu_{ij}) = row + beta0 + LV_i*theta_j + X_i*beta_j
boral.default <- function (y, X = NULL, traits = NULL, which.traits = NULL, family, trial.size = 1, num.lv = 0, row.eff = "none", row.ids = NULL, save.model = FALSE, calc.ics = TRUE, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123), prior.control = list(type = c("normal","normal","normal","uniform"), hypparams = c(100, 20, 100, 50), ssvs.index = -1, ssvs.g = 1e-6), do.fit = TRUE, model.name = NULL, ...) {
	new.format <- FALSE

	if(is.null(dim(y))) { 
		message("Converting y into a one column matrix."); 
		y <- matrix(y, ncol = 1) 
		}
	if(!is.null(X) && is.null(dim(X))) { 
		message("Converting X into a one column matrix."); 
		X <- matrix(X, ncol = 1) 
		}
	if(!is.null(traits) & is.null(dim(traits))) { 
		message("Converting traits into a one column matrix."); 
		traits <- matrix(traits, ncol = 1) 
		}
	if(!is.null(X)) { 
		if(!is.matrix(X)) X <- as.matrix(X)
		if(any(apply(X,2,function(x) all(x == 1)))) { stop("No intercept column should be included in X.") } 
		}
	
	
	if(!("type" %in% names(prior.control))) prior.control$type <- c("normal","normal","normal","uniform")
	if(!("hypparams" %in% names(prior.control))) prior.control$hypparams <- c(100, 20, 100, 50)
	if(!("ssvs.index" %in% names(prior.control))) prior.control$ssvs.index <- -1		
	if(!("ssvs.g" %in% names(prior.control))) prior.control$ssvs.g <- 1e-6
	
	if(length(prior.control$hypparams) != 4 || length(prior.control$type) != 4) 
		stop("prior.control$type and prior.control$hypparams must be a vector of four elements. Please see boral help file as to what the elements correspond to.\n")
	if(!all(prior.control$type[-4] %in% c("normal","uniform","cauchy")))	
		stop("At least one of the first three elements of prior.control$type is not supported in the current version of boral...sorry!")
	if(!(prior.control$type[4] %in% c("uniform","halfcauchy","halfnormal")))	
		stop("The fourth element of prior.control$type is not supported in the current version of boral...sorry!")
		
 	if(!is.null(which.traits) & any(prior.control$ssvs.index > -1)) 
		stop("The current version of boral only supports ssvs.index = -1 when traits are supplied...sorry!")
	if(!is.null(traits)) { 
		if(!is.matrix(traits)) traits <- as.matrix(traits) 
		if(any(apply(traits,2,function(x) all(x == 1)))) 
			stop("No intercept column should be included in traits. It will be included automatically.")  
		}
    
    
	if(num.lv == 1) 
		warning("We won't stop you, but one latent variable is unlikely to be successful in capturing between column correlation!")
	if(num.lv > 5) 
		warning("We won't stop you, but please consider if you really want more than five latent variables in the model!")
	
	
	if(length(family) != ncol(y) & length(family) != 1) 
		stop("Number of elements in family must either one or the # of columns in y") 
	if(length(family) == 1) family <- rep(family, ncol(y))
	if(!all(family %in% c("negative.binomial", "poisson", "binomial", "normal", "lnormal", "tweedie", "ordinal", "exponential", "gamma", "beta"))) 
		stop("At least one of the elements in family is not supported in current version of boral...sorry!")
	if(any(family == "ordinal")) {
		if(sum(y[, family == "ordinal"] == 0) > 0) 
			stop("For ordinal data, please shift minimum level to 1.")
 		if(!is.null(traits) & (sum(family == "ordinal") == 1)) 
			message("The intercept for the single ordinal response is set to zero and not regressed traits for parameter identifiability reasons.")
		}
	
	
	if(row.eff == FALSE) row.eff <- "none"; 
	if(row.eff == TRUE) row.eff <- "fixed"
	if(!(row.eff %in% c("none", "fixed", "random"))) 
		stop("row.eff must be one of none/fixed/random.")
	if(row.eff != "none" && is.null(row.ids)) {
		row.ids <- matrix(1:nrow(y), ncol = 1)
		message("row.ids assumed to be matrix with one column and elements 1,2,...nrow(y) i.e., a row-specific intercept.")
		}
	if(!is.null(row.ids)) {
		row.ids <- as.matrix(row.ids)
		if(nrow(row.ids) != nrow(y)) 
			stop("Number of rows in the matrix row.ids should be equal to number of rows in y.")
		if(is.null(colnames(row.ids))) colnames(row.ids) <- paste0("ID", 1:ncol(row.ids))
		}
		
		
	if(!is.null(X)) { num.X <- ncol(X) } else { num.X <- 0 }
	if(!is.null(traits)) { num.traits <- ncol(traits) } else { num.traits <- 0 }
	
	if(num.X == 0 & num.traits > 0) 
		stop("num.traits > 0 suggests traits are to be regressed against covariates X, so please supply X.") 
 	if(num.traits > 0 & is.null(which.traits)) 
		stop("If traits are supplied, then please also supply which.traits to inform what traits are regressed against which covariates.") 
 	if(!is.null(which.traits) & ((num.X+1) != length(which.traits))) 
		stop("which.traits should have equal to 1+length(ncol(X))") 
 	if(!is.null(which.traits)) if(any(sapply(which.traits,length) > num.traits)) 
		stop("Each element in the list which.traits should have at most ncol(traits) elements.") 
 	if(!is.null(which.traits)) if(any(sapply(which.traits, function(x) any(x > ncol(traits))))) 
		stop("The values contained in the list which.traits can be takes from 1 to ncol(traits).") 
	
	if(!(length(prior.control$ssvs.index) %in% c(1, ncol(X)))) 
		stop("Number of elements in prior.control$ssvs.index must either be one or the # of columns in X.")
	if(length(prior.control$ssvs.index) == 1 & num.X > 0) prior.control$ssvs.index <- rep(prior.control$ssvs.index, ncol(X))
	if(any(prior.control$ssvs.index < -1)) 
		stop("Elements of prior.control$ssvs.index can only take values in -1, 0, or any positive integer; please see help file for more information.")
	
	
	if(any(family == "binomial") & !(length(trial.size) %in% c(1, length(family)))) 
        stop("trial.size needs to be specified if any columns are binomially distributed; can either be a single element or a vector equal to the # of columns in y.")
	if(any(family == "binomial") & length(trial.size) == 1) {
		complete.trial.size <- rep(0, ncol(y))
		complete.trial.size[which(family == "binomial")] <- trial.size 
		}
	if(any(family == "binomial") & length(trial.size) == ncol(y)) complete.trial.size <- trial.size
	if(all(family != "binomial")) complete.trial.size <- rep(0, ncol(y))
	if(all(family == "binomial") & all(complete.trial.size == 1)) family <- rep("bernoulli",ncol(y))
	
	
	if(all(family != "ordinal")) num.ord.levels <- 0
	if(any(family == "ordinal")) num.ord.levels <- max(y[, family == "ordinal"])
	if(all(family != "multinom")) { num.multinom.levels <- 0; index.multinom.cols <- NULL }
# 	if(any(family == "multinom")) { 
# 		num.multinom.levels <- apply(y[, family == "multinom"], 2, max)
# 		index.multinom.cols <- which(family == "multinom") 
# 		}


	if(!("n.burnin" %in% names(mcmc.control))) mcmc.control$n.burnin <- 10000
	if(!("n.iteration" %in% names(mcmc.control))) mcmc.control$n.iteration <- 40000
	if(!("n.thin" %in% names(mcmc.control))) mcmc.control$n.thin <- 30
	if(!("seed" %in% names(mcmc.control))) mcmc.control$seed <- 123		

	
	#####
	n <- nrow(y); p <- ncol(y)
 	n.chains <- 1; ## Run one chain only to avoid arbitrary rotation problems
 	
 	if(num.lv > 0) 
		make.jagsboralmodel(family = family, num.X = num.X, num.traits = num.traits, which.traits = which.traits, num.lv = num.lv, row.eff = row.eff, row.ids = row.ids, trial.size = complete.trial.size, n = n, p = p, model.name = model.name, prior.control = prior.control)
	if(num.lv == 0)  
		make.jagsboralnullmodel(family = family, num.X = num.X, num.traits = num.traits, which.traits = which.traits, row.eff = row.eff, row.ids = row.ids, trial.size = complete.trial.size, n = n, p = p, model.name = model.name, prior.control = prior.control)
 	if(!do.fit) { 
		message("JAGS model file created only. Thank you, come again!")
		return() 
		}
	
		
	jags.data <- list("y", "n", "p", "num.lv", "num.X", "num.traits", "num.ord.levels", "num.multinom.levels")
	if(num.X > 0) jags.data <- c("X", jags.data)
	if(num.traits > 0) jags.data <- c("traits", jags.data)
	if(any(family == "ordinal")) { ones <- matrix(1, n, p); jags.data <- c(jags.data, "ones") }
	if(row.eff != "none") { 
		n.ID <- apply(row.ids, 2, function(x) length(unique(x)))
		jags.data <- c(jags.data, "row.ids", "n.ID") }
	
	
	jags.params <- c("all.params")
	if(num.lv > 0) jags.params <- c(jags.params, "lvs")
	if(row.eff != "none") jags.params <- c(jags.params, paste0("row.params.ID",1:ncol(row.ids)))
	if(row.eff == "random") jags.params <- c(jags.params, paste0("row.ranef.sigma.ID",1:ncol(row.ids)))
	if(num.X > 0 & any(family != "multinom")) jags.params <- c(jags.params, "X.params")
	#if(num.X > 0 & any(family == "multinom")) jags.params <- c(jags.params, "X.multinom.params")
	if(num.traits > 0) jags.params <- c(jags.params, "traits.int", "traits.params", "sigma.trait")
	if(any(family == "tweedie")) jags.params <- c(jags.params, "powerparam")
	if(any(family == "ordinal")) jags.params <- c(jags.params, "alpha", "ordinal.ranef.sigma")
	if(any(prior.control$ssvs.index == 0)) jags.params <- c(jags.params, paste0("probindX", which(prior.control$ssvs.index == 0)))
	if(any(prior.control$ssvs.index > 0)) jags.params <- c(jags.params, paste0("probGpX", unique(prior.control$ssvs.index[prior.control$ssvs.index > 0])))
	
	
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
			message("MCMC sampling through JAGS failed. This is likely due to the prior on the dispersion (size) parameter of the negative binomial distribution been too uninformative (see below). For instance, if the error message refers to all.params[25,4], then this means the MCMC sampling ran into issues for column (species) 25 in y.\n
			Please consider the following advice: 1) remove very rare species like singletons and doubletons, as they can potentially issus for MCMC sampling, and do not provide much information about the species community in general any way; 2) consider switching to a Poisson family for those response that don't appear to actually be overdispersed; 3) adopt a tougher prior for the overdispersion parameter e.g., keeping with a uniform prior but reducing hypparams[4], or using a half-cauchy or half-normal prior with small variance. ")
			print(jagsfit) }

		else {
			message("MCMC fitting through JAGS failed:")
			print(jagsfit) }

		message("boral fit failed...Exiting. Sorry!") 
		return()
		}
    
# 	return(jagsfit)
	## Format into big matrix; 
	fit.mcmcBase <- jagsfit$BUGSoutput
	fit.mcmc <- mcmc(fit.mcmcBase$sims.matrix, start = 1, thin = mcmc.control$n.thin) ## Thanks to Guilliaume Blanchet for this!
	if(n.chains == 1) {
		combined.fit.mcmc <- fit.mcmc
		}		
		

# 	if(n.chains > 1) {
# 		get.rhats <- process.rhats(sims.matrix = fit.mcmcBase$sims.matrix)

# 		exceed.rhatcutoff <- sum(sapply(get.rhats, function(x) sum(x > rhat.cutoff)))
# 		message("There were", exceed.rhatcutoff, "(", 100*exceed.rhatcutoff/sum(sapply(get.rhats,length)), "%) parameters whose Rhat exceeded the prespecified cutoff of", rhat.cutoff, "\n")		
# 		}	
	rm(fit.mcmc)
    
    
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


	lv.coefs.arr <- abind(
		matrix(apply(combined.fit.mcmc[, grep("all.params", colnames(combined.fit.mcmc))], 2, median), nrow = p),
		matrix(apply(combined.fit.mcmc[, grep("all.params", colnames(combined.fit.mcmc))], 2, mean), nrow = p),
		matrix(apply(combined.fit.mcmc[, grep("all.params", colnames(combined.fit.mcmc))], 2, IQR), nrow = p),
		matrix(apply(combined.fit.mcmc[, grep("all.params", colnames(combined.fit.mcmc))], 2, sd), nrow = p), 
		along = 3)
	
	out.fit <- list()

	if(num.lv > 0) {
		lv.arr <- abind(
			matrix(apply(combined.fit.mcmc[, grep("lvs", colnames(combined.fit.mcmc))], 2, median), nrow = n),
			matrix(apply(combined.fit.mcmc[, grep("lvs", colnames(combined.fit.mcmc))], 2, mean), nrow = n),
			matrix(apply(combined.fit.mcmc[, grep("lvs", colnames(combined.fit.mcmc))], 2, IQR), nrow = n),
			matrix(apply(combined.fit.mcmc[, grep("lvs", colnames(combined.fit.mcmc))], 2, sd), nrow = n),
			along = 3)
		dimnames(lv.arr) <- list(rows = rownames(y), lvs = paste0("LV", 1:num.lv), type = c("median","mean","iqr","sd"))
		
		if(new.format) { out.fit$lv <- lv.arr }
		if(!new.format) { 
			out.fit$lv.median <- lv.arr[,,1]; out.fit$lv.mean <- lv.arr[,,2]; 
			out.fit$lv.iqr <- lv.arr[,,3]; out.fit$lv.sd <- lv.arr[,,4]
			}
			
		if(dim(lv.coefs.arr)[2] == (num.lv+2)) 
			dimnames(lv.coefs.arr) <- list(cols = colnames(y), coefficients = c("beta0", paste0("theta", 1:num.lv), "Dispersion"), type = c("median","mean","iqr","sd"))
		if(dim(lv.coefs.arr)[2] == (num.lv+1)) 
			dimnames(lv.coefs.arr) <- list(cols = colnames(y), coefficients = c("beta0", paste0("theta", 1:num.lv)), type = c("median","mean","iqr","sd"))
		}
		
	if(num.lv == 0) {
		if(dim(lv.coefs.arr)[2] == 2) 
			dimnames(lv.coefs.arr) <- list(cols = colnames(y), coefficients = c("beta0", "Dispersion"), type = c("median","mean","iqr","sd"))
		if(dim(lv.coefs.arr)[2] == 1) 
			dimnames(lv.coefs.arr) <- list(cols = colnames(y), coefficients = c("beta0"), type = c("median","mean","iqr","sd"))
		}	
	if(new.format) { out.fit$lv.coefs <- lv.coefs.arr }
	if(!new.format) { 
		out.fit$lv.coefs.median <- lv.coefs.arr[,,1]; out.fit$lv.coefs.mean <- lv.coefs.arr[,,2]
		out.fit$lv.coefs.iqr <- lv.coefs.arr[,,3]; out.fit$lv.coefs.sd <- lv.coefs.arr[,,4]
		}
	
	if(row.eff != "none") {
		out.fit$row.coefs <- vector("list", ncol(row.ids))
		names(out.fit$row.coefs) <- colnames(row.ids)
		for(k in 1:ncol(row.ids)) {
			row.coefs.arr <- cbind(
				apply(combined.fit.mcmc[, grep(paste0("row.params.ID",k), colnames(combined.fit.mcmc))], 2, median), 
				apply(combined.fit.mcmc[, grep(paste0("row.params.ID",k), colnames(combined.fit.mcmc))], 2, mean), 
				apply(combined.fit.mcmc[, grep(paste0("row.params.ID",k), colnames(combined.fit.mcmc))], 2, IQR), 
				apply(combined.fit.mcmc[, grep(paste0("row.params.ID",k), colnames(combined.fit.mcmc))], 2, sd))
			rownames(row.coefs.arr) <- 1:n.ID[k]; colnames(row.coefs.arr) <- c("median","mean","iqr","sd")
				
			if(new.format) { out.fit$row.coefs[[k]] <- row.coefs.arr }
			if(!new.format) { 
				out.fit$row.coefs[[k]] <- list(median = row.coefs.arr[,1], mean = row.coefs.arr[,2], iqr = row.coefs.arr[,3], sd = row.coefs.arr[,4])
				}
			}
	
		if(row.eff == "random") {
			out.fit$row.sigma <- vector("list", ncol(row.ids))
			names(out.fit$row.sigma) <- colnames(row.ids)
			for(k in 1:ncol(row.ids)) {
				row.sigma.vec <- c(
					median(combined.fit.mcmc[, grep(paste0("row.ranef.sigma.ID",k), colnames(combined.fit.mcmc))]),
					mean(combined.fit.mcmc[, grep(paste0("row.ranef.sigma.ID",k), colnames(combined.fit.mcmc))]),
					IQR(combined.fit.mcmc[, grep(paste0("row.ranef.sigma.ID",k), colnames(combined.fit.mcmc))]),
					sd(combined.fit.mcmc[, grep(paste0("row.ranef.sigma.ID",k), colnames(combined.fit.mcmc))]))
				names(row.sigma.vec) <- c("median","mean","iqr","sd")
				
				if(new.format) { out.fit$row.sigma[[k]] <- row.sigma.vec }
				if(!new.format) { 
					out.fit$row.sigma[[k]] <- list(median = row.sigma.vec[1], mean = row.sigma.vec[2], iqr = row.sigma.vec[3], sd = row.sigma.vec[4])				
					}
				}
			}
		}
		
	if(num.X > 0) {
		X.coefs.arr <- abind(
			matrix(apply(combined.fit.mcmc[, grep("X.params", colnames(combined.fit.mcmc))], 2, median), nrow = p),
			matrix(apply(combined.fit.mcmc[, grep("X.params", colnames(combined.fit.mcmc))], 2, mean), nrow = p),
			matrix(apply(combined.fit.mcmc[, grep("X.params", colnames(combined.fit.mcmc))], 2, IQR), nrow = p),
			matrix(apply(combined.fit.mcmc[, grep("X.params", colnames(combined.fit.mcmc))], 2, sd), nrow = p),
			along = 3)
		dimnames(X.coefs.arr) <- list(cols = colnames(y), X = colnames(X), type = c("median","mean","iqr","sd"))

		if(new.format) { out.fit$X.coefs <- X.coefs.arr }
		if(!new.format) { 
			out.fit$X.coefs.median <- X.coefs.arr[,,1]; out.fit$X.coefs.mean <- X.coefs.arr[,,2]
			out.fit$X.coefs.iqr <- X.coefs.arr[,,3]; out.fit$X.coefs.sd <- X.coefs.arr[,,4]
			}
				
		if(any(prior.control$ssvs.index == 0) & num.traits == 0) {
			ssvs.indcoefs.arr <- array(0, dim = c(p, sum(prior.control$ssvs.index == 0), 2))
			ssvs.indcoefs.arr[,,1] <- matrix(apply(combined.fit.mcmc[, grep("probindX", colnames(combined.fit.mcmc))], 
                2, mean), nrow = p)
			ssvs.indcoefs.arr[,,2] <- matrix(apply(combined.fit.mcmc[, grep("probindX", colnames(combined.fit.mcmc))], 
                2, sd), nrow = p)
			dimnames(ssvs.indcoefs.arr) <- list(cols = colnames(y), X = colnames(X)[which(prior.control$ssvs.index == 0)], type = c("mean","sd"))
                
			if(new.format) { out.fit$ssvs.indcoefs <- ssvs.indcoefs.arr }
			if(!new.format) { 
				out.fit$ssvs.indcoefs.mean <- ssvs.indcoefs.arr[,,1]; out.fit$ssvs.indcoefs.sd <- ssvs.indcoefs.arr[,,2]
				}
			}
		if(any(prior.control$ssvs.index == 0) & num.traits > 0) { ## Currently not used since SSVS not permitted on traits
			ssvs.indcoefs.arr <- array(0, dim = c(p, num.traits, 2))
			ssvs.indcoefs.arr[,,1] <- matrix(apply(combined.fit.mcmc[, grep("probindX", colnames(combined.fit.mcmc))], 2, mean), ncol = ncol(traits), byrow = TRUE)
			ssvs.indcoefs.arr[,,1][ssvs.indcoefs.arr[,,1] == 2] <- NA 
			ssvs.indcoefs.arr[,,2] <- matrix(apply(combined.fit.mcmc[, grep("probindX", colnames(combined.fit.mcmc))], 2, sd), ncol = ncol(traits), byrow = TRUE)
			ssvs.indcoefs.arr[,,2][ssvs.indcoefs.arr[,,1] == 2] <- NA 
			dimnames(ssvs.indcoefs.arr) <- list(cols = colnames(y), traits = colnames(traits), type = c("mean","sd"))
				
			if(new.format) { out.fit$ssvs.indcoefs <- ssvs.indcoefs.arr }
			if(!new.format) { 
				out.fit$ssvs.indcoefs.mean <- ssvs.indcoefs.arr[,,1]; out.fit$ssvs.indcoefs.sd <- ssvs.indcoefs.arr[,,2]
				}
			}
		if(any(prior.control$ssvs.index > 0)) {
			ssvs.gpcoefs.arr <- cbind(
				apply(as.matrix(combined.fit.mcmc[, grep("probGpX", colnames(combined.fit.mcmc))]), 2, mean),
				apply(as.matrix(combined.fit.mcmc[, grep("probGpX", colnames(combined.fit.mcmc))]), 2, sd))
			rownames(ssvs.gpcoefs.arr) <- paste0("Gp", unique(prior.control$ssvs.index[prior.control$ssvs.index > 0]))	
			colnames(ssvs.gpcoefs.arr) <- c("mean","sd")

			if(new.format) { out.fit$ssvs.gpcoefs <- ssvs.gpcoefs.arr }
			if(!new.format) { 
				out.fit$ssvs.gpcoefs.mean <- ssvs.gpcoefs.arr[,1]; out.fit$ssvs.gpcoefs.sd <- ssvs.gpcoefs.arr[,2]
				}
			}
		}
		
		
	if(num.traits > 0) {
		traits.coefs.arr <- array(0, dim = c(num.X+1, num.traits + 2, 4))
		traits.coefs.arr[,,1] <- cbind(
			apply(combined.fit.mcmc[, grep("traits.int", colnames(combined.fit.mcmc))], 2, median), 
			matrix(apply(combined.fit.mcmc[, grep("traits.params", colnames(combined.fit.mcmc))], 2, median), nrow = num.X+1), 
			apply(combined.fit.mcmc[, grep("sigma.trait", colnames(combined.fit.mcmc))], 2, median))
		traits.coefs.arr[,,2] <- cbind(
			apply(combined.fit.mcmc[, grep("traits.int", colnames(combined.fit.mcmc))], 2, mean), matrix(apply(combined.fit.mcmc[, 
			grep("traits.params", colnames(combined.fit.mcmc))], 2, mean), nrow = num.X+1), apply(combined.fit.mcmc[, 
			grep("sigma.trait", colnames(combined.fit.mcmc))], 2, mean))
		traits.coefs.arr[,,3] <- cbind(
			apply(combined.fit.mcmc[, grep("traits.int", colnames(combined.fit.mcmc))], 2, IQR), matrix(apply(combined.fit.mcmc[, 
			grep("traits.params", colnames(combined.fit.mcmc))], 2, IQR), nrow = num.X+1), apply(combined.fit.mcmc[, 
			grep("sigma.trait", colnames(combined.fit.mcmc))], 2, IQR))
		traits.coefs.arr[,,4] <- cbind(
			apply(combined.fit.mcmc[, grep("traits.int", colnames(combined.fit.mcmc))], 2, sd), matrix(apply(combined.fit.mcmc[, 
			grep("traits.params", colnames(combined.fit.mcmc))], 2, sd), nrow = num.X+1), apply(combined.fit.mcmc[, 
			grep("sigma.trait", colnames(combined.fit.mcmc))], 2, sd))
		dimnames(traits.coefs.arr) <- list(X.coefficients = c("beta0",colnames(X)), traits.coefficients =  c("kappa0",colnames(traits),"sigma"), type = c("median","mean","iqr","sd"))
			
		if(new.format) { out.fit$traits.coefs <- traits.coefs.arr }
		if(!new.format) { 
			out.fit$traits.coefs.median <- traits.coefs.arr[,,1]; out.fit$traits.coefs.mean <- traits.coefs.arr[,,2]
			out.fit$traits.coefs.iqr <- traits.coefs.arr[,,3]; out.fit$traits.coefs.sd <- traits.coefs.arr[,,4]
			}
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
		cutoffs.arr <- cbind(
			apply(combined.fit.mcmc[, grep("alpha", colnames(combined.fit.mcmc))], 2, median),
			apply(combined.fit.mcmc[, grep("alpha", colnames(combined.fit.mcmc))], 2, mean),
			apply(combined.fit.mcmc[, grep("alpha", colnames(combined.fit.mcmc))], 2, IQR),
			apply(combined.fit.mcmc[, grep("alpha", colnames(combined.fit.mcmc))], 2, sd))
		rownames(cutoffs.arr) <- paste0(1:(num.ord.levels - 1), "|", 2:num.ord.levels) 
		colnames(cutoffs.arr) <- c("median","mean","iqr","sd")
			
		if(new.format) { out.fit$cutoffs <- cutoffs.arr }
		if(!new.format) { 
			out.fit$cutoffs.median <- cutoffs.arr[,1]; out.fit$cutoffs.mean <- cutoffs.arr[,2]
			out.fit$cutoffs.iqr <- cutoffs.arr[,3]; out.fit$cutoffs.sd <- cutoffs.arr[,4]
			}

		if(sum(family == "ordinal") > 1 & is.null(traits)) { ## If there are traits, then ordinal random intercept is either zero (if there is only 1 ordinal column, or has the sigma.trait (if there are >1 ordinal columns)
			ordinal.sigma.vec <- c(
				median(combined.fit.mcmc[, grep("ordinal.ranef.sigma", colnames(combined.fit.mcmc))]),
				mean(combined.fit.mcmc[, grep("ordinal.ranef.sigma", colnames(combined.fit.mcmc))]),
				IQR(combined.fit.mcmc[, grep("ordinal.ranef.sigma", colnames(combined.fit.mcmc))]),
				sd(combined.fit.mcmc[, grep("ordinal.ranef.sigma", colnames(combined.fit.mcmc))]))
			names(ordinal.sigma.vec) <- c("median","mean","iqr","sd")

			if(new.format) { out.fit$ordinal.sigma <- ordinal.sigma.vec }
			if(!new.format) { 
				out.fit$ordinal.sigma.median <- ordinal.sigma.vec[1]; out.fit$ordinal.sigma.mean <- ordinal.sigma.vec[2]
				out.fit$ordinal.sigma.iqr <- ordinal.sigma.vec[3]; out.fit$ordinal.sigma.sd <- ordinal.sigma.vec[4]            
				}
			}
		}

	if(any(family == "tweedie")) {
		powerparam.vec <- c(
			median(combined.fit.mcmc[, grep("powerparam", colnames(combined.fit.mcmc))]),
			mean(combined.fit.mcmc[, grep("powerparam", colnames(combined.fit.mcmc))]),
			IQR(combined.fit.mcmc[, grep("powerparam", colnames(combined.fit.mcmc))]),
			sd(combined.fit.mcmc[, grep("powerparam", colnames(combined.fit.mcmc))]))
		names(powerparam.vec) <- c("median","mean","iqr","sd")
			
		if(new.format) { out.fit$powerparam <- powerparam.vec }
		if(!new.format) { 
			out.fit$powerparam.median <- powerparam.vec[1]; out.fit$powerparam.mean <- powerparam.vec[2]
			out.fit$powerparam.iqr <- powerparam.vec[3]; out.fit$powerparam.sd <- powerparam.vec[4]
			}
		}
	rm(list = ls(pattern = ".arr"))	
		

	#print(out.fit$lv.coefs.mean)
	get.hpds <- get.hpdintervals(y, X, traits, row.ids = row.ids, combined.fit.mcmc, num.lv)
	out.fit$hpdintervals <- get.hpds
	if(calc.ics) {
		message("Calculating Information criteria...")
		get.ics <- get.measures(y = y, X = X, family = family, trial.size = complete.trial.size, row.eff = row.eff, row.ids = row.ids, num.lv = num.lv, fit.mcmc = combined.fit.mcmc)
		ics <- c(get.dic(jagsfit), get.ics$waic, get.ics$eaic, get.ics$ebic)
		names.ics <- c("Conditional DIC", "WAIC", "EAIC", "EBIC")
		if(get.ics$do.marglik.ics) {
			ics <- c(ics, get.ics$aic.median, get.ics$bic.median)
			names.ics  <- c(names.ics , "AIC at post. median", "BIC at post. median")
			}
		names(ics) <- names.ics; out.fit$ics <- ics
		}
			
	if(save.model) { out.fit$jags.model <- jagsfit }

	out.fit$call <- match.call()
	out.fit$n <- n; out.fit$p <- p
	out.fit$X <- X; out.fit$traits <- traits; out.fit$y <- y
	out.fit$row.eff <- row.eff; out.fit$row.ids <- row.ids

	out.fit$family <- family; if(all(family == "bernoulli")) out.fit$family <- rep("binomial",p) ## Switch it back to binomial for compatibility
	out.fit$num.lv <- num.lv
	out.fit$num.X <- num.X; out.fit$num.traits <- num.traits
	out.fit$which.traits <- which.traits
	out.fit$calc.ics <- calc.ics
	out.fit$trial.size <- complete.trial.size
	out.fit$prior.control <- prior.control	
	out.fit$num.ord.levels <- num.ord.levels
	out.fit$mcmc.control <- mcmc.control
	out.fit$format <- new.format 
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
	
	col.seq <- rep("black", length(x$hpdintervals$X.coefs[,covname,"lower"]))
	col.seq[x$hpdintervals$X.coefs[,covname,"lower"] < 0 & x$hpdintervals$X.coefs[,covname,"upper"] > 0] <- "grey"
	
	At.y <- rev(1:nrow(x$X.coefs.median)) ## So that the graphs plots in the same order as rownames of x$X.coefs.median

	if(est == "median")
		plot(x = x$X.coefs.median[,covname], y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = covname, xlim = c(min(x$hpdintervals$X.coefs[,covname,"lower"]), max(x$hpdintervals$X.coefs[,covname,"upper"])), pch = "x", ...)
	if(est == "mean")
		plot(x = x$X.coefs.mean[,covname], y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = covname, xlim = c(min(x$hpdintervals$X.coefs[,covname,"lower"]), max(x$hpdintervals$X.coefs[,covname,"upper"])), pch = "x", ...)

	segments(x0 = x$hpdintervals$X.coefs[,covname,"lower"], y0 = At.y, x1 = x$hpdintervals$X.coefs[,covname,"upper"], y1 = At.y, col = col.seq, ...)  
	abline(v=0, lty=3)
	
	if(is.null(labely)) { axis(2, at=At.y, labels = rownames(x$X.coefs.mean), las=1, ...) } 
	if(!is.null(labely)) {
		if(length(labely) == nrow(x$X.coefs.mean)) axis(2, at=At.y, labels=labely, las=1, ...) 
		if(length(labely) == 1) mtext(text = labely, side = 2, line = 3, las = 3, ...)
		} 

	if(exists("ssvs.indcoefs.mean", x)) {
		message("Posterior probabilities of inclusion for ", covname, ":")
		print(round(x$ssvs.indcoefs.mean[,covname],3))
		message()
		}
	}

		
lvsplot <- function(x, jitter = FALSE, biplot = TRUE, ind.spp = NULL, alpha = 0.5, main = NULL, est = "median", which.lvs = c(1,2), return.vals = FALSE, ...) {
 	if(x$num.lv == 0) stop("No latent variables to plot.")
	if(length(which.lvs) != 2) stop("which.lvs should be a vector of length 2, indicating which axes to plot. This argument is ignored if x$num.lv = 1.")
	if(x$num.lv > 2 & any(which.lvs > x$num.lv)) stop("Fewer latent variables than chosen by which.lvs.")
 
 	n <- nrow(x$lv.median); p <- nrow(x$lv.coefs.median)
 	if(!is.null(ind.spp)) { if(ind.spp > p) { ind.spp <- p } }
	if(biplot == TRUE & !is.null(ind.spp)) { 
		message("Only the first ", ind.spp, " ``most important'' latent variable coefficients included in biplot.") 
		}
	if(biplot == TRUE & is.null(ind.spp)) { 
		ind.spp <- p; message("All latent variable coefficients included in biplot.") 
		}

 	if(x$num.lv == 1) {
		choose.lvs <- x$lv.median; choose.lvs.coefs <- x$lv.coefs.median[,2]
		main <- "Plot of latent variable posterior medians"
		if(est == "mean") { choose.lvs <- x$lv.mean; choose.lvs.coefs <- x$lv.coefs.mean[,2] }

		if(!biplot) {
			if(is.null(main) & est == "median") { main = "Plot of latent variable posterior medians" }
			if(is.null(main) & est == "mean") { main = "Plot of latent variable posterior means" }
			plot(1:n, choose.lvs, xlab = "Row index", ylab = "Latent variable 1", main = main, type = "n", ...)
			if(!jitter) text(x = 1:n, y = x$lv.median, label = 1:n, ...)
			if(jitter) text(x = 1:n, y = jitter(x$lv.median), label = 1:n, ...)
			}

		if(biplot) {
			if(is.null(main) & est == "median") { main = "Biplot of latent variable posterior medians" }
			if(is.null(main) & est == "mean") { main = "Biplot of latent variable posterior means" }
			plot(1:n, choose.lvs, xlab = "Row index", ylab = "Latent variable 1", main = main, type = "n", ...)
			if(!jitter) text(x = 1:n, y = x$lv.median, label = 1:n, ...)
			if(jitter) text(x = 1:n, y = jitter(x$lv.median), label = 1:n, ...)
			text(choose.lv.coefs, label = rownames(x$lv.coefs.mean), col = "red", ...)	
			}
		}


 	if(x$num.lv > 1) {
   		testcov <- x$lv.median%*%t(x$lv.coefs.median[,2:(x$num.lv+1)])
		if(est == "mean") { testcov <- x$lv.mean%*%t(x$lv.coefs.mean[,2:(x$num.lv+1)]) }

		do.svd <- svd(testcov, x$num.lv, x$num.lv)   		
   		choose.lvs <- scale(do.svd$u*matrix(do.svd$d[1:x$num.lv]^alpha,nrow=x$n,ncol=x$num.lv,byrow=TRUE),center = TRUE, scale = FALSE)
   		choose.lv.coefs <- scale(do.svd$v*matrix(do.svd$d[1:x$num.lv]^(1-alpha),nrow=x$p,ncol=x$num.lv,byrow=TRUE), center = TRUE, scale = FALSE)
   		
		if(!biplot) {
			if(is.null(main) & est == "median") main <- "Plot of latent variable posterior medians"
			if(is.null(main) & est == "mean") main <- "Plot of latent variable posterior meas, axis"
			plot(choose.lvs, xlab = paste("Latent variable", which.lvs[1]), ylab = paste("Latent variable", which.lvs[2]), main = main, type = "n", ...)
			if(!jitter) text(choose.lvs[,which.lvs], label = 1:n, ...)
			if(jitter) text(jitter(choose.lvs[,which.lvs[1]]), jitter(choose.lvs[,which.lvs[2]]), label = 1:n, ...)
			}

		if(biplot) {
			if(is.null(main) & est == "median") main <- "Biplot of latent variable posterior medians"
			if(is.null(main) & est == "mean") main <- "Biplot of latent variable posterior means"
			largest.lnorms <- order(rowSums(choose.lv.coefs^2),decreasing=TRUE)[1:ind.spp] ## Open question as to whether largest.lnorms should be based on which.lvs only
			
			plot(rbind(choose.lvs,choose.lv.coefs)[,which.lvs], xlab = paste("Latent variable", which.lvs[1]), ylab = paste("Latent variable", which.lvs[2]), main = main, type = "n", xlim = 1.1*range(rbind(choose.lvs,choose.lv.coefs)[,which.lvs[1]]), ylim = 1.1*range(rbind(choose.lvs,choose.lv.coefs)[,which.lvs[2]]), ...)
			if(!jitter) text(choose.lvs[,which.lvs], label = 1:n, ...)
			if(jitter) text(jitter(choose.lvs[,which.lvs[1]]), jitter(choose.lvs[,which.lvs[2]]), label = 1:n, ...)
			text(choose.lv.coefs[largest.lnorms,which.lvs], label = rownames(x$lv.coefs.mean[largest.lnorms,]), col = "red", ...)	
			}
 		}	

 	out <- list(scaled.lvs = choose.lvs, scaled.lv.coefs = choose.lv.coefs)
 	if(return.vals) return(out)
 	}


print.boral <- function(x, ...) {
 	message("Call:")
 	print(x$call)
 	message()
 	message("Response matrix attributes\n \t# of rows: ", x$n, "\n\t# of columns: ", x$p) 
 	message("Model attributes\n \tColumn families used: ", unique(x$family), "\n\t# of latent variables: ", x$num.lv, "\n\tRow effects included (none/fixed/random)? ", x$row.eff, "\n") 
 	if(any(x$family == "binomial")) message("Trial sizes used (columns with binomial families): ", x$trial.size)
 	if(any(x$family == "ordinal")) message("Number of levels for ordinal data: ", x$num.ord.levels)
 	if(x$num.X > 0) message("Model matrix with ", x$num.X, " covariates also fitted\n")
 	if(x$num.traits > 0) message("Trait matrix with ", x$num.traits, " traits also included\n")
 	if(any(x$prior.control$ssvs.index > -1)) message("SSVS performed on covariates with indices ", x$prior.control$ssvs.index, "\n")
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
 	
 	if(any(x$family == "ordinal")) { message("Proportional odds (Cumulative probit) cutoffs"); print(x$cutoffs); message() }
 	if(any(x$family == "tweedie")) { message("Tweedie power parameter"); print(x$powerparam); message() }
 
 	if(x$calc.ics) {
 		message("Some Information Criteria:")
		print(x$ics)
 		}

	if(!is.null(x$ssvs.indcoefs.prob)) { message("SSVS probabilities on individual coefficients"); print(x$ssvs.indcoefs.prob); message() }
	if(!is.null(x$ssvs.gpcoefs.prob)) { message("SSVS probabilities on groups of coefficients"); print(x$ssvs.gpcoefs.prob); message() }			
	}	
		
		
summary.boral <- function(object, est = "median", ...) {
	if(est == "median") {
 		gather.output <- list(call = object$call, coefficients = round(object$lv.coefs.median,3))
 		if(object$row.eff != "none") {
			for(k in 1:ncol(object$row.ids)) gather.output$row.coefficients[[k]] = round(object$row.coefs[[k]]$median,3)
			}
 		if(object$num.X > 0) gather.output$X.coefficients = round(object$X.coefs.median,3)
 		if(object$num.traits > 0) gather.output$traits.coefficients = round(object$traits.coefs.median,3)
 		if(any(object$family == "ordinal")) gather.output$cutoffs = round(object$cutoffs.median,3)
 		if(any(object$family == "tweedie")) gather.output$powerparam = round(object$powerparam.median,3)
 		if(!is.null(object$X.multinom.coefs.median)) gather.output$X.multinom.coefficients = round(object$X.multinom.coefs.median,3) }
 
 	if(est == "mean") {
 		gather.output <- list(call = object$call, coefficients = round(object$lv.coefs.mean,3))
 		if(object$row.eff != "none") {
			for(k in 1:ncol(object$row.ids)) gather.output$row.coefficients[[k]] = round(object$row.coefs[[k]]$mean,3)
			}
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

	if(object$calc.ics) gather.output$ics <- round(object$ics,3)
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
 	qqline(y = get.ds.res2[is.finite(get.ds.res2)], ...)
 	
 	palette("default")
 	#par(.pardefault) 	
 	}

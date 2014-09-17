##########
## Auxilary functions
##########
## Calculate conditional logl
## loglik = sum_{i=1}^n sum_{j=1}^s \log( f(y_ij|b_i) )
#  lv.coefs = matrix(fit.mcmc[t,grep("all.params",colnames(fit.mcmc))],nrow=p)
#  X.coefs = get.X.coefs
#  row.coefs = get.row.coefs
#  lv = matrix(fit.mcmc[t,grep("lvs", colnames(fit.mcmc))],nrow=n)
#  cutoffs = get.cutoffs
#, X.multinom.coefs = NULL
calc.condlogLik <- function(y, X = NULL, family, trial.size = NULL, lv.coefs, X.coefs = NULL, row.coefs = NULL, lv, cutoffs = NULL, powerparam = NULL) {
	if(is.null(lv) | is.null(lv.coefs)) stop("lv and lv.coefs must be given. Please use calc.loglik.lv0 to calculate likelihood in boral models with no latent variables.")
	if(length(family) != ncol(y) & length(family) != 1) { stop("Number of elements in family is either 1 or equal to # of columns in y") }
	if(length(family) == 1) complete.family <- rep(family,ncol(y))
	if(length(family) > 1) complete.family <- family

	if(any(complete.family == "ordinal") & is.null(cutoffs)) stop("Ordinal data requires cutoffs to be supplied") 
	if(any(family == "tweedie") & (powerparam < 1 || powerparam > 2)) stop("Common power parameter for tweedie must be between 1 and 2.") 
	#if(any(complete.family == "multinom") & is.null(X.multinom.coefs)) stop("Multinomial data requires X.multinom.coefs to be supplied.") 

	n <- nrow(y); p <- ncol(y); num.lv <- ncol(lv)
	loglik <- 0; loglik.comp <- matrix(NA,n,p) 
	if(is.null(row.coefs)) row.coefs <- rep(0,n)

	if(any(family == "binomial") & is.null(trial.size)) 
		stop("trial.size needs to be specified if any columns are binomially distributed; can either be a single element or a vector equal to the # of columns in y. The latter will assume the specified trial size for all rows labelled binomial in the family argument.")
	if(length(trial.size) == 1) complete.trial.size <- rep(trial.size,ncol(y))
	if(length(trial.size) > 1) complete.trial.size <- trial.size

	all.etas <- cbind(1,lv)%*%t(as.matrix(lv.coefs[,1:(num.lv+1)]))
	if(!is.null(X.coefs)) all.etas <- all.etas + as.matrix(X)%*%t(X.coefs)
	## Assumes the columns in X.coefs corresponding to multinomial are set to 0

	index.multinom.cols <- which(complete.family == "multinom")
	for(j in 1:p) {
		species.etas <- all.etas[,j] + row.coefs
		
		if(complete.family[j] == "binomial") 
			loglik.comp[,j] <- dbinom(as.vector(unlist(y[,j])), complete.trial.size[j], exp(species.etas)/(1+exp(species.etas)), log = T)
		if(complete.family[j] == "poisson") 
			loglik.comp[,j] <- dpois(as.vector(unlist(y[,j])), exp(species.etas), log = T)
		if(complete.family[j] == "negative.binomial") 
			loglik.comp[,j] <- dnbinom(as.vector(unlist(y[,j])), size=1/(lv.coefs[j,ncol(lv.coefs)]+1e-6), mu=exp(species.etas), log = T)
		if(complete.family[j] == "exponential") 
			loglik.comp[,j] <- dexp(as.vector(unlist(y[,j])), 1/exp(species.etas), log = T)
		if(complete.family[j] == "gamma") 
			loglik.comp[,j] <- dgamma(as.vector(unlist(y[,j])), shape=exp(species.etas)*lv.coefs[j,ncol(lv.coefs)], rate = lv.coefs[j,ncol(lv.coefs)], log = T)
		if(complete.family[j] == "beta") 
			loglik.comp[,j] <- dbeta(as.vector(unlist(y[,j])), lv.coefs[j,ncol(lv.coefs)]*exp(species.etas)/(1+exp(species.etas)), lv.coefs[j,ncol(lv.coefs)]*(1-exp(species.etas)/(1+exp(species.etas))),log = T)
		if(complete.family[j] == "normal") 
			loglik.comp[,j] <- dnorm(as.vector(unlist(y[,j])), mean=species.etas, sd=sqrt(lv.coefs[j,ncol(lv.coefs)]+1e-6), log = T)
		if(complete.family[j] == "lnormal") 
			loglik.comp[,j] <- dlnorm(as.vector(unlist(y[,j])), meanlog=species.etas, sdlog=sqrt(lv.coefs[j,ncol(lv.coefs)]+1e-6), log = T)
			if(complete.family[j] == "tweedie") 
			loglik.comp[,j] <- dTweedie(as.vector(unlist(y[,j])), mu = exp(species.etas), phi = lv.coefs[j,ncol(lv.coefs)]+1e-6, p = powerparam, LOG = TRUE) 
		if(complete.family[j] == "ordinal") { 
			get.probs <- ordinal.conversion.spp(n, lv, lv.coefs[j,], num.lv, row.coefs, X, X.coefs[j,], cutoffs); 	
			for(i in 1:n) { loglik.comp[i,j] <- log(get.probs[i,as.vector(y[i,j])]+1e-5) } }	
# 		if(complete.family[j] == "multinom") { 
# 			if(!is.null(X.multinom.coefs)) spp.etas <- matrix(rep(species.etas,dim(X.multinom.coefs)[3]),nrow=n) + as.matrix(X)%*%X.multinom.coefs[which(index.multinom.cols == j),,]
# 			get.probs <- exp(spp.etas)/apply(exp(spp.etas),1,sum)
# 			for(i in 1:n) { loglik.comp[i,j] <- log(get.probs[as.vector(y[i,j])]+1e-5) } }	
		} 

	return(list(logLik = sum(loglik.comp), logLik.comp = apply(loglik.comp,1,sum))) }

## Calculate logl for models with no latent variables
## Conditional and marginal log-likelihood are the same in such case
## lv.coefs still need to be provided though at it contains the row effects, spp effects, and species-specific dispersion parameters
calc.logLik.lv0 <- function(y, X = NULL, family, trial.size = NULL, lv.coefs, X.coefs = NULL, row.coefs = NULL, cutoffs = NULL, powerparam = NULL) {
	if(is.null(lv.coefs)) stop("lv.coefs must be given, as it contains the column-specific intercepts.")
	if(length(family) != ncol(y) & length(family) != 1) { stop("Number of elements in family is either 1 or equal to # of columns in y") }
	if(length(family) == 1) complete.family <- rep(family,ncol(y))
	if(length(family) > 1) complete.family <- family

	n <- nrow(y); p <- ncol(y); 
	logl <- 0; logl.comp <- matrix(0,n,p)
	if(is.null(row.coefs)) row.coefs <- rep(0,n)
	if(any(complete.family == "ordinal") & is.null(cutoffs)) stop("Ordinal data requires cutoffs to be supplied") 
	if(any(family == "tweedie") & (powerparam < 1 || powerparam > 2)) stop("Common power parameter for tweedie must be between 1 and 2.") 

	if(any(family == "binomial") & is.null(trial.size)) 
		stop("trial.size needs to be specified if any columns are binomially distributed; can either be a single element or a vector equal to the # of columns in y. 
			The latter will assume the specified trial size for all rows labelled binomial in the family argument.")
	if(length(trial.size) == 1) complete.trial.size <- rep(trial.size,ncol(y))
	if(length(trial.size) > 1) complete.trial.size <- trial.size
	#if(any(complete.family == "multinom") & is.null(X.multinom.coefs)) stop("Multinomial data requires X.multinom.coefs to be supplied.") 

	index.multinom.cols <- which(complete.family == "multinom")
	for(j in 1:p) {
		eta <- row.coefs + lv.coefs[j,1]
		if(!is.null(X.coefs)) eta <- eta + as.matrix(X)%*%X.coefs[j,]
		## Assumes the columns in X.coefs corresponding to multinomial are set to 0

		if(complete.family[j] == "poisson") logl.comp[,j] <- (dpois(as.vector(unlist(y[,j])), lambda = exp(eta), log = T))
		if(complete.family[j] == "binomial") logl.comp[,j] <- (dbinom(as.vector(unlist(y[,j])), complete.trial.size[j], prob = exp(eta)/(1+exp(eta)), log = T))
		if(complete.family[j] == "negative.binomial") logl.comp[,j] <- (dnbinom(as.vector(unlist(y[,j])), mu = exp(eta), size = 1/lv.coefs[j,2], log = T)) 
		if(complete.family[j] == "exponential") logl.comp[,j] <- (dexp(as.vector(unlist(y[,j])), rate = 1/exp(eta), log = T))
		if(complete.family[j] == "dgamma") logl.comp[,j] <- (dgamma(as.vector(unlist(y[,j])), shape=exp(eta)*lv.coefs[j,2], rate = lv.coefs[j,2], log = T))
		if(complete.family[j] == "beta") logl.comp[,j] <- (dbeta(as.vector(unlist(y[,j])), lv.coefs[j,2]*exp(eta)/(1+exp(eta)),lv.coefs[j,2]*(1-exp(eta)/(1+exp(eta))), log = T))
		if(complete.family[j] == "normal") logl.comp[,j] <- (dnorm(as.vector(unlist(y[,j])), mean = eta, sd = sqrt(lv.coefs[j,2]), log = T)) 
		if(complete.family[j] == "lnormal") logl.comp[,j] <- (dlnorm(as.vector(unlist(y[,j])), meanlog = eta, sdlog = sqrt(lv.coefs[j,2]), log = T)) 
		if(complete.family[j] == "tweedie") logl.comp[,j] <- (dTweedie(as.vector(unlist(y[,j])), mu = exp(eta), phi = lv.coefs[j,2], p = powerparam, LOG = T)) 
		if(complete.family[j] == "ordinal") { 
			get.probs <- ordinal.conversion.spp(lv = NULL, lv.coefs[j,], num.lv = 0, row.coefs, X, X.coefs[j,], cutoffs); 
			for(i in 1:n) { loglik.comp[i,j] <- log(get.probs[i,as.vector(y[i,j])]+1e-5) } }	
# 		if(complete.family[j] == "multinom") {
# 			eta <- matrix(rep(eta,dim(X.multinom.coefs)[3]),nrow=n) + as.matrix(X)%*%X.multinom.coefs[which(index.multinom.cols == j),,]
# 			get.probs <- exp(eta)/apply(exp(eta),1,sum)
# 			for(i in 1:n) { logl.comp[i,j] <- log(get.probs[as.vector(y[i,j])]) } }
		}

	return(list(logLik = sum(logl.comp), logLik.row.comp = apply(logl.comp,1,sum), logLik.col.comp = apply(logl.comp,2,sum))) }

## Calculate marginal logl
## loglik = sum_{i=1}^n \log( \int \prod_{j=1}^s f(y_ij|b_i) f(b_i) db_i )
# lv.coefs = coef.mat; X.coefs = get.X.coefs; row.coefs = get.row.coefs; X.mc = NULL; cutoffs = get.cutoffs; X.multinom.coefs = NULL
calc.marglogLik <- function(y, X = NULL, family, trial.size = NULL, lv.coefs, X.coefs = NULL, row.coefs = NULL, num.lv, X.mc = NULL, cutoffs = NULL, powerparam = NULL) {
	if(num.lv == 0) stop("Please use calc.loglik.lv0 to calculate likelihood in borla models with no latent variables.")
	if(is.null(lv.coefs)) stop("lv.coefs must be given. Please use calc.loglik.lv0 to calculate likelihood in boral models with no latent variables.")
	
	n <- nrow(y); p <- ncol(y); 
	loglik <- 0; loglik.comp <- numeric(n)
	if(length(family) != ncol(y) & length(family) != 1) { stop("Number of elements in family is either 1 or equal to # of columns in y") }
	if(length(family) == 1) complete.family <- rep(family,ncol(y))
	if(length(family) > 1) complete.family <- family

	if(any(family == "binomial") & is.null(trial.size)) 
		stop("trial.size needs to be specified if any columns are binomially distributed; can either be a single element or a vector equal to the # of columns in y. 
			The latter will assume the specified trial size for all rows labelled binomial in the family argument.")
	if(length(trial.size) == 1) complete.trial.size <- rep(trial.size,ncol(y))
	if(length(trial.size) > 1) complete.trial.size <- trial.size

	if(is.null(X.mc)) { X.mc <- cbind(1, rmvnorm(5000, rep(0,num.lv))) }  	
	if(is.null(row.coefs)) row.coefs <- rep(0,n)
	if(any(complete.family == "ordinal") & is.null(cutoffs)) stop("Ordinal data requires cutoffs to be supplied") 
	#if(any(complete.family == "multinom") & is.null(X.multinom.coefs)) stop("Multinomial data requires X.multinom.coefs to be supplied.") 
	if(any(family == "tweedie") & (powerparam < 1 || powerparam > 2)) stop("Common power parameter for tweedie must be between 1 and 2.") 

	## Internal function - Given the coefficients and cutoffs, return the multinomial probabilities for proportional odds regression for element of y
	ordinal.conversion.special <- function(X.mc, lv.coefs.j, num.lv, row.coefs.i = NULL, X.i = NULL, X.coefs.j = NULL, cutoffs) {
		etas <- matrix(NA,nrow(X.mc),length(cutoffs)) ## dim(X.mc) x num.ord.levels - 1
		for(k in 1:length(cutoffs)) {
			etas[,k] <- X.mc%*%c(cutoffs[k],-lv.coefs.j[2:(num.lv+1)])-row.coefs.i-lv.coefs.j[1] ## Don't forget the negative sign!
			if(!is.null(X.coefs.j)) etas[,k] <- etas[,k]-t(as.matrix(X.i))%*%X.coefs.j
			}
		probs <- matrix(NA,nrow(X.mc),length(cutoffs)+1) ## dim(X.mc) x num.ord.levels
		probs[,1] <- exp(etas[,1])/(1+exp(etas[,1]))
		for(k in 2:ncol(etas)) { probs[,k] <- exp(etas[,k])/(1+exp(etas[,k])) - exp(etas[,k-1])/(1+exp(etas[,k-1])) }		
		probs[,length(cutoffs)+1] <- 1-exp(etas[,length(cutoffs)])/(1+exp(etas[,length(cutoffs)]))		
		rm(etas); probs
		}	
	
	index.multinom.cols <- which(complete.family == "multinom")
	for(i in 1:n) {
		spp.f <- matrix(NA,nrow=nrow(X.mc),ncol=p)
		spp.att.eta <- X.mc%*%t(lv.coefs[,1:(num.lv+1)]) + row.coefs[i]
		if(!is.null(X.coefs)) spp.att.eta <- spp.att.eta + matrix(t(as.matrix(X[i,]))%*%t(X.coefs),nrow=nrow(X.mc),ncol=p,byrow=T) 
		## Assumes the columns in X.coefs corresponding to multinomial are set to 0
				
		for(j in 1:p) {
			if(complete.family[j] == "binomial") {
				spp.f[,j] <- dbinom(rep(as.vector(y[i,j]),nrow(X.mc)), complete.trial.size[j], exp(spp.att.eta[,j])/(1+exp(spp.att.eta[,j]))) }
			if(complete.family[j] == "poisson") { 
				spp.f[,j] <- dpois(rep(as.vector(y[i,j]),nrow(X.mc)), exp(spp.att.eta[,j])) } 
			if(complete.family[j] == "negative.binomial") { 
				spp.f[,j] <- dnbinom(rep(as.vector(y[i,j]),nrow(X.mc)), size=1/(lv.coefs[j,ncol(lv.coefs)]+1e-6), mu=exp(spp.att.eta[,j])) }
			if(complete.family[j] == "exponential") { 
				spp.f[,j] <- dexp(rep(as.vector(y[i,j]),nrow(X.mc)), rate = 1/exp(spp.att.eta[,j])) } 
			if(complete.family[j] == "gamma") { 
				spp.f[,j] <- dgamma(rep(as.vector(y[i,j]),nrow(X.mc)), shape = exp(spp.att.eta[,j])*lv.coefs[j,ncol(lv.coefs)], rate = lv.coefs[j,ncol(lv.coefs)]) } 
			if(complete.family[j] == "beta") {
				spp.f[,j] <- dbeta(rep(as.vector(y[i,j]),nrow(X.mc)), lv.coefs[j,ncol(lv.coefs)]*exp(spp.att.eta[,j])/(1+exp(spp.att.eta[,j])), lv.coefs[j,ncol(lv.coefs)]*(1-exp(spp.att.eta[,j])/(1+exp(spp.att.eta[,j])))) }
			if(complete.family[j] == "normal") { 
				spp.f[,j] <- dnorm(rep(as.vector(y[i,j]),nrow(X.mc)), mean = spp.att.eta[,j], sd = sqrt(lv.coefs[j,ncol(lv.coefs)])) }
			if(complete.family[j] == "lnormal") { 
				spp.f[,j] <- dlnorm(rep(as.vector(y[i,j]),nrow(X.mc)), meanlog = spp.att.eta[,j], sdlog = sqrt(lv.coefs[j,ncol(lv.coefs)])) }
			if(complete.family[j] == "tweedie") { 
				spp.f[,j] <- dTweedie(rep(as.vector(y[i,j]),nrow(X.mc)), mu = exp(spp.att.eta[,j]), phi = lv.coefs[j,ncol(lv.coefs)]+1e-6, p = powerparam, LOG = FALSE); spp.f[,j][which(spp.f[,j] == 0)] <- 1 }
			if(complete.family[j] == "ordinal") { 
				get.probs <- ordinal.conversion.special(X.mc, lv.coefs.j = lv.coefs[j,], num.lv, row.coefs.i = row.coefs[i], X.i = X[i,], X.coefs.j = X.coefs[j,], cutoffs = cutoffs); spp.f[,j] <- get.probs[,as.vector(y[i,j])]+1e-5 }
# 			if(complete.family[j] == "multinom") {
# 				num.multinom.levels <- dim(X.multinom.coefs)[3]
# 				if(!is.null(X.multinom.coefs)) { 
# 					spp.att.eta2 <- spp.att.eta[,j] + matrix(t(as.matrix(X[i,]))%*%X.multinom.coefs[which(index.multinom.cols == j),,],nrow(X.mc),num.multinom.levels,byrow=T) }
# 				get.probs <- exp(spp.att.eta2)/apply(exp(spp.att.eta2),1,sum)
# 				spp.f[,j] <- get.probs[,as.vector(y[i,j])]+1e-5 }
			}

		spp.f[!is.finite(spp.f)] = 1; spp.f <- matrix(spp.f,nrow=nrow(X.mc),byrow=F)
		Q <- mean(apply(spp.f,1,prod))           
		loglik <- loglik + log(Q); loglik.comp[i] <- log(Q)
		}

	return(list(logLik = loglik, logLik.comp = loglik.comp)) }

create.life <- function(true.lv = NULL, lv.coefs, X = NULL, X.coefs = NULL, family, row.coefs = NULL, trial.size = NULL, cutoffs = NULL, powerparam = NULL, manual.dim = NULL) {
	num.lv <- max(ncol(true.lv),0)
	n <- max(nrow(true.lv), nrow(X), length(row.coefs)) 
	s <- max(nrow(lv.coefs),nrow(X.coefs),length(cutoffs))
	if(is.null(dim(lv.coefs))) { lv.coefs <- matrix(lv.coefs,ncol=1) }
	
	if((is.null(n) | is.null(s)) & is.null(manual.dim)) 
		stop("Sorry, but boral cannot determine the number of rows and columns for the response matrix. Please supply manual.dim as vector containing n and p.") 
	if((is.null(n) | is.null(s)) & !is.null(manual.dim)) { n <- manual.dim[1]; s <- manual.dim[2] }  

 	if(any(family != "multinom")) { if((is.null(X) & !is.null(X.coefs)) | (!is.null(X) & is.null(X.coefs))) 
 		stop("If there are covariates to be included, then both X and X.coefs must be specified.") }
	if(length(family) != s & length(family) != 1) 
		stop("Number of elements in family must be either 1 or equal to # of rows in lv.coefs or X.coefs or the second number in manual.dim") 
	if(length(family) == 1) family <- rep(family, s)
	if(!all(family %in% c("negative.binomial","poisson","binomial","normal","lnormal","tweedie","ordinal","exponential","beta"))) 
		stop("One of the elements in family is not compatible with current version of boral...sorry!") 
	
	if(any(family == "binomial") & is.null(trial.size)) 
		stop("trial.size needs to be specified if any columns are binomially distributed; can either be a single element or a vector equal to the # of elements in family. The latter will assume the specified trial size for all rows labelled binomial in the family argument.")
	if(length(trial.size) == 1) trial.size <- rep(trial.size,s)			
	
	if(any(family == "ordinal") & is.null(cutoffs)) 
		stop("cutoffs (an ascending vector of intercepts for proportional odds regression) must be supplied if any columns are ordinal data.")
	if(!is.null(cutoffs)) { 
		num.ord.levels <- length(cutoffs) + 1; cutoffs <- sort(cutoffs); print("Sorting cutoffs...just in case") }

	if(any(family == "tweedie") & is.null(powerparam)) 
		stop("Common powerparam must be supplied if any columns are tweedie data (Var = dispersion*mu^powerparam)")
	
# 	if(any(family == "multinom")) {
# 		if((is.null(X) & !is.null(X.multinom.coefs)) | (!is.null(X) & is.null(X.multinom.coefs))) 
# 		stop("If any columns are multinomially distributed and if there are covariates to be included, X.multinom.coefs must be specified") }
# 	if(!is.null(X.multinom.coefs)) { if(sum(family == "multinom") != dim(X.multinom.coefs)[1])
# 		stop("X.multinom.coefs must have dimensions (# of multinomial columns)*(# of columns in X)*(number of levels)") 
# 		if(!is.null(X.coefs)) print("For many columns set to multinomial, please ensure corresponding rows X.coefs are set entirely to zero. Thanks.") }
# 	if(!is.null(X.multinom.coefs)) { 
# 		X.multinom.coefs[,,1] <- 0; index.multinom.cols <- which(family == "multinom") 
# 		print("Setting X.multinom.coefs[,,1] = 0...just in case.") }
		
	if(num.lv > 3) warnings("We won't stop you, but please consider if you really want more than three latent variables in the model.")

	sim.y <- matrix(NA,n,s)
	if(!is.null(true.lv)) eta <- cbind(1,true.lv)%*%t(lv.coefs[,1:(num.lv+1)])
	if(is.null(true.lv)) eta <- rep(1,n)%*%t(lv.coefs[,1:(num.lv+1)])
	if(!is.null(X.coefs)) { eta <- eta + as.matrix(X)%*%t(X.coefs) } ## Assumes the X.coefs corresponding to multinomial columns are all 0
	for(i in 1:n) { if(!is.null(row.coefs)) eta[i,] <- eta[i,] + row.coefs[i] }
		
	for(j in 1:s) {
		if(family[j] == "binomial") sim.y[,j] <- rbinom(n,size=trial.size[j],prob=exp(eta[,j])/(1+exp(eta[,j])))
		if(family[j] == "poisson") sim.y[,j] <- rpois(n,lambda=exp(eta[,j]))
		if(family[j] == "negative.binomial") sim.y[,j] <- rnbinom(n,mu=exp(eta[,j]),size=1/(lv.coefs[j,ncol(lv.coefs)]+1e-6))
		if(family[j] == "exponential") sim.y[,j] <- rexp(n,rate=1/exp(eta[,j]))
		if(family[j] == "gamma") sim.y[,j] <- rgamma(n,shape=exp(eta[,j])*lv.coefs[j,ncol(lv.coefs)], rate=lv.coefs[j,ncol(lv.coefs)])
		if(family[j] == "beta") sim.y[,j] <- rbeta(n,shape1=lv.coefs[j,ncol(lv.coefs)]*exp(eta[,j])/(1+exp(eta[,j])),shape2=lv.coefs[j,ncol(lv.coefs)]*(1-exp(eta[,j])/(1+exp(eta[,j]))))
		if(family[j] == "normal") sim.y[,j] <- rnorm(n,mean=eta[,j],sd=sqrt(lv.coefs[j,ncol(lv.coefs)])) 
		if(family[j] == "lnormal") sim.y[,j] <- rlnorm(n,meanlog=eta[,j],sdlog=sqrt(lv.coefs[j,ncol(lv.coefs)])) 
		if(family[j] == "tweedie") sim.y[,j] <- rTweedie(n,mu=exp(eta[,j]),phi=lv.coefs[j,ncol(lv.coefs)], p = powerparam) 
		if(family[j] == "ordinal") {
			get.probs <- ordinal.conversion.spp(n, lv = true.lv, lv.coefs[j,], num.lv, row.coefs, X, X.coefs[j,], cutoffs); 	
			for(i in 1:n) { sim.y[i,j] <- sample(1:num.ord.levels,1,prob=get.probs[i,]) } }
# 		if(family[j] == "multinom") {
# 			if(!is.null(X.multinom.coefs)) { neweta <- eta[,j] + as.matrix(X)%*%X.multinom.coefs[which(index.multinom.cols == j),,] }
# 			neweta[is.na(neweta)] <- -Inf
# 			get.probs <- exp(neweta)/apply(exp(neweta),1,sum)
# 			for(i in 1:n) { sim.y[i,j] <- sample(1:dim(X.multinom.coefs)[3],1,prob=get.probs[i,]) } }
		}

	return(sim.y)	
	}
	
## Dunn-Smyth residuals
## Create a confusion matrix for ordinal and multinomial data
ds.residuals <- function(object, est = "median") {  
	n <- object$n; p <- object$p; 
	num.lv <- object$num.lv; num.ord.levels <- object$num.ord.levels; 
	X <- object$X; y <- object$y
	mus <- fitted.boral(object, X, est = est)
	if(est == "median") phis <- object$lv.coefs.median[,num.lv+2]
	if(est == "mean") phis <- object$lv.coefs.mean[,num.lv+2]

	if(any(object$family == "ordinal")) {
		print("One or more columns of y have ordinal responses. Constructing a single table of agreement for these.")
		true.resp <- as.matrix(y[,which(object$family == "ordinal")])
		pred.resp <- matrix(NA,n,ncol(true.resp)) }
# 	if(any(object$family == "multinom")) {
# 		print("One or more columns of y have multinomial responses. Constructing a single table of agreement for these.")
# 		true.multinom.resp <- as.matrix(y[,which(object$family == "multinom")])
# 		pred.multinom.resp <- matrix(NA,n,ncol(true.multinom.resp)) }
	if(any(object$family == "tweedie")) {
		if(est == "median") powerparam <- object$powerparam.median
		if(est == "mean") powerparam <- object$powerparam.mean }

	ds.res.out <- matrix(NA,n,p)
	rownames(ds.res.out) <- rownames(y); colnames(ds.res.out) <- colnames(y)
	for(i in 1:n) { for(j in 1:p) {
		if(object$family[j] == "poisson") { 
			a <- ppois(as.vector(unlist(y[i,j]))-1, mus$out[i,j]); b <- ppois(as.vector(unlist(y[i,j])), mus$out[i,j]); 		
			u <- runif(n = 1, min = a, max = b); ds.res.out[i,j] <- qnorm(u) }
		if(object$family[j] == "negative.binomial") { 
			a <- pnbinom(as.vector(unlist(y[i,j]))-1, mu=mus$out[i,j], size=1/phis[j]); b <- pnbinom(as.vector(unlist(y[i,j])), mu=mus$out[i,j], size=1/phis[j])
			u <- runif(n = 1, min = a, max = b); ds.res.out[i,j] <- qnorm(u) }
		if(object$family[j] == "binomial") { 
			a <- pbinom(as.vector(unlist(y[i,j]))-1, object$trial.size[j], mus$out[i,j]); b <- pbinom(as.vector(unlist(y[i,j])), object$trial.size[j], mus$out[i,j])
			u <- runif(n = 1, min = a, max = b); ds.res.out[i,j] <- qnorm(u) }
		if(object$family[j] == "exponential") { 
			a <- pexp(as.vector(unlist(y[i,j])), rate=1/mus$out[i,j]); ds.res.out[i,j] <- qnorm(a) }
		if(object$family[j] == "gamma") { 
			a <- pgamma(as.vector(unlist(y[i,j])), shape=mus$out[i,j]*phis[j], rate=phis[j]); ds.res.out[i,j] <- qnorm(a) }
		if(object$family[j] == "beta") { 
			a <- pbeta(as.vector(unlist(y[i,j])), shape1=phis[j]*mus$out[i,j], shape2=phis[j]*(1-mus$out[i,j])); ds.res.out[i,j] <- qnorm(a) }
		if(object$family[j] == "normal") { 
			a <- pnorm(as.vector(unlist(y[i,j])), mus$out[i,j], sqrt(phis[j])); ds.res.out[i,j] <- qnorm(a) }
# 			X2 <- cbind(1,X); hatmat <- X2%*%solve(t(X2)%*%X2)%*%t(X2)
# 			ds.res.out[i,j] <- (y[i,j]-mus$out[i,j])/(sqrt(phis[j])*sqrt(1-hatmat[i,i])) }
		if(object$family[j] == "lnormal") { 
			a <- plnorm(as.vector(unlist(y[i,j])), log(mus$out[i,j]), sqrt(phis[j])); ds.res.out[i,j] <- qnorm(a) }
		if(object$family[j] == "tweedie") { 
			a <- pTweedie(as.vector(unlist(y[i,j])), mu = mus$out[i,j], phi = phis[j], p = powerparam); ds.res.out[i,j] <- qnorm(a) }

		if(object$family[j] == "ordinal") { ## get max predicted probability
			pred.resp[i,which(object$family == "ordinal")==j] <- which.max(mus$ordinal.probs[i,j,]) }
# 		if(object$family[j] == "multinom") { ## get max predicted probability
# 			pred.resp[i,which(object$family == "multinom")==j] <- which.max(mus$multinom.probs[i,j,]) }
 		} }

	if(sum(object$family == "ordinal") > 0) { agree.tab <- table(as.vector(pred.resp), as.vector(true.resp)); } else { agree.tab <- NULL }
	#if(sum(object$family == "multinom") > 0) { agree.multinom.tab <- table(as.vector(pred.multinom.resp), as.vector(true.multinom.resp)); }	else { agree.multinom.tab <- NULL }
	
	return(list(agree.ordinal = agree.tab, residuals = ds.res.out))
	}

## Fitted values
## For ordinal and multinomial data, returns a matrix of probabilities for each vector of rows
fitted.boral <- function(object, est = "median",...) {
	n <- object$n; p <- object$p; num.lv <- object$num.lv; 
	X <- object$X; y <- object$y
	fitted.out <- matrix(NA,n,p)
	rownames(fitted.out) <- rownames(y); colnames(fitted.out) <- colnames(y)
	if(any(object$family == "ordinal")) { 
		fitted.ordinal.probs <- array(NA,dim=c(n,p,object$num.ord.levels)) 
		dimnames(fitted.ordinal.probs) <- list(r = rownames(y), c = colnames(y), levels = 1:object$num.ord.levels) } 
	else { fitted.ordinal.probs <- NULL }
# 	if(any(object$family == "multinom")) { 
# 		fitted.multinom.probs <- array(NA,dim=c(n,p,object$num.multinom.levels))
# 		dimnames(fitted.multinom.probs) <- list(r = rownames(y), c = colnames(y), levels = 1:object$num.multinom.levels) } 
# 	else { fitted.multinom.probs <- NULL }

	if(is.null(object$lv.median)) { eta <- matrix(1,n,1)%*%t(object$lv.coefs.median[,1:(num.lv+1)]) }
	if(!is.null(object$lv.median)) { eta <- cbind(1,object$lv.median)%*%t(object$lv.coefs.median[,1:(num.lv+1)])	}
	if(!is.null(object$X.coefs.median)) { eta <- eta + as.matrix(X)%*%t(object$X.coefs.median) }
	if(est == "mean") {
		if(is.null(object$lv.mean)) { eta <- matrix(1,n,1)%*%t(object$lv.coefs.mean[,1:(num.lv+1)]) }
		if(!is.null(object$lv.mean)) { eta <- cbind(1,object$lv.mean)%*%t(object$lv.coefs.mean[,1:(num.lv+1)]) } 
		if(!is.null(object$X.coefs.mean)) { eta <- eta + as.matrix(X)%*%t(object$X.coefs.mean) } }

	for(i in 1:n) {
		if(est == "median") { if(!is.null(object$row.coefs.median)) eta[i,] <- eta[i,] + object$row.coefs.median[i] }
		if(est == "mean") { if(!is.null(object$row.coefs.mean)) eta[i,] <- eta[i,] + object$row.coefs.mean[i] } }
	
	index.multinom.cols <- which(object$family == "multinom")
	for(j in 1:p) {
		if(object$family[j] %in% c("binomial","beta")) fitted.out[,j] <- exp(eta[,j])/(1+exp(eta[,j]))
		if(object$family[j] %in% c("poisson","lnormal","negative.binomial","tweedie","exponential","gamma")) fitted.out[,j] <- exp(eta[,j])
		if(object$family[j] == "normal") fitted.out[,j] <- (eta[,j]) 
# 		if(object$family[j] == "multinom") {
# 			if(est == "median") { if(!is.null(object$X.multinom.coefs.median)) eta2 <- eta[,j] + as.matrix(X)%*%object$X.multinom.coefs.median[which(index.multinom.cols == j),,] }
# 			if(est == "mean") { if(!is.null(object$X.multinom.coefs.mean)) eta2 <- eta[,j] + as.matrix(X)%*%object$X.multinom.coefs.mean[which(index.multinom.cols == j),,] }
# 			get.probs <- exp(eta2)/apply(exp(eta2),1,sum)	
# 			fitted.multinom.probs[,j,] <- get.probs
# 			}

		if(object$family[j] == "ordinal") {
			fitted.out[,j] <- NA
				if(est == "median")
					fitted.ordinal.probs[,j,] <- ordinal.conversion.spp(n, lv = object$lv.median, lv.coefs.j = object$lv.coefs.median[j,], num.lv = num.lv, row.coefs = object$row.coefs.median, X = X, X.coefs.j = object$X.coefs.median[j,], cutoffs = object$cutoffs.median)
				if(est == "mean")
					fitted.ordinal.probs[,j,] <- ordinal.conversion.spp(n, lv = object$lv.mean, lv.coefs.j = object$lv.coefs.mean[j,], num.lv = num.lv, row.coefs = object$row.coefs.mean, X = X, X.coefs.j = object$X.coefs.mean[j,], cutoffs = object$cutoffs.mean)
			}
		}	

	if(all(object$family == "ordinal") | all(object$family == "multinom")) fitted.out <- NULL
	return(list(ordinal.probs = fitted.ordinal.probs, out = fitted.out))
	}

## Calculates DIC based on the conditional log-likelihood
## May not work well given the results of Millar (2009)
get.dic <- function(jagsfit) { 
	jagsfit$BUGSoutput$DIC
	}
	
get.hpdintervals <- function(y, X = NULL, fit.mcmc, num.lv, prob = 0.95) {
	n <- nrow(y); p <- ncol(y)

	get.int <- HPDinterval(fit.mcmc, prob = prob); 
	hpd.lower <- get.int[,1]; hpd.upper <- get.int[,2]
	final.list <- list(lv.coefs.lower = matrix(hpd.lower[grep("all.params",names(hpd.lower))], nrow=p), lv.coefs.upper = matrix(hpd.upper[grep("all.params",names(hpd.upper))], nrow=p))
	rownames(final.list$lv.coefs.lower) <- rownames(final.list$lv.coefs.upper) <- colnames(y)

	if(num.lv > 0) {
		final.list$lv.lower = matrix(hpd.lower[grep("lvs", names(hpd.lower))], nrow=n)
		final.list$lv.upper = matrix(hpd.upper[grep("lvs", names(hpd.upper))], nrow=n)

		rownames(final.list$lv.lower) <- rownames(final.list$lv.upper) <- rownames(y)
		colnames(final.list$lv.lower) <- colnames(final.list$lv.upper) <- paste("lv",1:num.lv,sep="")
		colnames(final.list$lv.coefs.lower) <- colnames(final.list$lv.coefs.upper) <- c("theta0",paste("theta",1:num.lv,sep=""),"Dispersion")
		}		

	if(num.lv == 0) { colnames(final.list$lv.coefs.lower) <- colnames(final.list$lv.coefs.upper) <- c("theta0","Dispersion") }	

	if(length(grep("row.params", names(hpd.lower))) > 0) {
		final.list$row.lv.coefs.lower <- hpd.lower[grep("row.params", names(hpd.lower))]
		final.list$row.lv.coefs.upper <- hpd.upper[grep("row.params", names(hpd.upper))]
		names(final.list$row.lv.coefs.lower) <- names(final.list$row.lv.coefs.upper) <- rownames(y) }

	if(length(grep("X.params", names(hpd.lower))) > 0) {
		final.list$X.coefs.lower <- matrix(hpd.lower[grep("X.params", names(hpd.lower))],nrow=p)
		final.list$X.coefs.upper <- matrix(hpd.upper[grep("X.params", names(hpd.upper))],nrow=p)
		rownames(final.list$X.coefs.lower) <- rownames(final.list$X.coefs.upper) <- colnames(y)
		colnames(final.list$X.coefs.lower) <- colnames(final.list$X.coefs.upper) <- colnames(X) 
		
		}

# 	if(length(grep("X.multinom.params", names(hpd.lower))) > 0) {
# 		final.list$X.multinom.coefs.lower <- array(matrix(hpd.lower[grep("X.multinom.params", names(hpd.lower))],dim=c(length(index.multinom.cols),ncol(X),ncol(all.X.multinom.coefs.lower)/ncol(X))), dimnames = list("1" = index.multinom.cols, "2" = colnames(X), "level" = 1:(ncol(all.X.multinom.coefs.upper)/ncol(X))))
# 		final.list$X.multinom.coefs.upper <- array(matrix(hpd.lower[grep("X.multinom.params", names(hpd.upper))],dim=c(length(index.multinom.cols),ncol(X),ncol(all.X.multinom.coefs.lower)/ncol(X))), dimnames = list("1" = index.multinom.cols, "2" = colnames(X), "level" = 1:(ncol(all.X.multinom.coefs.upper)/ncol(X))))
# 		}

	if(length(grep("alpha", names(hpd.lower))) > 0) { ## If alpha exists, then cutoffs are there and some columns involved ordinal responses
		final.list$cutoffs.lower <- hpd.lower[grep("alpha", names(hpd.lower))]
		final.list$cutoffs.upper <- hpd.upper[grep("alpha", names(hpd.upper))]
		num.ord.levels <- length(final.list$cutoffs.lower) + 1
		names(final.list$cutoffs.lower) <- names(final.list$cutoffs.upper) <- paste(1:(num.ord.levels-1),"|",2:num.ord.levels,sep="") 
		}
				
	if(length(grep("powerparams", names(hpd.lower))) > 0) { ## If powerparam exists, then power parameters are there and some columns involved tweedie responses
		final.list$powerparam.lower <- hpd.lower[grep("powerparam", names(hpd.lower))]
		final.list$powerparam.upper <- hpd.upper[grep("powerparam", names(hpd.upper))]
		}

	return(final.list) }

## Calculates conditional WAIC, EAIC, EBIC. All of these may not work well, given the result from Millar (2009) on DIC 
## Calculate BF using harmonic mean (dangerous), returned on a log10 scale: Dangerous, since it works with the conditional and not the marginal likelihood!
## Pragmatism: Calculate the marginal likelihood at component medians, and base a AIC and BIC on that
#family = out.fit$family; trial.size = out.fit$trial.size; row.eff = out.fit$row.eff; num.lv = out.fit$num.lv; fit.mcmc = as.mcmc(out.fit$jags.model); more.measures = FALSE
get.measures <- function(y, X = NULL, family, trial.size = NULL, row.eff, num.lv, fit.mcmc, more.measures = FALSE) {
	if(length(family) != ncol(y) & length(family) != 1) { stop("Number of elements in family is either 1 or equal to # of columns in y") }
	if(length(family) == 1) complete.family <- rep(family,ncol(y))
	if(length(family) > 1) complete.family <- family
	if(any(family == "binomial") & is.null(trial.size)) 
		stop("trial.size needs to be specified if any columns are binomially distributed; can either be a single element or a vector equal to the # of columns in y. The latter will assume the specified trial size for all rows labelled binomial in the family argument.")
	if(num.lv == 0 & more.measures == TRUE) stop("For boral models with no latent variables, the marginal and conditional likelihoods are equivalent, and there is nothing to gain from setting more.measures = TRUE")

	n <- nrow(y); p <- ncol(y)
	all.lppd <- matrix(NA,nrow(fit.mcmc),n); 

	index.multinom.cols <- which(complete.family == "multinom")
	for(t in 1:nrow(fit.mcmc)) {
		if(row.eff) { get.row.coefs <- fit.mcmc[t,grep("row.params", colnames(fit.mcmc))] } else { get.row.coefs <- NULL }
		if(!is.null(X)) { get.X.coefs <- matrix(fit.mcmc[t,grep("X.params", colnames(fit.mcmc))],nrow=p) } else { get.X.coefs <- NULL }
		if(any(complete.family == "ordinal")) { get.cutoffs <- fit.mcmc[t,grep("alpha", colnames(fit.mcmc))] } else { get.cutoffs <- NULL }
		if(any(complete.family == "tweedie")) { get.powerparam <- fit.mcmc[t,grep("powerparam", colnames(fit.mcmc))] } else { get.powerparam <- NULL }
# 		if(any(complete.family == "multinom") & !is.null(X)) { 
# 			get.X.multinom.coefs <- array(matrix(fit.mcmc[t,grep("X.multinom.params", colnames(fit.mcmc))],dim = c(length(index.multinom.cols),ncol(X),ncol(get.X.multinom.coefs)/ncol(X))), dimnames=NULL) } else { get.X.multinom.coefs <- NULL }

		if(num.lv > 0) 
			get.out <- calc.condlogLik(y, X, complete.family, trial.size, lv.coefs = matrix(fit.mcmc[t,grep("all.params",colnames(fit.mcmc))],nrow=p), 
				X.coefs = get.X.coefs, row.coefs = get.row.coefs, lv = matrix(fit.mcmc[t,grep("lvs", colnames(fit.mcmc))],nrow=n), cutoffs = get.cutoffs, powerparam = get.powerparam)
		if(num.lv == 0) {
			get.out <- calc.logLik.lv0(y, X, complete.family, trial.size, lv.coefs = matrix(fit.mcmc[t,grep("all.params",colnames(fit.mcmc))],nrow=p), 
				X.coefs = get.X.coefs, row.coefs = get.row.coefs, cutoffs = get.cutoffs, powerparam = get.powerparam); get.out$logLik.comp <- get.out$logLik.row.comp }
				
		get.out$logLik.comp[!is.finite(get.out$logLik.comp)] <- NA
		all.lppd[t,] <- get.out$logLik.comp }	
		
	all.cond.logl <- apply(all.lppd,1,sum); 
	waic.out <- -2*sum(log(apply(exp(all.lppd),2,mean,na.rm=T))) + 2*sum(apply(all.lppd,2,var,na.rm=T))
	#bf.hm <- -log10(mean(1/exp(all.cond.logl[is.finite(all.cond.logl)])))
	
	## Calculate marginal logL at component medians
	coef.mat <- matrix(apply(fit.mcmc[,grep("all.params",colnames(fit.mcmc))],2,median),nrow=p)
	if(row.eff) { get.row.coefs <- apply(fit.mcmc[,grep("row.params", colnames(fit.mcmc))],2,median) } else { get.row.coefs <- NULL }
	if(!is.null(X)) { get.X.coefs <- matrix(apply(fit.mcmc[,grep("X.params",colnames(fit.mcmc))],2,median),nrow=p) } else { get.X.coefs <- NULL }
	if(any(complete.family == "ordinal")) { get.cutoffs <- apply(fit.mcmc[,grep("alpha", colnames(fit.mcmc))],2,median) } else { get.cutoffs <- NULL }
	if(any(complete.family == "tweedie")) { get.powerparam <- median(fit.mcmc[,grep("powerparam", colnames(fit.mcmc))]) } else { get.powerparam <- NULL }
# 	if(any(complete.family == "multinom") & !is.null(X)) { 
# 		get.X.multinom.coefs <- array(matrix(apply(fit.mcmc[,grep("X.multinom.params", colnames(fit.mcmc))],2,median),dim = c(length(index.multinom.cols),ncol(X),ncol(get.X.multinom.coefs)/ncol(X))), dimnames=NULL) } else { get.X.multinom.coefs <- NULL }

	if(num.lv > 0)
		median.marglogl <- calc.marglogLik(y, X, complete.family, trial.size, lv.coefs = coef.mat, X.coefs = get.X.coefs, row.coefs = get.row.coefs, num.lv, X.mc = NULL, cutoffs = get.cutoffs, powerparam = get.powerparam)
	if(num.lv == 0) 
		median.marglogl <- calc.logLik.lv0(y, X, complete.family, trial.size, lv.coefs = coef.mat, X.coefs = get.X.coefs, row.coefs = get.row.coefs, cutoffs = get.cutoffs, powerparam = get.powerparam)

	num.params <- sum(coef.mat!=0) + n*as.numeric(row.eff) + sum(get.X.coefs!=0)*as.numeric(!is.null(X)) + any(complete.family == "ordinal")*sum(get.cutoffs!=0) + (1-is.null(get.powerparam)) #+ any(complete.family == "multinom")*sum(get.X.multinom.coefs!=0)
	marg.aic <- -2*median.marglogl$logLik  + 2*num.params
	marg.bic <- -2*median.marglogl$logLik + log(n)*num.params

	eaic <- -2*mean(all.cond.logl, na.rm = T) + 2*num.params
	ebic <- -2*mean(all.cond.logl, na.rm = T) + log(n)*num.params
	
	## Compound Laplace-Metroplis estimator from Lewis and Raftery, at component medians
	get.bic2.det <- fit.mcmc[,grep("all.params",colnames(fit.mcmc))]
	if(row.eff) get.bic2.det <- cbind(get.bic2.det,fit.mcmc[,grep("row.params", colnames(fit.mcmc))])
	if(!is.null(X)) get.bic2.det <- cbind(get.bic2.det,fit.mcmc[,grep("X.params", colnames(fit.mcmc))])
	if(any(complete.family == "ordinal")) get.bic2.det <- cbind(get.bic2.det,fit.mcmc[,grep("alpha", colnames(fit.mcmc))])	
	if(any(complete.family == "tweedie")) get.bic2.det <- cbind(get.bic2.det,fit.mcmc[,grep("powerparam", colnames(fit.mcmc))])	
	#if(any(complete.family == "multinom")) get.bic2.det <- cbind(get.bic2.det,fit.mcmc[,grep("X.multinom.params", colnames(fit.mcmc))])	
	get.bic2.det <- get.bic2.det[,which(colSums(get.bic2.det)!=0)] ## Remove all zero columns -- non estimated parameters	
	marg.bic2 <- -2*median.marglogl$logLik - num.params*log(2*pi) - log(det(solve(cov(get.bic2.det)+1e-5))) 
	if(!is.finite(marg.bic2)) marg.bic2 <- -2*median.marglogl$logLik - num.params*log(2*pi) - log(det(ginv(cov(get.bic2.det)+1e-5))) 

	rm(get.bic2.det)
     
     out.list <- list(waic = waic.out, eaic = eaic, ebic = ebic, aic.median = marg.aic, bic.median = marg.bic, #bf.hm = bf.hm, 
		comp.lm = marg.bic2, all.cond.logLik = all.cond.logl, num.params = num.params)
		
	if(more.measures) {
		cat("Calculating additional information criteria...")
		more.measures <- get.more.measures(y = y, X = X, family = family, trial.size = trial.size, num.lv = num.lv, fit.mcmc = fit.mcmc, row.eff = row.eff, verbose=TRUE)
		out.list <- c(out.list, more.measures)
		}

	return(out.list)
	}
	

## Calculates marginal logl for all samples to produce a proper AIC and BIC
## Calculates WAIC based on the marginal likelihood
## Calculates DIC based on the marginal likelihood
get.more.measures <- function(y, X = NULL, family, trial.size = NULL, row.eff, num.lv, fit.mcmc, verbose=TRUE) {
	if(num.lv == 0) stop("For boral models with no latent variables, the marginal and conditional likelihoods are equivalent, and there is nothing to gain from using get.more.measures")
	if(length(family) != ncol(y) & length(family) != 1) { stop("Number of elements in family is either 1 or equal to # of columns in y") }
	if(length(family) == 1) complete.family <- rep(family,ncol(y))
	if(length(family) > 1) complete.family <- family
	if(any(family == "binomial") & is.null(trial.size)) 
		stop("trial.size needs to be specified if any columns are binomially distributed; can either be a single element or a vector equal to the # of columns in y. The latter will assume the specified trial size for all rows labelled binomial in the family argument.")

	n <- nrow(y); p <- ncol(y)
	big.X <- cbind(1,rmvnorm(2000,rep(0,num.lv))); all.marg.logl <- matrix(NA,nrow(fit.mcmc),n)

	index.multinom.cols <- which(complete.family == "multinom")
	## Calculate marginal likelihood at all iterations
	for(t in 1:nrow(fit.mcmc)) {
		if(verbose == TRUE & t %% 100 == 0) cat("Onto mcmc sample",t,"\n")

		coef.mat <- matrix(fit.mcmc[t,grep("all.params",colnames(fit.mcmc))],nrow=p)
		if(row.eff) { get.row.coefs <- fit.mcmc[t,grep("row.params", colnames(fit.mcmc))] } else { get.row.coefs <- NULL }
		if(!is.null(X)) { get.X.coefs <- matrix(fit.mcmc[t,grep("X.params",colnames(fit.mcmc))],nrow=p) } else { get.X.coefs <- NULL }
		if(any(complete.family == "ordinal")) { get.cutoffs <- fit.mcmc[t,grep("alpha", colnames(fit.mcmc))] } else { get.cutoffs <- NULL }
		if(any(complete.family == "tweedie")) { get.powerparam <- fit.mcmc[t,grep("powerparam", colnames(fit.mcmc))] } else { get.powerparam <- NULL }
# 		if(any(complete.family == "multinom") & !is.null(X)) { 
# 			get.X.multinom.coefs <- array(matrix(fit.mcmc[t,grep("X.multinom.params", colnames(fit.mcmc))],dim = c(length(index.multinom.cols),ncol(X),ncol(get.X.multinom.coefs)/ncol(X))), dimnames = NULL) } else { get.X.multinom.coefs <- NULL }
		
		get.mll <- calc.marglogLik(y, X, complete.family, trial.size, lv.coefs = coef.mat, X.coefs = get.X.coefs, row.coefs = get.row.coefs, num.lv, X.mc = big.X, cutoffs = get.cutoffs, powerparam = get.powerparam)
		all.marg.logl[t,] <- get.mll$logLik.comp }
     
	## Calculate WAIC based on marginal
	marg.waic <- -2*sum(log(apply(exp(all.marg.logl),2,mean,na.rm=T))) + 2*sum(apply(all.marg.logl,2,var,na.rm=T))

	## Calculate AIC, BIC, DIC at posterior mode	
	num.params <- sum(coef.mat!=0) + n*as.numeric(row.eff) + sum(get.X.coefs!=0)*as.numeric(!is.null(X)) + any(complete.family == "ordinal")*sum(get.cutoffs!=0) + (1-is.null(get.powerparam))#+ any(complete.family == "multinom")*sum(get.X.multinom.coefs!=0)
	bic1 <- -2*max(rowSums(all.marg.logl)) + log(n)*num.params
	aic1 <- -2*max(rowSums(all.marg.logl)) + 2*num.params

	## Calculate DIC based on marginal
	coef.mat <- matrix(apply(fit.mcmc[,grep("all.params",colnames(fit.mcmc))],2,mean),nrow=p)
	if(row.eff) { get.row.coefs <- apply(fit.mcmc[,grep("row.params", colnames(fit.mcmc))],2,mean) } else { get.row.coefs <- NULL }
	if(!is.null(X)) { get.X.coefs <- matrix(apply(fit.mcmc[,grep("X.params",colnames(fit.mcmc))],2,mean),nrow=p) } else { get.X.coefs <- NULL }
	if(any(complete.family == "ordinal")) { get.cutoffs <- apply(fit.mcmc[,grep("alpha", colnames(fit.mcmc))],2,mean) } else { get.cutoffs <- NULL }
	if(any(complete.family == "tweedie")) { get.powerparam <- mean(fit.mcmc[,grep("powerparam", colnames(fit.mcmc))]) } else { get.powerparam <- NULL }
# 	if(any(complete.family == "multinom") & !is.null(X)) { 
# 		get.X.multinom.coefs <- array(matrix(apply(fit.mcmc[,grep("X.multinom.params", colnames(fit.mcmc))],2,mean),dim = c(length(index.multinom.cols),ncol(X),ncol(get.X.multinom.cols)/ncol(X))), dimnames = NULL) } else { get.X.multinom.coefs <- NULL }
		
 	marg.dic <- -2*calc.marglogLik(y, X, complete.family, trial.size, lv.coefs = coef.mat, X.coefs = get.X.coefs, row.coefs = get.row.coefs, num.lv, X.mc = big.X, cutoffs = get.cutoffs, powerparam = get.powerparam)$logLik
 	marg.dic <- marg.dic + 2*(2*var(rowSums(all.marg.logl),na.rm=T))
		
	return(list(aic.mode = aic1, bic.mode = bic1, marg.dic = marg.dic, marg.waic = marg.waic, all.marg.logLik = rowSums(all.marg.logl), num.params = num.params)) }	
	
## Produce the correlation due to similarity of responses to X
get.enviro.cor <- function(y, X, fit.mcmc, est = "median") {
	if(length(grep("X.params", colnames(fit.mcmc))) == 0) stop("Cannot find MCMC sample corresponding to coefficients for X.")

	n <- nrow(y); p <- ncol(y)
	enviro.cor.mat <- matrix(0,p,p)
	if(is.null(colnames(y))) colnames(y) <- 1:ncol(y); 
	rownames(enviro.cor.mat) <- colnames(enviro.cor.mat) <- colnames(y)
	all.enviro.cor.mat <- array(0,dim=c(nrow(fit.mcmc),p,p))

	for(t in 1:nrow(fit.mcmc)) {
		cw.X.coefs <- matrix(fit.mcmc[t,grep("X.params", colnames(fit.mcmc))],nrow=p)
		enviro.linpreds <- X%*%t(as.matrix(cw.X.coefs))
		all.enviro.cor.mat[t,,] <- cor(enviro.linpreds) }

	for(j in 1:p) { for(j2 in 1:p) { ## Average/Median over the MCMC samples
		if(est == "median") enviro.cor.mat[j,j2] <- median(all.enviro.cor.mat[,j,j2])
		if(est == "mean") enviro.cor.mat[j,j2] <- mean(all.enviro.cor.mat[,j,j2]) } }
		
	#return(list(residual.correlations = enviro.cor.mat))
	#corrplot(enviro.cor.mat, title = "Environmental correlations", type = "lower")
	return(enviro.cor.mat)
	}

## Produce the residual correlation based on latent variables
get.residual.cor <- function(y, fit.mcmc, est = "median") {
	if(length(grep("lvs", colnames(fit.mcmc))) == 0) stop("Cannot find MCMC samples corresponding to latent variables.")

	n <- nrow(y); p <- ncol(y)
	res.cor.mat <- matrix(0,p,p)
	if(is.null(colnames(y))) colnames(y) <- 1:ncol(y); 
	rownames(res.cor.mat) <- colnames(res.cor.mat) <- colnames(y)
	all.res.cor.mat <- array(0,dim=c(nrow(fit.mcmc),p,p))

	for(t in 1:nrow(fit.mcmc)) {
		lvs <- matrix(fit.mcmc[t,grep("lvs", colnames(fit.mcmc))],nrow=n)
		num.lv <- ncol(lvs)
		lvs.coefs <- matrix(fit.mcmc[t,grep("all.params", colnames(fit.mcmc))],nrow=p)
		res.linpreds <- lvs%*%t(as.matrix(lvs.coefs[,2:(num.lv+1)]))
		all.res.cor.mat[t,,] <- cor(res.linpreds) }

	for(j in 1:p) { for(j2 in 1:p) { ## Average/Median over the MCMC samples
		if(est == "median") res.cor.mat[j,j2] <- median(all.res.cor.mat[,j,j2])
		if(est == "mean") res.cor.mat[j,j2] <- mean(all.res.cor.mat[,j,j2]) } }
		
	#return(list(residual.correlations = res.cor.mat))
	#corrplot(res.cor.mat, title = "Residual correlations", type = "lower")
	return(res.cor.mat)
	}
		
#, index.multinom.cols = NULL
make.jagsboralmodel <- function(family, num.X = 0, row.eff, n, p, hypparams = c(100,100,100,100), ssvs.index = -1) {
	X.eff <- ifelse(num.X == 0, FALSE, TRUE)
	if(length(family) == 1) complete.family <- rep(family,p)
	if(length(family) > 1) complete.family <- family
	if(length(ssvs.index) == 1 & X.eff) ssvs.index <- rep(ssvs.index, num.X)
	mod.general.lv <- paste("model { \n\t C <- 10000 \n\t ##Likelihood \n\t for(i in 1:n) {",sep="")
	
	index.ord.cols <- which(complete.family == "ordinal")
	for(j in 1:p) {
		if(complete.family[j] != "multinom") {
			if(!X.eff & row.eff) mod.general.lv <- c(mod.general.lv, paste("\t\t eta[i,",j,"] <- row.params[i] + inprod(all.params[",j,",2:(num.lv+1)],lvs[i,])",sep="")) 
			if(X.eff & row.eff) mod.general.lv <- c(mod.general.lv, paste("\t\t eta[i,",j,"] <- row.params[i] + inprod(all.params[",j,",2:(num.lv+1)],lvs[i,]) + inprod(X.params[",j,",],X[i,])",sep="")) 
			if(!X.eff & !row.eff) mod.general.lv <- c(mod.general.lv, paste("\t\t eta[i,",j,"] <- inprod(all.params[",j,",2:(num.lv+1)],lvs[i,])",sep="")) 
			if(X.eff & !row.eff) mod.general.lv <- c(mod.general.lv, paste("\t\t eta[i,",j,"] <- inprod(all.params[",j,",2:(num.lv+1)],lvs[i,]) + inprod(X.params[",j,",],X[i,])",sep="")) 
			}

		if(complete.family[j] == "negative.binomial") {
			mod.general.lv <- c(mod.general.lv, paste("\t\t u[i,",j,"] ~ dgamma(all.params[",j,",num.lv+2], all.params[",j,",num.lv+2])",sep="")) 
			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dpois(exp(all.params[",j,",1] + eta[i,",j,"])*u[i,",j,"]) ## Parameterizing the NB as a multiplicative random effect models\n",sep="")) 
			}
		if(complete.family[j] == "normal") {
			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dnorm(all.params[",j,",1] + eta[i,",j,"],all.params[",j,",num.lv+2]) \n",sep="")) 
			}
		if(complete.family[j] == "binomial") {
			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dbin(ilogit(all.params[",j,",1] + eta[i,",j,"]),complete.trial.size[",j,"])\n",sep="")) 
			}
		if(complete.family[j] == "exponential") {
			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dexp(pow(exp(all.params[",j,",1] + eta[i,",j,"]),-1))\n",sep="")) 
			}
		if(complete.family[j] == "gamma") {
			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dgamma(exp(all.params[",j,",1] + eta[i,",j,"])*all.params[",j,",num.lv+2], all.params[",j,",num.lv+2])\n",sep="")) 
			}
		if(complete.family[j] == "beta") {
			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dbeta(ilogit(all.params[",j,",1] + eta[i,",j,"])*all.params[",j,",num.lv+2],(1-ilogit(all.params[",j,",1] + eta[i,",j,"]))*all.params[",j,",num.lv+2])\n",sep="")) 
			}
		if(complete.family[j] == "poisson") {
			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dpois(exp(all.params[",j,",1] + eta[i,",j,"]))\n",sep="")) 
			}		
		if(complete.family[j] == "lnormal") {
			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dlnorm(all.params[",j,",1] + eta[i,",j,"],all.params[",j,",num.lv+2]) \n",sep="")) 
			}		
		if(complete.family[j] == "tweedie") {
			mod.general.lv <- c(mod.general.lv, paste("\t\t lambdanum[i,",j,"] <- pow(exp(all.params[",j,",1] + eta[i,",j,"]),2-powerparam)/(all.params[",j,",num.lv+2]*(2-powerparam))",sep=""))
			mod.general.lv <- c(mod.general.lv, paste("\t\t numfish[i,",j,"] ~ dpois(lambdanum[i,",j,"])",sep=""))
			mod.general.lv <- c(mod.general.lv, paste("\t\t choose.shape[i,",j,",1] <- numfish[i,",j,"]*(2-powerparam)/(powerparam-1)",sep="")) 
			mod.general.lv <- c(mod.general.lv, paste("\t\t choose.rate[i,",j,",1] <- 1/(all.params[",j,",num.lv+2]*(powerparam-1)*pow(exp(all.params[",j,",1] + eta[i,",j,"]),powerparam-1))",sep="")) 
			mod.general.lv <- c(mod.general.lv, paste("\t\t choose.shape[i,",j,",2] <- 1",sep="")) 			
			mod.general.lv <- c(mod.general.lv, paste("\t\t choose.rate[i,",j,",2] <- exp(-lambdanum[i,",j,"])",sep="")) 
			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dgamma(choose.shape[i,",j,",1+equals(y[i,",j,"],0)],choose.rate[i,",j,",1+equals(y[i,",j,"],0)]) \n",sep="")) 
			}		
		if(complete.family[j] == "ordinal") {
    			mod.general.lv <- c(mod.general.lv, paste("\t\t prob[i,",which(index.ord.cols == j),",1] <- ilogit(alpha[1]-eta[i,",j,"]-all.params[",j,",1])",sep=""))
    			mod.general.lv <- c(mod.general.lv, paste("\t\t for(k in 2:(num.ord.levels-1)) {",sep=""))			
    			mod.general.lv <- c(mod.general.lv, paste("\t\t\t prob[i,",which(index.ord.cols == j),",k] <- ilogit(alpha[k]-eta[i,",j,"]-all.params[",j,",1]) - ilogit(alpha[k-1]-eta[i,",j,"]-all.params[",j,",1]) }",sep=""))
    			mod.general.lv <- c(mod.general.lv, paste("\t\t prob[i,",which(index.ord.cols == j),",num.ord.levels] <- 1-ilogit(alpha[num.ord.levels-1]-eta[i,",j,"]-all.params[",j,",1])",sep=""))
    			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dcat(prob[i,",which(index.ord.cols == j),",]+0.001)\n",sep="")) 
			}

		if(complete.family[j] == "multinom") { ## Coefficients for LVs are constrained to be same for all levels! Otherwise identifiability constraints are hard!
			stop("You shouldn't have gotten here!")
#    			mod.general.lv <- c(mod.general.lv, paste("\t\t for(k in 1:num.multinom.levels[",j,"]) {",sep=""))
# 			if(!X.eff & row.eff) mod.general.lv <- c(mod.general.lv, paste("\t\t\t mu[i,",which(index.multinom.cols == j),",k] <- exp(row.params[i] + inprod(all.params[",j,",2:(num.lv+1)],lvs[i,]))",sep="")) 
# 			if(X.eff & row.eff) mod.general.lv <- c(mod.general.lv, paste("\t\t\t mu[i,",which(index.multinom.cols == j),",k] <- exp(row.params[i] + inprod(all.params[",j,",2:(num.lv+1)],lvs[i,]) + inprod(X.multinom.params[",which(index.multinom.cols == j),",,k],X[i,]))",sep="")) 
# 			if(!X.eff & !row.eff) mod.general.lv <- c(mod.general.lv, paste("\t\t\t mu[i,",which(index.multinom.cols == j),",k] <- exp(inprod(all.params[",j,",2:(num.lv+1)],lvs[i,]))",sep="")) 
# 			if(X.eff & !row.eff) mod.general.lv <- c(mod.general.lv, paste("\t\t\t mu[i,",which(index.multinom.cols == j),",k] <- exp(inprod(all.params[",j,",2:(num.lv+1)],lvs[i,]) + inprod(X.multinom.params[",which(index.multinom.cols == j),",,k],X[i,]))",sep="")) 			
# 			mod.general.lv <- c(mod.general.lv, paste("\t\t\t prob[i,",which(index.multinom.cols == j),",k] <- mu[i,",which(index.multinom.cols == j),",k]/sum(mu[i,",which(index.multinom.cols == j),",]) }",sep="")) 
# 			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dcat(prob[i,",which(index.multinom.cols == j),",]+0.001)\n",sep="")) 
			}		
		}

	mod.general.lv <- c(mod.general.lv, paste("\n\n\t } for(i in 1:n) { for(k in 1:num.lv) { lvs[i,k] ~ dnorm(0,1) } } \n\n\t ## Latent variables",sep=""))

	prior.string <- paste("dnorm(0,",1/hypparams[1],")",sep="")
	# 	if(prior.type[1] == "normal") prior.string <- paste("dnorm(0,",1/hypparams[1],")",sep="")
	# 	if(prior.type[1] == "t") prior.string <- paste("dt(0,",1/hypparams[1],",1)",sep="")
	# 	if(prior.type[1] == "uniform") prior.string <- paste("dunif(-,",hypparams[1],",",hypparams[1],")",sep="")

	if(any(complete.family == "ordinal")) { ## If all ordinal responses, impose a sum to zero constraint on the intercepts for those columns
		if(length(index.ord.cols) > 1) {
			for(j in index.ord.cols[-1]) {
				mod.general.lv <- c(mod.general.lv, paste("\t all.params[",j,",1] ~ ", prior.string," ## Ordinal species intercept",sep="")) }
			mod.general.lv <- c(mod.general.lv, paste("\t all.params[",index.ord.cols[1],",1] <- -1*(",paste("all.params[",index.ord.cols[-1],",1]",sep="",collapse="+"),")",sep="")) ## Sum to zero constraint
			for(j in (1:p)[-index.ord.cols]) mod.general.lv <- c(mod.general.lv, paste("\t all.params[",j,",1] ~ ",prior.string,sep="")) ## All other intercepts
			}  
		if(length(index.ord.cols) == 1) {
			mod.general.lv <- c(mod.general.lv, paste("\t all.params[",index.ord.cols,",1] <- 0 ## Ordinal species intercept",sep=""))
			for(j in (1:p)[-index.ord.cols]) mod.general.lv <- c(mod.general.lv, paste("\t all.params[",j,",1] ~ ",prior.string,sep="")) ## All other intercepts
			}  
		}		
	if(!any(complete.family == "ordinal")) { mod.general.lv <- c(mod.general.lv, paste("\t for(i in 1:p) { \n\t\t all.params[i,1] ~ ", prior.string," } ## Species intercept \n",sep="")) }

	if(row.eff) mod.general.lv <- c(mod.general.lv, paste("\t for(i in 1:n) { row.params[i] ~ ", prior.string, " }",sep=""))		

	prior.string <- paste("dnorm(0,",1/hypparams[2],")",sep="")
# 	if(prior.type[2] == "normal") prior.string <- paste("dnorm(0,",1/hypparams[2],")",sep="")
# 	if(prior.type[2] == "t") prior.string <- paste("dt(0,",1/hypparams[2],",1)",sep="")
# 	if(prior.type[2] == "uniform") prior.string <- paste("dunif(-,",hypparams[2],",",hypparams[2],")",sep="")
	mod.general.lv <- c(mod.general.lv, paste("\t for(i in 1:(num.lv-1)) { for(j in (i+2):(num.lv+1)) { \n\t\t all.params[i,j] <- 0 } } ## Constraints to 0 on upper diagonal",sep=""))
	mod.general.lv <- c(mod.general.lv, paste("\t for(i in 1:num.lv) { \n\t\t all.params[i,i+1] ~ dunif(0,",hypparams[2],") } ## Sign constraints on diagonal elements",sep=""))
	mod.general.lv <- c(mod.general.lv, paste("\t for(i in 2:num.lv) { for(j in 2:i) { \n\t\t all.params[i,j] ~ ", prior.string, " } } ## Free lower diagonals",sep=""))
	mod.general.lv <- c(mod.general.lv, paste("\t for(i in (num.lv+1):p) { for(j in 2:(num.lv+1)) { \n\t\t all.params[i,j] ~ ", prior.string, " } } ## All other elements",sep=""))

	prior.string <- paste("dunif(0,",hypparams[4],")",sep="")
# 	if(prior.type[4] == "uniform") prior.string <- paste("dunif(0,",hypparams[4],")",sep="")
# 	if(prior.type[4] == "gamma") prior.string <- paste("dgamma(",1/hypparams[4],",",1/hypparams[4],")",sep="")
# 	if(prior.type[4] == "halft") prior.string <- paste("dt(0,",1/hypparams[4],",1)I(0,)",sep="")
	mod.general.lv <- c(mod.general.lv, paste("\t for(j in 1:p) { all.params[j,num.lv+2] ~ ",prior.string," }",sep=""))

	if(X.eff) {
		mod.general.lv <- c(mod.general.lv, paste("\n"))
		prior.string <- paste("dnorm(0,",1/hypparams[3],")",sep="")
		#if(prior.type[3] == "t") prior.string <- paste("dt(0,",1/hypparams[3],",1)",sep="")
		#if(prior.type[3] == "unif") prior.string <- paste("dunif(-,",hypparams[3],",",hypparams[3],")",sep="")
		for(i in 1:length(ssvs.index)) {
			## Covariate remains in model
			if(ssvs.index[i] == -1) mod.general.lv <- c(mod.general.lv, paste("\t for(j in 1:p) { X.params[j,",i,"] ~ ",prior.string," } ",sep=""))
			## Individual coefficients tested
			ssvs.prior.string <- paste("dnorm(0,pow(",hypparams[3],"*((1-probindX",i,"[j])*0.0001+probindX",i,"[j]),-1))",sep="")
			if(ssvs.index[i] == 0) mod.general.lv <- c(mod.general.lv, paste("\t for(j in 1:p) { X.params[j,",i,"] ~ ",ssvs.prior.string,"; probindX",i,"[j] ~ dbern(0.5) }",sep=""))
			## Coefficients tested collectively, potentially along with other covariates
			ssvs.prior.string <- paste("dnorm(0,pow(",hypparams[3],"*((1-probGpX",ssvs.index[i],")*0.0001+probGpX",ssvs.index[i],"),-1))",sep="")
			if(ssvs.index[i] > 0) mod.general.lv <- c(mod.general.lv, paste("\t for(j in 1:p) { X.params[j,",i,"] ~ ",ssvs.prior.string," } ",sep=""))
			}
		if(any(ssvs.index > 0)) {
			for(i in unique(ssvs.index[ssvs.index>0])) { mod.general.lv <- c(mod.general.lv, paste("\t probGpX",i," ~ dbern(0.5)",sep="")) }
			}
		}	
# 	if(X.eff & any(family == "multinom")) {
# 		mod.general.lv <- c(mod.general.lv, paste("\t for(j in 1:",length(index.multinom.cols),") { for(i in 1:num.X) { X.multinom.params[j,i,1] <- 0 } } ",sep=""))
# 		mod.general.lv <- c(mod.general.lv, paste("\t for(j in 1:",length(index.multinom.cols),") { for(i in 1:num.X) { for(k in 2:num.multinom.levels) { X.multinom.params[j,i,k] ~ dnorm(0,",1/hypparams[1],") } } } ",sep=""))
# 		}

	if(any(complete.family == "tweedie")) 
		mod.general.lv <- c(mod.general.lv, paste("\t powerparam ~ dunif(1,2)"))
	if(any(complete.family == "ordinal")) {
		mod.general.lv <- c(mod.general.lv, paste("\t for(k in 1:(num.ord.levels-1)) { alpha0[k] ~ dnorm(0,",1/hypparams[1],") }"))
		mod.general.lv <- c(mod.general.lv, paste("\t alpha[1:(num.ord.levels-1)] <- sort(alpha0)"))	
		} 

	mod.general.lv <- c(mod.general.lv, "\n\t }")
	write(mod.general.lv, file = "jagsboralmodel.txt")	
	}
		
make.jagsboralnullmodel <- function(family, num.X = 0, row.eff, n, p, hypparams = c(100,100,100,100), ssvs.index = -1) {
	X.eff <- ifelse(num.X == 0, FALSE, TRUE)
	if(length(family) == 1) complete.family <- rep(family,p)
	if(length(family) > 1) complete.family <- family
	if(length(ssvs.index) == 1 & X.eff) ssvs.index <- rep(ssvs.index, num.X)
	mod.general.lv <- paste("model { \n\t C <- 10000 \n\t ##Likelihood \n\t for(i in 1:n) {",sep="")

	index.ord.cols <- which(complete.family == "ordinal")
	for(j in 1:p) {
		if(complete.family[j] != "multinom") {
			if(!X.eff & row.eff) mod.general.lv <- c(mod.general.lv, paste("\t\t eta[i,",j,"] <- row.params[i]",sep="")) 
			if(X.eff & row.eff) mod.general.lv <- c(mod.general.lv, paste("\t\t eta[i,",j,"] <- row.params[i] + inprod(X.params[",j,",],X[i,])",sep="")) 
			if(!X.eff & !row.eff) mod.general.lv <- c(mod.general.lv, paste("\t\t eta[i,",j,"] <- 0",sep="")) 
			if(X.eff & !row.eff) mod.general.lv <- c(mod.general.lv, paste("\t\t eta[i,",j,"] <- inprod(X.params[",j,",],X[i,])",sep="")) 
		}
		
		if(complete.family[j] == "negative.binomial") {
			mod.general.lv <- c(mod.general.lv, paste("\t\t prob[i,",j,"] <- all.params[",j,",2]/(all.params[",j,",2]+exp(all.params[",j,",1] + eta[i,",j,"]))",sep="")) 
			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dnegbin(prob[i,",j,"], all.params[",j,",2])## Parameterizing the NB as a multiplicative random effect models\n",sep="")) 
			}
		if(complete.family[j] == "normal") {
			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dnorm(all.params[",j,",1] + eta[i,",j,"],all.params[",j,",2]) \n",sep="")) 
			}
		if(complete.family[j] == "binomial") {
			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dbin(ilogit(all.params[",j,",1] + eta[i,",j,"]),complete.trial.size[",j,"])\n",sep="")) 
			}
		if(complete.family[j] == "exponential") {
			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dexp(pow(exp(all.params[",j,",1] + eta[i,",j,"]),-1))\n",sep="")) 
			}
		if(complete.family[j] == "gamma") {
			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dgamma(exp(all.params[",j,",1] + eta[i,",j,"])*all.params[",j,",num.lv+2], all.params[",j,",num.lv+2])\n",sep="")) 
			}
		if(complete.family[j] == "beta") {
			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dbeta(ilogit(all.params[",j,",1] + eta[i,",j,"])*all.params[",j,",num.lv+2],(1-ilogit(all.params[",j,",1] + eta[i,",j,"]))*all.params[",j,",num.lv+2])\n",sep="")) 
			}
		if(complete.family[j] == "poisson") {
				mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dpois(exp(all.params[",j,",1] + eta[i,",j,"]))\n",sep="")) 
			}		
		if(complete.family[j] == "lnormal") {
			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dlnorm(all.params[",j,",1] + eta[i,",j,"],all.params[",j,",2]) \n",sep="")) 
			}
		if(complete.family[j] == "tweedie") {
			mod.general.lv <- c(mod.general.lv, paste("\t\t lambdanum[i,",j,"] <- pow(exp(all.params[",j,",1] + eta[i,",j,"]),2-powerparam)/(all.params[",j,",2]*(2-powerparam))",sep=""))			
			mod.general.lv <- c(mod.general.lv, paste("\t\t numfish[i,",j,"] ~ dpois(lambdanum[i,",j,"])",sep=""))
			mod.general.lv <- c(mod.general.lv, paste("\t\t choose.shape[i,",j,",1] <- numfish[i,",j,"]*(2-powerparam)/(powerparam-1)",sep="")) 
			mod.general.lv <- c(mod.general.lv, paste("\t\t choose.rate[i,",j,",1] <- 1/(all.params[",j,",2]*(powerparam-1)*pow(exp(all.params[",j,",1] + eta[i,",j,"]),powerparam-1))",sep="")) 
			mod.general.lv <- c(mod.general.lv, paste("\t\t choose.shape[i,",j,",2] <- 1",sep="")) 			
			mod.general.lv <- c(mod.general.lv, paste("\t\t choose.rate[i,",j,",2] <- exp(-lambdanum[i,",j,"])",sep="")) 
			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dgamma(choose.shape[i,",j,",1+equals(y[i,",j,"],0)],choose.rate[i,",j,",1+equals(y[i,",j,"],0)]) \n",sep="")) 
			}		
		if(complete.family[j] == "ordinal") {
    			mod.general.lv <- c(mod.general.lv, paste("\t\t prob[i,",which(index.ord.cols == j),",1] <- ilogit(alpha[1]-eta[i,",j,"]-all.params[",j,",1])",sep=""))
    			mod.general.lv <- c(mod.general.lv, paste("\t\t for(k in 2:(num.ord.levels-1)) {",sep=""))			
    			mod.general.lv <- c(mod.general.lv, paste("\t\t\t prob[i,",which(index.ord.cols == j),",k] <- ilogit(alpha[k]-eta[i,",j,"]-all.params[",j,",1]) - ilogit(alpha[k-1]-eta[i,",j,"]-all.params[",j,",1]) }",sep=""))
    			mod.general.lv <- c(mod.general.lv, paste("\t\t prob[i,",which(index.ord.cols == j),",num.ord.levels] <- 1-ilogit(alpha[num.ord.levels-1]-eta[i,",j,"]-all.params[",j,",1])",sep=""))
    			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dcat(prob[i,",which(index.ord.cols == j),",]+0.001)\n",sep="")) 
			}
			
		if(complete.family[j] == "multinom") { ## Coefficients for LVs are constrained to be same for all levels! Otherwise identifiability constraints are hard!
			stop("You shouldn't have gotten here!")
#    			mod.general.lv <- c(mod.general.lv, paste("\t\t for(k in 1:num.multinom.levels) {",sep=""))
# 			if(!X.eff & row.eff) mod.general.lv <- c(mod.general.lv, paste("\t\t\t mu[i,",which(index.multinom.cols == j),",k] <- exp(row.params[i] + all.params[",j,",1])",sep="")) 
# 			if(X.eff & row.eff) mod.general.lv <- c(mod.general.lv, paste("\t\t\t mu[i,",which(index.multinom.cols == j),",k] <- exp(row.params[i] + all.params[",j,",1] + inprod(X.multinom.params[",which(index.multinom.cols == j),",,k],X[i,]))",sep="")) 
# 			if(!X.eff & !row.eff) mod.general.lv <- c(mod.general.lv, paste("\t\t\t mu[i,",which(index.multinom.cols == j),",k] <- exp(all.params[",j,",1])",sep="")) 
# 			if(X.eff & !row.eff) mod.general.lv <- c(mod.general.lv, paste("\t\t\t mu[i,",which(index.multinom.cols == j),",k] <- exp(all.params[",j,",1] + inprod(X.multinom.params[",which(index.multinom.cols == j),",,k],X[i,]))",sep="")) 			
# 			mod.general.lv <- c(mod.general.lv, paste("\t\t\t prob[i,",which(index.multinom.cols == j),",k] <- mu[i,",which(index.multinom.cols == j),",k]/sum(mu[i,",which(index.multinom.cols == j),",]) }",sep="")) 
# 
# 			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dcat(prob[i,",which(index.multinom.cols == j),",]+0.001)\n",sep="")) 
			}		

		}

	mod.general.lv <- c(mod.general.lv, paste("\t\t } \n\n\t ## Prior",sep=""))
	prior.string <- paste("dnorm(0,",1/hypparams[1],")",sep="")
	#if(prior.type[1] == "t") prior.string <- paste("dt(0,",1/hypparams[1],",1)",sep="")
	#if(prior.type[1] == "unif") prior.string <- paste("dunif(-,",hypparams[1],",",hypparams[1],")",sep="")

	if(any(complete.family == "ordinal")) { ## If all ordinal responses, impose a sum to zero constraint on the intercepts for those columns
		if(length(index.ord.cols) > 1) {
			for(j in index.ord.cols[-1]) {
				mod.general.lv <- c(mod.general.lv, paste("\t all.params[",j,",1] ~ ",prior.string," ## Ordinal species intercept",sep="")) }
			mod.general.lv <- c(mod.general.lv, paste("\t all.params[",index.ord.cols[1],",1] <- -1*(",paste("all.params[",index.ord.cols[-1],",1]",sep="",collapse="+"),")",sep="")) ## Sum to zero constraint
			for(j in (1:p)[-index.ord.cols]) mod.general.lv <- c(mod.general.lv, paste("\t all.params[",j,",1] ~ ",prior.string,sep="")) ## All other intercepts
			}  
		if(length(index.ord.cols) == 1) {
			mod.general.lv <- c(mod.general.lv, paste("\t all.params[",index.ord.cols,",1] <- 0 ## Ordinal species intercept",sep=""))
			for(j in (1:p)[-index.ord.cols]) mod.general.lv <- c(mod.general.lv, paste("\t all.params[",j,",1] ~ ",prior.string,sep="")) ## All other intercepts
			}  
		}		
	if(!any(complete.family == "ordinal")) { mod.general.lv <- c(mod.general.lv, paste("\t for(i in 1:p) { \n\t\t all.params[i,1] ~ ",prior.string," } ## Species intercept \n",sep="")) }

	if(row.eff) mod.general.lv <- c(mod.general.lv, paste("\t for(i in 1:n) { row.params[i] ~ ",prior.string," }",sep=""))		

	mod.general.lv <- c(mod.general.lv, paste("\t for(j in 1:p) { all.params[j,2] ~ dunif(0,",hypparams[4],") }",sep=""))
	#if(prior.type[4] == "gamma") mod.general.lv <- c(mod.general.lv, paste("\t for(j in 1:p) { all.params[j,2] ~ dgamma(",1/hypparams[4],",",1/hypparams[4],") }",sep=""))
	#if(prior.type[4] == "halft") mod.general.lv <- c(mod.general.lv, paste("\t for(j in 1:p) { all.params[j,2] ~ dt(0,",1/hypparams[4],",1)I(0,) }",sep=""))

	if(X.eff) {
		mod.general.lv <- c(mod.general.lv, paste("\n"))
		prior.string <- paste("dnorm(0,",1/hypparams[3],")",sep="")
		#if(prior.type[3] == "t") prior.string <- paste("dt(0,",1/hypparams[3],",1)",sep="")
		#if(prior.type[3] == "unif") prior.string <- paste("dunif(-,",hypparams[3],",",hypparams[3],")",sep="")
		for(i in 1:length(ssvs.index)) {
			## Covariate remains in model
			if(ssvs.index[i] == -1) mod.general.lv <- c(mod.general.lv, paste("\t for(j in 1:p) { X.params[j,",i,"] ~ ",prior.string," } ",sep=""))
			## Individual coefficients tested
			ssvs.prior.string <- paste("dnorm(0,pow(",hypparams[3],"*((1-probindX",i,"[j])*0.0001+probindX",i,"[j]),-1))",sep="")
			if(ssvs.index[i] == 0) mod.general.lv <- c(mod.general.lv, paste("\t for(j in 1:p) { X.params[j,",i,"] ~ ",ssvs.prior.string,"; probindX",i,"[j] ~ dbern(0.5) }",sep=""))
			## Coefficients tested collectively, potentially along with other covariates
			ssvs.prior.string <- paste("dnorm(0,pow(",hypparams[3],"*((1-probGpX",ssvs.index[i],")*0.0001+probGpX",ssvs.index[i],"),-1))",sep="")
			if(ssvs.index[i] > 0) mod.general.lv <- c(mod.general.lv, paste("\t for(j in 1:p) { X.params[j,",i,"] ~ ",ssvs.prior.string," } ",sep=""))
			}
		if(any(ssvs.index > 0)) {
			for(i in unique(ssvs.index[ssvs.index>0])) { mod.general.lv <- c(mod.general.lv, paste("\t probGpX",i," ~ dbern(0.5)",sep="")) }
			}
		}	
# 	if(X.eff & any(family == "multinom")) {
# 		mod.general.lv <- c(mod.general.lv, paste("\t for(j in 1:",length(index.multinom.cols),") { for(i in 1:num.X) { X.multinom.params[j,i,1] <- 0 } } ",sep=""))
# 		mod.general.lv <- c(mod.general.lv, paste("\t for(j in 1:",length(index.multinom.cols),") { for(i in 1:num.X) { for(k in 2:num.multinom.levels) { X.multinom.params[j,i,k] ~ dnorm(0,",1/hypparams[1],") } } } ",sep=""))
# 		}

	if(any(complete.family == "tweedie")) 
		mod.general.lv <- c(mod.general.lv, paste("\t powerparam ~ dunif(1,2)"))
	if(any(complete.family == "ordinal")) {
		mod.general.lv <- c(mod.general.lv, paste("\t for(k in 1:(num.ord.levels-1)) { alpha0[k] ~ dnorm(0,",1/hypparams[1],") }"))
		mod.general.lv <- c(mod.general.lv, paste("\t alpha[1:(num.ord.levels-1)] <- sort(alpha0)"))	
		} 
	
	mod.general.lv <- c(mod.general.lv, "\n\t }")
	write(mod.general.lv, file = "jagsboralmodel.txt")	
	}


## Hidden Internal function - Given the lvs, coefficients and cutoffs, return the multinomial probabilities for proportional odds regression for specific element of y
ordinal.conversion.spp <- function(n, lv = NULL, lv.coefs.j = NULL, num.lv = NULL, row.coefs = NULL, X = NULL, X.coefs.j = NULL, cutoffs) {
	if((!is.null(X) & is.null(X.coefs.j)) | (is.null(X) & !is.null(X.coefs.j))) stop("X and X.coefs must be supplied simultaneously")
	if(is.null(num.lv)) num.lv <- 0
	if(num.lv < ncol(lv)) stop("# of latent variables in num.lv is inconsistent with length of lv")
	
	etas <- matrix(0,n,length(cutoffs)) ## num.ord.levels - 1
	for(k in 1:length(cutoffs)) {
		etas[,k] <- rep(cutoffs[k],n) - lv.coefs.j[1]
		if(!is.null(lv)) etas[,k] <- etas[,k] - as.matrix(lv)%*%lv.coefs.j[2:(num.lv+1)]
		if(!is.null(row.coefs)) etas[,k] <- etas[,k] - row.coefs
		if(!is.null(X)) etas[,k] <- etas[,k] - as.matrix(X)%*%X.coefs.j ## Don't forget the negative sign!
		}
	probs <- matrix(0,n,length(cutoffs)+1) ## num.ord.levels
	probs[,1] <- exp(etas[,1])/(1+exp(etas[,1]))
	for(k in 2:ncol(etas)) { probs[,k] <- exp(etas[,k])/(1+exp(etas[,k])) - exp(etas[,k-1])/(1+exp(etas[,k-1])) }		
	probs[,length(cutoffs)+1] <- 1 - exp(etas[,length(cutoffs)])/(1+exp(etas[,length(cutoffs)]))

	rm(etas); return(probs)
	}	
	

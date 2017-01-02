##########
## Auxilary functions
##########
.onAttach <- function(...) {
	packageStartupMessage("If you recently updated boral, please check news(package = \"boral\") for the updates in the latest version.", appendLF=TRUE)
}


## Calculate conditional logl
## loglik = sum_{i=1}^n sum_{j=1}^s \log( f(y_ij|z_i) )
#  lv.coefs = matrix(fit.mcmc[t,grep("all.params",colnames(fit.mcmc))],nrow=p)
#  X.coefs = cw.X.coefs
#  row.coefs = cw.row.coefs
#  lv = matrix(fit.mcmc[t,grep("lvs", colnames(fit.mcmc))],nrow=n)
#  cutoffs = cw.cutoffs
#, X.multinom.coefs = NULL
calc.condlogLik <- function(y, X = NULL, family, trial.size = 1, lv.coefs, X.coefs = NULL, row.coefs = NULL, row.ids = NULL, lv = NULL, cutoffs = NULL, powerparam = NULL) {
	if(length(family) != ncol(y) & length(family) != 1) { 
		stop("Number of elements in family is either 1 or equal to # of columns in y") }
	if(length(family) == 1) complete.family <- rep(family,ncol(y))
	if(length(family) > 1) complete.family <- family

	if(any(complete.family == "ordinal") & is.null(cutoffs)) 
		stop("Ordinal data requires cutoffs to be supplied") 
	if(any(family == "tweedie") & (powerparam < 1 || powerparam > 2)) 
		stop("Common power parameter for tweedie must be between 1 and 2.") 
	#if(any(complete.family == "multinom") & is.null(X.multinom.coefs)) stop("Multinomial data requires X.multinom.coefs to be supplied.") 

	if(!is.null(row.coefs)) {
		if(is.null(row.ids)) row.ids <- matrix(1:nrow(y), ncol = 1)
		if(!is.list(row.coefs))
			stop("row.coefs should be a list with length equal to the number of columns in row.coefs.")
		if(length(row.coefs) != ncol(row.ids))
			stop("row.coefs should be a list with length equal to the number of columns in row.coefs.")
		}
	if(!is.null(row.ids)) {
		row.ids <- as.matrix(row.ids)
		if(nrow(row.ids) != nrow(y)) 
			stop("Number of rows in the matrix row.ids should be equal to number of rows in y.")
		if(is.null(colnames(row.ids))) colnames(row.ids) <- paste0("ID", 1:ncol(row.ids))
		}

		
	## Finish checks
	n <- nrow(y); p <- ncol(y); 
	num.lv <- 0; if(!is.null(lv)) num.lv <- ncol(lv)
	loglik <- 0; loglik.comp <- matrix(NA, nrow = n, ncol = p) 

	if(length(trial.size) == 1) complete.trial.size <- rep(trial.size,ncol(y))
	if(length(trial.size) > 1) complete.trial.size <- trial.size

	all.etas <- matrix(lv.coefs[,1], nrow = n, ncol = p, byrow = TRUE)
	if(!is.null(lv)) all.etas <- all.etas + lv%*%t(as.matrix(lv.coefs[,2:(num.lv+1)]))
	if(!is.null(X.coefs)) all.etas <- all.etas + as.matrix(X)%*%t(X.coefs)
	## Assumes the columns in X.coefs corresponding to multinomial are set to 0

	index.multinom.cols <- which(complete.family == "multinom")
	for(j in 1:p) {
		species.etas <- all.etas[,j]
		if(!is.null(row.coefs)) { for(k in 1:ncol(row.ids)) species.etas <- species.etas + row.coefs[[k]][row.ids[,k]] }
		
		if(complete.family[j] == "binomial") 
			loglik.comp[,j] <- dbinom(as.vector(unlist(y[,j])), complete.trial.size[j], pnorm(species.etas), log = TRUE)
		if(complete.family[j] == "poisson") 
			loglik.comp[,j] <- dpois(as.vector(unlist(y[,j])), exp(species.etas), log = TRUE)
		if(complete.family[j] == "negative.binomial") 
			loglik.comp[,j] <- dnbinom(as.vector(unlist(y[,j])), size=1/(lv.coefs[j,ncol(lv.coefs)]+1e-5), mu=exp(species.etas), log = TRUE)
		if(complete.family[j] == "exponential") 
			loglik.comp[,j] <- dexp(as.vector(unlist(y[,j])), 1/exp(species.etas), log = TRUE)
		if(complete.family[j] == "gamma") 
			loglik.comp[,j] <- dgamma(as.vector(unlist(y[,j])), shape=exp(species.etas)*lv.coefs[j,ncol(lv.coefs)], rate = lv.coefs[j,ncol(lv.coefs)], log = TRUE)
		if(complete.family[j] == "beta") 
			loglik.comp[,j] <- dbeta(as.vector(unlist(y[,j])), lv.coefs[j,ncol(lv.coefs)]*exp(species.etas)/(1+exp(species.etas)), lv.coefs[j,ncol(lv.coefs)]*(1-exp(species.etas)/(1+exp(species.etas))),log = TRUE)
		if(complete.family[j] == "normal") 
			loglik.comp[,j] <- dnorm(as.vector(unlist(y[,j])), mean=species.etas, sd=(lv.coefs[j,ncol(lv.coefs)]+1e-6), log = TRUE)
		if(complete.family[j] == "lnormal") 
			loglik.comp[,j] <- dlnorm(as.vector(unlist(y[,j])), meanlog=species.etas, sdlog=(lv.coefs[j,ncol(lv.coefs)]+1e-6), log = TRUE)
			if(complete.family[j] == "tweedie") 
			loglik.comp[,j] <- dTweedie(as.vector(unlist(y[,j])), mu = exp(species.etas), phi = lv.coefs[j,ncol(lv.coefs)]+1e-6, p = powerparam, LOG = TRUE) 
		if(complete.family[j] == "ordinal") { 
			get.probs <- ordinal.conversion.spp(n = n, lv = lv, lv.coefs.j = lv.coefs[j,], num.lv = num.lv, row.coefs = row.coefs, row.ids = row.ids, X = X, X.coefs.j = X.coefs[j,], cutoffs = cutoffs); 
			for(i in 1:n) { loglik.comp[i,j] <- log(get.probs[i,as.vector(y[i,j])]+1e-5) } }	
# 		if(complete.family[j] == "multinom") { 
# 			if(!is.null(X.multinom.coefs)) spp.etas <- matrix(rep(species.etas,dim(X.multinom.coefs)[3]),nrow=n) + as.matrix(X)%*%X.multinom.coefs[which(index.multinom.cols == j),,]
# 			get.probs <- exp(spp.etas)/apply(exp(spp.etas),1,sum)
# 			for(i in 1:n) { loglik.comp[i,j] <- log(get.probs[as.vector(y[i,j])]+1e-5) } }	
		} 

	return(list(logLik = sum(loglik.comp), logLik.comp = loglik.comp)) 
	}

	
## Calculate logl for models with no latent variables
## Conditional and marginal log-likelihood are the same in such case, except that in the case of row.eff = "random" there is marginalization over the random row effect
## lv.coefs still need to be provided though at it contains the spp effects, and species-specific dispersion parameters
## Furthermore, note also that since it takes it X.coefs and lv.coefs, then: 1) marginalization is not done over the spp coefs it traits is not NULL; 2) marginalization is not done over the spp intercepts if > 1 are ordinal and hence the beta_{0j} are random effects...TOO DAMN HARD!
calc.logLik.lv0 <- function (y, X = NULL, family, trial.size = 1, lv.coefs, X.coefs = NULL, row.eff = "none", row.params = NULL, row.ids = NULL, cutoffs = NULL, powerparam = NULL) {
	if(is.null(lv.coefs)) 
		stop("lv.coefs must be given, as it contains the column-specific intercepts.")
	if(length(family) != ncol(y) & length(family) != 1) { 
		stop("Number of elements in family is either 1 or equal to # of columns in y") }
	if(length(family) == 1) complete.family <- rep(family, ncol(y))
	if(length(family) > 1) complete.family <- family
    

	if(row.eff == FALSE) row.eff <- "none"; 
	if(row.eff == TRUE) row.eff <- "fixed"
	if(!is.null(row.params)) {
		if(is.null(row.ids)) row.ids <- matrix(1:nrow(y), ncol = 1)
		if(!is.list(row.params))
			stop("row.params should be a list with length equal to the number of columns in row.ids.")
		if(length(row.params) != ncol(row.ids))
			stop("row.params should be a list with length equal to the number of columns in row.ids.")
		}
	if(!is.null(row.ids)) {
		row.ids <- as.matrix(row.ids)
		if(row.eff == "random" && (ncol(row.ids) > 1 || length(unique(row.ids[,1])) != nrow(y))) 
			stop("calc.logLik.lv0 currently only permits a single, unique random intercept for each row...sorry!")
		if(nrow(row.ids) != nrow(y)) 
			stop("Number of rows in the matrix row.ids should be equal to number of rows in y.")
		if(is.null(colnames(row.ids))) colnames(row.ids) <- paste0("ID", 1:ncol(row.ids))
		}
	if(row.eff == "fixed") { row.coefs <- row.params }
	if(row.eff == "none") { row.coefs <- NULL }


	if(any(complete.family == "ordinal") & is.null(cutoffs)) 
		stop("Ordinal data requires cutoffs to be supplied")
	if(any(family == "tweedie") & (powerparam < 1 || powerparam > 2)) 
		stop("Common power parameter for tweedie must be between 1 and 2.")
	if(length(trial.size) == 1) complete.trial.size <- rep(trial.size, ncol(y))
	if(length(trial.size) > 1) complete.trial.size <- trial.size
    
		
	## Checks done
	n <- nrow(y); p <- ncol(y); logl <- 0; logl.comp <- matrix(0, n, p)
	index.multinom.cols <- which(complete.family == "multinom")

	if(row.eff != "random") {
		for(j in 1:p) {
			eta <- lv.coefs[j, 1]
			if(!is.null(X.coefs)) eta <- eta + as.matrix(X) %*% X.coefs[j, ]
			if(!is.null(row.coefs)) { for(k in 1:ncol(row.ids)) eta <- eta + row.coefs[[k]][row.ids[,k]] }

			if(complete.family[j] == "poisson") logl.comp[, j] <- (dpois(as.vector(unlist(y[, j])), lambda = exp(eta), log = TRUE))
			if(complete.family[j] == "binomial") logl.comp[, j] <- (dbinom(as.vector(unlist(y[, j])), complete.trial.size[j], prob = pnorm(eta), log = TRUE))
			if(complete.family[j] == "negative.binomial") logl.comp[, j] <- (dnbinom(as.vector(unlist(y[, j])), mu = exp(eta), size = 1/(lv.coefs[j, 2]+1e-5), log = TRUE))
			if(complete.family[j] == "exponential") logl.comp[, j] <- (dexp(as.vector(unlist(y[, j])), rate = 1/exp(eta), log = TRUE))
			if(complete.family[j] == "gamma") logl.comp[, j] <- (dgamma(as.vector(unlist(y[, j])), shape = exp(eta) * lv.coefs[j, 2], rate = lv.coefs[j, 2], log = TRUE))
			if(complete.family[j] == "beta") logl.comp[, j] <- (dbeta(as.vector(unlist(y[, j])), lv.coefs[j, 2] * exp(eta)/(1 + exp(eta)), lv.coefs[j, 2] * (1 - exp(eta)/(1 + exp(eta))), log = TRUE))
			if(complete.family[j] == "normal") logl.comp[, j] <- (dnorm(as.vector(unlist(y[,j])), mean = eta, sd = (lv.coefs[j, 2]), log = TRUE))
			if(complete.family[j] == "lnormal") logl.comp[, j] <- (dlnorm(as.vector(unlist(y[,j])), meanlog = eta, sdlog = (lv.coefs[j,2]), log = TRUE))
			if(complete.family[j] == "tweedie") logl.comp[, j] <- (dTweedie(as.vector(unlist(y[,j])), mu = exp(eta), phi = lv.coefs[j, 2], p = powerparam, LOG = TRUE))
			if(complete.family[j] == "ordinal") {
				get.probs <- ordinal.conversion.spp(n = n, lv = NULL, lv.coefs.j = lv.coefs[j, ], num.lv = 0, row.coefs = row.coefs, row.ids = row.ids, X = X, X.coefs.j = X.coefs[j, ], cutoffs = cutoffs, est = "ignore")
				for(i in 1:n) { logl.comp[i, j] <- log(get.probs[i, as.vector(y[i,j])] + 1e-05) }
				}
				
				
# 		if(complete.family[j] == "multinom") {
# 			eta <- matrix(rep(eta,dim(X.multinom.coefs)[3]),nrow=n) + as.matrix(X)%*%X.multinom.coefs[which(index.multinom.cols == j),,]
# 			get.probs <- exp(eta)/apply(exp(eta),1,sum)
# 			for(i in 1:n) { logl.comp[i,j] <- log(get.probs[as.vector(y[i,j])]) } }
			}

		return(list(logLik = sum(logl.comp), logLik.comp = rowSums(logl.comp)))
		}

		
	if(row.eff == "random") {
		loglik <- 0; loglik.comp <- numeric(n); mc.row.eff <- 1000

		ordinal.conversion.special <- function(lv.coefs.j, row.coefs.i, X.i = NULL, X.coefs.j = NULL, cutoffs) {
			mc.row.eff <- length(row.coefs.i)
			etas <- matrix(NA, mc.row.eff, length(cutoffs))
			for(k in 1:length(cutoffs)) {
				etas[, k] <- cutoffs[k] - row.coefs.i - lv.coefs.j[1]
				if(!is.null(X.coefs.j)) 
				etas[, k] <- etas[, k] - t(as.matrix(X.i)) %*% X.coefs.j }
			probs <- matrix(NA, mc.row.eff, length(cutoffs) + 1)
			probs[, 1] <- pnorm(etas[,1])
			for(k in 2:ncol(etas)) { probs[, k] <- pnorm(etas[, k]) - pnorm(etas[,k-1]) }
			probs[, length(cutoffs) + 1] <- 1 - pnorm(etas[,length(cutoffs)])
			rm(etas)
			probs 
			}
        
		for(i in 1:n) {
			spp.f <- eta <- matrix(lv.coefs[, 1], nrow = mc.row.eff, ncol = p, byrow = TRUE)
			eta <- eta + rnorm(mc.row.eff, 0, sd = row.params[[1]]) 

			if(!is.null(X.coefs)) 
				eta <- eta + matrix(t(as.matrix(X[i, ])) %*% t(X.coefs), nrow = mc.row.eff, ncol = p, byrow = TRUE)
				for(j in 1:p) {
					if(complete.family[j] == "binomial") { spp.f[, j] <- dbinom(rep(as.vector(y[i, j]), mc.row.eff), complete.trial.size[j], pnorm(eta[, j])) }
					if(complete.family[j] == "poisson") { spp.f[, j] <- dpois(rep(as.vector(y[i, j]), mc.row.eff), exp(eta[, j])) }
					if(complete.family[j] == "negative.binomial") { spp.f[, j] <- dnbinom(rep(as.vector(y[i, j]), mc.row.eff), size =1/(lv.coefs[j, ncol(lv.coefs)]+1e-5), mu = exp(eta[, j])) }
					if(complete.family[j] == "exponential") { spp.f[, j] <- dexp(rep(as.vector(y[i, j]), mc.row.eff), rate = 1/exp(eta[, j])) }
					if(complete.family[j] == "gamma") { spp.f[, j] <- dgamma(rep(as.vector(y[i, j]), mc.row.eff), shape = exp(eta[, j]) * lv.coefs[j, ncol(lv.coefs)], rate = lv.coefs[j, ncol(lv.coefs)]) }
					if(complete.family[j] == "beta") { spp.f[, j] <- dbeta(rep(as.vector(y[i, j]), mc.row.eff), lv.coefs[j, ncol(lv.coefs)] * exp(eta[, j])/(1 + exp(eta[, j])), lv.coefs[j, ncol(lv.coefs)] * (1 - exp(eta[, j])/(1 + exp(eta[, j])))) }
					if(complete.family[j] == "normal") { spp.f[, j] <- dnorm(rep(as.vector(y[i, j]), mc.row.eff), mean = eta[, j], sd = (lv.coefs[j, ncol(lv.coefs)])) }
					if(complete.family[j] == "lnormal") { spp.f[, j] <- dlnorm(rep(as.vector(y[i, j]), mc.row.eff), meanlog = eta[, j], sdlog = (lv.coefs[j,ncol(lv.coefs)])) }
					if(complete.family[j] == "tweedie") { 
						spp.f[, j] <- dTweedie(rep(as.vector(y[i, j]), mc.row.eff), mu = exp(eta[, j]), phi = lv.coefs[j, ncol(lv.coefs)] + 1e-06, p = powerparam, LOG = FALSE)
						spp.f[, j][which(spp.f[, j] == 0)] <- 1 }
					if(complete.family[j] == "ordinal") {
						get.probs <- ordinal.conversion.special(lv.coefs.j = lv.coefs[j, ], row.coefs.i = rnorm(mc.row.eff, 0, sd = row.params[[1]]), X.i = X[i, ], X.coefs.j = X.coefs[j, ], cutoffs = cutoffs)
						spp.f[, j] <- get.probs[, as.vector(y[i, j])] + 1e-05 }
					}
				
				spp.f[!is.finite(spp.f)] = 1
				spp.f <- matrix(spp.f, nrow = mc.row.eff, byrow = FALSE)
				Q <- mean(apply(spp.f, 1, prod))
				loglik <- loglik + log(Q)
				loglik.comp[i] <- log(Q)
			}
        
        return(list(logLik = loglik, logLik.comp = loglik.comp)) 
        }
	}

	
## Calculate marginal logl
## loglik = sum_{i=1}^n \log( \int \prod_{j=1}^s f(y_ij|z_i) f(z_i) dz_i )
## Furthermore, note also that since it takes it X.coefs and lv.coefs, then: 1) marginalization is not done over the spp coefs it traits is not NULL; 2) marginalization is not done over the spp intercepts if > 1 are ordinal and hence the beta_{0j} are random effects...TOO DAMN HARD!
# lv.coefs = lv.coefs.mat; X.coefs = cw.X.coefs; row.coefs = cw.row.coefs; lv.mc = NULL; cutoffs = cw.cutoffs; X.multinom.coefs = NULL
calc.marglogLik <- function (y, X = NULL, family, trial.size = 1, lv.coefs, X.coefs = NULL, row.eff = "none", row.params = NULL, row.ids = NULL, num.lv, lv.mc = NULL, cutoffs = NULL, powerparam = NULL) {
	if(num.lv == 0) stop("Please use calc.loglik.lv0 to calculate likelihood in boral models with no latent variables.")
	if(is.null(lv.coefs)) stop("lv.coefs must be given. Please use calc.loglik.lv0 to calculate likelihood in boral models with no latent variables.")
    
	if(length(family) != ncol(y) & length(family) != 1) { 
		stop("Number of elements in family is either 1 or equal to # of columns in y") }
	if(length(family) == 1) complete.family <- rep(family, ncol(y))
	if(length(family) > 1) complete.family <- family
	if(any(family == "binomial") & !(length(trial.size) %in% c(1, length(family)))) 
	if(length(trial.size) == 1) complete.trial.size <- rep(trial.size, ncol(y))
	if(length(trial.size) > 1) complete.trial.size <- trial.size

	if(row.eff == FALSE) row.eff <- "none"; 
	if(row.eff == TRUE) row.eff <- "fixed"
	if(!is.null(row.params)) {
		if(is.null(row.ids)) row.ids <- matrix(1:nrow(y), ncol = 1)
		if(!is.list(row.params))
			stop("row.params should be a list with length equal to the number of columns in row.ids.")
		if(length(row.params) != ncol(row.ids))
			stop("row.params should be a list with length equal to the number of columns in row.ids.")
		}
	if(!is.null(row.ids)) {
		row.ids <- as.matrix(row.ids)
		if(row.eff == "random" && (ncol(row.ids) > 1 || length(unique(row.ids[,1])) != nrow(y))) 
			stop("calc.marglogLik currently only permits a single, unique random intercept for each row...sorry!")
		if(nrow(row.ids) != nrow(y)) 
			stop("Number of rows in the matrix row.ids should be equal to number of rows in y.")
		if(is.null(colnames(row.ids))) colnames(row.ids) <- paste0("ID", 1:ncol(row.ids))
		}
	if(row.eff == "fixed") { row.coefs <- row.params }

	if(any(complete.family == "ordinal") & is.null(cutoffs)) 
		stop("Ordinal data requires cutoffs to be supplied")
	if(any(family == "tweedie") & (powerparam < 1 || powerparam > 2)) 
		stop("Common power parameter for tweedie must be between 1 and 2.")
	#if(any(complete.family == "multinom") & is.null(X.multinom.coefs)) stop("Multinomial data requires X.multinom.coefs to be supplied.") 

		
	## Checks done
	n <- nrow(y); p <- ncol(y)
	loglik <- 0; loglik.comp <- numeric(n)
	if(is.null(lv.mc)) { lv.mc <- cbind(1, rmvnorm(2000, rep(0, num.lv))) }

	## Internal function - Given the coefficients and cutoffs, return the multinomial probabilities for proportional odds regression for element of y
	ordinal.conversion.special <- function(lv.mc, lv.coefs.j, num.lv, row.coefs.i = NULL, X.i = NULL, X.coefs.j = NULL, cutoffs) {
		etas <- matrix(NA, nrow(lv.mc), length(cutoffs))
		for(k in 1:length(cutoffs)) {
			etas[, k] <- lv.mc %*% c(cutoffs[k], -lv.coefs.j[2:(num.lv + 1)]) - row.coefs.i - lv.coefs.j[1]
			if(!is.null(X.coefs.j)) etas[, k] <- etas[, k] - t(as.matrix(X.i)) %*% X.coefs.j }
		
		probs <- matrix(NA, nrow(lv.mc), length(cutoffs) + 1)
		probs[, 1] <- pnorm(etas[,1])
		for(k in 2:ncol(etas)) { probs[, k] <- pnorm(etas[,k]) - pnorm(etas[,k-1]) }
		probs[, length(cutoffs) + 1] <- 1 - pnorm(etas[,length(cutoffs)])
		rm(etas)
		probs
		}
		
	index.multinom.cols <- which(complete.family == "multinom")
	for(i in 1:n) {
		spp.f <- matrix(NA, nrow = nrow(lv.mc), ncol = p)
		row.coefs.i <- 0
		if(row.eff == "fixed") { row.coefs.i <- 0; for(k in 1:ncol(row.ids)) row.coefs.i <- row.coefs.i + row.params[[k]][row.ids[i,k]] }
		if(row.eff == "random") row.coefs.i <- rnorm(nrow(lv.mc), 0, sd = row.params[[1]])
		spp.att.eta <- lv.mc %*% t(lv.coefs[, 1:(num.lv + 1)]) + row.coefs.i
		if(!is.null(X.coefs)) 
			spp.att.eta <- spp.att.eta + matrix(t(as.matrix(X[i, ])) %*% t(X.coefs), nrow = nrow(lv.mc), ncol = p, byrow = TRUE)
		
		for(j in 1:p) {
			if(complete.family[j] == "binomial") { 
				spp.f[, j] <- dbinom(rep(as.vector(y[i, j]), nrow(lv.mc)), complete.trial.size[j], pnorm(spp.att.eta[, j])) }
			if(complete.family[j] == "poisson") { 	
				spp.f[, j] <- dpois(rep(as.vector(y[i, j]), nrow(lv.mc)), exp(spp.att.eta[, j])) }
			if(complete.family[j] == "negative.binomial") { 
				spp.f[, j] <- dnbinom(rep(as.vector(y[i, j]), nrow(lv.mc)), size = 1/(lv.coefs[j, ncol(lv.coefs)]+1e-5), mu = exp(spp.att.eta[, j])) }
			if(complete.family[j] == "exponential") { 
				spp.f[, j] <- dexp(rep(as.vector(y[i, j]), nrow(lv.mc)), rate = 1/exp(spp.att.eta[, j])) }
			if(complete.family[j] == "gamma") { 
				spp.f[, j] <- dgamma(rep(as.vector(y[i, j]), nrow(lv.mc)), shape = exp(spp.att.eta[, j]) * lv.coefs[j, ncol(lv.coefs)], rate = lv.coefs[j, ncol(lv.coefs)]) }
			if(complete.family[j] == "beta") {
				spp.f[, j] <- dbeta(rep(as.vector(y[i, j]), nrow(lv.mc)), lv.coefs[j, ncol(lv.coefs)] * exp(spp.att.eta[, j])/(1 + exp(spp.att.eta[, j])), lv.coefs[j, ncol(lv.coefs)] * (1 - exp(spp.att.eta[, j])/(1 + exp(spp.att.eta[, j])))) }
			if(complete.family[j] == "normal") {
				spp.f[, j] <- dnorm(rep(as.vector(y[i, j]), nrow(lv.mc)), mean = spp.att.eta[, j], sd = (lv.coefs[j, ncol(lv.coefs)])) }
			if(complete.family[j] == "lnormal") {
				spp.f[, j] <- dlnorm(rep(as.vector(y[i, j]), nrow(lv.mc)), meanlog = spp.att.eta[, j], sdlog = (lv.coefs[j, ncol(lv.coefs)])) }
			if(complete.family[j] == "tweedie") {
				spp.f[, j] <- dTweedie(rep(as.vector(y[i, j]), nrow(lv.mc)), mu = exp(spp.att.eta[, j]), phi = lv.coefs[j, ncol(lv.coefs)] + 1e-06, p = powerparam, LOG = FALSE)
				spp.f[, j][which(spp.f[, j] == 0)] <- 1 }
			if(complete.family[j] == "ordinal") {
				get.probs <- ordinal.conversion.special(lv.mc, lv.coefs.j = lv.coefs[j, ], num.lv, row.coefs.i = row.coefs.i, X.i = X[i, ], X.coefs.j = X.coefs[j, ], cutoffs = cutoffs)
				spp.f[, j] <- get.probs[, as.vector(y[i, j])] + 1e-05 }
# 			if(complete.family[j] == "multinom") {
# 				num.multinom.levels <- dim(X.multinom.coefs)[3]
# 				if(!is.null(X.multinom.coefs)) { 
# 					spp.att.eta2 <- spp.att.eta[,j] + matrix(t(as.matrix(X[i,]))%*%X.multinom.coefs[which(index.multinom.cols == j),,],nrow(lv.mc),num.multinom.levels,byrow=TRUE) }
# 				get.probs <- exp(spp.att.eta2)/apply(exp(spp.att.eta2),1,sum)
# 				spp.f[,j] <- get.probs[,as.vector(y[i,j])]+1e-5 }
			}
        
		spp.f[!is.finite(spp.f)] = 1
		spp.f <- matrix(spp.f, nrow = nrow(lv.mc), byrow = FALSE)
		Q <- mean(apply(spp.f, 1, prod))	
		loglik <- loglik + log(Q)
		loglik.comp[i] <- log(Q)
		}

	return(list(logLik = loglik, logLik.comp = loglik.comp))
	}
	
	
create.life <- function (true.lv = NULL, lv.coefs, X = NULL, X.coefs = NULL, traits = NULL, traits.coefs = NULL, family, row.eff = "none", row.params = NULL, row.ids = NULL, trial.size = 1, cutoffs = NULL, powerparam = NULL, manual.dim = NULL, save.params = FALSE) {
	num.lv <- max(ncol(true.lv), 0)
	n <- max(nrow(true.lv), nrow(X))
	s <- max(nrow(lv.coefs), nrow(X.coefs), length(cutoffs))
	if(is.null(dim(lv.coefs))) { lv.coefs <- as.matrix(lv.coefs) }
	if((is.null(n) | is.null(s)) & is.null(manual.dim)) 
		stop("Sorry, but boral cannot determine the number of rows and columns for the response matrix. Please supply manual.dim as vector containing n and p.")
	if((is.null(n) | is.null(s)) & !is.null(manual.dim)) { n <- manual.dim[1]; s <- manual.dim[2] }
	
	if(!is.null(X)) { if(is.null(X.coefs) & is.null(traits))
			stop("If X is supplied, then either X.coefs or traits and traits.coefs must be supplied.")
			}
	if(!is.null(X)) { 
		if(!is.matrix(X)) X <- as.matrix(X)
		if(any(apply(X,2,function(x) all(x == 1)))) { stop("No intercept column should be included in X") } 
		}

	if((is.null(traits) & !is.null(traits.coefs)) | (!is.null(traits) & is.null(traits.coefs))) 
		stop("If traits is supplied, then traits.coefs must also be supplied.")
	if(!is.null(traits.coefs)) 
		message("Since trait.coefs has been supplied, then X.coefs will be ignored (X.coefs will instead be drawn as random effects based off trait.coefs)")
	if(!is.null(traits)) { 
		if(!is.matrix(traits)) traits <- as.matrix(traits) 
		if(any(apply(traits,2,function(x) all(x == 1)))) { stop("No intercept column should be included in traits. It will be included automatically") } 
		}


	if(length(family) != s & length(family) != 1)
		stop("Number of elements in family must be either 1 or equal to # of rows in lv.coefs/X.coefs/second number in manual.dim.")
	if(length(family) == 1) family <- rep(family, s)
	if(!all(family %in% c("negative.binomial", "poisson", "binomial", "normal", "lnormal", "tweedie", "ordinal", "exponential", "beta"))) 
		stop("One of the elements in family is not compatible with current version of boral...sorry!")
	if(any(family == "binomial") & !(length(trial.size) %in% c(1, length(family)))) 
		stop("trial.size needs to be specified if any columns are binomially distributed; can either be a single element or a vector equal to the # of rows in lv.coefs/X.coefs/second number in manual dim. The latter will assume the specified trial size for all rows labelled binomial in the family argument.")
	if(length(trial.size) == 1) trial.size <- rep(trial.size, s)
		
	if(any(family == "ordinal") & is.null(cutoffs)) 
		stop("cutoffs (an ascending vector of intercepts for proportional odds regression) must be supplied if any columns are ordinal data.")
	if(any(family == "ordinal")) 
		index.ord.cols <- which(family == "ordinal")
	if(!is.null(cutoffs)) {
		num.ord.levels <- length(cutoffs) + 1
		cutoffs <- sort(cutoffs)
		message("Sorting cutoffs...just in case") 
		}
	if(any(family == "tweedie") & is.null(powerparam)) 
		stop("Common powerparam must be supplied if any columns are tweedie data (Var = dispersion*mu^powerparam)")

	
	row.coefs <- NULL
	if(row.eff == FALSE) row.eff <- "none"; 
	if(row.eff == TRUE) row.eff <- "fixed"
	if(!is.null(row.params)) {
		if(is.null(row.ids)) row.ids <- matrix(1:n, ncol = 1)
		if(!is.list(row.params))
			stop("row.params should be a list with length equal to the number of columns in row.ids.")
		if(length(row.params) != ncol(row.ids))
			stop("row.params should be a list with length equal to the number of columns in row.ids.")
		}
	if(!is.null(row.ids)) {
		row.ids <- as.matrix(row.ids)
		if(nrow(row.ids) != n) 
			stop("Number of rows in the matrix row.ids should be equal to number of rows in y.")
		if(is.null(colnames(row.ids))) colnames(row.ids) <- paste0("ID", 1:ncol(row.ids))
		}
	if(row.eff == "fixed") { row.coefs <- row.params }
	if(row.eff == "random") {
		row.coefs <- vector("list", ncol(row.ids))
		for(k in 1:ncol(row.ids)) row.coefs[[k]] <- rnorm(length(unique(row.ids[,k])), mean = 0, sd = row.params[[k]]) 
		}

		
	if(num.lv > 5) 
		warnings("We won't stop you, but please consider if you really want more than five latent variables in the model.")
    
    
	sim.y <- matrix(NA, n, s)
	if(is.null(true.lv)) eta <- matrix(0,n,s)
	if(!is.null(true.lv)) eta <- true.lv%*%t(lv.coefs[,2:(num.lv + 1)])
	if(is.null(traits.coefs)) { eta <- eta + rep(1,n)%*%t(as.matrix(lv.coefs[,1])) }
	if(!is.null(X.coefs) & is.null(traits.coefs)) { eta <- eta + as.matrix(X) %*% t(X.coefs) }
	if(!is.null(traits.coefs)) {
		X.coefs <- matrix(0, s, ncol(X)) ## overwrite X.coefs
		lv.coefs[,1] <- rnorm(s, cbind(1,traits)%*%traits.coefs[1,-ncol(traits.coefs)], sd = traits.coefs[1,ncol(traits.coefs)]) ## overwrite spp-specific intercepts
          if(any(family == "ordinal")) {
               if(length(index.ord.cols) == 1) lv.coefs[index.ord.cols,1] <- 0 ## If there is just one ordinal column, then the random intercept for this column is zero for identifiability reasons
               }
		for(k in 1:ncol(X)) 
			X.coefs[,k] <- rnorm(s, cbind(1,traits)%*%traits.coefs[k+1,-ncol(traits.coefs)], sd = traits.coefs[k+1,ncol(traits.coefs)])
		eta <- eta + cbind(1,as.matrix(X))%*%t(cbind(lv.coefs[,1],X.coefs)) 
		}
	
	if(!is.null(row.coefs)) { for(k in 1:ncol(row.ids)) eta <- eta + row.coefs[[k]][row.ids[,k]] }
	
	for(j in 1:s) {
		if(family[j] == "binomial") sim.y[,j] <- rbinom(n, size = trial.size[j], prob = pnorm(eta[,j]))
		if(family[j] == "poisson") sim.y[,j] <- rpois(n, lambda = exp(eta[,j]))
		if(family[j] == "negative.binomial") sim.y[,j] <- rnbinom(n, mu = exp(eta[,j]), size = 1/(lv.coefs[j,ncol(lv.coefs)]+1e-5))
		if(family[j] == "exponential") sim.y[,j] <- rexp(n, rate = 1/exp(eta[, j]))
		if(family[j] == "gamma") sim.y[,j] <- rgamma(n, shape = exp(eta[, j])*lv.coefs[j,ncol(lv.coefs)], rate = lv.coefs[j, ncol(lv.coefs)])
		if(family[j] == "beta") 
			sim.y[, j] <- rbeta(n, shape1 = lv.coefs[j,ncol(lv.coefs)]*exp(eta[,j])/(1 + exp(eta[,j])), shape2 = lv.coefs[j,ncol(lv.coefs)]*(1-exp(eta[,j])/(1 +exp(eta[,j]))))
		if(family[j] == "normal") sim.y[, j] <- rnorm(n, mean = eta[, j], sd = (lv.coefs[j,ncol(lv.coefs)]))
		if(family[j] == "lnormal") sim.y[, j] <- rlnorm(n, meanlog = eta[, j], sdlog = (lv.coefs[j,ncol(lv.coefs)]))
		if(family[j] == "tweedie") sim.y[, j] <- rTweedie(n, mu = exp(eta[, j]), phi = lv.coefs[j,ncol(lv.coefs)], p = powerparam)
		if(family[j] == "ordinal") {
			get.probs <- ordinal.conversion.spp(n = n, lv = true.lv, lv.coefs.j = lv.coefs[j, ], num.lv = num.lv, row.coefs = row.coefs, row.ids = row.ids, X = X, X.coefs.j = X.coefs[j,], cutoffs = cutoffs, est = "ignore")
			for(i in 1:n) { sim.y[i, j] <- sample(1:num.ord.levels, 1, prob = get.probs[i,]) }
			}
		}
	
	if(!save.params) out <- sim.y
	if(save.params) out <- list(resp = sim.y, true.lv = true.lv, lv.coefs = lv.coefs, X.coefs = X.coefs, traits.coefs = traits.coefs, row.params = row.params, row.coefs = row.coefs, cutoffs = cutoffs, powerparam = powerparam)
	
	return(out)
	}
	
	
## Dunn-Smyth residuals
## Also create a confusion matrix for ordinal and multinomial data
ds.residuals <- function(object, est = "median") {  
	n <- object$n; p <- object$p; 
	num.lv <- object$num.lv; num.ord.levels <- object$num.ord.levels; 
	X <- object$X; y <- object$y
	mus <- fitted.boral(object, X, est = est)

	if(any(object$family == "ordinal")) {
		message("One or more columns of y have ordinal responses. Constructing a single confusion matrix for these.")
		true.resp <- as.matrix(y[,which(object$family == "ordinal")])
		pred.resp <- matrix(NA,n,ncol(true.resp)) 
		}
# 	if(any(object$family == "multinom")) {
# 		print("One or more columns of y have multinomial responses. Constructing a single confusion matrix for these.")
# 		true.multinom.resp <- as.matrix(y[,which(object$family == "multinom")])
# 		pred.multinom.resp <- matrix(NA,n,ncol(true.multinom.resp)) }
	if(any(object$family == "tweedie")) {
		if(est == "median") powerparam <- object$powerparam.median
		if(est == "mean") powerparam <- object$powerparam.mean 
		}

	ds.res.out <- matrix(NA,n,p)
	rownames(ds.res.out) <- rownames(y); colnames(ds.res.out) <- colnames(y)
	for(i in 1:n) { for(j in 1:p) {
		if(object$family[j] == "poisson") { 
			a <- ppois(as.vector(unlist(y[i,j]))-1, mus$out[i,j]); 
			b <- ppois(as.vector(unlist(y[i,j])), mus$out[i,j]); 		
			u <- runif(n = 1, min = a, max = b); ds.res.out[i,j] <- qnorm(u) 
			}
		if(object$family[j] == "negative.binomial") {
			if(est == "median") phis <- object$lv.coefs.median[,num.lv+2]+1e-5
			if(est == "mean") phis <- object$lv.coefs.mean[,num.lv+2]+1e-5
			a <- pnbinom(as.vector(unlist(y[i,j]))-1, mu=mus$out[i,j], size=1/phis[j]); 
			b <- pnbinom(as.vector(unlist(y[i,j])), mu=mus$out[i,j], size=1/phis[j])
			u <- runif(n = 1, min = a, max = b); ds.res.out[i,j] <- qnorm(u) 
			}
		if(object$family[j] == "binomial") { 
			a <- pbinom(as.vector(unlist(y[i,j]))-1, object$trial.size[j], mus$out[i,j]); 
			b <- pbinom(as.vector(unlist(y[i,j])), object$trial.size[j], mus$out[i,j])
			u <- runif(n = 1, min = a, max = b); ds.res.out[i,j] <- qnorm(u) 
			}
		if(object$family[j] == "exponential") { 
			a <- pexp(as.vector(unlist(y[i,j])), rate=1/mus$out[i,j]); 
			ds.res.out[i,j] <- qnorm(a) 
			}
		if(object$family[j] == "gamma") { 
			if(est == "median") phis <- object$lv.coefs.median[,num.lv+2]
			if(est == "mean") phis <- object$lv.coefs.mean[,num.lv+2]
			a <- pgamma(as.vector(unlist(y[i,j])), shape=mus$out[i,j]*phis[j], rate=phis[j]); 
			ds.res.out[i,j] <- qnorm(a) 
			}
		if(object$family[j] == "beta") { 
			if(est == "median") phis <- object$lv.coefs.median[,num.lv+2]
			if(est == "mean") phis <- object$lv.coefs.mean[,num.lv+2]
			a <- pbeta(as.vector(unlist(y[i,j])), shape1=phis[j]*mus$out[i,j], shape2=phis[j]*(1-mus$out[i,j])); ds.res.out[i,j] <- qnorm(a) }
		if(object$family[j] == "normal") { 
			if(est == "median") phis <- object$lv.coefs.median[,num.lv+2]
			if(est == "mean") phis <- object$lv.coefs.mean[,num.lv+2]
			a <- pnorm(as.vector(unlist(y[i,j])), mus$out[i,j], sd = (phis[j])); 
			ds.res.out[i,j] <- qnorm(a) 
			}
# 			X2 <- cbind(1,X); hatmat <- X2%*%solve(t(X2)%*%X2)%*%t(X2)
# 			ds.res.out[i,j] <- (y[i,j]-mus$out[i,j])/(sqrt(phis[j])*sqrt(1-hatmat[i,i])) }
		if(object$family[j] == "lnormal") { 
			if(est == "median") phis <- object$lv.coefs.median[,num.lv+2]
			if(est == "mean") phis <- object$lv.coefs.mean[,num.lv+2]
			a <- plnorm(as.vector(unlist(y[i,j])), log(mus$out[i,j]), sdlog = (phis[j])); ds.res.out[i,j] <- qnorm(a) 
			}
		if(object$family[j] == "tweedie") { 
			if(est == "median") phis <- object$lv.coefs.median[,num.lv+2]
			if(est == "mean") phis <- object$lv.coefs.mean[,num.lv+2]
			a <- pTweedie(as.vector(unlist(y[i,j])), mu = mus$out[i,j], phi = phis[j], p = powerparam); ds.res.out[i,j] <- qnorm(a) 
			}
		if(object$family[j] == "ordinal") { 
			pred.resp[,which(object$family == "ordinal")==j] <- mus$out[,which(object$family == "ordinal")==j] ## get max predicted probability
			cumsum.b <- sum(mus$ordinal.probs[i,j,1:(y[i,j])])
			cumsum.a <- sum(mus$ordinal.probs[i,j,1:(y[i,j]-1)])
			u <- runif(n = 1, min = cumsum.a, max = cumsum.b); 
			if(abs(u-1) < 1e-5) u <- 1; if(abs(u-0) < 1e-5) u <- 0
			ds.res.out[i,j] <- qnorm(u) 
			}
# 		if(object$family[j] == "multinom") { ## get max predicted probability
# 			pred.resp[i,which(object$family == "multinom")==j] <- which.max(mus$multinom.probs[i,j,]) }
 		} }

	if(sum(object$family == "ordinal") > 0) { 
		agree.tab <- table(as.vector(pred.resp), as.vector(true.resp)); 
		} 
	else { agree.tab <- NULL }
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
		fitted.ordinal.probs <- array(NA, dim=c(n,p,object$num.ord.levels)) 
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

	if(!is.null(object$row.ids) && est == "median") {
		for(j in 1:p) for(k in 1:ncol(object$row.ids)) eta[,j] <- eta[,j] + object$row.coefs[[k]]$median[object$row.ids[,k]] 
		}
	if(!is.null(object$row.ids) && est == "mean") {
		for(j in 1:p) for(k in 1:ncol(object$row.ids)) eta[,j] <- eta[,j] + object$row.coefs[[k]]$mean[object$row.ids[,k]] 
		}
	
	index.multinom.cols <- which(object$family == "multinom")
	for(j in 1:p) {
		if(object$family[j] %in% c("binomial")) fitted.out[,j] <- pnorm(eta[,j])
		if(object$family[j] %in% c("beta")) fitted.out[,j] <- exp(eta[,j])/(1+exp(eta[,j]))
		if(object$family[j] %in% c("poisson","lnormal","negative.binomial","tweedie","exponential","gamma")) fitted.out[,j] <- exp(eta[,j])
		if(object$family[j] == "normal") fitted.out[,j] <- (eta[,j]) 
# 		if(object$family[j] == "multinom") {
# 			if(est == "median") { if(!is.null(object$X.multinom.coefs.median)) eta2 <- eta[,j] + as.matrix(X)%*%object$X.multinom.coefs.median[which(index.multinom.cols == j),,] }
# 			if(est == "mean") { if(!is.null(object$X.multinom.coefs.mean)) eta2 <- eta[,j] + as.matrix(X)%*%object$X.multinom.coefs.mean[which(index.multinom.cols == j),,] }
# 			get.probs <- exp(eta2)/apply(exp(eta2),1,sum)	
# 			fitted.multinom.probs[,j,] <- get.probs
# 			}

		if(object$family[j] == "ordinal") {
			if(est == "median")
				fitted.ordinal.probs[,j,] <- ordinal.conversion.spp(n = n, lv = object$lv.median, lv.coefs.j = object$lv.coefs.median[j,], num.lv = num.lv, row.coefs = object$row.coefs, row.ids = object$row.ids, X = X, X.coefs.j = object$X.coefs.median[j,], cutoffs = object$cutoffs.median, est = "median")
			if(est == "mean")
				fitted.ordinal.probs[,j,] <- ordinal.conversion.spp(n = n, lv = object$lv.mean, lv.coefs.j = object$lv.coefs.mean[j,], num.lv = num.lv, row.coefs = object$row.coefs, row.ids = object$row.ids, X = X, X.coefs.j = object$X.coefs.mean[j,], cutoffs = object$cutoffs.mean, est = "mean")
			fitted.out[,j] <- apply(fitted.ordinal.probs[,j,],1,which.max) ## get max predicted probability
			}
		}	

	return(list(ordinal.probs = fitted.ordinal.probs, out = fitted.out))
	}

	
## Calculates DIC based on the conditional log-likelihood
get.dic <- function(jagsfit) { 
	jagsfit$BUGSoutput$DIC
	}
	
	
get.hpdintervals <- function(y, X = NULL, traits = NULL, row.ids = NULL, fit.mcmc, num.lv, prob = 0.95) {
	n <- nrow(y); p <- ncol(y)
	
	get.int <- HPDinterval(fit.mcmc, prob = prob); 
	hpd.lower <- get.int[,1]; hpd.upper <- get.int[,2]

	lv.coefs.arr <- abind(matrix(hpd.lower[grep("all.params",names(hpd.lower))], nrow=p), matrix(hpd.upper[grep("all.params",names(hpd.upper))], nrow=p), along = 3)	

	final.list <- list()
	
	if(num.lv > 0) {
		lv.arr <- abind(matrix(hpd.lower[grep("lvs", names(hpd.lower))], nrow=n), matrix(hpd.upper[grep("lvs", names(hpd.upper))], nrow=n), along = 3)
		dimnames(lv.arr) <- list(rows = rownames(y), lvs = paste0("LV", 1:num.lv), type = c("lower","upper"))		
		final.list$lv <- lv.arr

		if(dim(lv.coefs.arr)[2] == (num.lv+2)) 
			dimnames(lv.coefs.arr) <- list(cols = colnames(y), coefficients = c("beta0",paste0("theta",1:num.lv),"Dispersion"), type = c("lower","upper"))
		if(dim(lv.coefs.arr)[2] == (num.lv+1)) 
			dimnames(lv.coefs.arr) <- list(cols = colnames(y), coefficients = c("beta0",paste0("theta",1:num.lv)), type = c("lower","upper"))
		}
	if(num.lv == 0) { 
		if(dim(lv.coefs.arr)[2] == 2) 
			dimnames(lv.coefs.arr) <- list(cols = colnames(y), coefficients = c("beta0","Dispersion"), type = c("lower","upper"))
		if(dim(lv.coefs.arr)[2] == 1) 
			dimnames(lv.coefs.arr) <- list(cols = colnames(y), coefficients = c("beta0"), type = c("lower","upper"))
		}
	final.list$lv.coefs <- lv.coefs.arr
	
	if(length(grep("row.params", names(hpd.lower))) > 0) {
		n.ID <- apply(row.ids, 2, function(x) length(unique(x)))
		final.list$row.coefs <- vector("list", ncol(row.ids))
		names(final.list$row.coefs) <- colnames(row.ids)
		for(k in 1:ncol(row.ids)) {
			row.coefs.arr <- cbind(
					hpd.lower[grep(paste0("row.params.ID",k), names(hpd.lower))],
					hpd.upper[grep(paste0("row.params.ID",k), names(hpd.upper))])
			rownames(row.coefs.arr) <- 1:n.ID[k]; colnames(row.coefs.arr) <- c("lower","upper")
			
			final.list$row.coefs[[k]] <- row.coefs.arr
			}

		if(length(grep("row.ranef.sigma", names(hpd.lower))) > 0) { 
			final.list$row.sigma <- vector("list", ncol(row.ids))
			names(final.list$row.sigma) <- colnames(row.ids)
			for(k in 1:ncol(row.ids)) {
				row.sigma.vec <- c(
					hpd.lower[grep(paste0("row.ranef.sigma.ID",k), names(hpd.lower))],
					hpd.upper[grep(paste0("row.ranef.sigma.ID",k), names(hpd.upper))])
				names(row.sigma.vec) <- c("lower","upper")
				
				final.list$row.sigma[[k]] <- row.sigma.vec
				}
			}
		}

	if(length(grep("X.params", names(hpd.lower))) > 0) {
		X.coefs.arr <- abind(matrix(hpd.lower[grep("X.params", names(hpd.lower))],nrow=p), matrix(hpd.upper[grep("X.params", names(hpd.upper))],nrow=p), along = 3)
		dimnames(X.coefs.arr) <- list(cols = colnames(y), X = colnames(X), type = c("lower","upper"))
	
		final.list$X.coefs <- X.coefs.arr
		}

		
	if(length(grep("traits.params", names(hpd.lower))) > 0) { ## If T.params exists, then X.params are regressed against traits
		traits.coefs.arr <- abind(cbind(hpd.lower[grep("traits.int", names(hpd.lower))], matrix(hpd.lower[grep("traits.params", names(hpd.lower))],nrow=ncol(X)+1), hpd.lower[grep("sigma.trait", names(hpd.lower))]), cbind(hpd.upper[grep("traits.int", names(hpd.upper))], matrix(hpd.upper[grep("traits.params", names(hpd.upper))],nrow=ncol(X)+1), hpd.upper[grep("sigma.trait", names(hpd.upper))]), along = 3)
		dimnames(traits.coefs.arr) <- list(X.coefficients = c("beta0",colnames(X)), traits.coefficients = c("kappa0",colnames(traits),"sigma"), type = c("lower","upper"))
					
		final.list$traits.coefs <- traits.coefs.arr
		}

# 	if(length(grep("X.multinom.params", names(hpd.lower))) > 0) {
# 		final.list$X.multinom.coefs.lower <- array(matrix(hpd.lower[grep("X.multinom.params", names(hpd.lower))],dim=c(length(index.multinom.cols),ncol(X),ncol(all.X.multinom.coefs.lower)/ncol(X))), dimnames = list("1" = index.multinom.cols, "2" = colnames(X), "level" = 1:(ncol(all.X.multinom.coefs.upper)/ncol(X))))
# 		final.list$X.multinom.coefs.upper <- array(matrix(hpd.lower[grep("X.multinom.params", names(hpd.upper))],dim=c(length(index.multinom.cols),ncol(X),ncol(all.X.multinom.coefs.lower)/ncol(X))), dimnames = list("1" = index.multinom.cols, "2" = colnames(X), "level" = 1:(ncol(all.X.multinom.coefs.upper)/ncol(X))))
# 		}

	if(length(grep("alpha", names(hpd.lower))) > 0) { ## If alpha exists, then cutoffs are there and some columns involved ordinal responses
		cutoffs.arr <- cbind(hpd.lower[grep("alpha", names(hpd.lower))], hpd.upper[grep("alpha", names(hpd.upper))])
		num.ord.levels <- nrow(cutoffs.arr) + 1
		rownames(cutoffs.arr) <- paste0(1:(num.ord.levels-1),"|",2:num.ord.levels)
		colnames(cutoffs.arr) <- c("lower","upper")
	
		final.list$cutoffs <- cutoffs.arr

		if(length(grep("ordinal.ranef.sigma", names(hpd.lower))) > 0) { 
			ordinal.sigma.vec <- c(hpd.lower[grep("ordinal.ranef.sigma", names(hpd.lower))], hpd.upper[grep("ordinal.ranef.sigma", names(hpd.upper))])
			names(ordinal.sigma.vec) <- c("lower","upper")
			final.list$ordinal.sigma <- ordinal.sigma.vec
			}
		}
				
	if(length(grep("powerparams", names(hpd.lower))) > 0) { ## If powerparam exists, then power parameters are there and some columns involved tweedie responses
		powerparam.vec <- c(hpd.lower[grep("powerparam", names(hpd.lower))], hpd.upper[grep("powerparam", names(hpd.upper))])
		names(powerparam.vec) <- c("lower","upper")
		final.list$powerparam <- powerparam.vec
		}

	rm(list = ls(pattern = ".arr"))
	return(final.list) 
	}

	
## Calculates conditional WAIC, EAIC, EBIC. 
## Also calculate the marginal likelihood at component medians, and bases a AIC and BIC on this. Note this in cases were calc.marglogl and calc.logLik.lv0 actually produce a sensible result
get.measures <- function(y, X = NULL, family, trial.size = 1, row.eff = "none", row.ids = NULL, num.lv, fit.mcmc) {
	do.marglik.ics <- TRUE

	if(length(family) != ncol(y) & length(family) != 1) stop("Number of elements in family is either 1 or equal to # of columns in y")
	if(length(family) == 1) complete.family <- rep(family, ncol(y))
	if(length(family) > 1) complete.family <- family
	index.ordinal.cols <- which(family == "ordinal")
	
	if(any(family == "binomial") & !(length(trial.size) %in% c(1, length(family)))) 
        stop("trial.size needs to be specified if any columns are binomially distributed; can either be a single element or a vector equal to the # of columns in y. The latter will assume the specified trial size for all rows labelled binomial in the family argument.")
	
	if(row.eff == FALSE) row.eff <- "none"; 
	if(row.eff == TRUE) row.eff <- "fixed"
	if(!is.null(row.ids)) {
		row.ids <- as.matrix(row.ids)
		if(nrow(row.ids) != nrow(y)) stop("Number of rows in the matrix row.ids should be equal to number of rows in y.")
		if(row.eff == "random" && (ncol(row.ids) > 1 || length(unique(row.ids[,1])) != nrow(y))) do.marglik.ics <- FALSE
		if(is.null(colnames(row.ids))) colnames(row.ids) <- paste0("ID", 1:ncol(row.ids))
		}

	if(length(grep("traits.params", colnames(fit.mcmc))) > 1) { do.marglik.ics <- FALSE }	
	if(length(index.ordinal.cols) > 1) { do.marglik.ics <- FALSE }	
		
	## Checks done	
	n <- nrow(y); p <- ncol(y)
	all.lppd <- matrix(NA, nrow = nrow(fit.mcmc), ncol = n*p)
	index.multinom.cols <- which(complete.family == "multinom")

	for(t in 1:nrow(fit.mcmc)) {
		cw.params <- list(lv.coefs = matrix(fit.mcmc[t, grep("all.params", colnames(fit.mcmc))], nrow = p))
		if(row.eff != "none") {
			cw.params$row.coefs <- vector("list", ncol(row.ids))
			for(k in 1:ncol(row.ids)) cw.params$row.coefs[[k]] <- fit.mcmc[t, grep(paste0("row.params.ID",k), colnames(fit.mcmc))] 
			}
		if(!is.null(X)) cw.params$X.coefs <- matrix(fit.mcmc[t, grep("X.params", colnames(fit.mcmc))], nrow = p) 
		if(any(complete.family == "ordinal")) cw.params$cutoffs <- fit.mcmc[t, grep("alpha", colnames(fit.mcmc))] 
		if(any(complete.family == "tweedie")) cw.params$powerparam <- fit.mcmc[t, grep("powerparam", colnames(fit.mcmc))] 
# 		if(any(complete.family == "multinom") & !is.null(X)) { 
# 			get.X.multinom.coefs <- array(matrix(fit.mcmc[t,grep("X.multinom.params", colnames(fit.mcmc))],dim = c(length(index.multinom.cols),ncol(X),ncol(get.X.multinom.coefs)/ncol(X))), dimnames=NULL) } else { get.X.multinom.coefs <- NULL }
		
		if(num.lv > 0) {
			get.out <- calc.condlogLik(y, X, complete.family, trial.size, lv.coefs = cw.params$lv.coefs, X.coefs = cw.params$X.coefs, row.coefs = cw.params$row.coefs, row.ids = row.ids, lv = matrix(fit.mcmc[t, grep("lvs", colnames(fit.mcmc))], nrow = n), cutoffs = cw.params$cutoffs, powerparam = cw.params$powerparam)
			}
		if(num.lv == 0) {
			get.out <- calc.condlogLik(y, X, complete.family, trial.size, lv.coefs = cw.params$lv.coefs, X.coefs = cw.params$X.coefs, row.coefs = cw.params$row.coefs, row.ids = row.ids, lv = NULL, cutoffs = cw.params$cutoffs, powerparam = cw.params$powerparam)
			}
			
		get.out$logLik.comp[!is.finite(get.out$logLik.comp)] <- NA
		all.lppd[t,] <- as.vector(get.out$logLik.comp)
		}
		
		
	all.cond.logl <- rowSums(all.lppd)
	waic.out <- -2 * sum(log(colMeans(exp(all.lppd), na.rm = TRUE))) + 2 * sum(apply(all.lppd, 2, var, na.rm = TRUE))
	cond.num.params <- sum(cw.params$lv.coefs != 0) + n*num.lv + ## LVs and loadings
		sum(cw.params$X.coefs != 0)*as.numeric(!is.null(X)) + ## X.coefs
		any(complete.family == "ordinal")*sum(cw.params$cutoffs != 0) + (1 - is.null(cw.params$powerparam)) ## other parameters
	if(row.eff != "none") { 
		cond.num.params <- cond.num.params + sum(sapply(cw.params$row.coefs,length)) 
		} ## row effects

	eaic <- -2*mean(all.cond.logl, na.rm = TRUE) + 2*cond.num.params
	ebic <- -2*mean(all.cond.logl, na.rm = TRUE) + log(n*p)*cond.num.params
	out.list <- list(waic = waic.out, eaic = eaic, ebic = ebic, all.cond.logLik = all.cond.logl, cond.num.params = cond.num.params, do.marglik.ics = do.marglik.ics)

	if(do.marglik.ics) {
		## Calculate marginal logL at component medians
		params.median <- list(lv.coefs = matrix(apply(fit.mcmc[, grep("all.params", colnames(fit.mcmc))], 2, median), nrow = p))
		if(row.eff == "fixed") {
			params.median$row.coefs <- vector("list", ncol(row.ids))
			for(k in 1:ncol(row.ids)) params.median$row.coefs[[k]] <- apply(fit.mcmc[, grep(paste0("row.params.ID",k), colnames(fit.mcmc))], 2, median)
			}
		if(row.eff == "random") {
			params.median$row.coefs <- vector("list", ncol(row.ids))
			for(k in 1:ncol(row.ids)) params.median$row.coefs[[k]] <- median(fit.mcmc[, grep(paste0("row.ranef.sigma.ID",k), colnames(fit.mcmc))])
			}
		if(!is.null(X)) 
			params.median$X.coefs <- matrix(apply(fit.mcmc[, grep("X.params", colnames(fit.mcmc))], 2, median), nrow = p) 
		if(any(complete.family == "ordinal")) 
			params.median$cutoffs <- apply(fit.mcmc[, grep("alpha", colnames(fit.mcmc))], 2, median)
		if(any(complete.family == "tweedie")) 
			params.median$powerparam <- median(fit.mcmc[, grep("powerparam", colnames(fit.mcmc))]) 
	# 	if(any(complete.family == "multinom") & !is.null(X)) { 
	# 		get.X.multinom.coefs <- array(matrix(apply(fit.mcmc[,grep("X.multinom.params", colnames(fit.mcmc))],2,median),dim = c(length(index.multinom.cols),ncol(X),ncol(get.X.multinom.coefs)/ncol(X))), dimnames=NULL) } else { get.X.multinom.coefs <- NULL }

		if(num.lv > 0) {
			median.marglogl <- calc.marglogLik(y, X, complete.family, trial.size, lv.coefs = params.median$lv.coefs, X.coefs = params.median$X.coefs, row.eff = row.eff, row.ids = row.ids, row.params = params.median$row.coefs, num.lv = num.lv, lv.mc = NULL, cutoffs = params.median$cutoffs, powerparam = params.median$powerparam) 
			}
		if(num.lv == 0) 
			median.marglogl <- calc.logLik.lv0(y, X, complete.family, trial.size, lv.coefs = params.median$lv.coefs, X.coefs = params.median$X.coefs, row.eff = row.eff, row.ids = row.ids, row.params = params.median$row.coefs, cutoffs = params.median$cutoffs, powerparam = params.median$cw.powerparam)	

		marg.num.params <- sum(params.median$lv.coefs != 0) + ## Loadings 
			sum(params.median$X.coefs != 0)*as.numeric(!is.null(X)) + ## X coefs
			any(complete.family == "ordinal")*sum(params.median$cutoffs != 0) + (1 - is.null(params.median$powerparam)) ## other parameters
		if(row.eff != "none") { 
			marg.num.params <- marg.num.params + sum(sapply(params.median$row.coefs,length)) 
			} ## row effects

		marg.aic <- -2 * median.marglogl$logLik + 2 * marg.num.params
		marg.bic <- -2 * median.marglogl$logLik + log(n*p) * marg.num.params
		
		out.list$marg.num.params <- marg.num.params; out.list$aic.median <- marg.aic; out.list$bic.median <- marg.bic
		}
		
		
	return(out.list)
	}
	

## Calculates marginal logl for all samples to produce a proper AIC and BIC. 
## Calculates WAIC based on the marginal likelihood; DIC based on the marginal likelihood
## All of this is only permitted when the fitted boral model has a simple enough structure to allow these calculations!
get.more.measures <- function(y, X = NULL, family, trial.size = 1, row.eff = "none", row.ids = NULL, num.lv, fit.mcmc, verbose = TRUE) {
	do.marglik.ics <- TRUE
	index.ordinal.cols <- which(family == "ordinal")

	if(num.lv == 0) 
		stop("For boral models with no latent variables, the marginal and conditional likelihoods are equivalent, and there is nothing to gain from using get.more.measures")
	if(length(family) != ncol(y) & length(family) != 1) 
		stop("Number of elements in family is either 1 or equal to # of columns in y")
	if(length(family) == 1) complete.family <- rep(family, ncol(y))
	if(length(family) > 1) complete.family <- family
	
	if(any(family == "binomial") & !(length(trial.size) %in% c(1, length(family)))) 
		stop("trial.size needs to be specified if any columns are binomially distributed; can either be a single element or a vector equal to the # of columns in y. The latter will assume the specified trial size for all rows labelled binomial in the family argument.")

	if(row.eff == FALSE) row.eff <- "none"
	if(row.eff == TRUE) row.eff <- "fixed"
	if(!is.null(row.ids)) {
		row.ids <- as.matrix(row.ids)
		if(nrow(row.ids) != nrow(y)) stop("Number of rows in the matrix row.ids should be equal to number of rows in y.")
		if(row.eff == "random" && (ncol(row.ids) > 1 || length(unique(row.ids[,1])) != nrow(y))) do.marglik.ics <- FALSE
		if(is.null(colnames(row.ids))) colnames(row.ids) <- paste0("ID", 1:ncol(row.ids))
		}

	if(length(grep("traits.params", colnames(fit.mcmc))) > 1) { do.marglik.ics <- FALSE }	
	if(length(index.ordinal.cols) > 1) { do.marglik.ics <- FALSE }	
        
	## Checks done
	if(!do.marglik.ics) {
		message("The current version of boral does not implement information criterion based on the marginal likelihood for the specified model, because the number of random effects included in the model and/or the random effects structure is too complicated...sorry!")
		return()
		}
	
	n <- nrow(y); p <- ncol(y)
	big.lv <- cbind(1, rmvnorm(2000, rep(0, num.lv)))
	all.marg.logl <- matrix(NA, nrow(fit.mcmc), n)
	index.multinom.cols <- which(complete.family == "multinom")
	
	## Calculate marginal likelihood at all iterations
	for(t in 1:nrow(fit.mcmc)) {
		if(verbose == TRUE & t%%100 == 0) message("Onto mcmc sample ", t)
		cw.params <- list(lv.coefs = matrix(fit.mcmc[t, grep("all.params", colnames(fit.mcmc))], nrow = p))		
		if(row.eff == "fixed") {
			cw.params$row.coefs <- vector("list", ncol(row.ids))
			for(k in 1:ncol(row.ids)) cw.params$row.coefs[[k]] <- fit.mcmc[t, grep(paste0("row.params.ID",k), colnames(fit.mcmc))]
			}
		if(row.eff == "random") {
			cw.params$row.coefs <- vector("list", ncol(row.ids))
			for(k in 1:ncol(row.ids)) cw.params$row.coefs[[k]] <- fit.mcmc[t, grep(paste0("row.ranef.sigma.ID",k), colnames(fit.mcmc))]
			}
		if(!is.null(X)) cw.params$X.coefs <- matrix(fit.mcmc[t, grep("X.params", colnames(fit.mcmc))], nrow = p) 
		if(any(complete.family == "ordinal")) cw.params$cutoffs <- fit.mcmc[t, grep("alpha", colnames(fit.mcmc))] 
		if(any(complete.family == "tweedie")) cw.params$powerparam <- fit.mcmc[t, grep("powerparam", colnames(fit.mcmc))] 
# 		if(any(complete.family == "multinom") & !is.null(X)) { 
# 		get.X.multinom.coefs <- array(matrix(fit.mcmc[t,grep("X.multinom.params", colnames(fit.mcmc))],dim = c(length(index.multinom.cols),ncol(X),ncol(get.X.multinom.coefs)/ncol(X))), dimnames = NULL) } else { get.X.multinom.coefs <- NULL }

		get.mll <- calc.marglogLik(y, X, complete.family, trial.size, lv.coefs = cw.params$lv.coefs, X.coefs = cw.params$X.coefs, row.eff = row.eff, row.ids = row.ids, row.params = cw.params$row.coefs, num.lv, lv.mc = big.lv, cutoffs = cw.params$cutoffs, powerparam = cw.params$powerparam)
		all.marg.logl[t,] <- get.mll$logLik.comp
		}
				
	## Calculate WAIC based on marginal
	marg.waic <- -2 * sum(log(apply(exp(all.marg.logl), 2, mean, na.rm = TRUE))) + 2 * sum(apply(all.marg.logl, 2, var, na.rm = TRUE))

	## Calculate AIC, BIC at posterior mode	
	marg.num.params <- sum(cw.params$lv.coefs != 0) + ## Loadings 
		sum(cw.params$X.coefs != 0)*as.numeric(!is.null(X)) +  ## X coefs
		any(complete.family == "ordinal")*sum(cw.params$cutoffs != 0) + (1 - is.null(cw.params$powerparam)) ## other parameters	
	if(row.eff != "none") { marg.num.params <- marg.num.params + sum(sapply(cw.params$row.coefs,length)) } ## row effects
	bic1 <- -2 * max(rowSums(all.marg.logl)) + log(n)*marg.num.params
	aic1 <- -2 * max(rowSums(all.marg.logl)) + 2*marg.num.params
	
	## Calculate DIC based on marginal
	params.mean <- list(lv.coefs = matrix(apply(fit.mcmc[, grep("all.params", colnames(fit.mcmc))], 2, mean), nrow = p))
	if(row.eff == "fixed") {
		params.mean$row.coefs <- vector("list", ncol(row.ids))
		for(k in 1:ncol(row.ids)) params.mean$row.coefs[[k]] <- fit.mcmc[t, grep(paste0("row.params.ID",k), colnames(fit.mcmc))]
		}
	if(row.eff == "random") {
		params.mean$row.coefs <- vector("list", ncol(row.ids))
		for(k in 1:ncol(row.ids)) params.mean$row.coefs[[k]] <- fit.mcmc[t, grep(paste0("row.ranef.sigma.ID",k), colnames(fit.mcmc))]
		}
	if(!is.null(X)) params.mean$X.coefs <- matrix(apply(fit.mcmc[, grep("X.params", colnames(fit.mcmc))], 2, mean), nrow = p) 
	if(any(complete.family == "ordinal")) params.mean$cutoffs <- apply(fit.mcmc[, grep("alpha", colnames(fit.mcmc))], 2, mean)
	if(any(complete.family == "tweedie")) params.mean$powerparam <- mean(fit.mcmc[, grep("powerparam", colnames(fit.mcmc))]) 
	
	marg.dic <- -2*calc.marglogLik(y, X, complete.family, trial.size, lv.coefs = params.mean$lv.coefs, X.coefs = params.mean$X.coefs, row.eff = row.eff, row.ids = row.ids, row.params = params.mean$row.coefs, num.lv, lv.mc = big.lv, cutoffs = params.mean$cutoffs, powerparam = params.mean$powerparam)$logLik
	marg.dic <- marg.dic + 2*(2*var(rowSums(all.marg.logl), na.rm = TRUE))

	return(list(aic.mode = aic1, bic.mode = bic1, marg.dic = marg.dic, marg.waic = marg.waic, all.marg.logLik = rowSums(all.marg.logl), marg.num.params = marg.num.params))
	}
	

## Produce the correlation due to similarity of responses to X
get.enviro.cor <- function(object, est = "median", prob = 0.95) {
	fit.mcmc <- object$jags.model$BUGSoutput
	if(is.null(fit.mcmc)) stop("MCMC samples not found")
	fit.mcmc <- mcmc(object$jags.model$BUGSoutput$sims.matrix, start = 1, thin = object$mcmc.control$n.thin)
	y <- object$y; X <- object$X
	
	if(length(grep("X.params", colnames(fit.mcmc))) == 0) stop("Cannot find MCMC sample corresponding to coefficients for X.")

	n <- nrow(y); p <- ncol(y)
	enviro.cor.mat <- enviro.cov.mat <- matrix(0,p,p)
	sig.enviro.cor.mat <- matrix(0,p,p)
	if(is.null(colnames(y))) colnames(y) <- 1:ncol(y); 
	rownames(enviro.cor.mat) <- rownames(enviro.cov.mat) <- rownames(sig.enviro.cor.mat) <- colnames(y)
	colnames(enviro.cor.mat) <- colnames(enviro.cov.mat) <- colnames(sig.enviro.cor.mat) <- colnames(y)
	all.enviro.cov.mat <- all.enviro.cor.mat <- array(0,dim=c(nrow(fit.mcmc),p,p))

	
	for(t in 1:nrow(fit.mcmc)) {
		cw.X.coefs <- matrix(fit.mcmc[t,grep("X.params", colnames(fit.mcmc))],nrow=p)
		enviro.linpreds <- X%*%t(as.matrix(cw.X.coefs))
		all.enviro.cov.mat[t,,] <- cov(enviro.linpreds)
		all.enviro.cor.mat[t,,] <- cor(enviro.linpreds) 
		}

	for(j in 1:p) { for(j2 in 1:p) { ## Average/Median over the MCMC samples
		if(est == "median") { 
			enviro.cov.mat[j,j2] <- median(all.enviro.cov.mat[,j,j2])
			enviro.cor.mat[j,j2] <- median(all.enviro.cor.mat[,j,j2]) }
		if(est == "mean") {
			enviro.cov.mat[j,j2] <- mean(all.enviro.cov.mat[,j,j2])
			enviro.cor.mat[j,j2] <- mean(all.enviro.cor.mat[,j,j2]) } 
		
		sig.enviro.cor.mat[j,j2] <- enviro.cor.mat[j,j2]
		get.hpd.cors <- HPDinterval(as.mcmc(all.enviro.cor.mat[,j,j2]), prob = prob)
		if(0 > get.hpd.cors[1] & 0 < get.hpd.cors[2]) sig.enviro.cor.mat[j,j2] <- 0
		} }
		
	#return(list(residual.correlations = enviro.cor.mat))
	#corrplot(enviro.cor.mat, title = "Environmental correlations", type = "lower")
	return(list(cor = enviro.cor.mat, sig.cor = sig.enviro.cor.mat, cov = enviro.cov.mat))
	}

	
## Produce the residual correlation based on latent variables
get.residual.cor <- function(object, est = "median", prob = 0.95) {
	fit.mcmc <- object$jags.model$BUGSoutput
	if(is.null(fit.mcmc)) stop("MCMC samples not found")
	fit.mcmc <- mcmc(object$jags.model$BUGSoutput$sims.matrix, start = 1, thin = object$mcmc.control$n.thin)
	y <- object$y; X <- object$X
	num.lv <- object$num.lv

	if(length(grep("lvs", colnames(fit.mcmc))) == 0) stop("Cannot find MCMC samples corresponding to latent variables.")

	n <- nrow(y); p <- ncol(y)
	sig.rescor.mat <- rescor.mat <- rescov.mat <- matrix(0,p,p)
	if(is.null(colnames(y))) colnames(y) <- 1:ncol(y); 
	rownames(rescor.mat) <- colnames(rescor.mat) <- colnames(y)
	rownames(sig.rescor.mat) <- colnames(sig.rescor.mat) <- colnames(y)
	rownames(rescov.mat) <- colnames(rescov.mat) <- colnames(y)
	all.rescor.mat <- all.rescov.mat <- array(0,dim=c(nrow(fit.mcmc),p,p))
	all.trace.rescor <- numeric(nrow(fit.mcmc))

	for(t in 1:nrow(fit.mcmc)) {
		#lvs <- matrix(fit.mcmc[t,grep("lvs", colnames(fit.mcmc))],nrow=n)
		lvs.coefs <- matrix(fit.mcmc[t,grep("all.params", colnames(fit.mcmc))],nrow=p)
# 		if(all(object$family == "binomial") & all(object$trial.size == 1)) 
# 			lvs.coefs[,2:(num.lv+1)] <- lvs.coefs[,2:(num.lv+1)]/matrix(sqrt(1-rowSums(lvs.coefs[,2:(num.lv+1)]^2)),nrow=p,ncol=num.lv,byrow=FALSE) ## If data is Bernoulli, then scale the coefficients to acocunt for constraints (see Knott and Bartholomew, Chapter 4)
		
		lambdalambdaT <- as.matrix(lvs.coefs[,2:(num.lv+1)])%*%t(as.matrix(lvs.coefs[,2:(num.lv+1)]))
		all.rescov.mat[t,,] <- (lambdalambdaT) 
		all.trace.rescor[t] <- sum(diag(lambdalambdaT))
		
 		if(all(object$family == "negative.binomial")) {
   			get.var.phis <- numeric(p); 
   			## Multiplicative Poisson gamma model implies a log gamma random effect on the linear predictors
   			for(j in 1:p) 
				get.var.phis[j] <- var(log(rgamma(2000,shape=1/lvs.coefs[j,ncol(lvs.coefs)],rate=1/lvs.coefs[j,ncol(lvs.coefs)])))
			all.rescov.mat[t,,] <- lambdalambdaT + diag(x=get.var.phis,nrow=p)
# #   			for(j in 1:p) { for(j2 in 1:p) { all.rescor.mat[t,j,j2] <- lambdalambdaT[j,j2]/sqrt((lambdalambdaT[j,j]+get.var.phis[j])*(lambdalambdaT[j2,j2]+get.var.phis[j2])) } }
# #   			all.trace.rescor[t] <- sum(diag(lambdalambdaT)+get.var.phis)
#   			}
			}
		all.rescor.mat[t,,] <- cov2cor(all.rescov.mat[t,,]) 
		}
		
	for(j in 1:p) { for(j2 in 1:p) { ## Average/Median over the MCMC samples
		if(est == "median") { rescor.mat[j,j2] <- median(all.rescor.mat[,j,j2]); rescov.mat[j,j2] <- median(all.rescov.mat[,j,j2]) }
		if(est == "mean") { rescor.mat[j,j2] <- mean(all.rescor.mat[,j,j2]); rescov.mat[j,j2] <- mean(all.rescov.mat[,j,j2]) }
		
		sig.rescor.mat[j,j2] <- rescor.mat[j,j2]
		get.hpd.cors <- HPDinterval(as.mcmc(all.rescor.mat[,j,j2]), prob = 0.95)
		if(0 > get.hpd.cors[1] & 0 < get.hpd.cors[2]) sig.rescor.mat[j,j2] <- 0
		} }

	if(est == "median") final.trace <- median(all.trace.rescor)
	if(est == "mean") final.trace <- mean(all.trace.rescor) 	
		
	#return(list(residual.correlations = rescor.mat))
	#corrplot(rescor.mat, title = "Residual correlations", type = "lower")
	return(list(cor = rescor.mat, sig.cor = sig.rescor.mat, cov = rescov.mat, trace = final.trace))
	}
		
		
## Wrapper for create.life: Takes a boral model and applies create.life to it if possible
simulate.boral <- function(object, nsim = 1, seed = NULL, est = "median", ...) {
	if(class(object) != "boral") { stop("object must be of class boral. Thanks!") }
	
	if(est == "mean") {
		true.mod <- list(lv.coefs = object$lv.coefs.mean, lvs = object$lv.mean, X.coefs = object$X.coefs.mean, traits = object$traits, traits.coefs = object$traits.coefs.mean, cutoffs = object$cutoffs.mean, powerparam = object$powerparam.mean) 
		if(object$row.eff == "fixed") { true.mod$row.params <- object$row.coefs.mean }
		if(object$row.eff == "random") { true.mod$row.params <- object$row.sigma.mean }
		}
	
	if(est == "median") {
		true.mod <- list(lv.coefs = object$lv.coefs.median, lvs = object$lv.median, X.coefs = object$X.coefs.median, traits = object$traits, traits.coefs = object$traits.coefs.median, cutoffs = object$cutoffs.median, powerparam = object$powerparam.median) 
		if(object$row.eff == "fixed") { true.mod$row.params <- object$row.coefs.mean }
		if(object$row.eff == "random") { true.mod$row.params <- object$row.sigma.mean }
		}

	if(!is.null(seed)) set.seed(seed)
	
	out <- replicate(nsim, create.life(true.lv = true.mod$lvs, lv.coefs = true.mod$lv.coefs, X = object$X, X.coefs = true.mod$X.coefs, traits = object$traits, traits.coefs = true.mod$traits.coefs, family = object$family, row.eff = object$row.eff, row.params = true.mod$row.params, trial.size = object$trial.size, cutoffs = true.mod$cutoffs, powerparam = true.mod$powerparam))
		
	return(out)
	}

	

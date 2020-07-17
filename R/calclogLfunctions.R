##########
## Auxilary functions to calculate log-likelihood values
##########

## Calculate conditional logl
## loglik = sum_{i=1}^n sum_{j=1}^s \log( f(y_ij|z_i) )
#  lv.coefs = matrix(fit.mcmc[t,grep("lv.coefs",colnames(fit.mcmc))],nrow=p)
#  X.coefs = cw.X.coefs
#  row.coefs = cw.row.coefs
#  lv = matrix(fit.mcmc[t,grep("lvs", colnames(fit.mcmc))],nrow=n)
#  cutoffs = cw.cutoffs
#, X.multinom.coefs = NULL
calc.condlogLik <- function(y, X = NULL, family, trial.size = 1, lv.coefs, X.coefs = NULL, 
     row.coefs = NULL, row.ids = NULL, offset = NULL, lv = NULL, cutoffs = NULL, powerparam = NULL) {

     deprecate_warn("1.9", "boral::calc.marglogLik()", details = "We will be phasing out all functions to calculate log-likelihoods of any sort (too hard to maintain)!")

     if(length(family) != ncol(y) & length(family) != 1) { 
        stop("Number of elements in family is either 1 or equal to # of columns in y") }
    if(length(family) == 1) 
        complete_family <- rep(family,ncol(y))
    if(length(family) > 1) 
        complete_family <- family

    if(any(complete_family == "ordinal") & is.null(cutoffs)) 
        stop("Ordinal data requires cutoffs to be supplied.") 
    if(any(family == "tweedie") & (powerparam < 1 || powerparam > 2)) 
        stop("Common power parameter for tweedie must be between 1 and 2.") 
    #if(any(complete_family == "multinom") & is.null(X.multinom.coefs)) stop("Multinomial data requires X.multinom.coefs to be supplied") 

    if(!is.null(row.coefs)) {
        if(is.null(row.ids)) {
            row.ids <- matrix(1:nrow(y), ncol = 1)
            colnames(row.ids) <- "ID1"
            }
        row.ids <- check_row_ids(row.ids = row.ids, y = y)	
            if(!is.list(row.coefs))
                stop("row.coefs should be a list with length equal to the number of columns in row.coefs")
            if(length(row.coefs) != ncol(row.ids))
                stop("row.coefs should be a list with length equal to the number of columns in row.coefs")
        }
    check_offset(offset = offset, y = y)	
            

    ## Finish checks
    n <- nrow(y); p <- ncol(y); 
    num.lv <- 0
    if(!is.null(lv)) 
        num.lv <- ncol(lv)
    loglik <- 0
    loglik_comp <- matrix(NA, nrow = n, ncol = p) 

    if(length(trial.size) == 1) 
        complete_trial_size <- rep(trial.size,ncol(y))
    if(length(trial.size) > 1) 
        complete_trial_size <- trial.size

    all_etas <- matrix(lv.coefs[,1], nrow = n, ncol = p, byrow = TRUE)
    if(!is.null(lv)) 
        all_etas <- all_etas + tcrossprod(lv, as.matrix(lv.coefs[,2:(num.lv+1)]))
    if(!is.null(X.coefs)) 
        all_etas <- all_etas + tcrossprod(as.matrix(X), X.coefs)
    if(!is.null(offset)) 
        all_etas <- all_etas + offset
    ## Assumes the columns in X.coefs corresponding to multinomial are set to 0

    index_multinom_cols <- which(complete_family == "multinom")
    for(j in 1:p) {
        species_etas <- all_etas[,j]
        if(!is.null(row.coefs)) { for(k in 1:ncol(row.ids)) species_etas <- species_etas + row.coefs[[k]][row.ids[,k]] }
            
        if(complete_family[j] == "binomial") 
            loglik_comp[,j] <- dbinom(as.vector(unlist(y[,j])), size = complete_trial_size[j], prob = pnorm(species_etas), log = TRUE)
        if(complete_family[j] == "poisson") 
            loglik_comp[,j] <- dpois(as.vector(unlist(y[,j])), lambda = exp(species_etas), log = TRUE)
        if(complete_family[j] == "ztpoisson") 
            loglik_comp[,j] <- dztpois(as.vector(unlist(y[,j])), lambda = exp(species_etas), log = TRUE)
        if(complete_family[j] == "negative.binomial") 
            loglik_comp[,j] <- dnbinom(as.vector(unlist(y[,j])), size=1/(lv.coefs[j,ncol(lv.coefs)]+1e-5), mu=exp(species_etas), log = TRUE)
        if(complete_family[j] == "ztnegative.binomial") 
            loglik_comp[,j] <- dztnbinom(as.vector(unlist(y[,j])), size=1/(lv.coefs[j,ncol(lv.coefs)]+1e-5), mu=exp(species_etas), log = TRUE)
        if(complete_family[j] == "exponential") 
            loglik_comp[,j] <- dexp(as.vector(unlist(y[,j])), 1/exp(species_etas), log = TRUE)
        if(complete_family[j] == "gamma") 
            loglik_comp[,j] <- dgamma(as.vector(unlist(y[,j])), shape=exp(species_etas)*lv.coefs[j,ncol(lv.coefs)], rate = lv.coefs[j,ncol(lv.coefs)], log = TRUE)
        if(complete_family[j] == "beta") 
            loglik_comp[,j] <- dbeta(as.vector(unlist(y[,j])), lv.coefs[j,ncol(lv.coefs)]*exp(species_etas)/(1+exp(species_etas)), lv.coefs[j,ncol(lv.coefs)]*(1-exp(species_etas)/(1+exp(species_etas))),log = TRUE)
        if(complete_family[j] == "normal") 
            loglik_comp[,j] <- dnorm(as.vector(unlist(y[,j])), mean=species_etas, sd=(lv.coefs[j,ncol(lv.coefs)]+1e-6), log = TRUE)
        if(complete_family[j] == "lnormal") 
            loglik_comp[,j] <- dlnorm(as.vector(unlist(y[,j])), meanlog=species_etas, sdlog=(lv.coefs[j,ncol(lv.coefs)]+1e-6), log = TRUE)
            if(complete_family[j] == "tweedie") 
                loglik_comp[,j] <- dTweedie(as.vector(unlist(y[,j])), mu = exp(species_etas), phi = lv.coefs[j,ncol(lv.coefs)]+1e-6, p = powerparam, LOG = TRUE) 
        if(complete_family[j] == "ordinal") { 
            get_probs <- ordinal_conversion(n = n, lv = lv, lv.coefs.j = lv.coefs[j,], num.lv = num.lv, row.coefs = row.coefs, row.ids = row.ids, X = X, X.coefs.j = X.coefs[j,], cutoffs = cutoffs, est = "ignore"); 
            for(i in 1:n) 
                loglik_comp[i,j] <- log(get_probs[i,as.vector(y[i,j])]+1e-5) 
            }	
# 		if(complete_family[j] == "multinom") { 
# 			if(!is.null(X.multinom.coefs)) spp.etas <- matrix(rep(species_etas,dim(X.multinom.coefs)[3]),nrow=n) + as.matrix(X)%*%X.multinom.coefs[which(index_multinom_cols == j),,]
# 			get_probs <- exp(spp.etas)/apply(exp(spp.etas),1,sum)
# 			for(i in 1:n) { loglik_comp[i,j] <- log(get_probs[as.vector(y[i,j])]+1e-5) } }	
        } 

    return(list(logLik = sum(loglik_comp), logLik.comp = loglik_comp)) 
    }

	
## Calculate logl for models with no latent variables
## Conditional and marginal log-likelihood are the same in such case, except that in the case of row.eff = "random" there is marginalization over the random row effect
## lv.coefs still need to be provided though as it contains the spp effects and column-specific dispersion parameters
## Furthermore, note also that since it takes it X.coefs and lv.coefs:
## 1) it traits is not NULL, since marginalization is not done over the spp coefs; 
## 2) if > 1 are ordinal and hence the beta_{0j} are random effects, then marginalization is not done over the spp intercepts

calc.logLik.lv0 <- function (y, X = NULL, family, trial.size = 1, lv.coefs, 
     X.coefs = NULL, row.eff = "none", row.params = NULL, row.ids = NULL, offset = NULL, 
     cutoffs = NULL, powerparam = NULL) {

     deprecate_warn("1.9", "boral::calc.marglogLik()", details = "We will be phasing out all functions to calculate log-likelihoods of any sort (too hard to maintain)!")

     if(length(family) != ncol(y) & length(family) != 1) 
        stop("Number of elements in family is either 1 or equal to the number of columns in y.") 
    if(length(family) == 1) 
        complete_family <- rep(family, ncol(y))
    if(length(family) > 1) 
        complete_family <- family

    if(row.eff != "none") {
        if(is.null(row.ids)) {
            row.ids <- matrix(1:nrow(y), ncol = 1)
            colnames(row.ids) <- "ID1"
            }
        row.ids <- check_row_ids(row.ids = row.ids, y = y)	
        check_row_params(row.params = row.params, y = y, row.ids = row.ids)
        }
    if(row.eff == "fixed") 
        row.coefs <- row.params
    if(row.eff == "none") 
        row.coefs <- NULL

    check_offset(offset = offset, y = y)

    if(any(complete_family == "ordinal") & is.null(cutoffs)) 
        stop("Ordinal data requires cutoffs to be supplied.")
    if(any(family == "tweedie") & (powerparam < 1 || powerparam > 2)) 
        stop("Common power parameter for tweedie must be between 1 and 2.")
    if(length(trial.size) == 1) 
        complete_trial_size <- rep(trial.size, ncol(y))
    if(length(trial.size) > 1) 
        complete_trial_size <- trial.size

            
    ## Checks done
    n <- nrow(y); p <- ncol(y)
    logl <- 0
    loglik_comp <- matrix(0, nrow = n, ncol = p)
    index_multinom_cols <- which(complete_family == "multinom")

    if(row.eff != "random") {
        for(j in 1:p) {
            eta <- lv.coefs[j, 1]
            if(!is.null(X.coefs)) 
                eta <- eta + as.matrix(X) %*% X.coefs[j, ]
            if(!is.null(row.coefs)) { 
                for(k in 1:ncol(row.ids)) 
                    eta <- eta + row.coefs[[k]][row.ids[,k]] 
                }
            if(!is.null(offset)) 
                eta <- eta + offset[,j]

            if(complete_family[j] == "poisson")
                loglik_comp[, j] <- (dpois(as.vector(unlist(y[, j])), lambda = exp(eta), log = TRUE))
            if(complete_family[j] == "ztpoisson")
                loglik_comp[, j] <- (dztpois(as.vector(unlist(y[, j])), lambda = exp(eta), log = TRUE))
            if(complete_family[j] == "binomial") 
                loglik_comp[, j] <- (dbinom(as.vector(unlist(y[, j])), complete_trial_size[j], prob = pnorm(eta), log = TRUE))
            if(complete_family[j] == "negative.binomial") 
                loglik_comp[, j] <- (dnbinom(as.vector(unlist(y[, j])), mu = exp(eta), size = 1/(lv.coefs[j, 2]+1e-5), log = TRUE))
            if(complete_family[j] == "ztnegative.binomial") 
                loglik_comp[, j] <- (dztnbinom(as.vector(unlist(y[, j])), mu = exp(eta), size = 1/(lv.coefs[j, 2]+1e-5), log = TRUE))
            if(complete_family[j] == "exponential") 
                loglik_comp[, j] <- (dexp(as.vector(unlist(y[, j])), rate = 1/exp(eta), log = TRUE))
            if(complete_family[j] == "gamma") 
                loglik_comp[, j] <- (dgamma(as.vector(unlist(y[, j])), shape = exp(eta) * lv.coefs[j, 2], rate = lv.coefs[j, 2], log = TRUE))
            if(complete_family[j] == "beta") 
                loglik_comp[, j] <- (dbeta(as.vector(unlist(y[, j])), lv.coefs[j, 2] * exp(eta)/(1 + exp(eta)), lv.coefs[j, 2] * (1 - exp(eta)/(1 + exp(eta))), log = TRUE))
            if(complete_family[j] == "normal") 
                loglik_comp[, j] <- (dnorm(as.vector(unlist(y[,j])), mean = eta, sd = (lv.coefs[j, 2]), log = TRUE))
            if(complete_family[j] == "lnormal") 
                loglik_comp[, j] <- (dlnorm(as.vector(unlist(y[,j])), meanlog = eta, sdlog = (lv.coefs[j,2]), log = TRUE))
            if(complete_family[j] == "tweedie") 
                loglik_comp[, j] <- (dTweedie(as.vector(unlist(y[,j])), mu = exp(eta), phi = lv.coefs[j, 2], p = powerparam, LOG = TRUE))
            if(complete_family[j] == "ordinal") {
                get_probs <- ordinal_conversion(n = n, lv = NULL, lv.coefs.j = lv.coefs[j, ], num.lv = 0, row.coefs = row.coefs, row.ids = row.ids, X = X, X.coefs.j = X.coefs[j, ], offset.j = offset[,j], cutoffs = cutoffs, est = "ignore")
                for(i in 1:n) 
                    loglik_comp[i, j] <- log(get_probs[i, as.vector(y[i,j])] + 1e-05)
                }
                            
                            
# 		if(complete_family[j] == "multinom") {
# 			eta <- matrix(rep(eta,dim(X.multinom.coefs)[3]),nrow=n) + as.matrix(X)%*%X.multinom.coefs[which(index_multinom_cols == j),,]
# 			get_probs <- exp(eta)/apply(exp(eta),1,sum)
# 			for(i in 1:n) { loglik_comp[i,j] <- log(get_probs[as.vector(y[i,j])]) } }
                    }

    return(list(logLik = sum(loglik_comp), logLik.comp = rowSums(loglik_comp)))
    }

            
    if(row.eff == "random") {
        loglik <- 0; loglik_comp <- numeric(n); 
        mc_row_eff <- 1000
        mc_row_coefs <- vector("list", ncol(row.ids))
        for(k in 1:ncol(row.ids)) 
        mc_row_coefs[[k]] <- matrix(rnorm(length(unique(row.ids[,k]))*mc_row_eff, mean = 0, sd = row.params[[k]]), nrow = mc_row_eff) 
    
        ordinal.conversion.special <- function(lv.coefs.j, row.coefs.i, X.i = NULL, X.coefs.j = NULL, offset.j = NULL, cutoffs) {
            mc_row_eff <- length(row.coefs.i)
            etas <- matrix(NA, mc_row_eff, length(cutoffs))
            for(k in 1:length(cutoffs)) {
                etas[, k] <- cutoffs[k] - row.coefs.i - lv.coefs.j[1]
                if(!is.null(X.coefs.j)) 
                    etas[, k] <- etas[, k] - t(as.matrix(X.i)) %*% X.coefs.j 
                if(!is.null(offset.j)) 
                    etas[, k] <- etas[, k] - matrix(offset.j, ncol=1)
                }

            probs <- matrix(NA, mc_row_eff, length(cutoffs) + 1)
            probs[, 1] <- pnorm(etas[,1])
            for(k in 2:ncol(etas)) 
                probs[, k] <- pnorm(etas[, k]) - pnorm(etas[,k-1])
            probs[, length(cutoffs) + 1] <- 1 - pnorm(etas[,length(cutoffs)])
            rm(etas)
            return(probs)
            }

        for(i in 1:n) {
            spp_f <- eta <- matrix(lv.coefs[, 1], nrow = mc_row_eff, ncol = p, byrow = TRUE)
          for(k in 1:ncol(row.ids)) 
               eta <- eta + mc_row_coefs[[k]][,row.ids[i,k]]

               if(!is.null(X.coefs)) 
                    eta <- eta + matrix(t(as.matrix(X[i, ])) %*% t(X.coefs), nrow = mc_row_eff, ncol = p, byrow = TRUE)
               if(!is.null(offset)) 
                    eta <- eta + matrix(offset[i,], nrow = mc_row_eff, ncol = p, byrow = TRUE)

                    for(j in 1:p) {
                         if(complete_family[j] == "binomial") 
                         spp_f[, j] <- dbinom(rep(as.vector(y[i, j]), mc_row_eff), size = complete_trial_size[j], prob = pnorm(eta[, j]))
                         if(complete_family[j] == "poisson") 
                         spp_f[, j] <- dpois(rep(as.vector(y[i, j]), mc_row_eff), lambda = exp(eta[, j]))
                         if(complete_family[j] == "ztpoisson") 
                         spp_f[, j] <- dztpois(rep(as.vector(y[i, j]), mc_row_eff), lambda = exp(eta[, j]))
                         if(complete_family[j] == "negative.binomial") 
                         spp_f[, j] <- dnbinom(rep(as.vector(y[i, j]), mc_row_eff), size =1/(lv.coefs[j, ncol(lv.coefs)]+1e-5), mu = exp(eta[, j]))
                         if(complete_family[j] == "ztnegative.binomial") 
                         spp_f[, j] <- dztnbinom(rep(as.vector(y[i, j]), mc_row_eff), size =1/(lv.coefs[j, ncol(lv.coefs)]+1e-5), mu = exp(eta[, j]))
                         if(complete_family[j] == "exponential") 
                         spp_f[, j] <- dexp(rep(as.vector(y[i, j]), mc_row_eff), rate = 1/exp(eta[, j]))
                         if(complete_family[j] == "gamma") 
                         spp_f[, j] <- dgamma(rep(as.vector(y[i, j]), mc_row_eff), shape = exp(eta[, j]) * lv.coefs[j, ncol(lv.coefs)], rate = lv.coefs[j, ncol(lv.coefs)])
                         if(complete_family[j] == "beta") 
                         spp_f[, j] <- dbeta(rep(as.vector(y[i, j]), mc_row_eff), lv.coefs[j, ncol(lv.coefs)] * exp(eta[, j])/(1 + exp(eta[, j])), lv.coefs[j, ncol(lv.coefs)] * (1 - exp(eta[, j])/(1 + exp(eta[, j]))))
                         if(complete_family[j] == "normal") 
                         spp_f[, j] <- dnorm(rep(as.vector(y[i, j]), mc_row_eff), mean = eta[, j], sd = (lv.coefs[j, ncol(lv.coefs)]))
                         if(complete_family[j] == "lnormal") 
                         spp_f[, j] <- dlnorm(rep(as.vector(y[i, j]), mc_row_eff), meanlog = eta[, j], sdlog = (lv.coefs[j,ncol(lv.coefs)]))
                         if(complete_family[j] == "tweedie") { 
                         spp_f[, j] <- dTweedie(rep(as.vector(y[i, j]), mc_row_eff), mu = exp(eta[, j]), phi = lv.coefs[j, ncol(lv.coefs)] + 1e-06, p = powerparam, LOG = FALSE)
                         spp_f[, j][which(spp_f[, j] == 0)] <- 1 
                         }
                         if(complete_family[j] == "ordinal") {
                         get_probs <- ordinal.conversion.special(lv.coefs.j = lv.coefs[j, ], row.coefs.i = rnorm(mc_row_eff, 0, sd = row.params[[1]]), X.i = X[i, ], X.coefs.j = X.coefs[j, ], cutoffs = cutoffs)
                         spp_f[, j] <- get_probs[, as.vector(y[i, j])] + 1e-05 
                         }
                         }
                         
                    spp_f[!is.finite(spp_f)] = 1
                    spp_f <- matrix(spp_f, nrow = mc_row_eff, byrow = FALSE)
                    Q <- mean(apply(spp_f, 1, prod))
                    loglik <- loglik + log(Q)
                    loglik_comp[i] <- log(Q)
               }

          return(list(logLik = loglik, logLik.comp = loglik_comp)) 
          }
    }

	
## Calculate marginal logl
## loglik = sum_{i=1}^n \log( \int \prod_{j=1}^s f(y_ij|z_i) f(z_i) dz_i ) where z_i is assumed to be independent
## Furthermore, note also that since it takes it X.coefs and lv.coefs, and because is assumed the z_i are independent:
## 1) it traits is not NULL, since marginalization is not done over the spp coefs; 
## 2) if > 1 are ordinal and hence the beta_{0j} are random effects, then marginalization is not done over the spp intercepts
## 3) The marginalization is incorrect if z_i is structured
# family <- spider.fit.nb$family; trial.size = 1; lv.coefs = spider.fit.nb$lv.coefs.mean; X.coefs = spider.fit.nb$X.coefs.mean; row.eff = "random"; row.params = list(ID1 = spider.fit.nb$row.sigma$ID1$mean, ID2 = spider.fit.nb$row.sigma$ID2$mean); row.ids = spider.fit.nb$row.ids; num.lv <- spider.fit.nb$num.lv; lv.mc = NULL; cutoffs = NULL; powerparam = NULL

calc.marglogLik <- function (y, X = NULL, family, trial.size = 1, lv.coefs, 
     X.coefs = NULL, row.eff = "none", row.params = NULL, row.ids = NULL, offset = NULL, num.lv, 
     lv.mc = NULL, cutoffs = NULL, powerparam = NULL) { 
	
     deprecate_warn("1.6", "boral::calc.marglogLik()")
     
    if(num.lv == 0) 
        stop("Please use calc.loglik.lv0 to calculate likelihood in boral models with no latent variables.")
    if(is.null(lv.coefs)) 
        stop("lv.coefs must be given. Please use calc.loglik.lv0 to calculate likelihood in boral models with no latent variables.")

    if(is.null(lv.mc))
        lv.mc <- rmvnorm(1000, mean = rep(0, num.lv))
    numrows_mclv <- 1000
    if(!is.null(lv.mc))
        numrows_mclv <- nrow(lv.mc)

    if(length(family) != ncol(y) & length(family) != 1) { 
        stop("Number of elements in family is either 1 or equal to # of columns in y.") }
    if(length(family) == 1) 
        complete_family <- rep(family, ncol(y))
    if(length(family) > 1) 
        complete_family <- family
    if(any(family == "binomial") & !(length(trial.size) %in% c(1, length(family)))) {
        if(length(trial.size) == 1) 
            complete_trial_size <- rep(trial.size, ncol(y))
        }
    if(length(trial.size) > 1) 
        complete_trial_size <- trial.size

    if(row.eff != "none") {
        if(is.null(row.ids)) {
            row.ids <- matrix(1:nrow(y), ncol = 1)
            colnames(row.ids) <- "ID1"
            }
        row.ids <- check_row_ids(row.ids = row.ids, y = y)	
        check_row_params(row.params = row.params, y = y, row.ids = row.ids)
        }
    if(row.eff == "fixed") 
        row.coefs <- row.params

    check_offset(offset = offset, y = y)
    
    if(any(complete_family == "ordinal") & is.null(cutoffs)) 
        stop("Ordinal data requires cutoffs to be supplied.")
    if(any(family == "tweedie") & (powerparam < 1 || powerparam > 2)) 
        stop("Common power parameter for tweedie must be between 1 and 2.")
    #if(any(complete_family == "multinom") & is.null(X.multinom.coefs)) stop("Multinomial data requires X.multinom.coefs to be supplied") 

            
    ## Checks done
    n <- nrow(y); p <- ncol(y)
    loglik <- 0
    loglik_comp <- numeric(n)

    ## Internal function - Given the coefficients and cutoffs, return the multinomial probabilities for proportional odds regression for element of y
    ordinal.conversion.special <- function(lv.mc, lv.coefs.j, num.lv, row.coefs.i = NULL, X.i = NULL, X.coefs.j = NULL, offset.j = NULL, cutoffs) {
        etas <- matrix(NA, numrows_mclv, length(cutoffs))
        for(k in 1:length(cutoffs)) {
            etas[, k] <- cbind(1,lv.mc) %*% c(cutoffs[k], -lv.coefs.j[2:(num.lv + 1)]) - row.coefs.i - lv.coefs.j[1]
            if(!is.null(X.coefs.j)) 
                etas[, k] <- etas[, k] - t(as.matrix(X.i)) %*% X.coefs.j 
            if(!is.null(offset.j)) 
                etas[, k] <- etas[, k] - matrix(offset.j, ncol = 1) 
            }
            
        probs <- matrix(NA, numrows_mclv, length(cutoffs) + 1)
        probs[, 1] <- pnorm(etas[,1])
        for(k in 2:ncol(etas)) 
            probs[, k] <- pnorm(etas[,k]) - pnorm(etas[,k-1])
        probs[, length(cutoffs) + 1] <- 1 - pnorm(etas[,length(cutoffs)])
        rm(etas)
        probs
        }
            
    index_multinom_cols <- which(complete_family == "multinom")
    if(row.eff == "random") {
        mc_row_coefs <- vector("list", ncol(row.ids))
        for(k in 1:ncol(row.ids)) 
            mc_row_coefs[[k]] <- matrix(rnorm(length(unique(row.ids[,k]))*numrows_mclv, mean = 0, sd = row.params[[k]]), nrow = numrows_mclv) 
        }

        
    for(i in 1:n) {
        spp_f <- matrix(NA, nrow = numrows_mclv, ncol = p)
        spp.att.eta <- tcrossprod(cbind(1,lv.mc), lv.coefs[, 1:(num.lv + 1)])
        if(row.eff == "fixed") { 
            row.coefs.i <- 0; 
            for(k in 1:ncol(row.ids)) 
                row.coefs.i <- row.coefs.i + row.params[[k]][row.ids[i,k]] 
            spp.att.eta <- spp.att.eta + row.coefs.i
            }
        if(row.eff == "random") {
            for(k in 1:ncol(row.ids)) 
                spp.att.eta <- spp.att.eta + mc_row_coefs[[k]][,row.ids[i,k]]
            }
        if(!is.null(X.coefs)) 
            spp.att.eta <- spp.att.eta + matrix(t(as.matrix(X[i, ])) %*% t(X.coefs), nrow = numrows_mclv, ncol = p, byrow = TRUE)
        if(!is.null(offset)) 
            spp.att.eta <- spp.att.eta + matrix(offset[i,], nrow = numrows_mclv, ncol = p, byrow = TRUE)
            
        for(j in 1:p) {
            if(complete_family[j] == "binomial") 
                spp_f[, j] <- dbinom(rep(as.vector(y[i, j]), numrows_mclv), size = complete_trial_size[j], prob = pnorm(spp.att.eta[, j])) 
            if(complete_family[j] == "poisson") 	
                spp_f[, j] <- dpois(rep(as.vector(y[i, j]), numrows_mclv), lambda = exp(spp.att.eta[, j])) 
            if(complete_family[j] == "poisson") 
                spp_f[, j] <- dztpois(rep(as.vector(y[i, j]), numrows_mclv), lambda = exp(spp.att.eta[, j])) 
            if(complete_family[j] == "negative.binomial") 
                spp_f[, j] <- dnbinom(rep(as.vector(y[i, j]), numrows_mclv), size = 1/(lv.coefs[j, ncol(lv.coefs)]+1e-5), mu = exp(spp.att.eta[, j])) 
            if(complete_family[j] == "ztnegative.binomial") 
                spp_f[, j] <- dztnbinom(rep(as.vector(y[i, j]), numrows_mclv), size = 1/(lv.coefs[j, ncol(lv.coefs)]+1e-5), mu = exp(spp.att.eta[, j])) 
            if(complete_family[j] == "exponential") 
                spp_f[, j] <- dexp(rep(as.vector(y[i, j]), numrows_mclv), rate = 1/exp(spp.att.eta[, j])) 
            if(complete_family[j] == "gamma")  
                spp_f[, j] <- dgamma(rep(as.vector(y[i, j]), numrows_mclv), shape = exp(spp.att.eta[, j]) * lv.coefs[j, ncol(lv.coefs)], rate = lv.coefs[j, ncol(lv.coefs)]) 
            if(complete_family[j] == "beta") 
                spp_f[, j] <- dbeta(rep(as.vector(y[i, j]), numrows_mclv), lv.coefs[j, ncol(lv.coefs)] * exp(spp.att.eta[, j])/(1 + exp(spp.att.eta[, j])), lv.coefs[j, ncol(lv.coefs)] * (1 - exp(spp.att.eta[, j])/(1 + exp(spp.att.eta[, j])))) 
            if(complete_family[j] == "normal") 
                spp_f[, j] <- dnorm(rep(as.vector(y[i, j]), numrows_mclv), mean = spp.att.eta[, j], sd = (lv.coefs[j, ncol(lv.coefs)])) 
            if(complete_family[j] == "lnormal") 
                spp_f[, j] <- dlnorm(rep(as.vector(y[i, j]), numrows_mclv), meanlog = spp.att.eta[, j], sdlog = (lv.coefs[j, ncol(lv.coefs)])) 
            if(complete_family[j] == "tweedie") {
                spp_f[, j] <- dTweedie(rep(as.vector(y[i, j]), numrows_mclv), mu = exp(spp.att.eta[, j]), phi = lv.coefs[j, ncol(lv.coefs)] + 1e-06, p = powerparam, LOG = FALSE)
                spp_f[, j][which(spp_f[, j] == 0)] <- 1 
                }
            if(complete_family[j] == "ordinal") {
                get_probs <- ordinal.conversion.special(lv.mc, lv.coefs.j = lv.coefs[j, ], num.lv, row.coefs.i = row.coefs.i, X.i = X[i, ], X.coefs.j = X.coefs[j, ], cutoffs = cutoffs)
                spp_f[, j] <- get_probs[, as.vector(y[i, j])] + 1e-05 
                }
# 			if(complete_family[j] == "multinom") {
# 				num.multinom.levels <- dim(X.multinom.coefs)[3]
# 				if(!is.null(X.multinom.coefs)) { 
# 					spp.att.eta2 <- spp.att.eta[,j] + matrix(t(as.matrix(X[i,]))%*%X.multinom.coefs[which(index_multinom_cols == j),,],numrows_mclv,num.multinom.levels,byrow=TRUE) }
# 				get_probs <- exp(spp.att.eta2)/apply(exp(spp.att.eta2),1,sum)
# 				spp_f[,j] <- get_probs[,as.vector(y[i,j])]+1e-5 }
            }

        spp_f[!is.finite(spp_f)] = 1
        spp_f <- matrix(spp_f, nrow = numrows_mclv, byrow = FALSE)
        Q <- mean(apply(spp_f, 1, prod))	
        loglik <- loglik + log(Q)
        loglik_comp[i] <- log(Q)
        }

    return(list(logLik = loglik, logLik.comp = loglik_comp))
    }
	
	

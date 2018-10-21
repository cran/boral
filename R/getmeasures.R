## Calculates conditional WAIC, EAIC, EBIC. 
## Also calculate the marginal likelihood at component medians, and bases a AIC and BIC on this. Note this only done in cases where calc.marglogLik and calc.logLik.lv0 actually produce a sensible result
get.measures <- function(y, X = NULL, family, trial.size = 1, row.eff = "none", row.ids = NULL, 
     offset = NULL, num.lv, fit.mcmc) {
        
    warning("Please note that as of version 1.6, functions to calculate information criteria will no longer be updated. Use at your peril!")
    
    
    if(length(family) != ncol(y) & length(family) != 1) 
            stop("Number of elements in family is either 1 or equal to # of columns in y")
    if(length(family) == 1) 
            complete_family <- rep(family, ncol(y))
    if(length(family) > 1) 
            complete_family <- family
    index_ordinal_cols <- which(family == "ordinal")
    
    if(any(family == "binomial") & !(length(trial.size) %in% c(1, length(family)))) 
    stop("trial.size needs to be specified if any columns are binomially distributed; can either be a single element or a vector equal to the # of columns in y. The latter will assume the specified trial size for all rows labelled binomial in the family argument")
    
# 	lv.control <- check.lv.control(num.lv, lv.control)
# 	num.lv <- lv.control$num.lv
    
    if(row.eff != "none") {
            if(is.null(row.ids)) {
                    row.ids <- matrix(1:nrow(y), ncol = 1)
                    colnames(row.ids) <- "ID1"
                    }
        row.ids <- check.row.ids(row.ids = row.ids, y = y)	
        }
    check.offset(offset = offset, y = y) 

    if(length(grep("ssvs",colnames(fit.mcmc))) > 0)
        warnings("Calculation of information criterion in the presence of any SSVS is problematic at best!")

    do.marglik.ics <- check.domarglik.ics(fit.mcmc.names = colnames(fit.mcmc), index.ordinal.cols = index_ordinal_cols)
        
    ##################
    ## Checks done 
    ##################

    n <- nrow(y); p <- ncol(y)
    all.lppd <- matrix(NA, nrow = nrow(fit.mcmc), ncol = n*p)
    index_multinom_cols <- which(complete_family == "multinom")

    for(k0 in 1:nrow(fit.mcmc)) {
            cw.params <- list(lv.coefs = matrix(fit.mcmc[k0, grep("lv.coefs", colnames(fit.mcmc))], nrow = p))
            if(row.eff != "none") {
                    cw.params$row.coefs <- vector("list", ncol(row.ids))
                    for(k in 1:ncol(row.ids)) 
                            cw.params$row.coefs[[k]] <- fit.mcmc[k0, grep(paste0("row.coefs.ID",k,"\\["), colnames(fit.mcmc))] 
                    }
            if(!is.null(X)) 
                    cw.params$X.coefs <- matrix(fit.mcmc[k0, grep("X.coefs", colnames(fit.mcmc))], nrow = p) 
            if(any(complete_family == "ordinal")) 
                    cw.params$cutoffs <- fit.mcmc[k0, grep("cutoffs", colnames(fit.mcmc))] 
            if(any(complete_family == "tweedie")) 
                    cw.params$powerparam <- fit.mcmc[k0, grep("powerparam", colnames(fit.mcmc))] 
# 		if(any(complete_family == "multinom") & !is.null(X)) { 
# 			get.X.multinom.coefs <- array(matrix(fit.mcmc[k0, grep("X.multinom.params", colnames(fit.mcmc))],dim = c(length(index_multinom_cols),ncol(X),ncol(get.X.multinom.coefs)/ncol(X))), dimnames=NULL) } else { get.X.multinom.coefs <- NULL }
            
            if(num.lv > 0) {
                    get_out <- calc.condlogLik(y, X = X, family = complete_family, trial.size, lv.coefs = cw.params$lv.coefs, X.coefs = cw.params$X.coefs, row.coefs = cw.params$row.coefs, row.ids = row.ids, offset = offset, lv = matrix(fit.mcmc[k0, grep("lvs", colnames(fit.mcmc))], nrow = n), cutoffs = cw.params$cutoffs, powerparam = cw.params$powerparam)
                    }
            if(num.lv == 0) {
                    get_out <- calc.condlogLik(y, X = X, family = complete_family, trial.size, lv.coefs = cw.params$lv.coefs, X.coefs = cw.params$X.coefs, row.coefs = cw.params$row.coefs, row.ids = row.ids, offset = offset, lv = NULL, cutoffs = cw.params$cutoffs, powerparam = cw.params$powerparam)
                    }
                    
            get_out$logLik.comp[!is.finite(get_out$logLik.comp)] <- NA
            all.lppd[k0,] <- as.vector(get_out$logLik.comp)
            }
            
            
    all.cond.logl <- rowSums(all.lppd)
    waic.out <- -2 * sum(log(colMeans(exp(all.lppd), na.rm = TRUE))) + 2 * sum(apply(all.lppd, 2, var, na.rm = TRUE))
    cond.num.params <- sum(cw.params$lv.coefs != 0) + n*num.lv + ## lv and loadings
            sum(cw.params$X.coefs != 0)*as.numeric(!is.null(X)) + ## X.coefs
            any(complete_family == "ordinal")*sum(cw.params$cutoffs != 0) + (1 - is.null(cw.params$powerparam)) ## other parameters
    if(row.eff != "none") { 
            cond.num.params <- cond.num.params + sum(sapply(cw.params$row.coefs,length)) 
            } ## row effects

    eaic <- -2*mean(all.cond.logl, na.rm = TRUE) + 2*cond.num.params
    ebic <- -2*mean(all.cond.logl, na.rm = TRUE) + log(n*p)*cond.num.params
    out.list <- list(waic = waic.out, eaic = eaic, ebic = ebic, all.cond.logLik = all.cond.logl, cond.num.params = cond.num.params, do.marglik.ics = do.marglik.ics)

    if(do.marglik.ics) {
            ## Calculate marginal logL at component medians
            params.median <- list(lv.coefs = matrix(apply(fit.mcmc[, grep("lv.coefs", colnames(fit.mcmc))], 2, median), nrow = p))
# 		if(type != "independent") {
# 			params.median$lv.covparams <- apply(as.matrix(fit.mcmc[, grep("lv.covparams", colnames(fit.mcmc))]), 2, median)
# 			}
            if(row.eff == "fixed") {
                    params.median$row.coefs <- vector("list", ncol(row.ids))
                    for(k in 1:ncol(row.ids)) 
                params.median$row.coefs[[k]] <- apply(fit.mcmc[, grep(paste0("row.coefs.ID",k,"\\["), colnames(fit.mcmc))], 2, median)
                    }
            if(row.eff == "random") {
                    params.median$row.coefs <- vector("list", ncol(row.ids))
                    for(k in 1:ncol(row.ids)) 
                params.median$row.coefs[[k]] <- median(fit.mcmc[, grep(paste0("row.sigma.ID",k,"$"), colnames(fit.mcmc))])
                    }
            if(!is.null(X)) 
                    params.median$X.coefs <- matrix(apply(fit.mcmc[, grep("X.coefs", colnames(fit.mcmc))], 2, median), nrow = p) 
            if(any(complete_family == "ordinal")) 
                    params.median$cutoffs <- apply(fit.mcmc[, grep("cutoffs", colnames(fit.mcmc))], 2, median)
            if(any(complete_family == "tweedie")) 
                    params.median$powerparam <- median(fit.mcmc[, grep("powerparam", colnames(fit.mcmc))]) 
    # 	if(any(complete_family == "multinom") & !is.null(X)) { 
    # 		get.X.multinom.coefs <- array(matrix(apply(fit.mcmc[,grep("X.multinom.params", colnames(fit.mcmc))],2,median),dim = c(length(index_multinom_cols),ncol(X),ncol(get.X.multinom.coefs)/ncol(X))), dimnames=NULL) } else { get.X.multinom.coefs <- NULL }

            if(num.lv > 0) {
                    median.marglogl <- calc.marglogLik(y, X = X, family = complete_family, trial.size, lv.coefs = params.median$lv.coefs, X.coefs = params.median$X.coefs, row.eff = row.eff, row.ids = row.ids, row.params = params.median$row.coefs, offset = offset, num.lv = num.lv, lv.mc = NULL, cutoffs = params.median$cutoffs, powerparam = params.median$powerparam) 
                    }
            if(num.lv == 0) 
                    median.marglogl <- calc.logLik.lv0(y, X = X, family = complete_family, trial.size, lv.coefs = params.median$lv.coefs, X.coefs = params.median$X.coefs, row.eff = row.eff, row.ids = row.ids, row.params = params.median$row.coefs, offset = offset, cutoffs = params.median$cutoffs, powerparam = params.median$cw.powerparam)	

            marg.num.params <- sum(params.median$lv.coefs != 0) + ## Loadings 
                    sum(params.median$X.coefs != 0)*as.numeric(!is.null(X)) + ## X coefs
                    any(complete_family == "ordinal")*sum(params.median$cutoffs != 0) + (1 - is.null(params.median$powerparam)) ## other parameters 
            if(row.eff != "none") { 
                    marg.num.params <- marg.num.params + sum(sapply(params.median$row.coefs,length)) 
                    } ## row effects

            marg.aic <- -2 * median.marglogl$logLik + 2 * marg.num.params
            marg.bic <- -2 * median.marglogl$logLik + log(n*p) * marg.num.params
            
            out.list$median.logLik <- median.marglogl$logLik
            out.list$marg.num.params <- marg.num.params
            out.list$aic.median <- marg.aic; out.list$bic.median <- marg.bic
            }
            
            
    return(out.list)
    }
	

## Calculates marginal logl for all samples to produce a proper AIC and BIC. 
## Calculates WAIC based on the marginal likelihood; DIC based on the marginal likelihood
## All of this is only permitted when the fitted boral model has a simple enough structure to allow these calculations!
get.more.measures <- function(y, X = NULL, family, trial.size = 1, row.eff = "none", row.ids = NULL, offset = NULL, num.lv, fit.mcmc, verbose = TRUE) { #lv.control

    warning("Please note that as of version 1.6, functions to calculate information criteria will no longer be updated. Use at your peril!")

    index_ordinal_cols <- which(family == "ordinal")

    if(num.lv == 0) 
            stop("For boral models with no latent variables, the marginal and conditional likelihoods are equivalent, and there is nothing to gain from using get.more.measures")
    if(length(family) != ncol(y) & length(family) != 1) 
            stop("Number of elements in family is either 1 or equal to # of columns in y")
    if(length(family) == 1) 
            complete_family <- rep(family, ncol(y))
    if(length(family) > 1) 
            complete_family <- family
    
    if(any(family == "binomial") & !(length(trial.size) %in% c(1, length(family)))) 
            stop("trial.size needs to be specified if any columns are binomially distributed; can either be a single element or a vector equal to the # of columns in y. The latter will assume the specified trial size for all rows labelled binomial in the family argument")

# 	lv.control <- check.lv.control(num.lv, lv.control)
# 	num.lv <- lv.control$num.lv

if(row.eff != "none") {
            if(is.null(row.ids)) {
                    row.ids <- matrix(1:nrow(y), ncol = 1)
                    colnames(row.ids) <- "ID1"
                    }
    row.ids <- check.row.ids(row.ids = row.ids, y = y)	
    }
    check.offset(offset = offset, y = y) 

if(length(grep("ssvs",colnames(fit.mcmc))) > 0)
    warnings("Calculation of information criterion in the presence of any SSVS is problematic at best!")

do.marglik.ics <- check.domarglik.ics(fit.mcmc.names = colnames(fit.mcmc), index.ordinal.cols = index_ordinal_cols)
    
            
    ## Checks done ##
    
    
    if(!do.marglik.ics) {
            message("The current version of boral does not implement information criterion based on the marginal likelihood for the specified model, because the number of random effects included in the model and/or the random effects structure is too complicated...sorry!")
            return()
            }
    
    n <- nrow(y); p <- ncol(y)
    big_lv <- rmvnorm(2000, rep(0, num.lv))
    all.marg.logl <- matrix(NA, nrow(fit.mcmc), n)
    index_multinom_cols <- which(complete_family == "multinom")
    
    ## Calculate marginal likelihood at all iterations
    for(k0 in 1:nrow(fit.mcmc)) {
            if(verbose == TRUE & k0%%100 == 0) 
                    message("Onto mcmc sample ", k0)
            cw.params <- list(lv.coefs = matrix(fit.mcmc[k0, grep("lv.coefs", colnames(fit.mcmc))], nrow = p))
#           if(lv.control$type != "independent") {
#                cw.params$lv.covparams <- fit.mcmc[k0, grep("lv.covparams", colnames(fit.mcmc))]
#                }
            if(row.eff == "fixed") {
                    cw.params$row.coefs <- vector("list", ncol(row.ids))
                    for(k in 1:ncol(row.ids)) 
                            cw.params$row.coefs[[k]] <- fit.mcmc[k0, grep(paste0("row.coefs.ID",k,"\\["), colnames(fit.mcmc))]
                    }
            if(row.eff == "random") {
                    cw.params$row.coefs <- vector("list", ncol(row.ids))
                    for(k in 1:ncol(row.ids)) 
                            cw.params$row.coefs[[k]] <- fit.mcmc[k0, grep(paste0("row.sigma.ID",k,"$"), colnames(fit.mcmc))]
                    }
            if(!is.null(X)) 
                    cw.params$X.coefs <- matrix(fit.mcmc[k0, grep("X.coefs", colnames(fit.mcmc))], nrow = p) 
            if(any(complete_family == "ordinal")) 
                    cw.params$cutoffs <- fit.mcmc[k0, grep("cutoffs", colnames(fit.mcmc))] 
            if(any(complete_family == "tweedie")) 
                    cw.params$powerparam <- fit.mcmc[k0, grep("powerparam", colnames(fit.mcmc))] 
# 		if(any(complete_family == "multinom") & !is.null(X)) { 
# 		get.X.multinom.coefs <- array(matrix(fit.mcmc[k0,grep("X.multinom.params", colnames(fit.mcmc))],dim = c(length(index_multinom_cols),ncol(X),ncol(get.X.multinom.coefs)/ncol(X))), dimnames = NULL) } else { get.X.multinom.coefs <- NULL }

            get.mll <- calc.marglogLik(y, X = X, family = complete_family, trial.size, lv.coefs = cw.params$lv.coefs, X.coefs = cw.params$X.coefs, row.eff = row.eff, row.ids = row.ids, row.params = cw.params$row.coefs, offset = offset, num.lv = num.lv, lv.mc = big_lv, cutoffs = cw.params$cutoffs, powerparam = cw.params$powerparam) 
            all.marg.logl[k0,] <- get.mll$logLik.comp
            }
                            
    ## Calculate WAIC based on marginal
    marg.waic <- -2 * sum(log(apply(exp(all.marg.logl), 2, mean, na.rm = TRUE))) + 2 * sum(apply(all.marg.logl, 2, var, na.rm = TRUE))

    ## Calculate AIC, BIC at posterior mode	
    marg.num.params <- sum(cw.params$lv.coefs != 0) + ## Loadings 
            sum(cw.params$X.coefs != 0)*as.numeric(!is.null(X)) +  ## X coefs
            any(complete_family == "ordinal")*sum(cw.params$cutoffs != 0) + (1 - is.null(cw.params$powerparam)) ## other parameters	
    if(row.eff != "none") { ## row effects
    marg.num.params <- marg.num.params + sum(sapply(cw.params$row.coefs,length)) 
    } 
    bic1 <- -2 * max(rowSums(all.marg.logl)) + log(n)*marg.num.params
    aic1 <- -2 * max(rowSums(all.marg.logl)) + 2*marg.num.params
    
    ## Calculate DIC based on marginal
    params.mean <- list(lv.coefs = matrix(apply(fit.mcmc[, grep("lv.coefs", colnames(fit.mcmc))], 2, mean), nrow = p))
# 	if(lv.control$type != "independent") {
# 		params.mean$lv.covparams <- apply(as.matrix(fit.mcmc[, grep("lv.covparams", colnames(fit.mcmc))],2,mean))
# 		}
    if(row.eff == "fixed") {
            params.mean$row.coefs <- vector("list", ncol(row.ids))
            for(k in 1:ncol(row.ids)) 
                    params.mean$row.coefs[[k]] <- apply(fit.mcmc[, grep(paste0("row.coefs.ID",k,"\\["), colnames(fit.mcmc))],2,mean)
            }
    if(row.eff == "random") {
            params.mean$row.coefs <- vector("list", ncol(row.ids))
            for(k in 1:ncol(row.ids)) 
                    params.mean$row.coefs[[k]] <- apply(fit.mcmc[, grep(paste0("row.sigma.ID",k,"$"), colnames(fit.mcmc))],2,mean)
            }
    if(!is.null(X)) 
            params.mean$X.coefs <- matrix(apply(fit.mcmc[, grep("X.coefs", colnames(fit.mcmc))], 2, mean), nrow = p) 
    if(any(complete_family == "ordinal")) 
            params.mean$cutoffs <- apply(fit.mcmc[, grep("cutoffs", colnames(fit.mcmc))], 2, mean)
    if(any(complete_family == "tweedie")) 
            params.mean$powerparam <- mean(fit.mcmc[, grep("powerparam", colnames(fit.mcmc))]) 
    
    marg.dic <- -2*calc.marglogLik(y, X = X, family = complete_family, trial.size, lv.coefs = params.mean$lv.coefs, X.coefs = params.mean$X.coefs, row.eff = row.eff, row.ids = row.ids, row.params = params.mean$row.coefs, offset = offset, num.lv = num.lv, lv.mc = big_lv, cutoffs = params.mean$cutoffs, powerparam = params.mean$powerparam)$logLik 
    marg.dic <- marg.dic + 2*(2*var(rowSums(all.marg.logl), na.rm = TRUE))

    return(list(aic.mode = aic1, bic.mode = bic1, marg.dic = marg.dic, marg.waic = marg.waic, all.marg.logLik = rowSums(all.marg.logl), marg.num.params = marg.num.params))
    }
	


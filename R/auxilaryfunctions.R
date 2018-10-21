##########
## Auxilary functions
##########
.onAttach <- function(...) {
    packageStartupMessage("If you recently updated boral, please check news(package = \"boral\") for the updates in the latest version", appendLF=TRUE)
}



## Variance partitioning on the linear predictor scale; the variance is calculated and and averaged wrt to the posterior distribution 
calc.varpart <- function(object, groupX = NULL) {
    fit.mcmc <- get.mcmcsamples(object)
    num.lv <- object$num.lv
    if(is.null(object$X))
        stop("Variance partitioning is not particularly useful when there are no covariates X included in the model")
    if(is.null(object$jags.model)) 
        stop("MCMC samples not found")
    if(!is.null(groupX)) { if(length(groupX) != (object$num.X+1)) 
        stop("If groupX is supplied, then it must be a vector with the length (object$num.X + 1)") }

    row_var <- lv_var <- X_var <- matrix(0, nrow(fit.mcmc), object$p)
    if(!is.null(groupX))
        groupX_var <- array(0, dim = c(length(unique(groupX)), nrow(fit.mcmc), object$p))          
    if(!is.null(object$traits)) {
        trait.linpred <- array(0, dim = c(nrow(fit.mcmc), object$n, object$p))
        all_cors_spp <- matrix(0, nrow(fit.mcmc), object$p)
        }

        
    for(k in 1:nrow(fit.mcmc)) {
        cw_X_coefs <- matrix(fit.mcmc[k, grep("X.coefs", colnames(fit.mcmc))], nrow = object$p)
        cw_lv_coefs <- matrix(fit.mcmc[k, grep("lv.coefs", colnames(fit.mcmc))], nrow = object$p) ## Need for spp intercept

        fullX <- cbind(1,object$X)
        full.coefs <- cbind(cw_lv_coefs[,1],as.matrix(cw_X_coefs))
        cw.X.linpred <- tcrossprod(fullX, full.coefs)
        if(!is.null(groupX)) {
            cw.groupX.linpred <- vector("list", length(unique(groupX)))
            for(k2 in 1:length(unique(groupX))) 
                cw.groupX.linpred[[k2]] <- tcrossprod(fullX[,which(groupX==k2)], full.coefs[,which(groupX==k2)])
            }
        if(object$row.eff == "random") {
            for(k2 in 1:ncol(object$row.ids)) 
                row_var[k,] <- row_var[k,] + rep(fit.mcmc[k, grep(paste0("row.sigma.ID",k2,"$"), colnames(fit.mcmc))]^2,object$p) 
            }
        if(object$row.eff == "fixed") {
            cw.row.linpred <- matrix(0, object$n, object$p)
            for(k2 in 1:ncol(object$row.ids)) 
                cw.row.linpred <- cw.row.linpred + matrix(fit.mcmc[k, grep(paste0("row.coefs.ID",k2,"\\["), colnames(fit.mcmc))][object$row.ids[,k2]], nrow = object$n, ncol = object$p, byrow = FALSE)
            row_var[k,] <- apply(cw.row.linpred, 2, var)
            }

        X_var[k,] <- apply(cw.X.linpred,2,var)
        if(!is.null(groupX)) {
            for(k2 in 1:length(unique(groupX))) 
                groupX_var[k2,k,] <- apply(cw.groupX.linpred[[k2]],2,var)
            }
        if(num.lv > 0) 
            lv_var[k,] <- rowSums(cw_lv_coefs[,2:(1+object$num.lv)]^2) 
        if(!is.null(object$traits)) {
            cw.traits.coefs <- cbind(fit.mcmc[k, grep("traits.int",colnames(fit.mcmc))], matrix(fit.mcmc[k, grep("traits.coefs",colnames(fit.mcmc))], nrow = ncol(object$X)+1))
            rownames(cw.traits.coefs) <- c("beta0", colnames(object$X))
            trait.X.coefs <- tcrossprod(cbind(1,object$traits), cw.traits.coefs) ## beta = intercept + trait %*% trait.coefs
            cw.trait.linpred <- tcrossprod(cbind(1,object$X), trait.X.coefs)
            all_cors_spp[k,] <-  sapply(1:object$p, function(i) cor(cw.X.linpred[,i], cw.trait.linpred[,i])^2)
            }
        }
    

    total.var <- X_var + row_var + lv_var
    var.X <- colMeans(X_var/total.var); names(var.X) <- colnames(object$y)
    var_lv <- NULL
    if(num.lv > 0) { 
        var_lv <- colMeans(lv_var/total.var)
        names(var_lv) <- colnames(object$y) 
        }
    var_row <- NULL
    if(object$row.eff != "none") { 
        var_row <- colMeans(row_var/total.var)
        names(var_row) <- colnames(object$y) 
        }

    ## As soon as groupX is supplied, change the way variance decomposition is done
    if(!is.null(groupX)) {
        total.var <- apply(groupX_var, c(2,3), sum) + row_var + lv_var ## Note this is not equal to total.var
        var.X <- matrix(0, length(unique(groupX)), object$p) 
        for(k2 in 1:length(unique(groupX))) 
            var.X[k2,] <- colMeans(groupX_var[k2,,]/total.var)
        rownames(var.X) <- unique(groupX)
        colnames(var.X) <- colnames(object$y)
        var_lv <- NULL
        if(num.lv > 0) { 
            var_lv <- colMeans(lv_var/total.var)
            names(var_lv) <- colnames(object$y) 
            }
        var_row <- NULL
        if(object$row.eff != "none") { 
            var_row <- colMeans(row_var/total.var)
            names(var_row) <- colnames(object$y) 
            }          
        }
    
    out <- list(varpart.X = var.X, varpart.lv = var_lv, varpart.row = var_row)
    if(!is.null(object$traits)) {
        out$R2.traits <- colMeans(all_cors_spp)
        names(all_cors_spp) <- colnames(object$y) 
        }
    return(out)
    }


## Dunn-Smyth residuals
## Also create a confusion matrix for ordinal and multinomial data
ds.residuals <- function(object, est = "median") {  
    n <- object$n; p <- object$p; 
    num.lv <- object$num.lv
    num.ord.levels <- object$num.ord.levels; 
    X <- object$X; y <- object$y
    mus <- fitted.boral(object, est = est)

    if(any(object$family == "ordinal")) {
        message("One or more columns of y have ordinal responses. Constructing a single confusion matrix for these")
        true_resp <- as.matrix(y[,which(object$family == "ordinal")])
        pred_resp <- matrix(NA,n,ncol(true_resp)) 
        }
# 	if(any(object$family == "multinom")) {
# 		print("One or more columns of y have multinomial responses. Constructing a single confusion matrix for these")
# 		true.multinom.resp <- as.matrix(y[,which(object$family == "multinom")])
# 		pred.multinom.resp <- matrix(NA,n,ncol(true.multinom.resp)) }
    if(any(object$family == "tweedie")) {
        if(est == "median") 
            powerparam <- object$powerparam.median
        if(est == "mean") 
            powerparam <- object$powerparam.mean 
        }

    dsres_out <- matrix(NA,n,p)
    rownames(dsres_out) <- rownames(y)
    colnames(dsres_out) <- colnames(y)
    for(i in 1:n) { for(j in 1:p) {
        if(object$family[j] == "poisson") { 
            a <- ppois(as.vector(unlist(y[i,j]))-1, mus$out[i,j]); 
            b <- ppois(as.vector(unlist(y[i,j])), mus$out[i,j]); 		
            u <- runif(n = 1, min = a, max = b); dsres_out[i,j] <- qnorm(u) 
            }
        if(object$family[j] == "negative.binomial") {
            if(est == "median") 
                phis <- object$lv.coefs.median[,num.lv+2]+1e-5
            if(est == "mean") 
                phis <- object$lv.coefs.mean[,num.lv+2]+1e-5
            a <- pnbinom(as.vector(unlist(y[i,j]))-1, mu=mus$out[i,j], size=1/phis[j]); 
            b <- pnbinom(as.vector(unlist(y[i,j])), mu=mus$out[i,j], size=1/phis[j])
            u <- runif(n = 1, min = a, max = b); dsres_out[i,j] <- qnorm(u) 
            }
        if(object$family[j] == "binomial") { 
            a <- pbinom(as.vector(unlist(y[i,j]))-1, object$trial.size[j], prob = mus$out[i,j]); 
            b <- pbinom(as.vector(unlist(y[i,j])), object$trial.size[j], prob = mus$out[i,j])
            u <- runif(n = 1, min = a, max = b); dsres_out[i,j] <- qnorm(u) 
            }
        if(object$family[j] == "exponential") { 
            a <- pexp(as.vector(unlist(y[i,j])), rate=1/mus$out[i,j]); 
            dsres_out[i,j] <- qnorm(a) 
            }
        if(object$family[j] == "gamma") { 
            if(est == "median") 
                phis <- object$lv.coefs.median[,num.lv+2]
            if(est == "mean") 
                phis <- object$lv.coefs.mean[,num.lv+2]
            a <- pgamma(as.vector(unlist(y[i,j])), shape=mus$out[i,j]*phis[j], rate=phis[j]); 
            dsres_out[i,j] <- qnorm(a) 
            }
        if(object$family[j] == "beta") { 
            if(est == "median") 
                phis <- object$lv.coefs.median[,num.lv+2]
            if(est == "mean") 
                phis <- object$lv.coefs.mean[,num.lv+2]
            a <- pbeta(as.vector(unlist(y[i,j])), shape1=phis[j]*mus$out[i,j], shape2=phis[j]*(1-mus$out[i,j]))
            dsres_out[i,j] <- qnorm(a) 
            }
        if(object$family[j] == "normal") { 
            if(est == "median") 
                phis <- object$lv.coefs.median[,num.lv+2]
            if(est == "mean") 
                phis <- object$lv.coefs.mean[,num.lv+2]
            a <- pnorm(as.vector(unlist(y[i,j])), mus$out[i,j], sd = (phis[j])); 
            dsres_out[i,j] <- qnorm(a) 
            }
# 			X2 <- cbind(1,X); hatmat <- X2%*%solve(crossprod(X2))%*%t(X2)
# 			dsres_out[i,j] <- (y[i,j]-mus$out[i,j])/(sqrt(phis[j])*sqrt(1-hatmat[i,i])) }
        if(object$family[j] == "lnormal") { 
            if(est == "median") 
                phis <- object$lv.coefs.median[,num.lv+2]
            if(est == "mean") 
                phis <- object$lv.coefs.mean[,num.lv+2]
            a <- plnorm(as.vector(unlist(y[i,j])), log(mus$out[i,j]), sdlog = (phis[j])); dsres_out[i,j] <- qnorm(a) 
            }
        if(object$family[j] == "tweedie") { 
            if(est == "median") 
                phis <- object$lv.coefs.median[,num.lv+2]
            if(est == "mean") 
                phis <- object$lv.coefs.mean[,num.lv+2]
            a <- pTweedie(as.vector(unlist(y[i,j])), mu = mus$out[i,j], phi = phis[j], p = powerparam); dsres_out[i,j] <- qnorm(a) 
            }
        if(object$family[j] == "ordinal") { 
            pred_resp[,which(object$family == "ordinal")==j] <- mus$out[,which(object$family == "ordinal")==j] ## get max predicted probability
            cumsum.b <- sum(mus$ordinal.probs[i,j,1:(y[i,j])])
            cumsum.a <- sum(mus$ordinal.probs[i,j,1:(y[i,j]-1)])
            u <- runif(n = 1, min = cumsum.a, max = cumsum.b); 
            if(abs(u-1) < 1e-5) 
                u <- 1
            if(abs(u-0) < 1e-5) 
                u <- 0
            dsres_out[i,j] <- qnorm(u) 
            }
# 		if(object$family[j] == "multinom") { ## get max predicted probability
# 			pred_resp[i,which(object$family == "multinom")==j] <- which.max(mus$multinom.probs[i,j,]) }
        } }

    if(sum(object$family == "ordinal") > 0) {
        agree_tab <- table(as.vector(pred_resp), as.vector(true_resp))
        }
    else { 
        agree_tab <- NULL 
        }
    #if(sum(object$family == "multinom") > 0) { agree.multinom.tab <- table(as.vector(pred.multinom.resp), as.vector(true.multinom.resp)); }	else { agree.multinom.tab <- NULL }
    
    return(list(agree.ordinal = agree_tab, residuals = dsres_out))
    }

	
## Fitted values
## For ordinal and multinomial data, returns a matrix of probabilities for each vector of rows
fitted.boral <- function(object, est = "median",...) {
    n <- object$n; p <- object$p
    num.lv <- object$num.lv; 
    X <- object$X
    y <- object$y
    fitted_out <- matrix(NA,n,p)
    rownames(fitted_out) <- rownames(y)
    colnames(fitted_out) <- colnames(y)

    if(any(object$family == "ordinal")) { 
        fitted_ordinal_probs <- array(NA, dim=c(n,p,object$num.ord.levels)) 
        dimnames(fitted_ordinal_probs) <- list(r = rownames(y), c = colnames(y), levels = 1:object$num.ord.levels) 
        } 
    else { 
        fitted_ordinal_probs <- NULL 
        }
# 	if(any(object$family == "multinom")) { 
# 		fitted_multinom_probs <- array(NA,dim=c(n,p,object$num.multinom.levels))
# 		dimnames(fitted_multinom_probs) <- list(r = rownames(y), c = colnames(y), levels = 1:object$num.multinom.levels) } 
# 	else { fitted_multinom_probs <- NULL }
    if(length(object$lv.coefs.median) == p) {
        object$lv.coefs.median <- matrix(object$lv.coefs.median, ncol = 1)
        object$lv.coefs.mean <- matrix(object$lv.coefs.mean, ncol = 1)
        }
    if(length(object$X.coefs.median) == p) {
        object$X.coefs.median <- matrix(object$X.coefs.median, ncol = 1)
        object$X.coefs.mean <- matrix(object$X.coefs.mean, ncol = 1)
        }

        
    if(is.null(object$lv.median)) 
        eta <- tcrossprod(matrix(1,n,1), matrix(object$lv.coefs.median[,1],p,1))
    if(!is.null(object$lv.median)) 
        eta <- tcrossprod(cbind(1,object$lv.median), object$lv.coefs.median[,1:(num.lv+1)])
    if(!is.null(object$X.coefs.median)) 
        eta <- eta + tcrossprod(as.matrix(X), object$X.coefs.median)
    if(!is.null(object$offset)) 
        eta <- eta + object$offset
    
    if(est == "mean") {
        if(is.null(object$lv.mean)) 
            eta <- tcrossprod(matrix(1,n,1), object$lv.coefs.mean[,1:(num.lv+1)])
        if(!is.null(object$lv.mean)) 
            eta <- tcrossprod(cbind(1,object$lv.mean), object$lv.coefs.mean[,1:(num.lv+1)]) 
        if(!is.null(object$X.coefs.mean)) 
            eta <- eta + tcrossprod(as.matrix(X), object$X.coefs.mean) 
            }

    if(!is.null(object$row.ids) && est == "median") {
        for(j in 1:p) 
            for(k in 1:ncol(object$row.ids)) 
                eta[,j] <- eta[,j] + object$row.coefs[[k]]$median[object$row.ids[,k]] 
        }
    if(!is.null(object$row.ids) && est == "mean") {
        for(j in 1:p) 
            for(k in 1:ncol(object$row.ids)) 
                eta[,j] <- eta[,j] + object$row.coefs[[k]]$mean[object$row.ids[,k]] 
        }
    
    index_multinom_cols <- which(object$family == "multinom")
    for(j in 1:p) {
        if(object$family[j] %in% c("binomial")) 
            fitted_out[,j] <- pnorm(eta[,j])
        if(object$family[j] %in% c("beta")) 
            fitted_out[,j] <- exp(eta[,j])/(1+exp(eta[,j]))
        if(object$family[j] %in% c("poisson","lnormal","negative.binomial","tweedie","exponential","gamma")) 
            fitted_out[,j] <- exp(eta[,j])
        if(object$family[j] == "normal") 
            fitted_out[,j] <- (eta[,j]) 
# 		if(object$family[j] == "multinom") {
# 			if(est == "median") { if(!is.null(object$X.multinom.coefs.median)) eta2 <- eta[,j] + as.matrix(X)%*%object$X.multinom.coefs.median[which(index_multinom_cols == j),,] }
# 			if(est == "mean") { if(!is.null(object$X.multinom.coefs.mean)) eta2 <- eta[,j] + as.matrix(X)%*%object$X.multinom.coefs.mean[which(index_multinom_cols == j),,] }
# 			get_probs <- exp(eta2)/apply(exp(eta2),1,sum)	
# 			fitted_multinom_probs[,j,] <- get_probs
# 			}

        if(object$family[j] == "ordinal") {
            if(est == "median")
                fitted_ordinal_probs[,j,] <- ordinal.conversion.spp(n = n, lv = object$lv.median, lv.coefs.j = object$lv.coefs.median[j,], num.lv = num.lv, row.coefs = object$row.coefs, row.ids = object$row.ids, X = X, X.coefs.j = object$X.coefs.median[j,], cutoffs = object$cutoffs.median, est = "median")
            if(est == "mean")
                fitted_ordinal_probs[,j,] <- ordinal.conversion.spp(n = n, lv = object$lv.mean, lv.coefs.j = object$lv.coefs.mean[j,], num.lv = num.lv, row.coefs = object$row.coefs, row.ids = object$row.ids, X = X, X.coefs.j = object$X.coefs.mean[j,], cutoffs = object$cutoffs.mean, est = "mean")
            fitted_out[,j] <- apply(fitted_ordinal_probs[,j,],1,which.max) ## get max predicted probability
            }
        }	

    return(list(ordinal.probs = fitted_ordinal_probs, out = fitted_out))
    }

	
## Calculates DIC based on the conditional log-likelihood
get.dic <- function(jagsfit) { 
    jagsfit$BUGSoutput$DIC
    }
	
	
get.hpdintervals <- function(y, X = NULL, traits = NULL, row.ids = NULL, fit.mcmc, lv.control, prob = 0.95, num.lv = NULL) {
	
    n <- nrow(y); p <- ncol(y)
    lv.control <- check.lv.control(num.lv = num.lv, lv.control = lv.control, need.distmat = FALSE)
    num.lv <- lv.control$num.lv
            
    intervals_out <- HPDinterval(fit.mcmc, prob = prob); 
    hpd_lower <- intervals_out[,1]
    hpd_upper <- intervals_out[,2]

    lv_coefs_arr <- abind(matrix(hpd_lower[grep("lv.coefs",names(hpd_lower))], nrow=p), matrix(hpd_upper[grep("lv.coefs",names(hpd_upper))], nrow=p), along = 3)	

    final_list <- list()
    
    if(num.lv > 0) {
        lv_arr <- abind(matrix(hpd_lower[grep("lvs", names(hpd_lower))], nrow=n), matrix(hpd_upper[grep("lvs", names(hpd_upper))], nrow=n), along = 3)
        dimnames(lv_arr) <- list(rows = rownames(y), lv = paste0("lv", 1:num.lv), type = c("lower","upper"))		
        final_list$lv <- lv_arr

        if(dim(lv_coefs_arr)[2] == (num.lv+2)) 
            dimnames(lv_coefs_arr) <- list(cols = colnames(y), coefficients = c("beta0",paste0("theta",1:num.lv),"Dispersion"), type = c("lower","upper"))
        if(dim(lv_coefs_arr)[2] == (num.lv+1)) 
            dimnames(lv_coefs_arr) <- list(cols = colnames(y), coefficients = c("beta0",paste0("theta",1:num.lv)), type = c("lower","upper"))
            
        if(lv.control$type != "independent") {
            lv_covparams_arr <- cbind(hpd_lower[grep("lv.covparams", names(hpd_lower))], hpd_upper[grep("lv.covparams", names(hpd_upper))])
            if(nrow(lv_covparams_arr) == 1)
                rownames(lv_covparams_arr) <- c("spatialscale (tau1)")
            if(nrow(lv_covparams_arr) == 2)
                rownames(lv_covparams_arr) <- c("spatialscale (tau1)", "spatialpower (tau2)")
            colnames(lv_covparams_arr) <- c("lower","upper")
            final_list$lv_covparams <- lv_covparams_arr
            }
        }
    if(num.lv == 0) { 
        if(dim(lv_coefs_arr)[2] == 2) 
            dimnames(lv_coefs_arr) <- list(cols = colnames(y), coefficients = c("beta0","Dispersion"), type = c("lower","upper"))
        if(dim(lv_coefs_arr)[2] == 1) 
            dimnames(lv_coefs_arr) <- list(cols = colnames(y), coefficients = c("beta0"), type = c("lower","upper"))
        }
    final_list$lv.coefs <- lv_coefs_arr
    
    if(length(grep("row.coefs", names(hpd_lower))) > 0) {
        n.ID <- apply(row.ids, 2, function(x) length(unique(x)))
        final_list$row.coefs <- vector("list", ncol(row.ids))
        names(final_list$row.coefs) <- colnames(row.ids)
        for(k in 1:ncol(row.ids)) {
            row_coefs_arr <- cbind(
                hpd_lower[grep(paste0("row.coefs.ID",k,"\\["), names(hpd_lower))],
                hpd_upper[grep(paste0("row.coefs.ID",k,"\\["), names(hpd_upper))])
            rownames(row_coefs_arr) <- 1:n.ID[k]; colnames(row_coefs_arr) <- c("lower","upper")
                    
            final_list$row.coefs[[k]] <- row_coefs_arr
            }

        if(length(grep("row.sigma.ID", names(hpd_lower))) > 0) { 
            final_list$row.sigma <- vector("list", ncol(row.ids))
            names(final_list$row.sigma) <- colnames(row.ids)
            for(k in 1:ncol(row.ids)) {
                row_sigma_vec <- c(
                    hpd_lower[grep(paste0("row.sigma.ID",k,"$"), names(hpd_lower))],
                    hpd_upper[grep(paste0("row.sigma.ID",k,"$"), names(hpd_upper))])
                names(row_sigma_vec) <- c("lower","upper")
                            
                final_list$row.sigma[[k]] <- row_sigma_vec
                }
            }
        }

    if(length(grep("X.coefs", names(hpd_lower))) > 0) {
        X_coefs_arr <- abind(matrix(hpd_lower[grep("X.coefs", names(hpd_lower))],nrow=p), matrix(hpd_upper[grep("X.coefs", names(hpd_upper))],nrow=p), along = 3)
        dimnames(X_coefs_arr) <- list(cols = colnames(y), X = colnames(X), type = c("lower","upper"))

        final_list$X.coefs <- X_coefs_arr
        }

            
    if(length(grep("traits.coefs", names(hpd_lower))) > 0) { ## If T.params exists, then X.coefs are regressed against traits
        traitscoefs_arr <- abind(
            cbind(hpd_lower[grep("traits.int", names(hpd_lower))], matrix(hpd_lower[grep("traits.coefs", names(hpd_lower))],nrow=ncol(X)+1), hpd_lower[grep("trait.sigma", names(hpd_lower))]), 
            cbind(hpd_upper[grep("traits.int", names(hpd_upper))], matrix(hpd_upper[grep("traits.coefs", names(hpd_upper))],nrow=ncol(X)+1), hpd_upper[grep("trait.sigma", names(hpd_upper))]), 
            along = 3)
        dimnames(traitscoefs_arr) <- list(X.coefficients = c("beta0",colnames(X)), traits.coefficients = c("kappa0",colnames(traits),"sigma"), type = c("lower","upper"))
                                    
        final_list$traits.coefs <- traitscoefs_arr
        }

# 	if(length(grep("X.multinom.params", names(hpd_lower))) > 0) {
# 		final_list$X.multinom.coefs.lower <- array(matrix(hpd_lower[grep("X.multinom.params", names(hpd_lower))],dim=c(length(index_multinom_cols),ncol(X),ncol(all.X.multinom.coefs.lower)/ncol(X))), dimnames = list("1" = index_multinom_cols, "2" = colnames(X), "level" = 1:(ncol(all.X.multinom.coefs.upper)/ncol(X))))
# 		final_list$X.multinom.coefs.upper <- array(matrix(hpd_lower[grep("X.multinom.params", names(hpd_upper))],dim=c(length(index_multinom_cols),ncol(X),ncol(all.X.multinom.coefs.lower)/ncol(X))), dimnames = list("1" = index_multinom_cols, "2" = colnames(X), "level" = 1:(ncol(all.X.multinom.coefs.upper)/ncol(X))))
# 		}

    if(length(grep("cutoffs", names(hpd_lower))) > 0) { ## If cutoffs exists, then cutoffs are there and some columns involved ordinal responses
        cutoffs_arr <- cbind(hpd_lower[grep("cutoffs", names(hpd_lower))], hpd_upper[grep("cutoffs", names(hpd_upper))])
        num.ord.levels <- nrow(cutoffs_arr) + 1
        rownames(cutoffs_arr) <- paste0(1:(num.ord.levels-1),"|",2:num.ord.levels)
        colnames(cutoffs_arr) <- c("lower","upper")
    
        final_list$cutoffs <- cutoffs_arr

        if(length(grep("ordinal.sigma", names(hpd_lower))) > 0) { 
            ordinal.sigma.vec <- c(hpd_lower[grep("ordinal.sigma", names(hpd_lower))], hpd_upper[grep("ordinal.sigma", names(hpd_upper))])
            names(ordinal.sigma.vec) <- c("lower","upper")
            final_list$ordinal.sigma <- ordinal.sigma.vec
            }
        }
                            
    if(length(grep("powerparam", names(hpd_lower))) > 0) { ## If powerparam exists, then power parameters are there and some columns involved tweedie responses
        powerparam_vec <- c(hpd_lower[grep("powerparam", names(hpd_lower))], hpd_upper[grep("powerparam", names(hpd_upper))])
        names(powerparam_vec) <- c("lower","upper")
        final_list$powerparam <- powerparam_vec
        }

    rm(list = ls(pattern = ".arr"))
    return(final_list) 
    }

	
		
## Simple extraction of MCMC samples from fitted boral object
get.mcmcsamples <- function(object) {
    fit.mcmc <- object$jags.model$BUGSoutput
    if(is.null(fit.mcmc)) 
    stop("MCMC samples not found. Please use save.model = TRUE to save MCMC samples when using boral")
    fit.mcmc <- mcmc(fit.mcmc$sims.matrix, start = 1, thin = object$mcmc.control$n.thin) ## Thanks to Guilliaume Blanchet for the original formatting!

    return(fit.mcmc)
    }

     
## Intervals from marginal predictions will typically be wider than those from conditional predictions, because the former takes into account the additional uncertainty due to the lv being unknown. 
## predict.type = "conditonal": Predictions conditional on the latent variables observed and everything else e.g., random row effects. 
## predict.type = "marginal": marginalizes from the linear predictor down i.e., LVs and random row effects. Does not marginalize over species random effects as it is not clear that you need to take into account their uncertainty i.e., you are not predicting to new species!
## In both cases, if newX is supplied, then checks are made to ensure things are compatible e.g., row effects and lv match nrow(newX), 
# object = fit.traits; newX = NULL; newrow.ids = NULL; predict.type = "marginal"; est = "median"; prob = 0.95; lv.mc = 1000

predict.boral <- function(object, newX = NULL, newrow.ids = NULL, distmat = NULL, predict.type = "conditional", 
    est = "median", prob = 0.95, lv.mc = 1000, return.alllinpred = FALSE, ...) {  
     
    num.lv <- object$num.lv
    predict.type <- match.arg(predict.type, choices = c("conditional","marginal"))
    if(predict.type == "marginal") {
        message("Marginal predictions take a long time, because there is a lot of (Monte-Carlo) integration involved. Apologies in advance!")
        if(num.lv == 0) {
            message("Please note if there are no latent variables in the model, then marginal and conditional predictions are equivalent")
            predict.type <- "conditional"
            }
        }
    if(predict.type == "conditional" & !is.null(newrow.ids)) { 
        message("For conditional predictions, newrow.ids is ignored since predictions are made conditional on the set of row effects i.e., on the same set of sites")
        newrow.ids <- NULL
        }

        
    ## Check properties of X
    if(!is.null(newX)) {
        X <- as.matrix(newX) 
        if(is.null(object$X.coefs.mean)) 
            stop("Cannot find coefficients for X in object, even though you supplied newX")
        if(ncol(object$X.coefs.mean) != ncol(newX)) 
            stop("Number of columns in newX does not match number of columns in object$X.coefs.mean")
        if(predict.type == "conditional") {
            if(object$row.eff != "none") { 
                if(nrow(object$row.ids) != nrow(X)) 
                    stop("For conditional predictions, the number of rows in newX must be equal to number of rows in object$row.ids")
                }
            if(num.lv > 0) { 
                if(nrow(object$lv.mean) != nrow(X)) 
                    stop("For conditional predictions, the number of rows in newX must be equal to number of rows in object$lv.mean") 
                }
            }
        }
    if(is.null(newX)) 
        X <- object$X 
    n <- nrow(X) 

    ## Check properties of X
    if(object$lv.control$type != "independent" & is.null(distmat))
        stop("distmat needs to be supplied when object$lv.control$type is not independent")
    if(!is.null(distmat)) {
        if(nrow(distmat) != n || ncol(distmat) != n)
            stop("The dimensions of distmat do not match the dimensions of object$X or newX")
        }
    
    ## Check properties of newrow.ids; this should only be activated once the predictions are marginal 
    if(!is.null(newrow.ids)) { 
        if(object$row.eff == "none") 
            stop("Cannot find row effects parameters in object, even though you supplied in newrow.ids")
    
        newrow.ids <- as.matrix(newrow.ids)
        if(is.null(colnames(newrow.ids))) 
            colnames(newrow.ids) <- paste0("ID", 1:ncol(newrow.ids))
        if(ncol(object$row.ids) != ncol(newrow.ids)) 
            stop("The number of columns in newrow.ids must be equal to number of columns in object$row.ids")
        for(k in 1:ncol(newrow.ids)) {
            if(!all(unique(newrow.ids[,k]) %in% unique(object$row.ids[,k])))
                stop(paste0("Not all levels of newrow.ids[,",k,"] can be found in object$row.ids[,",k,"]. This is a problem as then thre will be some IDs which are unknown (as based on object)"))
            }
        }     
    if(is.null(newrow.ids)) { 
        newrow.ids <- object$row.ids 
        }     
    if(!is.null(X) & !is.null(newrow.ids)) {
        if(n != nrow(newrow.ids)) 
            stop("Number of rows in X does not match number of rows in newrow.ids")
        }

    if(is.null(object$jags.model)) 
        stop("MCMC samples not found")

#      if(control$type != "independent" & predict.type == "marginal") {
#           if(is.null(distmat))
#                stop("For marginal predictions where the latent variables are structured, please supply distmat. Note the fitted boral model does not save distmat, in order to save memory space")
#            if(nrow(distmat) != n || ncol(distmat) != n)
#                 stop("distmat should be a symmetric matrix with the same number of rows as X")
#           }
        
    #################
    ## Checks done
    #################
    
    combined_fit_mcmc <- get.mcmcsamples(object) 
    mcmc_names <- colnames(combined_fit_mcmc)
    all_linpred <- array(NA, dim = c(n, object$p, nrow(combined_fit_mcmc)))
    pt_pred <- lower_linpred <- upper_linpred <- matrix(NA, nrow = n, ncol = object$p)

    if(predict.type == "conditional") {
        for(k0 in 1:nrow(combined_fit_mcmc)) {
            if(object$row.eff != "none") 
                cw_row_coefs <- vector("list", ncol(newrow.ids))

            cw_lv_coefs <- matrix(combined_fit_mcmc[k0, grep("lv.coefs", mcmc_names)], nrow = object$p)
            cw_eta <- matrix(cw_lv_coefs[,1], nrow = n, ncol = object$p, byrow = TRUE) 
            if(num.lv > 0) {
                cw.lv <- matrix(combined_fit_mcmc[k0, grep("lvs", mcmc_names)], nrow = n)
                cw_eta <- cw_eta + tcrossprod(cw.lv, cw_lv_coefs[,2:(num.lv+1)])
                }
            if(!is.null(X)) {
                cw_X_coefs <- matrix(combined_fit_mcmc[k0, grep("X.coefs", mcmc_names)], nrow = object$p) 
                cw_eta <- cw_eta + tcrossprod(X, cw_X_coefs)
                }
            if(!is.null(object$offset)) 
                cw_eta <- cw_eta + object$offset
            if(object$row.eff != "none") { 
                for(k in 1:ncol(newrow.ids)) {
                    cw_row_coefs[[k]] <- combined_fit_mcmc[k0, grep(paste0("row.coefs.ID",k,"\\["), mcmc_names)] 
                    cw_eta <- cw_eta + cw_row_coefs[[k]][newrow.ids[,k]] 
                    }
                }

            all_linpred[,,k0] <- cw_eta
            }
        }


    if(predict.type == "marginal") {
        mc_lv <- rmvnorm(n*lv.mc, mean = rep(0,num.lv))
        mc_lv <- array(c(mc_lv), dim = c(lv.mc, n, num.lv))
        
        for(k0 in 1:nrow(combined_fit_mcmc)) {
        
            if(k0 %% 100 == 0) 
                message("Onto MCMC sample ", k0)              
            if(!is.null(X))# & is.null(object$traits)) 
                cw_X_coefs <- matrix(combined_fit_mcmc[k0, grep("X.coefs", mcmc_names)], nrow = object$p) 
#                if(!is.null(object$traits)) { ## This marginalizes over the X coefs wrt traits as well
#                     cw_trait_params <- cbind(combined_fit_mcmc[k0, grep("traits.int", mcmc_names)], 
#                          matrix(combined_fit_mcmc[k0, grep("traits.coefs", mcmc_names)], nrow = ncol(X)+1), 
#                          combined_fit_mcmc[k0, grep("trait.sigma", mcmc_names)]) 
#                     mean_X_coefs <- tcrossprod(cbind(1,object$traits), cw_trait_params[,-ncol(cw_trait_params)])
#                     cw_X_coefs <- array(0, dim = c(lv.mc, object$p, ncol(X)+1)) ## Includes intercept
#                     for(l0 in 1:(ncol(X)+1)) 
#                          cw_X_coefs[,,l0] <- rmvnorm(lv.mc, mean = mean_X_coefs[,l0], sigma = diag(x = cw_trait_params[l0,ncol(cw_trait_params)]^2, nrow=object$p))
#                     rm(mean_X_coefs, cw_trait_params)
#                     }
#                if(object$lv.control$type != "independent") { 
#                     cw_lv_covparams <- combined_fit_mcmc[k0, grep("lv.covparams", mcmc_names)] 
#                     } 
            if(object$row.eff == "fixed") { 
                cw_row_coefs <- vector("list", ncol(newrow.ids))
                for(k in 1:ncol(newrow.ids)) 
                    cw_row_coefs[[k]] <- combined_fit_mcmc[k0, grep(paste0("row.coefs.ID",k,"\\["), mcmc_names)] 
                } 
            if(object$row.eff == "random") {
                cw_row_coefs <- vector("list", ncol(newrow.ids))
                ## Need to generate from length(unique(object$row.ids[,k])) to account for fact that we may have less levels in newrow.ids (which is OK)
                for(k in 1:ncol(newrow.ids)) 
                    cw_row_coefs[[k]] <- matrix(rnorm(length(unique(object$row.ids[,k]))*lv.mc, mean = 0, sd = combined_fit_mcmc[k0, grep(paste0("row.sigma.ID",k,"$"), mcmc_names)]), ncol = lv.mc) 
                }
            cw_lv_coefs <- matrix(combined_fit_mcmc[k0, grep("lv.coefs", mcmc_names)], nrow = object$p)
            all_linpredmc <- array(NA, dim = c(n, object$p, lv.mc))

            cw_lv_covparams <- combined_fit_mcmc[k0, grep("lv.covparams", mcmc_names)]
            if(object$lv.control$type == "exponential")
                covmat_chol <- chol(exp(-distmat/cw_lv_covparams[1]))
            if(object$lv.control$type == "squared.exponential")
                covmat_chol <- (chol(exp(-(distmat/cw_lv_covparams[1])^2)))
           if(object$lv.control$type == "powered.exponential")
                covmat_chol <- (chol(exp(-(distmat/cw_lv_covparams[1])^cw_lv_covparams[2])))
           if(object$lv.control$type == "spherical")
                covmat_chol <- (chol((distmat < cw_lv_covparams[1])*(1 - 1.5*distmat/cw_lv_covparams[1] + 0.5*(distmat/cw_lv_covparams[1])^3)))


            for(b in 1:lv.mc) {
                if(object$lv.control$type != "independent")
                    mc_lv[b,,] <- crossprod(covmat_chol, mc_lv[b,,])                    
                
                cw_eta <- tcrossprod(mc_lv[b,,], cw_lv_coefs[,2:(num.lv+1)]) ## LVs to marginalize over
                if(!is.null(X))# & is.null(object$traits)) ## Spp intercepts, X_coefs, and no traits
                    cw_eta <- cw_eta + tcrossprod(cbind(1,X),cbind(cw_lv_coefs[,1],cw_X_coefs))
#                     if(!is.null(object$traits)) ## Spp intercepts, X_coefs to marginalize over due to traits
#                          cw_eta <- cw_eta + tcrossprod(cbind(1,X),cw_X_coefs[b,,])

                if(!is.null(object$offset)) 
                    cw_eta <- cw_eta + object$offset
                if(object$row.eff == "fixed") { 
                    for(k in 1:ncol(newrow.ids)) 
                        cw_eta <- cw_eta + cw_row_coefs[[k]][newrow.ids[,k]] 
                        }
                if(object$row.eff == "random") { 
                    for(k in 1:ncol(newrow.ids)) 
                        cw_eta <- cw_eta + cw_row_coefs[[k]][newrow.ids[,k],b] 
                    }
                all_linpredmc[,,b] <- cw_eta
                rm(cw_eta)
                }
                
            all_linpred[,,k0] <- apply(all_linpredmc, c(1,2), mean)
            rm(all_linpredmc)
            }
        }

        
    for(i in 1:n) { for(j in 1:object$p) {
        if(est == "mean") 
            pt_pred[i,j] <- mean(all_linpred[i,j,]) ## Posterior mean
        if(est == "median") 
            pt_pred[i,j] <- median(all_linpred[i,j,]) ## Posterior median
        lower_linpred[i,j] <- quantile(all_linpred[i,j,], probs = (1-prob)/2)
        upper_linpred[i,j] <- quantile(all_linpred[i,j,], probs = 1-(1-prob)/2)
        } }
    
    out <- list(linpred = pt_pred, lower = lower_linpred, upper = upper_linpred)
    if(return.alllinpred) 
        out$all.linpred <- all_linpred

    return(out)
    }
    


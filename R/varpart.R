## Variance partitioning on the linear predictor scale; the variance is calculated and and averaged wrt to the posterior distribution 
calc.varpart <- function(object, groupX = NULL) {
    fit.mcmc <- get.mcmcsamples(object)
    num.lv <- object$num.lv
    if(is.null(object$X))
        stop("Variance partitioning is not particularly useful when there are no covariates X included in the model.")
    if(is.null(object$jags.model)) 
        stop("MCMC samples not found.")
    if(!is.null(groupX)) { if(length(groupX) != (object$num.X+1)) 
        stop("If groupX is supplied, then it must be a vector with the length (object$num.X + 1).") }

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



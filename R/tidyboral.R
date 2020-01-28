tidyboral <- function(object) 
    {
    message("This function is purely to produce \"tidier\" output for boral. Please do not overwrite the original boral object!")
    tidyout_fit <- list()
    
    lv_coefs_arr <- abind(
          object$lv.coefs.median,
          object$lv.coefs.mean,
          object$lv.coefs.iqr,
          object$lv.coefs.sd,
        along = 3)

        
     if(object$num.lv > 0)
          {
          lv_arr <- abind(
               object$lv.median,
               object$lv.mean,
               object$lv.iqr,
               object$lv.sd,
               along = 3)
          dimnames(lv_arr) <- list(rows = rownames(object$y), lv = paste0("lv", 1:object$num.lv), type = c("median","mean","iqr","sd"))
          lv_arr <- melt(lv_arr, value.name = "estimate", varnames = names(dimnames(lv_arr)))
    
          tidyout_fit$lv <- lv_arr
          rm(lv_arr)
          
          if(dim(lv_coefs_arr)[2] == (object$num.lv+2)) 
               dimnames(lv_coefs_arr) <- list(cols = colnames(object$y), 
                    coefficients = c("beta0", paste0("theta", 1:object$num.lv), "Dispersion"), type = c("median","mean","iqr","sd"))
          if(dim(lv_coefs_arr)[2] == (object$num.lv+1)) 
               dimnames(lv_coefs_arr) <- list(cols = colnames(object$y), 
                    coefficients = c("beta0", paste0("theta", 1:object$num.lv)), type = c("median","mean","iqr","sd"))

          if(object$lv.control$type != "independent") 
               {
               lv_params_arr <- cbind(object$lv.covparams.median, object$lv.covparams.mean, object$lv.covparams.iqr, object$lv.covparams.sd) 
               if(nrow(lv_params_arr) == 1) 
                    rownames(lv_params_arr) <- c("spatialscale (tau1)")
               if(nrow(lv_params_arr) == 2) 
                    rownames(lv_params_arr) <- c("spatialscale (tau1)", "spatialpower (tau2)")
               colnames(lv_params_arr) <- c("median","mean","iqr","sd")
               lv_params_arr <- melt(lv_params_arr, value.name = "estimate", varnames = c("parameters", "type"))
               tidyout_fit$lv.covparams <- lv_params_arr
               rm(lv_params_arr)
               }
          }
    
    if(object$num.lv == 0)
          {
          if(dim(lv_coefs_arr)[2] == 2) 
            dimnames(lv_coefs_arr) <- list(cols = colnames(object$y), 
               coefficients = c("beta0", "Dispersion"), type = c("median","mean","iqr","sd"))
        if(dim(lv_coefs_arr)[2] == 1) 
            dimnames(lv_coefs_arr) <- list(cols = colnames(object$y), 
               coefficients = c("beta0"), type = c("median","mean","iqr","sd"))
          }
     
     lv_coefs_arr <- melt(lv_coefs_arr, value.name = "estimate", varnames = names(dimnames(lv_coefs_arr)))
     tidyout_fit$lv.coefs <- lv_coefs_arr
     rm(lv_coefs_arr)

     
     if(object$row.eff != "none") 
          {
          tidyout_fit$row.coefs <- vector("list", ncol(object$row.ids))
          names(tidyout_fit$row.coefs) <- colnames(object$row.ids)
          for(k in 1:ncol(object$row.ids)) 
               {
               row_coefs_arr <- cbind(object$row.coefs[[k]]$median, object$row.coefs[[k]]$mean, object$row.coefs[[k]]$iqr, object$row.coefs[[k]]$sd)
               dimnames(row_coefs_arr) <- list(rowID = 1:length(unique(object$row.ids[,k])), type = c("median","mean","iqr","sd"))
               row_coefs_arr <- melt(row_coefs_arr, value.name = "estimate", varnames = names(dimnames(row_coefs_arr)))
               tidyout_fit$row.coefs[[k]] <- row_coefs_arr
               rm(row_coefs_arr)
               }
    
          if(object$row.eff == "random") 
               {
               tidyout_fit$row.sigma <- vector("list", ncol(object$row.ids))
               names(tidyout_fit$row.sigma) <- colnames(object$row.ids)
               for(k in 1:ncol(object$row.ids)) 
                    {
                    tidyout_fit$row.sigma[[k]] <- melt(object$row.sigma[[k]], value.name = "estimate")
                    }
               }
          }

                    
     if(object$num.X > 0) 
          {
          X_coefs_arr <- abind(
            object$X.coefs.median,
            object$X.coefs.mean,
            object$X.coefs.iqr,
            object$X.coefs.sd,
            along = 3)
          dimnames(X_coefs_arr) <- list(cols = colnames(object$y), coefficients = colnames(object$X), type = c("median","mean","iqr","sd"))
          X_coefs_arr <- melt(X_coefs_arr, value.name = "estimate", varnames = names(dimnames(X_coefs_arr)))
          tidyout_fit$X.coefs <- X_coefs_arr
          rm(X_coefs_arr)

          if(any(object$prior.control$ssvs.index == 0)) 
               { ## You should not be able to enter this loop if num.traits > 0!
               ssvs_indcoefs_arr <- abind(
                    object$ssvs.indcoefs.mean,
                    object$ssvs.indcoefs.sd,
                    along = 3)
               dimnames(ssvs_indcoefs_arr) <- list(cols = colnames(object$y), coefficients = colnames(object$X), type = c("mean","sd"))
               ssvs_indcoefs_arr <- melt(ssvs_indcoefs_arr, value.name = "estimate", varnames = names(dimnames(ssvs_indcoefs_arr)))
               tidyout_fit$ssvs.indcoefs <- ssvs_indcoefs_arr
               rm(ssvs_indcoefs_arr)
               }
          if(any(object$prior.control$ssvs.index > 0))   
               {
               ssvs_gpcoefs_arr <- cbind(object$ssvs.gpcoefs.mean, object$ssvs.gpcoefs.sd)
               dimnames(ssvs_gpcoefs_arr) <- list(Gp = rownames(object$ssvs.gpcoefs.mean), type = c("mean","sd"))
               ssvs_gpcoefs_arr <- melt(ssvs_gpcoefs_arr, value.name = "estimate", varnames = names(dimnames(ssvs_gpcoefs_arr)))
               tidyout_fit$ssvs.gpcoefs <- ssvs_gpcoefs_arr
               rm(ssvs_gpcoefs_arr)
               }
        
          if(any(unlist(object$prior.control$ssvs.traitsindex) == 0)) 
               { 
               ssvs_traitcoefs_arr <- abind(
                    object$ssvs.traitscoefs.mean,
                    object$ssvs.traitscoefs.sd,
                    along = 3)               
               dimnames(ssvs_traitcoefs_arr) <- list(X.coefficients = c("beta0",colnames(object$X)), 
                    traits.coefficients =  colnames(object$traits), type = c("mean","sd"))
               ssvs_traitcoefs_arr <- melt(ssvs_traitcoefs_arr, value.name = "estimate", varnames = names(dimnames(ssvs_traitcoefs_arr)))
               tidyout_fit$ssvs.traitscoefs <- ssvs_traitcoefs_arr
               rm(ssvs_traitcoefs_arr)               
               }
          }
          
        
        if(object$num.traits > 0) 
            {
            traitcoefs_arr <- abind(
                object$traits.coefs.median,
                object$traits.coefs.mean,
                object$traits.coefs.iqr,
                object$traits.coefs.sd,
                along = 3)
            dimnames(traitcoefs_arr) <- list(X.coefficients = c("beta0",colnames(object$X)), traits.coefficients =  c("kappa0",colnames(object$traits),"sigma"), type = c("median","mean","iqr","sd"))
            traitcoefs_arr <- melt(traitcoefs_arr, value.name = "estimate", varnames = names(dimnames(traitcoefs_arr)))
            tidyout_fit$traits.coefs <- traitcoefs_arr
            rm(traitcoefs_arr)
            }
            
         
        if(any(object$family == "ordinal")) 
            {
            cutoffs_arr <- cbind(object$cutoffs.median, object$cutoffs.mean, object$cutoffs.iqr, object$cutoffs.sd)
            dimnames(cutoffs_arr) <- list(parameters = paste0(1:(object$num.ord.levels - 1), "|", 2:object$num.ord.levels), type = c("median","mean","iqr","sd"))
            cutoffs_arr <- melt(cutoffs_arr, value.name = "estimate", varnames = names(dimnames(cutoffs_arr)))
            tidyout_fit$cutoffs <- cutoffs_arr
            rm(cutoffs_arr)
              
            ## If there are traits, then ordinal random intercept is either zero (if there is only 1 ordinal column, or has the trait.sigma (if there are >1 ordinal columns)
            if(sum(object$family == "ordinal") > 1 & is.null(object$traits)) 
                { 
                ordinal_sigma_vec <- c(object$ordinal.sigma.median, object$ordinal.sigma.mean, object$ordinal.sigma.iqr, object$ordinal.sigma.sd)
                names(ordinal_sigma_vec) <- c("median","mean","iqr","sd")
                ordinal_sigma_vec <- melt(ordinal_sigma_vec)
                colnames(ordinal_sigma_vec) <- c("estimate","type")
                tidyout_fit$ordinal.sigma <- ordinal_sigma_vec
                rm(ordinal_sigma_vec)
                }
            }

        if(any(object$family == "tweedie")) 
            {
            powerparam_vec <- c(object$powerparam.median, object$powerparam.mean, object$powerparam.iqr, object$powerparam.sd)
            names(powerparam_vec) <- c("median","mean","iqr","sd")
            powerparam_vec <- melt(powerparam_vec)
            colnames(powerparam_vec) <- c("estimate","type")
            tidyout_fit$powerparam <- powerparam_vec
            rm(powerparam_vec)                
            }

        
    tidyout_fit$hpdintervals <- tidyboral.hpdintervals(object)
    
    return(tidyout_fit)
    }
    

    
tidyboral.hpdintervals <- function(object) 
    {
    n <- nrow(object$y)
    p <- ncol(object$y)
    num.lv <- object$lv.control$num.lv
    tidyout_fit <- list()

    
    if(num.lv > 0) 
        {
        tidyout_fit$lv <- melt(object$hpdintervals$lv, value.name = "estimate", varnames = names(dimnames(object$hpdintervals$lv)))
               
        if(object$lv.control$type != "independent") 
            {
            lv_covparams_arr <- melt(object$hpdintervals$lv.covparams, value.name = "estimate")
            colnames(lv_covparams_arr) <- c("parameters", "type", "estimate")
            tidyout_fit$lv.covparams <- lv_covparams_arr
            rm(lv_covparams_arr)
            }
          }
    tidyout_fit$lv.coefs <- melt(object$hpdintervals$lv.coefs, value.name = "estimate", varnames = names(dimnames(object$hpdintervals$lv.coefs)))
    
    if(!is.null(object$hpdintervals$row.coefs)) 
        {
        tidyout_fit$row.coefs <- vector("list", ncol(object$row.ids))
        names(tidyout_fit$row.coefs) <- colnames(object$row.ids)
        for(k in 1:ncol(object$row.ids)) 
            {
            tidyout_fit$row.coefs[[k]] <- melt(object$hpdintervals$row.coefs[[k]], value.name = "estimate")
            colnames(tidyout_fit$row.coefs[[k]]) <- c("rowID", "type", "estimate")
            }
    
        if(!is.null(object$hpdintervals$row.sigma)) 
            { 
            tidyout_fit$row.sigma <- vector("list", ncol(object$row.ids))
            names(tidyout_fit$row.sigma) <- colnames(object$row.ids)
            for(k in 1:ncol(object$row.ids)) 
                {
                tidyout_fit$row.sigma[[k]] <- melt(object$hpdintervals$row.sigma[[k]], value.name = "estimate")
                }
            }
        }
    
    if(!is.null(object$hpdintervals$X.coefs)) 
        {
        tidyout_fit$X.coefs <- melt(object$hpdintervals$X.coefs, value.name = "estimate", varnames = names(dimnames(object$hpdintervals$X.coefs)))    
        }
    
    
    if(!is.null(object$hpdintervals$traits.coefs))
        { ## If T.params exists, then X.coefs are regressed against traits
        tidyout_fit$traits.coefs <- melt(object$hpdintervals$traits.coefs, value.name = "estimate", varnames = names(dimnames(object$hpdintervals$traits.coefs)))    
        }
    
    if(!is.null(object$hpdintervals$cutoffs))
        { ## If cutoffs exists, then cutoffs are there and some columns involved ordinal responses
        dimnames(object$hpdintervals$cutoffs) <- list(parameters = paste0(1:(object$num.ord.levels - 1), "|", 2:object$num.ord.levels), type = c("median","mean","iqr","sd"))
        tidyout_fit$cutoffs <- melt(object$hpdintervals$cutoffs, value.name = "estimate", varnames = names(dimnames(object$hpdintervals$cutoffs)))
     
        if(!is.null(object$hpdintervals$ordinal.sigma))
            { 
            tidyout_fit$ordinal.sigma <- melt(object$hpdintervals$ordinal.sigma, value.name = "estimate")
            }
        }
                              
    if(!is.null(object$hpdintervals$powerparam))
          { ## If powerparam exists, then power parameters are there and some columns involved tweedie responses
          tidyout_fit$powerparam <- melt(object$hpdintervals$powerparam, value.name = "estimate")
          }

          
    return(tidyout_fit)
    }
    
    
    
    

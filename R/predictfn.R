## Intervals from marginal predictions will typically be wider than those from conditional predictions, because the former takes into account the additional uncertainty due to the lv being unknown. 
## predict.type = "conditonal": Predictions conditional on the latent variables observed and everything else e.g., random row effects. 
## predict.type = "marginal": marginalizes from the linear predictor down i.e., LVs and random row effects. Does not marginalize over species random effects as it is not clear that you need to take into account their uncertainty i.e., you are not predicting to new species!
## In both cases, if newX is supplied, then checks are made to ensure things are compatible e.g., row effects and lv match nrow(newX), 
# object = fit.traits; newX = NULL; newrow.ids = NULL; predict.type = "marginal"; est = "median"; prob = 0.95; lv.mc = 1000

predict.boral <- function(object, newX = NULL, newrow.ids = NULL, distmat = NULL, predict.type = "conditional", 
     scale = "link", est = "median", prob = 0.95, lv.mc = 1000, return.alllinpred = FALSE, ...) 
     {  
     
     num.lv <- object$num.lv
     scale <- match.arg(scale, choices = c("link","response"))
     predict.type <- match.arg(predict.type, choices = c("conditional","marginal"))
     if(predict.type == "marginal") 
          {
          message("Marginal predictions take a long time, because there is a lot of (Monte-Carlo) integration involved. Apologies in advance!")
          if(num.lv == 0) 
               {
               message("Please note if there are no latent variables in the model, then marginal and conditional predictions are equivalent.")
               predict.type <- "conditional"
               }
          }
     if(predict.type == "conditional" & !is.null(newrow.ids)) 
          { 
          message("For conditional predictions, newrow.ids is ignored since predictions are made conditional on the set of row effects i.e., on the same set of sites.")
          newrow.ids <- NULL
          }

        
     ## Check properties of X
     if(!is.null(newX)) 
          {
          X <- as.matrix(newX) 
          if(is.null(object$X.coefs.mean)) 
               stop("Cannot find coefficients for X in object, even though you supplied newX")
          if(ncol(object$X.coefs.mean) != ncol(newX)) 
               stop("Number of columns in newX does not match number of columns in object$X.coefs.mean")
          if(predict.type == "conditional") 
               {
               if(object$row.eff != "none") 
                    { 
                    if(nrow(object$row.ids) != nrow(X)) 
                         stop("For conditional predictions, the number of rows in newX must be equal to number of rows in object$row.ids")
                    }
               if(num.lv > 0) 
                    { 
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
     if(!is.null(newrow.ids)) 
          { 
          if(object$row.eff == "none") 
               stop("Cannot find row effects parameters in object, even though you supplied in newrow.ids")
     
          newrow.ids <- as.matrix(newrow.ids)
          if(is.null(colnames(newrow.ids))) 
               colnames(newrow.ids) <- paste0("ID", 1:ncol(newrow.ids))
          if(ncol(object$row.ids) != ncol(newrow.ids)) 
               stop("The number of columns in newrow.ids must be equal to number of columns in object$row.ids")
          if(object$row.eff == "fixed")
               {
               for(k1 in 1:ncol(newrow.ids)) 
                    {
                    if(!all(unique(newrow.ids[,k1]) %in% unique(object$row.ids[,k1])))
                         stop(paste0("Not all levels of newrow.ids[,",k1,"] can be found in object$row.ids[,",k1,"]. When the 
                         row effects are fixed, this is a problem as then thre will be some IDs which are unknown (as based on object)."))
                    }
               }
          }     
     if(is.null(newrow.ids)) 
          { 
          newrow.ids <- object$row.ids 
          }        
     if(!is.null(X) & !is.null(newrow.ids)) 
          {
          if(n != nrow(newrow.ids)) 
               stop("Number of rows in X does not match number of rows in newrow.ids")
          }

     if(is.null(object$jags.model)) 
          stop("MCMC samples not found")

    
     ##-----------------------------
     ## Checks done
     ##-----------------------------
     
     combined_fit_mcmc <- get.mcmcsamples(object) 
     mcmc_names <- colnames(combined_fit_mcmc)
     all_linpred <- array(NA, dim = c(n, object$p, nrow(combined_fit_mcmc)))
     pt_pred <- lower_linpred <- upper_linpred <- matrix(NA, nrow = n, ncol = object$p)

     if(predict.type == "conditional") 
          {
          for(k0 in 1:nrow(combined_fit_mcmc)) 
               {
               if(object$row.eff != "none") 
                    cw_row_coefs <- vector("list", ncol(newrow.ids))

               cw_lv_coefs <- matrix(combined_fit_mcmc[k0, grep("lv.coefs", mcmc_names)], nrow = object$p)
               cw_eta <- matrix(cw_lv_coefs[,1], nrow = n, ncol = object$p, byrow = TRUE) 
               if(num.lv > 0) 
                    {
                    cw.lv <- matrix(combined_fit_mcmc[k0, grep("lvs", mcmc_names)], nrow = n)
                    cw_eta <- cw_eta + tcrossprod(cw.lv, cw_lv_coefs[,2:(num.lv+1)])
                    }
               if(!is.null(X)) 
                    {
                    cw_X_coefs <- matrix(combined_fit_mcmc[k0, grep("X.coefs", mcmc_names)], nrow = object$p) 
                    cw_eta <- cw_eta + tcrossprod(X, cw_X_coefs)
                    }
               if(!is.null(object$offset)) 
                    cw_eta <- cw_eta + object$offset
               if(object$row.eff != "none") 
                    { 
                    for(k1 in 1:ncol(newrow.ids)) 
                         {
                         cw_row_coefs[[k1]] <- combined_fit_mcmc[k0, grep(paste0("row.coefs.ID",k1,"\\["), mcmc_names)] 
                         cw_eta <- cw_eta + cw_row_coefs[[k1]][newrow.ids[,k1]] 
                         }
                    }

               all_linpred[,,k0] <- cw_eta
               }
          }


     if(predict.type == "marginal") 
          {
          mc_lv <- rmvnorm(n*lv.mc, mean = rep(0,num.lv))
          mc_lv <- array(c(mc_lv), dim = c(lv.mc, n, num.lv))
          
          for(k0 in 1:nrow(combined_fit_mcmc)) 
               {
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
               if(object$row.eff == "fixed") 
                    { 
                    cw_row_coefs <- vector("list", ncol(newrow.ids))
                    for(k1 in 1:ncol(newrow.ids)) 
                         cw_row_coefs[[k1]] <- combined_fit_mcmc[k0, grep(paste0("row.coefs.ID",k1,"\\["), mcmc_names)] 
                    } 
               if(object$row.eff == "random") 
                    {
                    cw_row_coefs <- vector("list", ncol(newrow.ids))
                    for(k1 in 1:ncol(newrow.ids)) 
                         cw_row_coefs[[k1]] <- matrix(rnorm(length(unique(newrow.ids[,k1]))*lv.mc, mean = 0, sd = combined_fit_mcmc[k0, grep(paste0("row.sigma.ID",k1,"$"), mcmc_names)]), ncol = lv.mc) 
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


               for(b0 in 1:lv.mc) 
                    {
                    if(object$lv.control$type != "independent")
                         mc_lv[b0,,] <- crossprod(covmat_chol, mc_lv[b0,,])                    
                    
                    cw_eta <- tcrossprod(mc_lv[b0,,], cw_lv_coefs[,2:(num.lv+1)]) ## LVs to marginalize over
                    if(!is.null(X))# & is.null(object$traits)) ## Spp intercepts, X_coefs, and no traits
                         cw_eta <- cw_eta + tcrossprod(cbind(1,X),cbind(cw_lv_coefs[,1],cw_X_coefs))
     #                     if(!is.null(object$traits)) ## Spp intercepts, X_coefs to marginalize over due to traits
     #                          cw_eta <- cw_eta + tcrossprod(cbind(1,X),cw_X_coefs[b0,,])

                    if(!is.null(object$offset)) 
                         cw_eta <- cw_eta + object$offset
                    if(object$row.eff == "fixed") 
                         { 
                         for(k in 1:ncol(newrow.ids)) 
                         cw_eta <- cw_eta + cw_row_coefs[[k]][newrow.ids[,k]] 
                         }
                    if(object$row.eff == "random") 
                         { 
                         for(k in 1:ncol(newrow.ids)) 
                         cw_eta <- cw_eta + cw_row_coefs[[k]][newrow.ids[,k],b0] 
                         }
                    all_linpredmc[,,b0] <- cw_eta
                    rm(cw_eta)
                    }
                    
               all_linpred[,,k0] <- apply(all_linpredmc, c(1,2), mean)
               rm(all_linpredmc)
               }
          }
        
        
     if(scale == "response") 
          {
          for(j in 1:object$p) 
               {
               if(object$family[j] == "binomial") 
                    all_linpred[,j,] <- pnorm(all_linpred[,j,])
               if(object$family[j] %in% c("poisson","negative.binomial","exponential","gamma","lnormal","tweedie")) 
                    all_linpred[,j,] <- exp(all_linpred[,j,])
               if(object$family[j] == "beta") 
                    all_linpred[,j,] <- binomial()$linkinv(all_linpred[,j,])
               if(object$family[j] == "ordinal") {
                    stop("Predictions on the response scale are currently not automatically available for ordinal responses, due to the added complexity in its calculation...sorry!")
                    }
               }
          }
        
        
     ## Construct numerical measures for predictions
     for(i in 1:n) { for(j in 1:object$p) 
          {
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


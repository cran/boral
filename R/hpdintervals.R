get.hpdintervals <- function(y, X = NULL, traits = NULL, row.ids = NULL, fit.mcmc, lv.control, prob = 0.95, num.lv = NULL) 
     {
     n <- nrow(y); p <- ncol(y)
     lv.control <- check_lv_control(num.lv = num.lv, lv.control = lv.control, need.distmat = FALSE)
     num.lv <- lv.control$num.lv
               
     intervals_out <- HPDinterval(fit.mcmc, prob = prob); 
     hpd_lower <- intervals_out[,1]
     hpd_upper <- intervals_out[,2]

     lv_coefs_arr <- abind(matrix(hpd_lower[grep("lv.coefs",names(hpd_lower))], nrow=p), 
          matrix(hpd_upper[grep("lv.coefs",names(hpd_upper))], nrow=p), 
          along = 3)
     final_list <- list()
     
     if(num.lv > 0) 
          {
          lv_arr <- abind(matrix(hpd_lower[grep("lvs", names(hpd_lower))], nrow=n), matrix(hpd_upper[grep("lvs", names(hpd_upper))], nrow=n), along = 3)
          dimnames(lv_arr) <- list(rows = rownames(y), lv = paste0("lv", 1:num.lv), type = c("lower","upper"))		
          final_list$lv <- lv_arr

          if(dim(lv_coefs_arr)[2] == (num.lv+2)) 
               dimnames(lv_coefs_arr) <- list(cols = colnames(y), coefficients = c("beta0",paste0("theta",1:num.lv),"Dispersion"), type = c("lower","upper"))
          if(dim(lv_coefs_arr)[2] == (num.lv+1)) 
               dimnames(lv_coefs_arr) <- list(cols = colnames(y), coefficients = c("beta0",paste0("theta",1:num.lv)), type = c("lower","upper"))
               
          if(lv.control$type != "independent") 
               {
               lv_covparams_arr <- cbind(hpd_lower[grep("lv.covparams", names(hpd_lower))], hpd_upper[grep("lv.covparams", names(hpd_upper))])
               if(nrow(lv_covparams_arr) == 1)
                    rownames(lv_covparams_arr) <- c("spatialscale (tau1)")
               if(nrow(lv_covparams_arr) == 2)
                    rownames(lv_covparams_arr) <- c("spatialscale (tau1)", "spatialpower (tau2)")
               colnames(lv_covparams_arr) <- c("lower","upper")
               final_list$lv.covparams <- lv_covparams_arr
               }
          }
     if(num.lv == 0)     
          { 
          if(dim(lv_coefs_arr)[2] == 2) 
               dimnames(lv_coefs_arr) <- list(cols = colnames(y), coefficients = c("beta0","Dispersion"), type = c("lower","upper"))
          if(dim(lv_coefs_arr)[2] == 1) 
               dimnames(lv_coefs_arr) <- list(cols = colnames(y), coefficients = c("beta0"), type = c("lower","upper"))
          }
     final_list$lv.coefs <- lv_coefs_arr
     
     if(length(grep("row.coefs", names(hpd_lower))) > 0) 
          {
          n.ID <- apply(row.ids, 2, function(x) length(unique(x)))
          final_list$row.coefs <- vector("list", ncol(row.ids))
          names(final_list$row.coefs) <- colnames(row.ids)
          for(k in 1:ncol(row.ids)) 
               {
               row_coefs_arr <- cbind(
                    hpd_lower[grep(paste0("row.coefs.ID",k,"\\["), names(hpd_lower))],
                    hpd_upper[grep(paste0("row.coefs.ID",k,"\\["), names(hpd_upper))])
               rownames(row_coefs_arr) <- 1:n.ID[k]; colnames(row_coefs_arr) <- c("lower","upper")
                         
               final_list$row.coefs[[k]] <- row_coefs_arr
               }

          if(length(grep("row.sigma.ID", names(hpd_lower))) > 0) 
               { 
               final_list$row.sigma <- vector("list", ncol(row.ids))
               names(final_list$row.sigma) <- colnames(row.ids)
               for(k in 1:ncol(row.ids)) 
                    {
                    row_sigma_vec <- c(
                         hpd_lower[grep(paste0("row.sigma.ID",k,"$"), names(hpd_lower))],
                         hpd_upper[grep(paste0("row.sigma.ID",k,"$"), names(hpd_upper))])
                    names(row_sigma_vec) <- c("lower","upper")
                              
                    final_list$row.sigma[[k]] <- row_sigma_vec
                    }
               }
          }

     if(length(grep("X.coefs", names(hpd_lower))) > 0) 
          {
          X_coefs_arr <- abind(matrix(hpd_lower[grep("X.coefs", names(hpd_lower))],nrow=p), matrix(hpd_upper[grep("X.coefs", names(hpd_upper))],nrow=p), along = 3)
          dimnames(X_coefs_arr) <- list(cols = colnames(y), coefficients = colnames(X), type = c("lower","upper"))

          final_list$X.coefs <- X_coefs_arr
          }

               
     if(length(grep("traits.coefs", names(hpd_lower))) > 0) 
          { ## If T.params exists, then X.coefs are regressed against traits
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

     if(length(grep("cutoffs", names(hpd_lower))) > 0) 
          { ## If cutoffs exists, then cutoffs are there and some columns involved ordinal responses
          cutoffs_arr <- cbind(hpd_lower[grep("cutoffs", names(hpd_lower))], hpd_upper[grep("cutoffs", names(hpd_upper))])
          num.ord.levels <- nrow(cutoffs_arr) + 1
          rownames(cutoffs_arr) <- paste0(1:(num.ord.levels-1),"|",2:num.ord.levels)
          colnames(cutoffs_arr) <- c("lower","upper")
     
          final_list$cutoffs <- cutoffs_arr

          if(length(grep("ordinal.sigma", names(hpd_lower))) > 0) 
               { 
               ordinal.sigma.vec <- c(hpd_lower[grep("ordinal.sigma", names(hpd_lower))], hpd_upper[grep("ordinal.sigma", names(hpd_upper))])
               names(ordinal.sigma.vec) <- c("lower","upper")
               final_list$ordinal.sigma <- ordinal.sigma.vec
               }
          }
                              
     if(length(grep("powerparam", names(hpd_lower))) > 0) 
          { ## If powerparam exists, then power parameters are there and some columns involved tweedie responses
          powerparam_vec <- c(hpd_lower[grep("powerparam", names(hpd_lower))], hpd_upper[grep("powerparam", names(hpd_upper))])
          names(powerparam_vec) <- c("lower","upper")
          final_list$powerparam <- powerparam_vec
          }

     rm(list = ls(pattern = ".arr"))
     return(final_list) 
     }

	

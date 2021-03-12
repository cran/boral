##############
## Bayesian Ordination and Regression AnaLysis 
## Ordinal data handled as propotional odds regression, with common cutoff points for all spp, but spp intercepts included as random effects; This should make sense if the same ordinal scale is applied to all responses
## Multinomial regression is written, but CURRENTLY NOT ALLOWED DUE TO COMPUTATION TIME
################

## Changes from v1.9 (news files to be updated, and remember to run runningexamples.R)
## - Quite a bit of editing done to the help files and some names changes to returned objects. For help files, in particular note the notation for latent variables has been changed from z_i to u_i, and random intercept covariates represented by z_i
## - coefsplot now does fourth corner coefficient plotting
## - Attempting to gradually add response-specific random intercepts. This is done with ranef.ids. 
##   + R scripts done [NC for no change -- ALL DONE]: 
##   + Help files done [NC for no change]: 
## - Including an argument formula.X in the boral fitting function for users to more readily specifiy a covariate matrix
## - ICs and log-likelihood are now formally defunct and no longer maintained!

## TODO:
## -1) lvsplot version 2.0; after some consultation with David and inspired by a problem of Hao Ran -- see playingaround_lvsplot2
## 0) Detection probabilities?
## 1) How to allow for multiple sets of LVs, 
## 2) Allow for phylogenetic correlations on response coefficients between responses
## 3) Reduce rank response coefficients to get constrained ordination? HARD!!! 
## 5) allow for multiple chains, but don't check convergence on the LVs and their loadings. Also cannot combined chains for LV and loadings unless you post process them, which is annoying.
## 6) Include OFAL priors as a way to select and make sparse loadings? See massaging-sparsepriors.R, although it currently does not appear to work very well?! Too harsh after talking to Florian Hartig about this at vISEC2020
## 8) Scrap ICs formally and do not allow their calculation at all by setting calc.ics = FALSE everywhere. Same thing with all forms of log-likelihood?
##############
# library(R2jags); 
# library(mvtnorm); 
# library(mvabund); 
# library(coda); 
# library(reshape2); 
# library(MASS)
# library(fishMod)
# library(abind)


boral <- function(y, ...)
     UseMethod("boral")


## Model is g(mu_{ij}) = row + beta0 + LV_i*theta_j + X_i*beta_j
boral.default <- function(y, X = NULL, formula.X = NULL, X.ind = NULL, traits = NULL, which.traits = NULL, family, trial.size = 1,  
     lv.control = list(num.lv = 0, type = "independent", distmat = NULL), row.eff = "none", row.ids = NULL, ranef.ids = NULL, 
     offset = NULL, save.model = FALSE, calc.ics = FALSE, 
     mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = NULL), 
     prior.control = list(type = c("normal","normal","normal","uniform"), hypparams = c(10, 10, 10, 30), ssvs.index = -1, ssvs.g = 1e-6, ssvs.traitsindex = -1), 
     do.fit = TRUE, model.name = NULL, num.lv = NULL, ...) {
    
     ##-----------------------------
     ## DO CHECKS 
     ##-----------------------------
     if(is.null(dim(y))) 
          y <- as.matrix(y) 
     if(!is.null(traits) & is.null(dim(traits)))
          traits <- as.matrix(traits) 
     if(!is.null(formula.X)) {
          formula.X <- check_X_formula(formula.X = formula.X, X = X)          
          X <- model.matrix(object = formula.X, data = as.data.frame(X))[,-1] ## Remove intercept
          message("Note X has been overrided based on the form X = model.matrix(formula.X, X).")
          }
     if(!is.null(X)) { 
          if(!is.matrix(X)) 
               X <- as.matrix(X)
          if(any(apply(X,2,function(x) all(x == 1)))) 
               stop("No intercept column should be included in X.") 
          }
     num.X <- ifelse(!is.null(X), ncol(X), 0) 
     num.traits <- ifelse(!is.null(traits), ncol(traits), 0)
     
     prior.control <- fillin_prior_control(x = prior.control)
     check_prior_control(prior.control = prior.control)
     
     check_traits(traits = traits, y = y) 

     lv.control <- check_lv_control(num.lv = num.lv, lv.control = lv.control)
     num.lv <- lv.control$num.lv

     family <- check_family(family = family, y = y)
     
     row.eff <- match.arg(row.eff, choices = c("none", "fixed", "random"))
     if(row.eff != "none" & is.null(row.ids)) {
          row.ids <- matrix(1:nrow(y), ncol = 1)
          colnames(row.ids) <- "ID1"
          message("row.ids assumed to be a matrix with one column and elements 1,2,...nrow(y) i.e., a row-specific intercept.")
          row.ids <- check_row_ids(row.ids = row.ids, y = y)	
          }

     ranef.ids <- check_ranef_ids(ranef.ids = ranef.ids, y = y)

     check_offset(offset = offset, y = y) 

     check_which_traits(num.traits = num.traits, which.traits = which.traits, traits = traits, y = y, num.X = num.X)

     X.ind <- check_X_ind(X.ind = X.ind, p = ncol(y), num.X = num.X, prior.control = prior.control)
                    

     if(num.traits > 0 & any(prior.control$ssvs.index != -1)) {
          message("If traits are supplied, then prior.control$ssvs.index is ignored and prior.control$ssvs.traitsindex is looked at. That is, boral assumes a fourth corner model is being fitted and so SSVS cannot be be applied to X.") 
          prior.control$ssvs.index <- -1
          } 	
     if(!(length(prior.control$ssvs.index) %in% c(1, ncol(X)))) 
          stop("Number of elements in prior.control$ssvs.index must either be one or the # of columns in X.")
     if(length(prior.control$ssvs.index) == 1 & num.X > 0) 
          prior.control$ssvs.index <- rep(prior.control$ssvs.index, ncol(X))
     if(any(prior.control$ssvs.index < -1)) 
          stop("Elements of prior.control$ssvs.index can only take values in -1, 0, or any positive integer; please see help file for more information.")		
     if(num.traits > 0) {
          if(!is.list(prior.control$ssvs.traitsindex)) {
               prior.control$ssvs.traitsindex <- vector("list",num.X+1)
               for(k in 1:(num.X+1)) 
                    prior.control$ssvs.traitsindex[[k]] <- rep(-1,length(which.traits[[k]])) 
               }
          if(is.list(prior.control$ssvs.traitsindex)) {
               check_ssvstraits(prior.control$ssvs.traitsindex, which.traits)
               }
          }
     
     complete_trial_size <- check_trial_size(family = family, trial.size = trial.size, create.complete.trial.size = TRUE, y = y)
     
     if(all(family != "ordinal")) 
          num.ord.levels <- 0
     if(any(family == "ordinal")) 
          num.ord.levels <- max(y[, family == "ordinal"])
#     if(all(family != "multinom")) 
#           { 
#           num.multinom.levels <- 0
#           index_multinom_cols <- NULL 
#           }
# 	if(any(family == "multinom")) { 
# 		num.multinom.levels <- apply(y[, family == "multinom"], 2, max)
# 		index_multinom_cols <- which(family == "multinom") 
# 		}


     mcmc.control <- fillin_mcmc_control(x = mcmc.control)

    
     ##-----------------------------
     ## MAKE JAGS SCRIPT
     ##-----------------------------
     n <- nrow(y)
     p <- ncol(y)
     n.chains <- 1; ## Run one chain only to avoid arbitrary rotation problems
     if(num.lv > 0) 
          make.jagsboralmodel(family = family, num.X = num.X, X.ind = X.ind, num.traits = num.traits, which.traits = which.traits, lv.control = lv.control, 
               row.eff = row.eff, row.ids = row.ids, ranef.ids = ranef.ids, offset = offset, trial.size = complete_trial_size, 
               n = n, p = p, model.name = model.name, prior.control = prior.control)
     if(num.lv == 0)  
          make.jagsboralnullmodel(family = family, num.X = num.X, X.ind = X.ind, num.traits = num.traits, which.traits = which.traits, 
               row.eff = row.eff, row.ids = row.ids, ranef.ids = ranef.ids, offset = offset, trial.size = complete_trial_size, 
               n = n, p = p, model.name = model.name, prior.control = prior.control)
     if(!do.fit) { 
          message("JAGS model file created only. Thank you, come again!")
          return() 
          }
     
               
     ##-----------------------------
     ## FORM DATA
     ##-----------------------------
     jags_data <- list("y", "n", "p", "num.lv", "num.X", "num.traits", "num.ord.levels") #num.multinom.levels
     if(lv.control$type != "independent") {
          zero.lvs <- rep(1,n)        
          distmat <- lv.control$distmat
          lv.control$distmat <- NULL
          jags_data <- c(jags_data, "zero.lvs", "distmat")
          }
     if(num.X > 0) {
          if(is.null(X.ind))
               jags_data <- c(jags_data, "X")
          if(!is.null(X.ind)) {
               X.ind[X.ind == 1] <- 1e6 ## Modifying this so as to make use of indicator functions in JAGS
               jags_data <- c(jags_data, "X", "X.ind")
               }
          }
     if(num.traits > 0) 
          jags_data <- c(jags_data, "traits")
     if(any(family == "ordinal")) { 
          ones <- matrix(1, n, p)
          jags_data <- c(jags_data, "ones") 
          }
     if(row.eff != "none") { 
          n.ID <- apply(row.ids, 2, function(x) length(unique(x)))
          jags_data <- c(jags_data, "row.ids", "n.ID") 
          }
     if(!is.null(ranef.ids)) { 
          n.ranefID <- apply(ranef.ids, 2, function(x) length(unique(x)))
          jags_data <- c(jags_data, "ranef.ids", "n.ranefID") 
          }
     if(!is.null(offset)) 
          jags_data <- c(jags_data, "offset")
            
    
     ##-----------------------------
     ## FORM PARAMETERS
     ##-----------------------------
     jags_params <- c("lv.coefs")
     if(num.lv > 0) 
          jags_params <- c(jags_params, "lvs")
     if(lv.control$type != "independent") 
          jags_params <- c(jags_params, "lv.covparams")
     if(row.eff != "none") 
          jags_params <- c(jags_params, paste0("row.coefs.ID",1:ncol(row.ids)))
     if(row.eff == "random") 
          jags_params <- c(jags_params, paste0("row.sigma.ID",1:ncol(row.ids)))
     if(!is.null(ranef.ids)) 
          jags_params <- c(jags_params, paste0("ranef.coefs.ID",1:ncol(ranef.ids)), paste0("ranef.sigma.ID",1:ncol(ranef.ids)))
     if(num.X > 0 & any(family != "multinom")) 
          jags_params <- c(jags_params, "X.coefs")
     #if(num.X > 0 & any(family == "multinom")) jags_params <- c(jags_params, "X.multinom.params")
     if(num.traits > 0) 
          jags_params <- c(jags_params, "traits.int", "traits.coefs", "trait.sigma")
     if(any(family == "tweedie")) 
          jags_params <- c(jags_params, "powerparam")
     if(any(family == "ordinal")) 
          jags_params <- c(jags_params, "cutoffs", "ordinal.sigma")
     if(any(prior.control$ssvs.index == 0)) 
          jags_params <- c(jags_params, paste0("ssvs.indX", which(prior.control$ssvs.index == 0)))
     if(any(prior.control$ssvs.index > 0)) 
          jags_params <- c(jags_params, paste0("ssvs.gp", unique(prior.control$ssvs.index[prior.control$ssvs.index > 0])))
     if(any(unlist(prior.control$ssvs.traitsindex) == 0)) 
          jags_params <- c(jags_params, paste0("ssvs.traitscoefs", rep(1:(ncol(X)+1), times = sapply(prior.control$ssvs.traitsindex, function(x) sum(x == 0))), unlist(sapply(prior.control$ssvs.traitsindex, function(x) which(x == 0)))))

          jags.inits <- function() {
          initial.list <- list()
          if(any(family %in% "tweedie")) 
               initial.list$numfish = matrix(1, n, sum(family=="tweedie"))
          if(any(family %in% "ordinal")) 
               initial.list$cutoffs0 <- seq(-1, 1, length = num.ord.levels - 1)
          
     #         if(all(family %in% "bernoulli")) {
     #             Tau <- rWishart(1,p+1,diag(p))[,,1]
     #             Sigma <- solve(Tau)
     #             Z <- abs(t(rmvnorm(n,rep(0,p),Sigma)))
     #             Z <- ifelse(as.matrix(y), Z, -1 * Z)
     #             initial.list$Z <- Z 
     #             }
                    
          return(initial.list)
          }

     if(!is.null(mcmc.control$seed))
          set.seed(mcmc.control$seed)
     actual.filename <- model.name
     if(is.null(actual.filename)) 
          actual.filename <- "jagsboralmodel.txt"


     ##-----------------------------
     ## THE FIT
     ##-----------------------------
     jagsfit <- try(suppressWarnings(
          jags(data = jags_data, inits = jags.inits, parameters.to.save = jags_params, model.file = actual.filename, n.iter = mcmc.control$n.iteration, n.burnin = mcmc.control$n.burnin, n.chains = 1, n.thin = mcmc.control$n.thin)), 
          silent=TRUE)

     #print(jagsfit)
     if(inherits(jagsfit,"try-error")) {
          lookfornegbinerror <- grep("Slicer stuck at value with infinite density", jagsfit[[1]])
          
          if(any(family %in% c("ztnegative.binomial","negative.binomial")) & length(lookfornegbinerror) == 1) { 
               message("MCMC sampling through JAGS failed. This is likely due to the prior on the dispersion (size) parameter of the negative binomial distribution been too uninformative (see below). For instance, if the error message refers to lv.coefs[25,4], then this means the MCMC sampling ran into issues for column (response) 25 in y.\n
               Please consider the following solutions: 1) remove very rare responses like singletons and doubletons, as they can potentially issus for MCMC sampling, and do not provide much information about the community in general, 2) adopt a tougher prior for the overdispersion parameter e.g., keeping with a uniform prior but reducing hypparams[4], or using a half-cauchy prior, 3) consider switching to a Poisson family for responses which don't appear to be overdispersed;")
               print(jagsfit) 
               }
          else {
               message("MCMC fitting through JAGS failed:")
               print(jagsfit) 
               }

          message("boral fit failed...Exiting. Sorry!") 
          return()
          }
    
# 	return(jagsfit)
     ##-----------------------------
     ## FORMAT INTO BIG MATRIX
     ##-----------------------------
     fit.mcmcBase <- jagsfit$BUGSoutput
     fit.mcmc <- mcmc(fit.mcmcBase$sims.matrix, start = 1, thin = mcmc.control$n.thin) 
     if(n.chains == 1)
          combined_fit_mcmc <- fit.mcmc
     # 	if(n.chains > 1) {
     # 		get.rhats <- process.rhats(sims.matrix = fit.mcmcBase$sims.matrix)

     # 		exceed.rhatcutoff <- sum(sapply(get.rhats, function(x) sum(x > rhat.cutoff)))
     # 		message("There were", exceed.rhatcutoff, "(", 100*exceed.rhatcutoff/sum(sapply(get.rhats,length)), "%) parameters whose Rhat exceeded the prespecified cutoff of", rhat.cutoff, "\n")		
     # 		}	
     rm(fit.mcmc, fit.mcmcBase)

#   	## For any multinomial columns, set the corresponding rows in X.coefs to zero
# 	if(any(family == "multinom") & num.X > 0) {
# 		for(k in index_multinom_cols) {
# 			sel.multinom.col <- grep(paste("X.coefs\\[", k, ",+", sep = ""), colnames(combined_fit_mcmc))
# 			combined_fit_mcmc[, sel.multinom.col] <- 0 }
# 		}

    
     ##-----------------------------
     ## BLING THE OUTPUT
     ##-----------------------------
     mcmc_names <- colnames(combined_fit_mcmc)
     
     if(is.null(colnames(y))) 
          colnames(y) <- 1:ncol(y)
     if(is.null(rownames(y))) 
          rownames(y) <- 1:nrow(y)
     if(num.X > 0) { 
          if(is.null(colnames(X))) colnames(X) <- 1:ncol(X)
          if(is.null(rownames(X))) rownames(X) <- 1:nrow(X) 
          }
     if(num.traits > 0) { 
          if(is.null(colnames(traits))) colnames(traits) <- 1:ncol(traits)
          if(is.null(rownames(traits))) rownames(traits) <- 1:nrow(traits) 
          }


     lv_coefs_arr <- abind(
          matrix(apply(combined_fit_mcmc[, grep("lv.coefs", mcmc_names)], 2, median), nrow = p),
          matrix(apply(combined_fit_mcmc[, grep("lv.coefs", mcmc_names)], 2, mean), nrow = p),
          matrix(apply(combined_fit_mcmc[, grep("lv.coefs", mcmc_names)], 2, IQR), nrow = p),
          matrix(apply(combined_fit_mcmc[, grep("lv.coefs", mcmc_names)], 2, sd), nrow = p), 
          along = 3)
     
     out_fit <- list()

     if(num.lv > 0) {
          lv_arr <- abind(
                    matrix(apply(combined_fit_mcmc[, grep("lvs", mcmc_names)], 2, median), nrow = n),
                    matrix(apply(combined_fit_mcmc[, grep("lvs", mcmc_names)], 2, mean), nrow = n),
                    matrix(apply(combined_fit_mcmc[, grep("lvs", mcmc_names)], 2, IQR), nrow = n),
                    matrix(apply(combined_fit_mcmc[, grep("lvs", mcmc_names)], 2, sd), nrow = n),
                    along = 3)
          dimnames(lv_arr) <- list(rows = rownames(y), lv = paste0("lv", 1:num.lv), type = c("median","mean","iqr","sd"))
          
          out_fit$lv.median <- as.matrix(lv_arr[,,1]); out_fit$lv.mean <- as.matrix(lv_arr[,,2]); 
          out_fit$lv.iqr <- as.matrix(lv_arr[,,3]); out_fit$lv.sd <- as.matrix(lv_arr[,,4])
          rm(lv_arr)
                    
          if(dim(lv_coefs_arr)[2] == (num.lv+2)) 
                    dimnames(lv_coefs_arr) <- list(resp = colnames(y), coefficients = c("beta0", paste0("theta", 1:num.lv), "Dispersion"), type = c("median","mean","iqr","sd"))
          if(dim(lv_coefs_arr)[2] == (num.lv+1)) 
                    dimnames(lv_coefs_arr) <- list(resp = colnames(y), coefficients = c("beta0", paste0("theta", 1:num.lv)), type = c("median","mean","iqr","sd"))
          if(lv.control$type != "independent") {
               lv_params_arr <- cbind(
                    apply(combined_fit_mcmc[, grep("lv.covparams", mcmc_names),drop=FALSE], 2, median), 
                    apply(combined_fit_mcmc[, grep("lv.covparams", mcmc_names),drop=FALSE], 2, mean),
                    apply(combined_fit_mcmc[, grep("lv.covparams", mcmc_names),drop=FALSE], 2, IQR),
                    apply(combined_fit_mcmc[, grep("lv.covparams", mcmc_names),drop=FALSE], 2, sd))
               if(nrow(lv_params_arr) == 1) 
                    rownames(lv_params_arr) <- c("spatialscale (tau1)")
               if(nrow(lv_params_arr) == 2) 
                    rownames(lv_params_arr) <- c("spatialscale (tau1)", "spatialpower (tau2)")
               colnames(lv_params_arr) <- c("median","mean","iqr","sd")

               out_fit$lv.covparams.median <- lv_params_arr[,1]
               out_fit$lv.covparams.mean <- lv_params_arr[,2]
               out_fit$lv.covparams.iqr <- lv_params_arr[,3]
               out_fit$lv.covparams.sd <- lv_params_arr[,4]
               rm(lv_params_arr)
               }
          }
            
            
     if(num.lv == 0) {
          if(dim(lv_coefs_arr)[2] == 2) 
               dimnames(lv_coefs_arr) <- list(resp = colnames(y), coefficients = c("beta0", "Dispersion"), type = c("median","mean","iqr","sd"))
          if(dim(lv_coefs_arr)[2] == 1) 
               dimnames(lv_coefs_arr) <- list(resp = colnames(y), coefficients = c("beta0"), type = c("median","mean","iqr","sd"))
          }	
     out_fit$lv.coefs.median <- lv_coefs_arr[,,1]
     out_fit$lv.coefs.mean <- lv_coefs_arr[,,2]
     out_fit$lv.coefs.iqr <- lv_coefs_arr[,,3]
     out_fit$lv.coefs.sd <- lv_coefs_arr[,,4]
     rm(lv_coefs_arr)
     if(length(out_fit$lv.coefs.median) == p) {
          out_fit$lv.coefs.median <- matrix(out_fit$lv.coefs.median, ncol = 1)
          out_fit$lv.coefs.mean <- matrix(out_fit$lv.coefs.mean, ncol = 1)
          out_fit$lv.coefs.iqr <- matrix(out_fit$lv.coefs.iqr, ncol = 1)
          out_fit$lv.coefs.sd <- matrix(out_fit$lv.coefs.sd, ncol = 1)
          rownames(out_fit$lv.coefs.median) <- rownames(out_fit$lv.coefs.mean) <- rownames(out_fit$lv.coefs.iqr) <- rownames(out_fit$lv.coefs.sd) <- colnames(y)
          colnames(out_fit$lv.coefs.median) <- colnames(out_fit$lv.coefs.mean) <- colnames(out_fit$lv.coefs.iqr) <- colnames(out_fit$lv.coefs.sd) <- "beta0"
          }

     if(row.eff != "none") {
          out_fit$row.coefs <- vector("list", ncol(row.ids))
          names(out_fit$row.coefs) <- colnames(row.ids)
          for(k in 1:ncol(row.ids)) {
               row_coefs_arr <- cbind(
                    apply(combined_fit_mcmc[, grep(paste0("row.coefs.ID",k,"\\["), mcmc_names),drop=FALSE], 2, median), 
                    apply(combined_fit_mcmc[, grep(paste0("row.coefs.ID",k,"\\["), mcmc_names),drop=FALSE], 2, mean), 
                    apply(combined_fit_mcmc[, grep(paste0("row.coefs.ID",k,"\\["), mcmc_names),drop=FALSE], 2, IQR), 
                    apply(combined_fit_mcmc[, grep(paste0("row.coefs.ID",k,"\\["), mcmc_names),drop=FALSE], 2, sd))
               rownames(row_coefs_arr) <- 1:n.ID[k]
               colnames(row_coefs_arr) <- c("median","mean","iqr","sd")
                         
               out_fit$row.coefs[[k]] <- list(median = row_coefs_arr[,1], mean = row_coefs_arr[,2], iqr = row_coefs_arr[,3], sd = row_coefs_arr[,4])
               rm(row_coefs_arr)
               }
     
               if(row.eff == "random") {
                    out_fit$row.sigma <- vector("list", ncol(row.ids))
                    names(out_fit$row.sigma) <- colnames(row.ids)
                    for(k in 1:ncol(row.ids)) {
                         row_sigma_vec <- c(
                         median(combined_fit_mcmc[, grep(paste0("row.sigma.ID",k,"$"), mcmc_names)]),
                         mean(combined_fit_mcmc[, grep(paste0("row.sigma.ID",k,"$"), mcmc_names)]),
                         IQR(combined_fit_mcmc[, grep(paste0("row.sigma.ID",k,"$"), mcmc_names)]),
                         sd(combined_fit_mcmc[, grep(paste0("row.sigma.ID",k,"$"), mcmc_names)]))
                         names(row_sigma_vec) <- c("median","mean","iqr","sd")
                              
                         out_fit$row.sigma[[k]] <- row_sigma_vec
                         }
                    }
          }
          
          
     if(!is.null(ranef.ids)) {
          out_fit$ranef.coefs.median <- out_fit$ranef.coefs.mean <- out_fit$ranef.coefs.iqr <- out_fit$ranef.coefs.sd <- vector("list", ncol(ranef.ids))
          names(out_fit$ranef.coefs.median) <- names(out_fit$ranef.coefs.mean) <- 
               names(out_fit$ranef.coefs.iqr)  <- names(out_fit$ranef.coefs.sd) <- colnames(ranef.ids)
          for(k0 in 1:ncol(ranef.ids)) {
               out_fit$ranef.coefs.median[[k0]] <- matrix(apply(combined_fit_mcmc[, grep(paste0("ranef.coefs.ID",k0,"\\["), mcmc_names),drop=FALSE], 2, median), nrow = p)
               out_fit$ranef.coefs.mean[[k0]] <- matrix(apply(combined_fit_mcmc[, grep(paste0("ranef.coefs.ID",k0,"\\["), mcmc_names),drop=FALSE], 2, mean), nrow = p)
               out_fit$ranef.coefs.iqr[[k0]] <- matrix(apply(combined_fit_mcmc[, grep(paste0("ranef.coefs.ID",k0,"\\["), mcmc_names),drop=FALSE], 2, IQR), nrow = p)
               out_fit$ranef.coefs.sd[[k0]] <- matrix(apply(combined_fit_mcmc[, grep(paste0("ranef.coefs.ID",k0,"\\["), mcmc_names),drop=FALSE], 2, sd), nrow = p)      
               colnames(out_fit$ranef.coefs.median[[k0]]) <- colnames(out_fit$ranef.coefs.mean[[k0]]) <- 
                    colnames(out_fit$ranef.coefs.iqr[[k0]]) <- colnames(out_fit$ranef.coefs.sd[[k0]]) <- 1:n.ranefID[k0]
               rownames(out_fit$ranef.coefs.median[[k0]]) <- rownames(out_fit$ranef.coefs.mean[[k0]]) <- 
                    rownames(out_fit$ranef.coefs.iqr[[k0]]) <- rownames(out_fit$ranef.coefs.sd[[k0]]) <- colnames(y)                         
               }
          
          out_fit$ranef.sigma.median <- out_fit$ranef.sigma.mean <- out_fit$ranef.sigma.iqr <- out_fit$ranef.sigma.sd <- matrix(NA, nrow = p, ncol = ncol(ranef.ids))
          for(k0 in 1:ncol(ranef.ids)) {
               out_fit$ranef.sigma.median[,k0] <- apply(combined_fit_mcmc[, grep(paste0("ranef.sigma.ID",k0,"\\["), mcmc_names)], 2, median)
               out_fit$ranef.sigma.mean[,k0] <- apply(combined_fit_mcmc[, grep(paste0("ranef.sigma.ID",k0,"\\["), mcmc_names)], 2, mean)
               out_fit$ranef.sigma.iqr[,k0] <- apply(combined_fit_mcmc[, grep(paste0("ranef.sigma.ID",k0,"\\["), mcmc_names)], 2, IQR)
               out_fit$ranef.sigma.sd[,k0] <- apply(combined_fit_mcmc[, grep(paste0("ranef.sigma.ID",k0,"\\["), mcmc_names)], 2, sd)
               }
          colnames(out_fit$ranef.sigma.median) <- colnames(out_fit$ranef.sigma.mean) <- 
               colnames(out_fit$ranef.sigma.iqr) <- colnames(out_fit$ranef.sigma.sd) <- colnames(ranef.ids)
          rownames(out_fit$ranef.sigma.median) <- rownames(out_fit$ranef.sigma.mean) <- 
               rownames(out_fit$ranef.sigma.iqr) <- rownames(out_fit$ranef.sigma.sd) <- colnames(y)
               
          }
     
            
     if(num.X > 0) {
          X_coefs_arr <- abind(
               matrix(apply(combined_fit_mcmc[, grep("X.coefs", mcmc_names),drop=FALSE], 2, median), nrow = p),
               matrix(apply(combined_fit_mcmc[, grep("X.coefs", mcmc_names),drop=FALSE], 2, mean), nrow = p),
               matrix(apply(combined_fit_mcmc[, grep("X.coefs", mcmc_names),drop=FALSE], 2, IQR), nrow = p),
               matrix(apply(combined_fit_mcmc[, grep("X.coefs", mcmc_names),drop=FALSE], 2, sd), nrow = p),
               along = 3)
          dimnames(X_coefs_arr) <- list(resp = colnames(y), coefficients = colnames(X), type = c("median","mean","iqr","sd"))

          out_fit$X.coefs.median <- X_coefs_arr[,,1]
          out_fit$X.coefs.mean <- X_coefs_arr[,,2]
          out_fit$X.coefs.iqr <- X_coefs_arr[,,3]
          out_fit$X.coefs.sd <- X_coefs_arr[,,4]
          rm(X_coefs_arr)
          if(length(out_fit$X.coefs.median) == p) {
               out_fit$X.coefs.median <- matrix(out_fit$X.coefs.median, ncol = 1)
               out_fit$X.coefs.mean <- matrix(out_fit$X.coefs.mean, ncol = 1)
               out_fit$X.coefs.iqr <- matrix(out_fit$X.coefs.iqr, ncol = 1)
               out_fit$X.coefs.sd <- matrix(out_fit$X.coefs.sd, ncol = 1)
               rownames(out_fit$X.coefs.median) <- rownames(out_fit$X.coefs.mean) <- rownames(out_fit$X.coefs.iqr) <- rownames(out_fit$X.coefs.sd) <- colnames(y)
               colnames(out_fit$X.coefs.median) <- colnames(out_fit$X.coefs.mean) <- colnames(out_fit$X.coefs.iqr) <- colnames(out_fit$X.coefs.sd) <- colnames(X)
               }		
                              
          if(any(prior.control$ssvs.index == 0)) { ## You should not be able to enter this loop if num.traits > 0!
               ssvs_indcoefs_arr <- array(NA, dim = c(p, num.X, 2))
               for(k1 in 1:num.X) {
                    find.Xvars <- grep(paste0("ssvs.indX",k1,"\\["), mcmc_names)
                    if(length(find.Xvars) > 0) {
                         ssvs_indcoefs_arr[,k1,1] <- colMeans(combined_fit_mcmc[,find.Xvars])
                         ssvs_indcoefs_arr[,k1,2] <- apply(combined_fit_mcmc[,find.Xvars], 2, sd)
                         }
                    }
               dimnames(ssvs_indcoefs_arr) <- list(resp = colnames(y), coefficients = colnames(X), type = c("mean","sd"))
               
               out_fit$ssvs.indcoefs.mean <- ssvs_indcoefs_arr[,,1]
               out_fit$ssvs.indcoefs.sd <- ssvs_indcoefs_arr[,,2]
               rm(ssvs_indcoefs_arr)
               }
          if(any(prior.control$ssvs.index > 0)) {
               ssvs_gpcoefs_arr <- cbind(
                    apply(combined_fit_mcmc[, grep("ssvs.gp", mcmc_names),drop=FALSE], 2, mean),
                    apply(combined_fit_mcmc[, grep("ssvs.gp", mcmc_names),drop=FALSE], 2, sd))
                    rownames(ssvs_gpcoefs_arr) <- paste0("Gp", unique(prior.control$ssvs.index[prior.control$ssvs.index > 0]))	
                    colnames(ssvs_gpcoefs_arr) <- c("mean","sd")

                    out_fit$ssvs.gpcoefs.mean <- ssvs_gpcoefs_arr[,1]
                    out_fit$ssvs.gpcoefs.sd <- ssvs_gpcoefs_arr[,2]
                    rm(ssvs_gpcoefs_arr)
               }
          
          ## Convert to matrix form
          if(any(unlist(prior.control$ssvs.traitsindex) == 0)) { 
               ssvs_traitscoefs_arr <- array(NA, dim = c(num.X+1, num.traits, 2))
               dimnames(ssvs_traitscoefs_arr) <- list(X.coefficients = c("beta0",colnames(X)), traits.coefficients =  colnames(traits), type = c("mean","sd"))
               for(k1 in 1:(num.X+1)) { for(k2 in 1:num.traits) {
                    find.Xvars <- grep(paste0("ssvs.traitscoefs",k1,k2,"$"), mcmc_names)
                    if(length(find.Xvars) == 1) {
                         ssvs_traitscoefs_arr[k1,k2,1] <- mean(combined_fit_mcmc[,find.Xvars])
                         ssvs_traitscoefs_arr[k1,k2,2] <- sd(combined_fit_mcmc[,find.Xvars])
                         }
                    } }
               
               out_fit$ssvs.traitscoefs.mean <- ssvs_traitscoefs_arr[,,1]
               out_fit$ssvs.traitscoefs.sd <- ssvs_traitscoefs_arr[,,2]
               rm(ssvs_traitscoefs_arr)
               }
          }
            
            
     if(num.traits > 0) {
          traitscoefs_arr <- array(0, dim = c(num.X+1, num.traits + 2, 4))
          traitscoefs_arr[,,1] <- cbind(
               apply(combined_fit_mcmc[, grep("traits.int", mcmc_names),drop=FALSE], 2, median), 
               matrix(apply(combined_fit_mcmc[, grep("traits.coefs", mcmc_names),drop=FALSE], 2, median), nrow = num.X+1), 
               apply(combined_fit_mcmc[, grep("trait.sigma", mcmc_names),drop=FALSE], 2, median))
          traitscoefs_arr[,,2] <- cbind(
               apply(combined_fit_mcmc[, grep("traits.int", mcmc_names),drop=FALSE], 2, mean), 
               matrix(apply(combined_fit_mcmc[, grep("traits.coefs", mcmc_names),drop=FALSE], 2, mean), nrow = num.X+1), 
               apply(combined_fit_mcmc[, grep("trait.sigma", mcmc_names),drop=FALSE], 2, mean))
          traitscoefs_arr[,,3] <- cbind(
               apply(combined_fit_mcmc[, grep("traits.int", mcmc_names),drop=FALSE], 2, IQR), 
               matrix(apply(combined_fit_mcmc[, grep("traits.coefs", mcmc_names),drop=FALSE], 2, IQR), nrow = num.X+1), 
               apply(combined_fit_mcmc[, grep("trait.sigma", mcmc_names),drop=FALSE], 2, IQR))
          traitscoefs_arr[,,4] <- cbind(
               apply(combined_fit_mcmc[, grep("traits.int", mcmc_names),drop=FALSE], 2, sd), 
               matrix(apply(combined_fit_mcmc[, grep("traits.coefs", mcmc_names),drop=FALSE], 2, sd), nrow = num.X+1), 
               apply(combined_fit_mcmc[, grep("trait.sigma", mcmc_names),drop=FALSE], 2, sd))
          dimnames(traitscoefs_arr) <- list(X.coefficients = c("beta0",colnames(X)), traits.coefficients =  c("kappa0",colnames(traits),"sigma"), type = c("median","mean","iqr","sd"))
                    
          out_fit$traits.coefs.median <- traitscoefs_arr[,,1]
          out_fit$traits.coefs.mean <- traitscoefs_arr[,,2]
          out_fit$traits.coefs.iqr <- traitscoefs_arr[,,3]
          out_fit$traits.coefs.sd <- traitscoefs_arr[,,4]
          rm(traitscoefs_arr)
          }

#   	if(num.X > 0 & any(family == "multinom")) {
#   		out_fit$X.multinom.coefs.median <- array(apply(combined_fit_mcmc[,grep("X.multinom.params", mcmc_names)],2,median),dim=c(length(index_multinom_cols),num.X,num.multinom.levels))
#   		out_fit$X.multinom.coefs.iqr <- array(apply(combined_fit_mcmc[,grep("X.multinom.params", mcmc_names)],2,IQR),dim=c(length(index_multinom_cols),num.X,num.multinom.levels))
#   		out_fit$X.multinom.coefs.mean <- array(apply(combined_fit_mcmc[,grep("X.multinom.params", mcmc_names)],2,mean),dim=c(length(index_multinom_cols),num.X,num.multinom.levels))
#   		out_fit$X.multinom.coefs.sd <- array(apply(combined_fit_mcmc[,grep("X.multinom.params", mcmc_names)],2,sd),dim=c(length(index_multinom_cols),num.X,num.multinom.levels))
#   
#   		dimnames(out_fit$X.multinom.coefs.median) <- dimnames(out_fit$X.multinom.coefs.iqr) <- dimnames(out_fit$X.multinom.coefs.mean) <- dimnames(out_fit$X.multinom.coefs.sd) <- list("1" = index_multinom_cols, "2" = colnames(X), "level" = 1:num.multinom.levels)
#   		}

     if(any(family == "ordinal")) {
          cutoffs_arr <- cbind(
               apply(combined_fit_mcmc[, grep("cutoffs", mcmc_names),drop=FALSE], 2, median),
               apply(combined_fit_mcmc[, grep("cutoffs", mcmc_names),drop=FALSE], 2, mean),
               apply(combined_fit_mcmc[, grep("cutoffs", mcmc_names),drop=FALSE], 2, IQR),
               apply(combined_fit_mcmc[, grep("cutoffs", mcmc_names),drop=FALSE], 2, sd))
          rownames(cutoffs_arr) <- paste0(1:(num.ord.levels - 1), "|", 2:num.ord.levels) 
          colnames(cutoffs_arr) <- c("median","mean","iqr","sd")
                    
          out_fit$cutoffs.median <- cutoffs_arr[,1]
          out_fit$cutoffs.mean <- cutoffs_arr[,2]
          out_fit$cutoffs.iqr <- cutoffs_arr[,3]
          out_fit$cutoffs.sd <- cutoffs_arr[,4]
          rm(cutoffs_arr)

          ## If there are traits, then ordinal random intercept is either zero (if there is only 1 ordinal column, or has the trait.sigma (if there are >1 ordinal columns)
          if(sum(family == "ordinal") > 1 & is.null(traits)) { 
               ordinal_sigma_vec <- c(
                    median(combined_fit_mcmc[, grep("ordinal.sigma", mcmc_names)]),
                    mean(combined_fit_mcmc[, grep("ordinal.sigma", mcmc_names)]),
                    IQR(combined_fit_mcmc[, grep("ordinal.sigma", mcmc_names)]),
                    sd(combined_fit_mcmc[, grep("ordinal.sigma", mcmc_names)]))
               names(ordinal_sigma_vec) <- c("median","mean","iqr","sd")

               out_fit$ordinal.sigma.median <- ordinal_sigma_vec[1]
               out_fit$ordinal.sigma.mean <- ordinal_sigma_vec[2]
               out_fit$ordinal.sigma.iqr <- ordinal_sigma_vec[3]
               out_fit$ordinal.sigma.sd <- ordinal_sigma_vec[4]            
               }
          }

     if(any(family == "tweedie")) {
          powerparam_vec <- c(
               median(combined_fit_mcmc[, grep("powerparam", mcmc_names)]),
               mean(combined_fit_mcmc[, grep("powerparam", mcmc_names)]),
               IQR(combined_fit_mcmc[, grep("powerparam", mcmc_names)]),
               sd(combined_fit_mcmc[, grep("powerparam", mcmc_names)]))
          names(powerparam_vec) <- c("median","mean","iqr","sd")
                    
          out_fit$powerparam.median <- powerparam_vec[1]
          out_fit$powerparam.mean <- powerparam_vec[2]
          out_fit$powerparam.iqr <- powerparam_vec[3]
          out_fit$powerparam.sd <- powerparam_vec[4]
          }
     rm(list = ls(pattern = ".arr"))	
            
            
            
     ##-----------------------------
     ## HPD INTERVALS, ICS, ETC...
     ##-----------------------------
     get.hpds <- get.hpdintervals(y, X = X, traits = traits, row.ids = row.ids, ranef.ids = ranef.ids, fit.mcmc = combined_fit_mcmc, lv.control = lv.control)
     out_fit$hpdintervals <- get.hpds
     if(calc.ics) {
          deprecate_warn("1.6", "boral::get.measures()", details = "All functions to calculate information criteria are being phased out (too hard to maintain)!")
          get_ics <- try(get.measures(y = y, X = X, family = family, trial.size = complete_trial_size, row.eff = row.eff, row.ids = row.ids, offset = offset, num.lv = num.lv, fit.mcmc = combined_fit_mcmc), silent = TRUE)
          if(!inherits(get_ics, "try-error")) {
               ics <- c(get.dic(jagsfit), get_ics$waic, get_ics$eaic, get_ics$ebic)
               names_ics <- c("Conditional DIC", "WAIC", "EAIC", "EBIC")
               if(get_ics$do.marglik.ics) 
                    {
                    ics <- c(ics, get_ics$aic.median, get_ics$bic.median, get_ics$median.logLik)
                    names_ics  <- c(names_ics , "AIC at post. median", "BIC at post. median", "Marginal log-likelihood at post. median")
                    }
               names(ics) <- names_ics; out_fit$ics <- ics
               }
          }
                         
     if(save.model) 
          out_fit$jags.model <- jagsfit

     out_fit$call <- match.call()
     out_fit$n <- n
     out_fit$p <- p
     out_fit$X <- X
     out_fit$formula.X <- formula.X
     X.ind[X.ind == 1e6] <- 1 ## Reset this to original format
     out_fit$X.ind <- X.ind
     out_fit$traits <- traits; 
     out_fit$y <- y
     out_fit$offset <- offset
     out_fit$row.eff <- row.eff
     out_fit$row.ids <- row.ids
     out_fit$ranef.ids <- ranef.ids

     out_fit$geweke.diag <- process_geweke(fit.mcmc = combined_fit_mcmc, y = y, X = X, traits = traits, family = family, 
          num.lv = num.lv, row.eff = row.eff, row.ids = row.ids, ranef.ids = ranef.ids, num.ord.levels = num.ord.levels, prior.control = prior.control)
     out_fit$family <- family
     out_fit$num.lv <- num.lv
     out_fit$lv.control <- lv.control
     out_fit$lv.control$distmat <- NULL ## Do not save distance matrix to save space
     out_fit$num.X <- num.X; out_fit$num.traits <- num.traits
     out_fit$which.traits <- which.traits
     out_fit$calc.ics <- calc.ics
     out_fit$trial.size <- complete_trial_size
     out_fit$prior.control <- prior.control	
     out_fit$num.ord.levels <- num.ord.levels
     out_fit$mcmc.control <- mcmc.control
     
     #out_fit$n.chain <- out_fit$n.chains; 
     class(out_fit) <- "boral"
     if(!save.model) { 
          if(file.exists(actual.filename)) 
               file.remove(actual.filename) 
          }
          
     set.seed(NULL)
     return(out_fit) 
     }
 	

print.boral <- function(x, ...) {
     message("Call:")
     print(x$call)
     message()
     message("Response matrix attributes\n \t# of rows: ", x$n, "\n\t# of columns: ", x$p) 
     message("Model attributes\n \tColumn families used: ", unique(x$family), "\n\t# of latent variables: ", x$num.lv, 
          "\n\tLatent variable covariance structure", x$lv.control$type, 
          "\n\tRow effects included (none/fixed/random)? ", x$row.eff, 
          "\n\tResponse-specific random intercepts included (TRUE/FALSE)? ", !is.null(x$ranef.ids),".\n") 

     if(any(x$family == "binomial")) 
          message("Trial sizes used (columns with binomial families): ", x$trial.size)
     
     if(any(x$family == "ordinal")) 
          message("Number of levels for ordinal data: ", x$num.ord.levels)
     
     if(x$num.X > 0) 
          message("Covariate matrix with ", x$num.X, " covariates also fitted.\n")
     
     if(x$num.traits > 0) 
          message("Trait matrix with ", x$num.traits, " traits also included.\n")
     
     if(any(x$prior.control$ssvs.index > -1) || any(unlist(x$prior.control$ssvs.traitsindex) > -1)) 
          message("SSVS has been performed on covariates.\n")
     }

 	

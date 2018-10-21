##############
## Latent Variable Ordination and Regression using MCMC 
## Ordinal data handled as propotional odds regression, with common cutoff points for all spp, but spp intercepts included as random effects; This should make sense if the same ordinal scale is applied to all species
## Multinomial regression is available, but CURRENTLY NOT ALLOWED DUE TO COMPUTATION TIME
################

## Changes from v1.6 (news files to be updated!)
## - A X.ind argument is now implement which allows users to manually force particular covariates to be excluded for particular species. The argument takes the form of a matrix with 1s and 0s. Defaults to NULL, in which case all covariates are included for all species
## - predict.boral now has the option of returning all the full array of linear predictors (rows by cols by num_mcmc)
## - Allow structured LVs in a limited manner


## 1) How to allow for multiple sets of LVs, 
## 3) Reformat output for boral so that instead of a separate element for mean,median,IQR,sd; present it as a 4 column array or vector. Without forcing all users to adopt this, the only simple way of doing this is to create a function that reformats boral output?
## 4) Allow for phylogenetic correlations on species coefficients between species
## 5) Reduce rank species coefficients to get constrained ordination? HARD!!! Phase out functions to calculate likelihoods for no LV model???
## 6) Species specific random effects! In principle, with row.ids in place this should not be too difficult!
## 7) allow for multiple chains, but don't check convergence on the LVs and their loadings. Also cannot combined chains for LV and loadings unless you post process them, which is annoying.
## 8) Include OFAL priors as a way to select and make sparse loadings? See massaging-sparsepriors.R
## 10) Scrap group SSVS?

##############
# library(R2jags); 
# library(mvtnorm); 
# library(mvabund); 
# library(coda); 
# library(MASS)
# library(fishMod)
# library(abind)
# source("boral17/R/auxilaryfunctions.R")
# source("boral17/R/calclogLfunctions.R")
# source("boral17/R/makejagsboralmodel.R")
# source("boral17/R/makejagsboralnullmodel.R")
# source("boral17/R/simdatafunctions.R")
# source("boral17/R/unseenfunctions.R")
# source("boral17/R/getmeasures.R")
#  
# 
# data(spider)
# y <- spider$abun
# X <- scale(spider$x)
# n <- nrow(y)
# p <- ncol(y)
# example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, n.thin = 1)
# fakedistmat <- as.matrix(dist(1:n))
# family = "negative.binomial"
# lv.control = list(num.lv = 2, type = "exponential", lv.covparams = 5, distmat = fakedistmat)
# row.eff = "none"
# n.burnin = 10000
# n.iteration = 40000
# X.ind = NULL; traits = NULL; which.traits = NULL; trial.size = 1
# row.eff = "none"; row.ids = NULL; offset = NULL; save.model = FALSE; calc.ics = FALSE
# mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123)
# prior.control = list(type = c("normal","normal","normal","uniform"), hypparams = c(10, 10, 10, 30), ssvs.index = -1, ssvs.g = 1e-6, ssvs.traitsindex = -1)
# do.fit = TRUE; model.name = NULL; num.lv = 0


boral <- function(y, ...) UseMethod("boral")


## Model is g(mu_{ij}) = row + beta0 + LV_i*theta_j + X_i*beta_j
boral.default <- function(y, X = NULL, X.ind = NULL, traits = NULL, which.traits = NULL, family, trial.size = 1,  
     lv.control = list(num.lv = 0, type = "independent", distmat = NULL),
     row.eff = "none", row.ids = NULL, offset = NULL, save.model = FALSE, calc.ics = FALSE, 
     mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123), 
     prior.control = list(type = c("normal","normal","normal","uniform"), hypparams = c(10, 10, 10, 30), ssvs.index = -1, ssvs.g = 1e-6, ssvs.traitsindex = -1), 
     do.fit = TRUE, model.name = NULL, num.lv = NULL, ...) {
     
    new.format <- FALSE

    
    ###############
    ## DO CHECKS 
    ###############
    if(is.null(dim(y))) 
        y <- as.matrix(y) 
    if(!is.null(traits) & is.null(dim(traits)))
        traits <- as.matrix(traits) 
    if(!is.null(X)) { 
        if(!is.matrix(X)) 
            X <- as.matrix(X)
        if(any(apply(X,2,function(x) all(x == 1)))) { stop("No intercept column should be included in X") } 
        }
	
    
    prior.control <- fillin.prior.control(x = prior.control)
    
    if(length(prior.control$hypparams) != 4 || length(prior.control$type) != 4) 
            stop("prior.control$type and prior.control$hypparams must be a vector of four elements. Please see boral help file as to what the elements correspond to.\n")
    if(!all(prior.control$type[-4] %in% c("normal","uniform","cauchy")))	
            stop("At least one of the first three elements of prior.control$type is not supported in the current version of boral...sorry!")
    if(!(prior.control$type[4] %in% c("uniform","halfcauchy","halfnormal")))	
            stop("The fourth element of prior.control$type is not supported in the current version of boral...sorry!")
            
    if(!is.null(traits)) { 
            if(!is.matrix(traits)) 
                    traits <- as.matrix(traits) 
            if(any(apply(traits,2,function(x) all(x == 1)))) 
                    stop("No intercept column should be included in traits. It will be included automatically")  
            }


    lv.control <- check.lv.control(num.lv = num.lv, lv.control = lv.control)
    num.lv <- lv.control$num.lv
    if(num.lv > 5) 
        warning("We won't stop you, but please consider if you really want more than five latent variables in the model!", immediate. = TRUE)
    
    
    if(length(family) != ncol(y) & length(family) != 1) 
        stop("Number of elements in family must either one or the # of columns in y") 
    if(length(family) == 1) 
        family <- rep(family, ncol(y))
    family <- match.arg(family, choices = c("negative.binomial", "poisson", "binomial", "normal", "lnormal", "tweedie", "ordinal", "exponential", "gamma", "beta"), several.ok = TRUE)
    if(length(family) != ncol(y))
        stop("At least one of the elements in family is not supported in current version of boral...sorry!")		
    if(any(family == "ordinal")) {
        if(sum(y[, family == "ordinal"] == 0) > 0) 
            stop("For ordinal data, please shift minimum level to 1")
        if(!is.null(traits) & (sum(family == "ordinal") == 1)) 
            message("The intercept for the single ordinal response is set to zero and not regressed traits for parameter identifiability reasons")
        }
    
    
    row.eff <- match.arg(row.eff, choices = c("none", "fixed", "random"))
    if(row.eff != "none" && is.null(row.ids)) {
        row.ids <- matrix(1:nrow(y), ncol = 1)
        colnames(row.ids) <- "ID1"
        message("row.ids assumed to be a matrix with one column and elements 1,2,...nrow(y) i.e., a row-specific intercept")
        row.ids <- check.row.ids(row.ids = row.ids, y = y)	
        }

    check.offset(offset = offset, y = y) 
    num.X <- ifelse(!is.null(X), ncol(X), 0) 
    num.traits <- ifelse(!is.null(traits), ncol(traits), 0)

    if(num.X == 0 & num.traits > 0) 
        stop("num.traits > 0 suggests traits are to be regressed against covariates X, so please supply X") 
    X.ind <- check.X.ind(X.ind = X.ind, p = ncol(y), num.X = num.X, prior.control = prior.control)
            
    if(num.traits > 0) {
        if(is.null(which.traits)) 
        stop("If traits are supplied, then please also supply which.traits to inform what traits are regressed against which covariates") 
        if((num.X+1) != length(which.traits))
            stop("which.traits should have equal to 1+length(ncol(X))") 
        if(any(sapply(which.traits,length) > num.traits)) 
            stop("Each element in the list which.traits should have at most ncol(traits) elements") 
        if(any(sapply(which.traits, function(x) any(x > ncol(traits))))) 
            stop("The values contained in the list which.traits can be takes from 1 to ncol(traits)") 
        }

    if(num.traits > 0 & any(prior.control$ssvs.index != -1)) {
            message("If traits are supplied, then prior.control$ssvs.index is ignored and prior.control$ssvs.traitsindex is looked at. That is, boral assumes a fourth corner model is being fitted and so SSVS cannot be be applied to X") 
            prior.control$ssvs.index <- -1
        } 	
    if(!(length(prior.control$ssvs.index) %in% c(1, ncol(X)))) 
            stop("Number of elements in prior.control$ssvs.index must either be one or the # of columns in X")
    if(length(prior.control$ssvs.index) == 1 & num.X > 0) 
            prior.control$ssvs.index <- rep(prior.control$ssvs.index, ncol(X))
    if(any(prior.control$ssvs.index < -1)) 
            stop("Elements of prior.control$ssvs.index can only take values in -1, 0, or any positive integer; please see help file for more information")		
    if(num.traits > 0) {
        if(!is.list(prior.control$ssvs.traitsindex)) {
            prior.control$ssvs.traitsindex <- vector("list",num.X+1)
            for(k in 1:(num.X+1)) 
                prior.control$ssvs.traitsindex[[k]] <- rep(-1,length(which.traits[[k]])) 
            }
        if(is.list(prior.control$ssvs.traitsindex)) {
            check.ssvstraits(prior.control$ssvs.traitsindex, which.traits)
            }
        }
        
    if(any(family == "binomial") & !(length(trial.size) %in% c(1, length(family)))) 
    stop("trial.size needs to be specified if any columns are binomially distributed; can either be a single element or a vector equal to the # of columns in y")
    if(any(family == "binomial") & length(trial.size) == 1) {
            complete_trial_size <- rep(0, ncol(y))
            complete_trial_size[which(family == "binomial")] <- trial.size 
            }
    if(any(family == "binomial") & length(trial.size) == ncol(y)) 
            complete_trial_size <- trial.size
    if(all(family != "binomial")) 
            complete_trial_size <- rep(0, ncol(y))
    if(all(family == "binomial") & all(complete_trial_size == 1)) 
            family <- rep("bernoulli",ncol(y))
	
    
    if(all(family != "ordinal")) 
            num.ord.levels <- 0
    if(any(family == "ordinal")) 
            num.ord.levels <- max(y[, family == "ordinal"])
    if(all(family != "multinom")) { 
            num.multinom.levels <- 0
            index_multinom_cols <- NULL 
            }
# 	if(any(family == "multinom")) { 
# 		num.multinom.levels <- apply(y[, family == "multinom"], 2, max)
# 		index_multinom_cols <- which(family == "multinom") 
# 		}

    
    mcmc.control <- fillin.mcmc.control(x = mcmc.control)

    
    ###############
    ## MAKE JAGS SCRIPT
    ###############    
    n <- nrow(y)
    p <- ncol(y)
    n.chains <- 1; ## Run one chain only to avoid arbitrary rotation problems
    if(num.lv > 0) 
            make.jagsboralmodel(family = family, num.X = num.X, X.ind = X.ind, num.traits = num.traits, which.traits = which.traits, lv.control = lv.control, row.eff = row.eff, row.ids = row.ids, offset = offset, trial.size = complete_trial_size, n = n, p = p, model.name = model.name, prior.control = prior.control)
    if(num.lv == 0)  
            make.jagsboralnullmodel(family = family, num.X = num.X, X.ind = X.ind, num.traits = num.traits, which.traits = which.traits, row.eff = row.eff, row.ids = row.ids, offset = offset, trial.size = complete_trial_size, n = n, p = p, model.name = model.name, prior.control = prior.control)
    if(!do.fit) { 
            message("JAGS model file created only. Thank you, come again!")
            return() 
            }
    
            
    ###############
    ## FORM DATA
    ###############
    jags_data <- list("y", "n", "p", "num.lv", "num.X", "num.traits", "num.ord.levels", "num.multinom.levels")
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
            X.ind[X.ind == 1] <- 1e6
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
    if(!is.null(offset)) 
        jags_data <- c(jags_data, "offset")
            
    
    ###############
    ## FORM PARAMETERS
    ###############
    jags_params <- c("lv.coefs")
    if(num.lv > 0) 
        jags_params <- c(jags_params, "lvs")
    if(lv.control$type != "independent") 
    	jags_params <- c(jags_params, "lv.covparams")
    if(row.eff != "none") 
        jags_params <- c(jags_params, paste0("row.coefs.ID",1:ncol(row.ids)))
    if(row.eff == "random") 
        jags_params <- c(jags_params, paste0("row.sigma.ID",1:ncol(row.ids)))
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
        
        if(all(family %in% "bernoulli")) {
            Tau <- rWishart(1,p+1,diag(p))[,,1]
            Sigma <- solve(Tau)
            Z <- abs(t(rmvnorm(n,rep(0,p),Sigma)))
            Z <- ifelse(as.matrix(y), Z, -1 * Z)
            initial.list$Z <- Z }
                
        return(initial.list)
        }

    set.seed(mcmc.control$seed)
    actual.filename <- model.name
    if(is.null(actual.filename)) 
        actual.filename <- "jagsboralmodel.txt"


    ############ 
    ## THE FIT
    ###########
    jagsfit <- try(suppressWarnings(jags(data = jags_data, inits = jags.inits, parameters.to.save = jags_params, model.file = actual.filename, n.iter = mcmc.control$n.iteration, n.burnin = mcmc.control$n.burnin, n.chains = 1, n.thin = mcmc.control$n.thin)),silent=TRUE)

    #print(jagsfit)
    if(inherits(jagsfit,"try-error")) {
        lookfornegbinerror <- grep("Slicer stuck at value with infinite density", jagsfit[[1]])
        if(any(family == "negative.binomial") & length(lookfornegbinerror) == 1) { 
            message("MCMC sampling through JAGS failed. This is likely due to the prior on the dispersion (size) parameter of the negative binomial distribution been too uninformative (see below). For instance, if the error message refers to lv.coefs[25,4], then this means the MCMC sampling ran into issues for column (species) 25 in y.\n
            Please consider the following solutions: 1) remove very rare species like singletons and doubletons, as they can potentially issus for MCMC sampling, and do not provide much information about the species community in general any way, 2) adopt a tougher prior for the overdispersion parameter e.g., keeping with a uniform prior but reducing hypparams[4], or using a half-cauchy prior, 3) consider switching to a Poisson family for responses which don't appear to be overdispersed; ")
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
    ##################
    ## FORMAT INTO BIG MATRIX
    ##################
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

    
    ###########################
    ## BLING THE OUTPUT
    ##########################
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
    
    out.fit <- list()

    if(num.lv > 0) {
        lv_arr <- abind(
                matrix(apply(combined_fit_mcmc[, grep("lvs", mcmc_names)], 2, median), nrow = n),
                matrix(apply(combined_fit_mcmc[, grep("lvs", mcmc_names)], 2, mean), nrow = n),
                matrix(apply(combined_fit_mcmc[, grep("lvs", mcmc_names)], 2, IQR), nrow = n),
                matrix(apply(combined_fit_mcmc[, grep("lvs", mcmc_names)], 2, sd), nrow = n),
                along = 3)
        dimnames(lv_arr) <- list(rows = rownames(y), lv = paste0("lv", 1:num.lv), type = c("median","mean","iqr","sd"))
        
        if(new.format) 
                out.fit$lv <- lv_arr
        if(!new.format) { 
                out.fit$lv.median <- as.matrix(lv_arr[,,1]); out.fit$lv.mean <- as.matrix(lv_arr[,,2]); 
                out.fit$lv.iqr <- as.matrix(lv_arr[,,3]); out.fit$lv.sd <- as.matrix(lv_arr[,,4])
                }
                
        if(dim(lv_coefs_arr)[2] == (num.lv+2)) 
                dimnames(lv_coefs_arr) <- list(cols = colnames(y), coefficients = c("beta0", paste0("theta", 1:num.lv), "Dispersion"), type = c("median","mean","iqr","sd"))
        if(dim(lv_coefs_arr)[2] == (num.lv+1)) 
                dimnames(lv_coefs_arr) <- list(cols = colnames(y), coefficients = c("beta0", paste0("theta", 1:num.lv)), type = c("median","mean","iqr","sd"))
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
            if(new.format) 
                out.fit$lv.covparams.arr <- lv_params_arr
            if(!new.format) {
                out.fit$lv.covparams.mean <- lv_params_arr[,1]
                out.fit$lv.covparams.median <- lv_params_arr[,2]
                out.fit$lv.covparams.iqr <- lv_params_arr[,3]
                out.fit$lv.covparams.sd <- lv_params_arr[,4]
                }
            }
        }
            
            
    if(num.lv == 0) {
        if(dim(lv_coefs_arr)[2] == 2) 
            dimnames(lv_coefs_arr) <- list(cols = colnames(y), coefficients = c("beta0", "Dispersion"), type = c("median","mean","iqr","sd"))
        if(dim(lv_coefs_arr)[2] == 1) 
            dimnames(lv_coefs_arr) <- list(cols = colnames(y), coefficients = c("beta0"), type = c("median","mean","iqr","sd"))
        }	
    if(new.format) { 
        out.fit$lv.coefs <- lv_coefs_arr 
        }
    if(!new.format) { 
        out.fit$lv.coefs.median <- lv_coefs_arr[,,1]
        out.fit$lv.coefs.mean <- lv_coefs_arr[,,2]
        out.fit$lv.coefs.iqr <- lv_coefs_arr[,,3]
        out.fit$lv.coefs.sd <- lv_coefs_arr[,,4]
        if(length(out.fit$lv.coefs.median) == p) {
            out.fit$lv.coefs.median <- matrix(out.fit$lv.coefs.median, ncol = 1)
            out.fit$lv.coefs.mean <- matrix(out.fit$lv.coefs.mean, ncol = 1)
            out.fit$lv.coefs.iqr <- matrix(out.fit$lv.coefs.iqr, ncol = 1)
            out.fit$lv.coefs.sd <- matrix(out.fit$lv.coefs.sd, ncol = 1)
            rownames(out.fit$lv.coefs.median) <- rownames(out.fit$lv.coefs.mean) <- rownames(out.fit$lv.coefs.iqr) <- rownames(out.fit$lv.coefs.sd) <- colnames(y)
            colnames(out.fit$lv.coefs.median) <- colnames(out.fit$lv.coefs.mean) <- colnames(out.fit$lv.coefs.iqr) <- colnames(out.fit$lv.coefs.sd) <- "beta0"
            }
        }

    if(row.eff != "none") {
        out.fit$row.coefs <- vector("list", ncol(row.ids))
        names(out.fit$row.coefs) <- colnames(row.ids)
        for(k in 1:ncol(row.ids)) {
            row_coefs_arr <- cbind(
                apply(combined_fit_mcmc[, grep(paste0("row.coefs.ID",k,"\\["), mcmc_names),drop=FALSE], 2, median), 
                apply(combined_fit_mcmc[, grep(paste0("row.coefs.ID",k,"\\["), mcmc_names),drop=FALSE], 2, mean), 
                apply(combined_fit_mcmc[, grep(paste0("row.coefs.ID",k,"\\["), mcmc_names),drop=FALSE], 2, IQR), 
                apply(combined_fit_mcmc[, grep(paste0("row.coefs.ID",k,"\\["), mcmc_names),drop=FALSE], 2, sd))
            rownames(row_coefs_arr) <- 1:n.ID[k]
            colnames(row_coefs_arr) <- c("median","mean","iqr","sd")
                    
            if(new.format) 
                out.fit$row.coefs[[k]] <- row_coefs_arr
            if(!new.format) 
                out.fit$row.coefs[[k]] <- list(median = row_coefs_arr[,1], mean = row_coefs_arr[,2], iqr = row_coefs_arr[,3], sd = row_coefs_arr[,4])
            }
    
            if(row.eff == "random") {
                out.fit$row.sigma <- vector("list", ncol(row.ids))
                names(out.fit$row.sigma) <- colnames(row.ids)
                for(k in 1:ncol(row.ids)) {
                    row_sigma_vec <- c(
                        median(combined_fit_mcmc[, grep(paste0("row.sigma.ID",k,"$"), mcmc_names)]),
                        mean(combined_fit_mcmc[, grep(paste0("row.sigma.ID",k,"$"), mcmc_names)]),
                        IQR(combined_fit_mcmc[, grep(paste0("row.sigma.ID",k,"$"), mcmc_names)]),
                        sd(combined_fit_mcmc[, grep(paste0("row.sigma.ID",k,"$"), mcmc_names)]))
                    names(row_sigma_vec) <- c("median","mean","iqr","sd")
                            
                    if(new.format) 
                        out.fit$row.sigma[[k]] <- row_sigma_vec
                    if(!new.format) 
                        out.fit$row.sigma[[k]] <- list(median = row_sigma_vec[1], mean = row_sigma_vec[2], iqr = row_sigma_vec[3], sd = row_sigma_vec[4])	
                    }
                }
        }
            
    if(num.X > 0) {
        X_coefs_arr <- abind(
            matrix(apply(combined_fit_mcmc[, grep("X.coefs", mcmc_names),drop=FALSE], 2, median), nrow = p),
            matrix(apply(combined_fit_mcmc[, grep("X.coefs", mcmc_names),drop=FALSE], 2, mean), nrow = p),
            matrix(apply(combined_fit_mcmc[, grep("X.coefs", mcmc_names),drop=FALSE], 2, IQR), nrow = p),
            matrix(apply(combined_fit_mcmc[, grep("X.coefs", mcmc_names),drop=FALSE], 2, sd), nrow = p),
            along = 3)
        dimnames(X_coefs_arr) <- list(cols = colnames(y), X = colnames(X), type = c("median","mean","iqr","sd"))

        if(new.format) 
            out.fit$X.coefs <- X_coefs_arr
        if(!new.format) { 
            out.fit$X.coefs.median <- X_coefs_arr[,,1]
            out.fit$X.coefs.mean <- X_coefs_arr[,,2]
            out.fit$X.coefs.iqr <- X_coefs_arr[,,3]
            out.fit$X.coefs.sd <- X_coefs_arr[,,4]
            if(length(out.fit$X.coefs.median) == p) {
                out.fit$X.coefs.median <- matrix(out.fit$X.coefs.median, ncol = 1)
                out.fit$X.coefs.mean <- matrix(out.fit$X.coefs.mean, ncol = 1)
                out.fit$X.coefs.iqr <- matrix(out.fit$X.coefs.iqr, ncol = 1)
                out.fit$X.coefs.sd <- matrix(out.fit$X.coefs.sd, ncol = 1)
                rownames(out.fit$X.coefs.median) <- rownames(out.fit$X.coefs.mean) <- rownames(out.fit$X.coefs.iqr) <- rownames(out.fit$X.coefs.sd) <- colnames(y)
                colnames(out.fit$X.coefs.median) <- colnames(out.fit$X.coefs.mean) <- colnames(out.fit$X.coefs.iqr) <- colnames(out.fit$X.coefs.sd) <- colnames(X)
                }		
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
            dimnames(ssvs_indcoefs_arr) <- list(cols = colnames(y), X = colnames(X), type = c("mean","sd"))
            
            if(new.format) 
                out.fit$ssvs.indcoefs <- ssvs_indcoefs_arr
                if(!new.format) { 
                    out.fit$ssvs.indcoefs.mean <- ssvs_indcoefs_arr[,,1]
                    out.fit$ssvs.indcoefs.sd <- ssvs_indcoefs_arr[,,2]
                    }
            }
        if(any(prior.control$ssvs.index > 0)) {
            ssvs_gpcoefs_arr <- cbind(
                apply(combined_fit_mcmc[, grep("ssvs.gp", mcmc_names),drop=FALSE], 2, mean),
                apply(combined_fit_mcmc[, grep("ssvs.gp", mcmc_names),drop=FALSE], 2, sd))
                rownames(ssvs_gpcoefs_arr) <- paste0("Gp", unique(prior.control$ssvs.index[prior.control$ssvs.index > 0]))	
                colnames(ssvs_gpcoefs_arr) <- c("mean","sd")

                if(new.format) 
                    out.fit$ssvs.gpcoefs <- ssvs_gpcoefs_arr
                if(!new.format) { 
                    out.fit$ssvs.gpcoefs.mean <- ssvs_gpcoefs_arr[,1]
                    out.fit$ssvs.gpcoefs.sd <- ssvs_gpcoefs_arr[,2]
                    }
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
            
            if(new.format) 
                out.fit$ssvs.traitscoefs <- ssvs_traitscoefs_arr
            if(!new.format) { 
                out.fit$ssvs.traitscoefs.mean <- ssvs_traitscoefs_arr[,,1]
                out.fit$ssvs.traitscoefs.sd <- ssvs_traitscoefs_arr[,,2]
                }
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
                
        if(new.format) 
            out.fit$traits.coefs <- traitscoefs_arr
        if(!new.format) { 
            out.fit$traits.coefs.median <- traitscoefs_arr[,,1]
            out.fit$traits.coefs.mean <- traitscoefs_arr[,,2]
            out.fit$traits.coefs.iqr <- traitscoefs_arr[,,3]
            out.fit$traits.coefs.sd <- traitscoefs_arr[,,4]
            }
        }

#   	if(num.X > 0 & any(family == "multinom")) {
#   		out.fit$X.multinom.coefs.median <- array(apply(combined_fit_mcmc[,grep("X.multinom.params", mcmc_names)],2,median),dim=c(length(index_multinom_cols),num.X,num.multinom.levels))
#   		out.fit$X.multinom.coefs.iqr <- array(apply(combined_fit_mcmc[,grep("X.multinom.params", mcmc_names)],2,IQR),dim=c(length(index_multinom_cols),num.X,num.multinom.levels))
#   		out.fit$X.multinom.coefs.mean <- array(apply(combined_fit_mcmc[,grep("X.multinom.params", mcmc_names)],2,mean),dim=c(length(index_multinom_cols),num.X,num.multinom.levels))
#   		out.fit$X.multinom.coefs.sd <- array(apply(combined_fit_mcmc[,grep("X.multinom.params", mcmc_names)],2,sd),dim=c(length(index_multinom_cols),num.X,num.multinom.levels))
#   
#   		dimnames(out.fit$X.multinom.coefs.median) <- dimnames(out.fit$X.multinom.coefs.iqr) <- dimnames(out.fit$X.multinom.coefs.mean) <- dimnames(out.fit$X.multinom.coefs.sd) <- list("1" = index_multinom_cols, "2" = colnames(X), "level" = 1:num.multinom.levels)
#   		}

    if(any(family == "ordinal")) {
        cutoffs_arr <- cbind(
            apply(combined_fit_mcmc[, grep("cutoffs", mcmc_names),drop=FALSE], 2, median),
            apply(combined_fit_mcmc[, grep("cutoffs", mcmc_names),drop=FALSE], 2, mean),
            apply(combined_fit_mcmc[, grep("cutoffs", mcmc_names),drop=FALSE], 2, IQR),
            apply(combined_fit_mcmc[, grep("cutoffs", mcmc_names),drop=FALSE], 2, sd))
        rownames(cutoffs_arr) <- paste0(1:(num.ord.levels - 1), "|", 2:num.ord.levels) 
        colnames(cutoffs_arr) <- c("median","mean","iqr","sd")
                
        if(new.format) 
            out.fit$cutoffs <- cutoffs_arr
        if(!new.format) { 
            out.fit$cutoffs.median <- cutoffs_arr[,1]; out.fit$cutoffs.mean <- cutoffs_arr[,2]
            out.fit$cutoffs.iqr <- cutoffs_arr[,3]; out.fit$cutoffs.sd <- cutoffs_arr[,4]
            }

        if(sum(family == "ordinal") > 1 & is.null(traits)) { ## If there are traits, then ordinal random intercept is either zero (if there is only 1 ordinal column, or has the trait.sigma (if there are >1 ordinal columns)
            ordinal.sigma.vec <- c(
                median(combined_fit_mcmc[, grep("ordinal.sigma", mcmc_names)]),
                mean(combined_fit_mcmc[, grep("ordinal.sigma", mcmc_names)]),
                IQR(combined_fit_mcmc[, grep("ordinal.sigma", mcmc_names)]),
                sd(combined_fit_mcmc[, grep("ordinal.sigma", mcmc_names)]))
            names(ordinal.sigma.vec) <- c("median","mean","iqr","sd")

            if(new.format) 
                out.fit$ordinal.sigma <- ordinal.sigma.vec
            if(!new.format) { 
                out.fit$ordinal.sigma.median <- ordinal.sigma.vec[1]; out.fit$ordinal.sigma.mean <- ordinal.sigma.vec[2]
                out.fit$ordinal.sigma.iqr <- ordinal.sigma.vec[3]; out.fit$ordinal.sigma.sd <- ordinal.sigma.vec[4]            
                }
            }
        }

    if(any(family == "tweedie")) {
        powerparam_vec <- c(
            median(combined_fit_mcmc[, grep("powerparam", mcmc_names)]),
            mean(combined_fit_mcmc[, grep("powerparam", mcmc_names)]),
            IQR(combined_fit_mcmc[, grep("powerparam", mcmc_names)]),
            sd(combined_fit_mcmc[, grep("powerparam", mcmc_names)]))
        names(powerparam_vec) <- c("median","mean","iqr","sd")
                
        if(new.format) 
            out.fit$powerparam <- powerparam_vec
        if(!new.format) { 
            out.fit$powerparam.median <- powerparam_vec[1]; out.fit$powerparam.mean <- powerparam_vec[2]
            out.fit$powerparam.iqr <- powerparam_vec[3]; out.fit$powerparam.sd <- powerparam_vec[4]
            }
        }
    rm(list = ls(pattern = ".arr"))	
            
            
            
    #############
    ## HPD INTERVALS, ICS, ETC...
    #############
    #print(out.fit$lv.coefs.mean)
    get.hpds <- get.hpdintervals(y, X = X, traits = traits, row.ids = row.ids, fit.mcmc = combined_fit_mcmc, lv.control = lv.control)
    out.fit$hpdintervals <- get.hpds
    if(calc.ics) {
        warning("Please note that as of version 1.6, functions to calculate information criteria will no longer be updated. Use at your peril!")
        get_ics <- get.measures(y = y, X = X, family = family, trial.size = complete_trial_size, row.eff = row.eff, row.ids = row.ids, offset = offset, num.lv = num.lv, fit.mcmc = combined_fit_mcmc)
        ics <- c(get.dic(jagsfit), get_ics$waic, get_ics$eaic, get_ics$ebic)
        names_ics <- c("Conditional DIC", "WAIC", "EAIC", "EBIC")
        if(get_ics$do.marglik.ics) {
            ics <- c(ics, get_ics$aic.median, get_ics$bic.median, get_ics$median.logLik)
            names_ics  <- c(names_ics , "AIC at post. median", "BIC at post. median", "Marginal log-likelihood at post. median")
            }
        names(ics) <- names_ics; out.fit$ics <- ics
        }
                    
    if(save.model) 
        out.fit$jags.model <- jagsfit

    out.fit$call <- match.call()
    out.fit$n <- n
    out.fit$p <- p
    out.fit$X <- X
    X.ind[X.ind == 1e6] <- 1
    out.fit$X.ind <- X.ind
    out.fit$traits <- traits; 
    out.fit$y <- y
    out.fit$offset <- offset
    out.fit$row.eff <- row.eff; out.fit$row.ids <- row.ids

    out.fit$geweke.diag <- process.geweke(fit.mcmc = combined_fit_mcmc, y = y, X = X, traits = traits, family = family, num.lv = num.lv, row.eff = row.eff, row.ids = row.ids, num.ord.levels = num.ord.levels, prior.control = prior.control)
    out.fit$family <- family
    if(all(family == "bernoulli")) 
        out.fit$family <- rep("binomial",p) ## Switch it back to binomial for compatibility
    out.fit$num.lv <- num.lv
    out.fit$lv.control <- lv.control
    out.fit$lv.control$distmat <- NULL ## Do not save distance matrix to save space
    out.fit$num.X <- num.X; out.fit$num.traits <- num.traits
    out.fit$which.traits <- which.traits
    out.fit$calc.ics <- calc.ics
    out.fit$trial.size <- complete_trial_size
    out.fit$prior.control <- prior.control	
    out.fit$num.ord.levels <- num.ord.levels
    out.fit$mcmc.control <- mcmc.control
    
    out.fit$format <- new.format 
    #out.fit$n.chain <- out.fit$n.chains; 
    class(out.fit) <- "boral"
    if(!save.model) { 
        if(file.exists(actual.filename)) 
            file.remove(actual.filename) 
        }

    return(out.fit) 
    }
 	

 	
coefsplot <- function(covname, x, labely = NULL, est = "median", ...) {
    if(!is.null(labely)) if(!(length(labely) == nrow(x$X.coefs.median) || length(labely)==1)) 
        stop("If labely is not NULL, then it must be either of length one or a vector as long as the number of rows in x$X.coefs.median (number of species). Thanks!")
    if(!(covname %in% colnames(x$X.coefs.mean))) 
        stop("covname not found among the covariates in the boral object x")
    
    col.seq <- rep("black", length(x$hpdintervals$X.coefs[,covname,"lower"]))
    col.seq[x$hpdintervals$X.coefs[,covname,"lower"] < 0 & x$hpdintervals$X.coefs[,covname,"upper"] > 0] <- "grey"
    
    At.y <- rev(1:nrow(x$X.coefs.median)) ## So that the graphs plots in the same order as rownames of x$X.coefs.median

    if(est == "median")
        plot(x = x$X.coefs.median[,covname], y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = covname, xlim = c(min(x$hpdintervals$X.coefs[,covname,"lower"]), max(x$hpdintervals$X.coefs[,covname,"upper"])), pch = "x", ...)
    if(est == "mean")
        plot(x = x$X.coefs.mean[,covname], y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = covname, xlim = c(min(x$hpdintervals$X.coefs[,covname,"lower"]), max(x$hpdintervals$X.coefs[,covname,"upper"])), pch = "x", ...)

    segments(x0 = x$hpdintervals$X.coefs[,covname,"lower"], y0 = At.y, x1 = x$hpdintervals$X.coefs[,covname,"upper"], y1 = At.y, col = col.seq, ...)  
    abline(v=0, lty=3)
    
    if(is.null(labely)) { axis(2, at=At.y, labels = rownames(x$X.coefs.mean), las=1, ...) } 
    if(!is.null(labely)) {
        if(length(labely) == nrow(x$X.coefs.mean)) axis(2, at=At.y, labels=labely, las=1, ...) 
        if(length(labely) == 1) mtext(text = labely, side = 2, line = 3, las = 3, ...)
        } 

    if(exists("ssvs.indcoefs.mean", x)) {
        message("Posterior probabilities of inclusion for ", covname, ":")
        print(round(x$ssvs.indcoefs.mean[,covname],3))
        message()
        }
    }

		

print.boral <- function(x, ...) {
    message("Call:")
    print(x$call)
    message()
    message("Response matrix attributes\n \t# of rows: ", x$n, "\n\t# of columns: ", x$p) 
    message("Model attributes\n \tColumn families used: ", unique(x$family), "\n\t# of latent variables: ", x$num.lv, "\n\tLatent variable covariance structure", x$lv.control$type, "\n\tRow effects included (none/fixed/random)? ", x$row.eff, "\n") 

    if(any(x$family == "binomial")) 
        message("Trial sizes used (columns with binomial families): ", x$trial.size)
    
    if(any(x$family == "ordinal")) 
        message("Number of levels for ordinal data: ", x$num.ord.levels)
    
    if(x$num.X > 0) 
        message("Model matrix with ", x$num.X, " covariates also fitted\n")
    
    if(x$num.traits > 0) 
        message("Trait matrix with ", x$num.traits, " traits also included\n")
    
    if(any(x$prior.control$ssvs.index > -1) || any(unlist(x$prior.control$ssvs.traitsindex) > -1)) 
        message("SSVS has been performed on covariates\n")
    }

 	
print.summary.boral <- function(x, ...) {
    message("Call:\n")
    print(x$call)
    message()
    
    if(x$est == "median")
        message("Median point estimates\n\n LV coefficients (thetas) and dispersion parameter (if applicable)"); 
    if(x$est == "mean")
        message("Mean point estimates\n\n LV coefficients (thetas) and dispersion parameter (if applicable)")
    print(x$coefficients)
    message() 
 	
    if(!is.null(x$lv.covparams)) { 
        message("LV covariance structure parameters\n")
        print(x$lv.covparams) 
        message() 
        }
    if(!is.null(x$row.coefficients)) { 
        message("Row coefficients\n")
        print(x$row.coefficients) 
        message() 
        }
    if(!is.null(x$X.coefficients)) { 
        message("X coefficients (betas)\n")
        print(x$X.coefficients)
        message() 
        }
    if(!is.null(x$X.multinom.coefficients)) { 
        message("There are also coefficients corresponding to multinomial columns which have not been printed")
        }
    if(!is.null(x$traits.coefficients)) { 
        message("Trait coefficients")
        print(x$traits.coefficients)
        message() 
        }
    
    if(any(x$family == "ordinal")) { 
        message("Proportional odds (Cumulative probit) cutoffs") 
        print(x$cutoffs)
        message() 
        }
    if(any(x$family == "tweedie")) { 
        message("Tweedie power parameter")
        print(x$powerparam) 
        message() 
        }
    }	
		
		
summary.boral <- function(object, est = "median", ...) {
    if(est == "median") {
        gather_output <- list(call = object$call, coefficients = round(object$lv.coefs.median,3))
        if(object$num.lv > 0)
            gather_output$lvs <- round(object$lv.median,3)
        if(object$lv.control$type != "independent")
            gather_output$lv.covparams <- round(object$lv.covparams.median,3)
        if(object$row.eff != "none") {
            for(k in 1:ncol(object$row.ids)) 
                gather_output$row.coefficients[[k]] = round(object$row.coefs[[k]]$median,3)
            }
        if(object$num.X > 0) 
            gather_output$X.coefficients <- round(object$X.coefs.median,3)
        if(object$num.traits > 0) 
            gather_output$traits.coefficients <- round(object$traits.coefs.median,3)
        if(any(object$family == "ordinal")) 
            gather_output$cutoffs <- round(object$cutoffs.median,3)
        if(any(object$family == "tweedie")) 
            gather_output$powerparam <- round(object$powerparam.median,3)
        if(!is.null(object$X.multinom.coefs.median)) 
            gather_output$X.multinom.coefficients <- round(object$X.multinom.coefs.median,3) 
        }

    if(est == "mean") {
        gather_output <- list(call = object$call, coefficients = round(object$lv.coefs.mean,3))
        if(object$num.lv > 0)
            gather_output$lvs <- round(object$lv.mean,3)
        if(object$lv.control$type != "independent")
            gather_output$lv.covparams <- round(object$lv.covparams.mean,3)
        if(object$row.eff != "none") {
            for(k in 1:ncol(object$row.ids)) 
                gather_output$row.coefficients[[k]] = round(object$row.coefs[[k]]$mean,3)
            }
        if(object$num.X > 0) 
            gather_output$X.coefficients <- round(object$X.coefs.mean,3)
        if(object$num.traits > 0) 
            gather_output$traits.coefficients <- round(object$traits.coefs.mean,3)
        if(any(object$family == "ordinal")) 
            gather_output$cutoffs <- round(object$cutoffs.mean,3)
        if(any(object$family == "tweedie")) 
            gather_output$powerparam <- round(object$powerparam.mean,3)
        if(!is.null(object$X.multinom.coefs.mean)) 
            gather_output$X.multinom.coefficients <- round(object$X.multinom.coefs.mean,3) 
        }


    gather_output$est <- est
    gather_output$calc.ics <- object$calc.ics
    gather_output$trial.size <- object$trial.size
    gather_output$num.ord.levels <- object$num.ord.levels
    gather_output$prior.control$ssvs.index <- object$prior.control$ssvs.index 


    if(any(object$prior.control$ssvs.index == 0)) 
        gather_output$ssvs.indcoefs.prob <- round(object$ssvs.indcoefs.mean,3)
    if(any(object$prior.control$ssvs.index > 0)) 
        gather_output$ssvs.gpcoefs.prob <- round(object$ssvs.gpcoefs.mean,3) 
    if(any(unlist(object$prior.control$ssvs.traitsindex) == 0)) 
        gather_output$ssvs.traitscoefs.prob <- round(object$ssvs.traitscoefs.mean,3) 

    class(gather_output) <- "summary.boral"
    gather_output 
    }
 			
 			
plot.boral <- function(x, est = "median", jitter = FALSE, ...) {
    #if(all(x$family %in% c("ordinal","multinom"))) stop("Residuals are not defined, and therefore residual analysis cannot be performed, if all columns of y are ordinal")
    if(any(x$family %in% c("ordinal","multinom"))) 
        warning("For all columns of y that are ordinal, the first plot constructed is of Dunn-Smyth residuals against fitted values (i.e., the level with the highest predicted probability). Note this can make things very confusing to interpret if only SOME of the columns in the response matrix are ordinal")
    get_mus <- fitted.boral(x, est = est)$out
    get_etas <- get_mus
    get_ds_res <- ds.residuals(object = x, est = est)
    print(get_ds_res$agree.ordinal)
    get_ds_res <- get_ds_res$residuals
 	
    for(j in 1:ncol(x$y)) {
        if(x$family[j] %in% c("beta")) 
            get_etas[,j] <- log((get_mus[,j]+1e-5)/(1-get_mus[,j]+1e-5))
        if(x$family[j] %in% c("binomial")) 
            get_etas[,j] <- qnorm(get_mus[,j]+1e-5)
        if(x$family[j] %in% c("poisson","lnormal","negative.binomial","tweedie","exponential")) 
            get_etas[,j] <- log(get_mus[,j]+1e-5)
        if(x$family[j] == "normal") 
            get_etas[,j] <- (get_mus[,j]) 
        if(x$family[j] == "ordinal") { } ## Fitted values are the class with highest probability, which is already attained from fitted.boral
        }


    #.pardefault <- par(no.readonly = TRUE)	
    #par(ask = TRUE, cex = a, mar = c(5,5,2,1), cex.lab = 0.8*a, cex.main = a, las = 1, ...) 
    palette(rainbow(ncol(get_etas)))

    matplot(get_etas, get_ds_res, ylab = "Dunn-Smyth Residuals", xlab = "Linear Predictors", type="n", ...)
    for(i in 1:ncol(get_etas)) 
        points(get_etas[,i], get_ds_res[,i], col=palette()[i], ...)
    abline(h=0, lty = 2, lwd = 2)
    
# 	matplot(get_mus, get_ds_res, ylab = "Dunn-Smyth Residuals", xlab = "Fitted Values", type="n")
# 	for(i in 1:ncol(get_mus)) { points(get_mus[,i], get_ds_res[,i], col=palette()[i]) }
# 	abline(h=0, lty = 2, lwd = 2)father

    matplot(get_ds_res, ylab = "Dunn-Smyth Residuals", xlab = "Row index",type="n", xaxt = "n", ...)
    axis(side = 1, at = 1:nrow(x$y), labels = rownames(x$lv.mean), ...)
    for (i in 1:ncol(get_mus)) { points(seq(1,nrow(x$y)),get_ds_res[,i], col=palette()[i], ...) }
    abline(0,0,lty=2)

    matplot(t(get_ds_res), ylab = "Dunn-Smyth Residuals", xlab = "Column index", type="n", xaxt = "n", ...)
    axis(side = 1, at = 1:ncol(x$y), labels = rownames(x$coefs.mean), ...)
    for(i in 1:ncol(get_mus)) { points(rep(i,nrow(get_etas)), get_ds_res[,i], col=palette()[i], ...) }
    abline(h=0, lty = 2, lwd = 2)

    get_ds_res2 <- as.vector(unlist(get_ds_res))
    qqnorm(get_ds_res2[is.finite(get_ds_res2)], main = "Normal Quantile Plot", ...)
    #qqline(y = get_ds_res2[is.finite(get_ds_res2)], ...)
    
    palette("default")
    #par(.pardefault) 	
    }

###############################
## Unseen check functions
###############################

check_domarglik_ics <- function(fit.mcmc.names, index.ordinal.cols) {	
     out <- TRUE
     
     if(length(grep("traits.coefs", fit.mcmc.names)) > 1) 
          out <- FALSE 
     if(length(index.ordinal.cols) > 1) 
          out <- FALSE
     if(length(grep("lv.covparams", fit.mcmc.names)) > 1) 
          out <- FALSE 
    
    return(out)
    }
	
	
check_family <- function(family, y, traits = NULL) {
     if(length(family) != ncol(y) & length(family) != 1) 
          stop("Number of elements in family must either one or the number of columns in y.") 
     if(length(family) == 1) 
          complete_family <- rep(family, ncol(y))
     if(length(family) == ncol(y)) 
          complete_family <- family
     complete_family <- match.arg(complete_family, choices = c("negative.binomial", "ztnegative.binomial", "poisson", "ztpoisson", "binomial", 
          "normal", "lnormal", "tweedie", "ordinal", "exponential", "gamma", "beta"), several.ok = TRUE)
     if(length(complete_family) != ncol(y))
          stop("At least one of the elements in family is not supported in current version of boral...sorry!")

     if(any(complete_family == "ordinal")) {
          if(sum(y[, complete_family == "ordinal", drop = FALSE] == 0) > 0) 
               stop("For ordinal data, please shift minimum level to 1.")
          if(!is.null(traits) & (sum(complete_family == "ordinal") == 1)) 
               message("The intercept for the single ordinal response is set to zero and not regressed traits for parameter identifiability reasons.")
          }

     if(any(complete_family %in% c("ztpoisson", "ztnegative.binomial"))) {
          if(sum(y[, complete_family == "ztpoisson", drop=FALSE] < 1) > 0)
                    stop("For zero truncated count data, all values have to be greater than or equal to 1.")
          if(sum(y[, complete_family == "ztnegative.binomial", drop=FALSE] < 1) > 0)
                    stop("For zero truncated count data, all values have to be greater than or equal to 1.")
          }
          
     return(complete_family)
     }
     
     
check_lv_control <- function(num.lv, lv.control, need.distmat = TRUE) {	
     if(is.null(lv.control$type))
          lv.control$type <- "independent"
     lv.control$type <- match.arg(lv.control$type, choices = c("independent","exponential","squared.exponential","powered.exponential","spherical"))

     if(!is.null(num.lv)) {
          warning("num.lv is now a redundant argument and replaced by lv.control$num.lv. Please set num.lv = NULL.", immediate. = TRUE)
          #lv.control$num.lv <- num.lv
          }
     if(is.null(lv.control$num.lv))
          lv.control$num.lv <- 0
        
     if(need.distmat) {
          if(lv.control$type != "independent" & is.null(lv.control$distmat))
               stop("If structured latent variables are used, then please supply a distance matrix to lv.control$distmat.")
          }
            
    if(lv.control$num.lv > 5) 
        warning("We won't stop you, but please consider if you really want more than five latent variables in the model!", immediate. = TRUE)

     return(lv.control)
    }

    
check_offset <- function(offset = NULL, y) {
     if(!is.null(offset)) { 
          if(!is.matrix(offset)) 
               stop("offset should be a matrix with the same dimensions as y.")
          if(nrow(offset) != nrow(y)) 
               stop("offset should be a matrix with the same dimensions as y.")
          if(ncol(offset) != ncol(y)) 
               stop("offset should be a matrix with the same dimensions as y.")
          } 
    }

    
check_prior_control <- function(prior.control) {
    if(length(prior.control$hypparams) != 4 || length(prior.control$type) != 4) 
            stop("prior.control$type and prior.control$hypparams must be a vector of four elements. Please see boral help file as to what the elements correspond to.\n")
    if(!all(prior.control$type[-4] %in% c("normal","uniform","cauchy")))	
            stop("At least one of the first three elements of prior.control$type is not supported in the current version of boral...sorry!")
    if(!(prior.control$type[4] %in% c("uniform","halfcauchy","halfnormal")))	
            stop("The fourth element of prior.control$type is not supported in the current version of boral...sorry!")            
     }

    
check_row_ids <- function(row.ids = NULL, y) {
     if(!is.null(row.ids)) {
          row.ids <- as.matrix(row.ids)
          if(nrow(row.ids) != nrow(y)) 
                    stop("Number of rows in the matrix row.ids should be equal to number of rows in y.")
          if(is.null(colnames(row.ids))) 
                    colnames(row.ids) <- paste0("ID", 1:ncol(row.ids))
          }
        
    return(row.ids)
    }
	
	
check_row_params <- function(row.params = NULL, y, row.ids = NULL) {
     if(!is.null(row.params)) {
          if(is.null(row.ids)) {
               row.ids <- matrix(1:nrow(y), ncol = 1)
               colnames(row.ids) <- "ID1"
               }
          if(!is.list(row.params))
                    stop("row.params should be a list with length equal to the number of columns in row.ids")
          if(length(row.params) != ncol(row.ids))
                    stop("row.params should be a list with length equal to the number of columns in row.ids")
          }
    }

check_ssvstraits <- function(ssvs.traitsindex, which.traits) {
     if(!is.null(which.traits)) {
          if(length(ssvs.traitsindex) != length(which.traits))
               stop("Both prior.control$ssvs.traitsindex and which.traits should be lists of equal length.")
          if(!all(unlist(which.traits) >= 0))
               stop("All elements of which.traits must be non-negative.")
          if(!all(unlist(ssvs.traitsindex) %in% c(-1,0)))
               stop("All elements in the list prior.control$ssvs.traitsindex should be either equal to -1 (no SSVS) or 0 (SSVS applied).")
          for(k in 1:length(which.traits)) {
               if(which.traits[[k]][1] == 0 & !all(ssvs.traitsindex[[k]] == -1))
                    stop(paste0("If which.traits[[",k,"]][1] == 0, then all the elements of ssvs.traitsindex[[",k,"]] must equal to -1. That is, if traits are not used then for covariate ",k," and no SSVS can be done on this."))
               if(any(which.traits[[k]] > 0) & length(ssvs.traitsindex[[k]]) != length(which.traits[[k]]))
                    stop(paste0("If the elements of which.traits[[",k,"]] are positive, then the length of ssvs.traitsindex[[",k,"]] must match the length of which.traits[[",k,"]]. That is, if traits are used then for covariate ",k, " then a corresponding index needs to be supplied to determine if SSVS needs to be done for this."))
               }
          }
    }

    
check_traits <- function(traits, y) {
     if(!is.null(traits)) { 
          if(!is.matrix(traits)) 
               traits <- as.matrix(traits) 
          if(nrow(traits) != ncol(y))
               stop("If traits are supplied, then please ensure the number of rows in traits i.e., number of species, is equal to the number of columns in the response matrix.") 
          if(any(apply(traits,2,function(x) all(x == 1)))) 
               stop("No intercept column should be included in traits. It will be included automatically.")
          }
     }

     
check_trial_size <- function(family, trial.size, create.complete.trial.size = FALSE, y = NULL) {
     if(any(family == "binomial") & !(length(trial.size) %in% c(1, length(family)))) 
          stop("trial.size needs to be specified if any columns are binomially distributed; can either be a single element or a vector equal to the # of columns in y.")
     
     if(create.complete.trial.size) {
          if(any(family == "binomial") & length(trial.size) == 1) {
               complete_trial_size <- rep(0, ncol(y))
               complete_trial_size[which(family == "binomial")] <- trial.size 
               }
          if(any(family == "binomial") & length(trial.size) == ncol(y)) 
               complete_trial_size <- trial.size
          if(all(family != "binomial")) 
               complete_trial_size <- rep(0, ncol(y))
          return(complete_trial_size)
          }
     }

     
check_X_ind <- function(X.ind = NULL, p, num.X, prior.control) {
     if(!is.null(X.ind)) {
          X.ind <- as.matrix(X.ind)
          if(nrow(X.ind) != p || ncol(X.ind) != num.X)
               stop("X.ind must be a matrix with the number of rows equal to the # of columns in y and the number of columns in X.")
          if(!all(X.ind %in% c(0,1)))
               stop("All elements of X.ind must either equal to 1 or 0, corresponding to whether a particular covariate is included or excluded for a particular column response, respectively.")
          }
     if(any(prior.control$ssvs.index != -1))
          message("X.ind is ignored for any columns on which SSVS is used.")

     return(X.ind)
     }

     
check_which_traits <- function(num.traits, which.traits, traits = NULL, y = NULL, num.X, makejagsboralfile_messages = FALSE) {

     if(num.traits > 0 & makejagsboralfile_messages == FALSE) {
          if(num.X == 0 & num.traits > 0) 
               stop("num.traits > 0 suggests traits are to be regressed against covariates X, so please supply X.") 
          if(is.null(which.traits)) 
               stop("If traits are supplied, then please also supply which.traits to inform what traits are regressed against which covariates.") 
          if(nrow(traits) != ncol(y))
               stop("If traits are supplied, then please ensure the number of rows in traits i.e., number of species, is equal to the number of columns in y.") 
          if((num.X+1) != length(which.traits))
               stop("which.traits should have equal to 1+ncol(X).") 
          if(any(sapply(which.traits, length) > num.traits)) 
               stop("Each element in the list which.traits should have at most ncol(traits) elements.") 
          if(any(sapply(which.traits, function(x) any(x > ncol(traits))))) 
               stop("The values contained in the list which.traits can be takes from 1 to ncol(traits).") 
          }

          
     if(num.traits > 0 & makejagsboralfile_messages == TRUE) {
          if(num.X == 0) 
               stop("num.traits > 0 suggests traits are to be regressed against covariates X, so please set num.X > 0.") 
          if(is.null(which.traits)) 
               stop("If num.traits > 0, then please supply which.traits to inform what traits are regressed against which covariates.") 
          if(!is.null(which.traits) & ((num.X+1) != length(which.traits))) 
               stop("which.traits should be a list with length 1+num.X.") 
          if(!is.null(which.traits) & any(sapply(which.traits,length) > num.traits)) 
               stop("Each element in the list which.traits should have at most num.traits elements.") 
          }
     }


     
     

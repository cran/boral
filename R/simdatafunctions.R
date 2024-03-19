##########
## Auxilary functions to simulate data
##########

## If true.lv is supplied, then create.life conditions on this and everthing in distmat is overwritten
## If true.lv is not supplied, then create.life will use the lv.control to simulate new lvs 
## Similarly, if true.ranef is supplied, then create.life conditions on this. But if true.ranef is not supplied, then create.life will use the ranef.params to to simulate response-specific random intercepts
create.life <- function(true.lv = NULL, lv.coefs, 
    lv.control = list(num.lv = 0, type = "independent", lv.covparams = NULL, distmat = NULL),
    X = NULL, X.coefs = NULL, traits = NULL, traits.coefs = NULL, 
    family, row.eff = "none", row.params = NULL, row.ids = NULL, 
    true.ranef = NULL, ranef.params = NULL, ranef.ids = NULL, 
    offset = NULL, trial.size = 1, cutoffs = NULL, powerparam = NULL, 
    manual.dim = NULL, save.params = FALSE) { 

     n <- max(nrow(true.lv), nrow(X))
     p <- max(nrow(lv.coefs), nrow(X.coefs), length(cutoffs))
     if(is.null(dim(lv.coefs))) 
          lv.coefs <- as.matrix(lv.coefs)
     
     if((is.null(n) | is.null(p)) & is.null(manual.dim)) 
          stop("Sorry, but boral cannot determine the number of rows and columns for the response matrix. Please supply manual.dim as vector containing n and p.")
     if((is.null(n) | is.null(p)) & !is.null(manual.dim)) { 
          n <- manual.dim[1]
          p <- manual.dim[2] 
          }
               
     ## Assumes simulation is based off true.lv, so reset lv.control
     if(!is.null(true.lv)) {
          lv.control <- list(num.lv = ncol(true.lv), type = "independent")
          }
     ## Assumes simulation requires generation of true.lv depending on lv.control$num.lv
     if(is.null(true.lv)) {        
          if(lv.control$num.lv > 0) {
               if(lv.control$type != "independent" & (is.null(lv.control$distmat) || is.null(lv.control$lv.covparams)))
                    stop("For structured latent variables are to be generated, please supply both lv.control$distmat and lv.control$lv.covparams")
               #if(lv.control$type != "independent" & (nrow(distmat) != n || ncol(distmat) != n))
               #    stop("lv.control$distmat should be a symmetric matrix with the same number of rows as X")

               if(lv.control$type == "exponential")
                    covmat_chol <- (chol(exp(-lv.control$distmat/lv.control$lv.covparams[1])))
               if(lv.control$type == "squared.exponential")
                    covmat_chol <- (chol(exp(-(lv.control$distmat/lv.control$lv.covparams[1])^2)))
               if(lv.control$type == "powered.exponential")
                    covmat_chol <- (chol(exp(-(lv.control$distmat/lv.control$lv.covparams[1])^lv.control$lv.covparams[2])))
               if(lv.control$type == "spherical")
                    covmat_chol <- (chol((lv.control$distmat < lv.control$lv.covparams[1])*(1 - 1.5*lv.control$distmat/lv.control$lv.covparams[1] + 0.5*(lv.control$distmat/lv.control$lv.covparams[1])^3)))
                    
               true.lv <- rmvnorm(n, mean = rep(0,lv.control$num.lv))
               if(lv.control$type != "independent")
                    true.lv <- crossprod(covmat_chol, true.lv)
               }
          }
     num.lv <- lv.control$num.lv
     #      
     if(!is.null(X)) { 
          if(is.null(X.coefs) & is.null(traits))
               stop("If X is supplied, then either X.coefs or traits and traits.coefs must be supplied.")
          }
     if(!is.null(X)) { 
          if(!is.matrix(X)) 
               X <- as.matrix(X)
          if(any(apply(X,2,function(x) all(x == 1)))) 
               stop("No intercept column should be included in X.") 
          }

     if((is.null(traits) & !is.null(traits.coefs)) | (!is.null(traits) & is.null(traits.coefs))) 
          stop("If traits is supplied, then traits.coefs must also be supplied.")
     if(!is.null(traits.coefs)) 
          message("Since trait.coefs has been supplied, then X.coefs will be ignored (it will instead be drawn as random effects based off trait.coefs).")
          check_traits <- function(traits = traits, y = matrix(NA,nrow=1,ncol=p)) 


     if(length(family) != p & length(family) != 1)
          stop("Number of elements in family must be either 1 or equal to # of rows in lv.coefs/X.coefs/second number in manual.dim")
     if(length(family) == 1) 
          family <- rep(family, p)
     if(!all(family %in% c("negative.binomial", "ztnegative.binomial", "poisson", "ztpoisson", "binomial", "normal", "lnormal", 
               "tweedie", "ordinal", "exponential", "beta"))) 
          stop("One of the elements in family is not compatible with current version of boral...sorry!")
          check_trial_size(family = family, trial.size = trial.size)
     if(length(trial.size) == 1) 
          trial.size <- rep(trial.size, p)
               
     if(any(family == "ordinal") & is.null(cutoffs)) 
          stop("cutoffs (an ascending vector of intercepts for proportional odds regression) must be supplied if any columns are ordinal data.")
     if(any(family == "ordinal")) 
          index.ord.cols <- which(family == "ordinal")
     if(!is.null(cutoffs)) {
          num.ord.levels <- length(cutoffs) + 1
          cutoffs <- sort(cutoffs)
          message("Sorting cutoffs...just in case!") 
          }
     if(any(family == "tweedie") & is.null(powerparam)) 
          stop("Common powerparam must be supplied if any columns are tweedie data (Var = dispersion*mu^powerparam).")

          
     sim_y <- matrix(NA, nrow = n, ncol = p)

     check_offset(offset = offset, y = sim_y)
     
     row.coefs <- NULL
     if(row.eff != "none") {
          if(is.null(row.ids)) {
               row.ids <- matrix(1:nrow(sim_y), ncol = 1)
               colnames(row.ids) <- "ID1"
               }
          row.ids <- check_row_ids(row.ids = row.ids, y = sim_y)
          check_row_params(row.params = row.params, y = sim_y, row.ids = row.ids)
          }
     if(row.eff == "fixed") 
          row.coefs <- row.params
     if(row.eff == "random") {
          row.coefs <- vector("list", ncol(row.ids))
          for(k in 1:ncol(row.ids)) 
               row.coefs[[k]] <- rnorm(length(unique(row.ids[,k])), mean = 0, sd = row.params[[k]]) 
          }
          
          
     ranef_coefs <- NULL     
     if(!is.null(ranef.ids) & (is.null(true.ranef) & is.null(ranef.params))) 
          stop("If you want response-specific random intercepts in the model, then you need supply either true.ranef or ranef.params, in conjunction with ranef.ids")
     ranef.ids <- check_ranef_ids(ranef.ids = ranef.ids, y = sim_y)
     check_ranef_params(ranef.params = ranef.params, y = sim_y, ranef.ids = ranef.ids)
     if(!is.null(true.ranef) & !is.null(ranef.ids)) {
          if(is.null(ranef.ids))
               stop("ranef.ids needs to be supplied if you want response-specific random intercepts in the model (as determined by true.ranef).")
          if(!is.list(true.ranef))
               stop("true.ranef must be a list of length ncol(ranef.ids).")
          for(k0 in 1:length(true.ranef)) {
               if(nrow(true.ranef[[k0]]) != p)
                    stop("The k-th element in the list true.ranef should be matrix of dimension p by length(unique(ranef.ids[,k])).")
               if(ncol(true.ranef[[k0]]) != length(unique(ranef.ids[,k0])))
                    stop("The k-th element in the list true.ranef should be matrix of dimension p by length(unique(ranef.ids[,k])).")
               }

          ranef_coefs <- true.ranef
          }
     if(is.null(true.ranef) & !is.null(ranef.ids)) {
          ranef_coefs <- vector("list", ncol(ranef.ids))
               for(k in 1:ncol(ranef.ids)) 
                    ranef_coefs[[k]] <- t(sapply(1:p, function(j) rnorm(length(unique(ranef.ids[,k])), mean = 0, sd = ranef.params[j,k])))
          }
          
     
     if(num.lv > 5) 
          warnings("We won't stop you, but please consider if you really want more than five latent variables in the model", immediate. = TRUE)
     
     
     if(num.lv == 0)
          eta <- matrix(0, nrow = n, ncol = p)
     if(num.lv > 0) 
          eta <- tcrossprod(true.lv, lv.coefs[,2:(num.lv + 1)])
     if(is.null(traits.coefs)) 
          eta <- eta + tcrossprod(rep(1,n), as.matrix(lv.coefs[,1]))
     if(!is.null(X.coefs) & is.null(traits.coefs)) 
          eta <- eta + tcrossprod(as.matrix(X), X.coefs)
     if(!is.null(offset)) 
          eta <- eta + offset
     if(!is.null(traits.coefs)) {
          X.coefs <- matrix(0, p, ncol(X)) ## overwrite X.coefs
          lv.coefs[,1] <- rnorm(p, cbind(1,traits)%*%traits.coefs[1,-ncol(traits.coefs)], sd = traits.coefs[1,ncol(traits.coefs)]) ## overwrite spp-specific intercepts
     
          if(any(family == "ordinal")) {
               if(length(index.ord.cols) == 1) 
                    lv.coefs[index.ord.cols,1] <- 0 ## If there is just one ordinal column, then the random intercept for this column is zero for identifiability reasons
                    }
          for(k in 1:ncol(X)) 
               X.coefs[,k] <- rnorm(p, cbind(1,traits)%*%traits.coefs[k+1,-ncol(traits.coefs)], sd = traits.coefs[k+1,ncol(traits.coefs)])
               eta <- eta + tcrossprod(cbind(1,as.matrix(X)), cbind(lv.coefs[,1],X.coefs)) 
               }    
     if(!is.null(row.coefs)) { 
          for(k in 1:ncol(row.ids)) 
               eta <- eta + row.coefs[[k]][row.ids[,k]] 
          }
     if(!is.null(ranef_coefs)) { 
          for(k in 1:ncol(ranef.ids)) 
               eta <- eta + t(ranef_coefs[[k]])[ranef.ids[,k],] 
          }
          
     for(j in 1:p) {
          if(family[j] == "binomial") 
               sim_y[,j] <- rbinom(n, size = trial.size[j], prob = pnorm(eta[,j]))
          if(family[j] == "poisson") 
               sim_y[,j] <- rpois(n, lambda = exp(eta[,j]))
          if(family[j] == "ztpoisson") 
               sim_y[,j] <- rztpois(n, lambda = exp(eta[,j]))
          if(family[j] == "negative.binomial") 
               sim_y[,j] <- rnbinom(n, mu = exp(eta[,j]), size = 1/(lv.coefs[j,ncol(lv.coefs)]+1e-5))
          if(family[j] == "negative.binomial") 
               sim_y[,j] <- rztnbinom(n, mu = exp(eta[,j]), size = 1/(lv.coefs[j,ncol(lv.coefs)]+1e-5))
          if(family[j] == "exponential") 
               sim_y[,j] <- rexp(n, rate = 1/exp(eta[, j]))
          if(family[j] == "gamma") 
               sim_y[,j] <- rgamma(n, shape = exp(eta[, j])*lv.coefs[j,ncol(lv.coefs)], rate = lv.coefs[j, ncol(lv.coefs)])
          if(family[j] == "beta") 
               sim_y[, j] <- rbeta(n, shape1 = lv.coefs[j,ncol(lv.coefs)]*exp(eta[,j])/(1 + exp(eta[,j])), shape2 = lv.coefs[j,ncol(lv.coefs)]*(1-exp(eta[,j])/(1 +exp(eta[,j]))))
          if(family[j] == "normal") 
               sim_y[, j] <- rnorm(n, mean = eta[, j], sd = (lv.coefs[j,ncol(lv.coefs)]))
          if(family[j] == "lnormal") 
               sim_y[, j] <- rlnorm(n, meanlog = eta[, j], sdlog = (lv.coefs[j,ncol(lv.coefs)]))
          if(family[j] == "tweedie") 
               sim_y[, j] <- rTweedie(n, mu = exp(eta[, j]), phi = lv.coefs[j,ncol(lv.coefs)], p = powerparam)
          if(family[j] == "ordinal") {
               get_probs <- ordinal_conversion(n = n, lv = true.lv, lv.coefs.j = lv.coefs[j, ], num.lv = num.lv, row.coefs = row.coefs, row.ids = row.ids, X = X, X.coefs.j = X.coefs[j,], cutoffs = cutoffs, est = "ignore")
               for(i in 1:n) 
                    sim_y[i, j] <- sample(1:num.ord.levels, 1, prob = get_probs[i,])
               }
          }
          
     if(!save.params) 
          out <- sim_y
     if(save.params) 
          out <- list(resp = sim_y, true.lv = true.lv, lv.coefs = lv.coefs, lv.covparams = lv.control$lv.covparams, X.coefs = X.coefs, traits.coefs = traits.coefs, row.params = row.params, row.coefs = row.coefs, row.ids = row.ids,
          ranef.params = ranef.params, ranef.coefs = ranef_coefs, ranef.ids = ranef.ids, 
          cutoffs = cutoffs, powerparam = powerparam)
     
     return(out)
     }
	


      
## Wrapper for create.life: Takes a boral model and applies create.life to it, if possible
## Note always new row effects if row.eff = "random", this is probably not ideal?
function() {
     object <- spiderfit_nb
     nsims <- 1
     seed = NULL
     new.lvs = FALSE
     new.ranefs = TRUE
     distmat = NULL
     est = "median"
     }
simulate.boral <- function(object, nsim = 1, seed = NULL, new.lvs = FALSE, new.ranefs = FALSE, distmat = NULL, est = "median", ...) { 
     if(!inherits(object, "boral"))
          stop("object must be of class boral.")
     
     if(est == "mean") {
          true_mod <- list(lv.coefs = object$lv.coefs.mean, lv = object$lv.mean, X.coefs = object$X.coefs.mean, 
               ranef.coefs = object$ranef.coefs.mean, ranef.params = object$ranef.params.mean, ranef.ids = object$ranef.ids,
               traits = object$traits, traits.coefs = object$traits.coefs.mean, cutoffs = object$cutoffs.mean, powerparam = object$powerparam.mean, 
               lv.covparams = object$lv.covparams.mean) 
          if(object$row.eff == "fixed") 
               true_mod$row.params <- lapply(object$row.coefs, function(x) x$mean)
          if(object$row.eff == "random") 
               true_mod$row.params <- lapply(object$row.sigma, function(x) x$mean)
          if(!is.null(object$ranef.ids)) 
               true_mod$ranef.params <- object$ranef.sigma.mean
          }
     
     if(est == "median") {
          true_mod <- list(lv.coefs = object$lv.coefs.median, lv = object$lv.median, X.coefs = object$X.coefs.median, traits = object$traits, 
               ranef.coefs = object$ranef.coefs.median, ranef.params = object$ranef.params.median, ranef.ids = object$ranef.ids,
               traits.coefs = object$traits.coefs.median, cutoffs = object$cutoffs.median, powerparam = object$powerparam.median, 
               lv.covparams = object$lv.covparams.median) 
          if(object$row.eff == "fixed") 
               true_mod$row.params <- lapply(object$row.coefs, function(x) x$median)
          if(object$row.eff == "random") 
               true_mod$row.params <- lapply(object$row.sigma, function(x) x$median)
          if(!is.null(object$ranef.ids)) 
               true_mod$ranef.params <- object$ranef.sigma.median
          }

     if(!is.null(seed)) 
          set.seed(seed)

          
     if(!new.lvs) {
          message("Same latent variables used across simulated datasets (if the model involves them).")          
          if(!new.ranefs) {
               message("Same response-specific random intercepts used across simulated datasets (if the model involves them).")          
               out <- replicate(nsim, 
                    create.life(true.lv = true_mod$lv, lv.coefs = true_mod$lv.coefs, X = object$X, X.coefs = true_mod$X.coefs, traits = object$traits, 
                         traits.coefs = true_mod$traits.coefs, family = object$family, 
                         row.eff = object$row.eff, row.params = true_mod$row.params, row.ids = object$row.ids, 
                         true.ranef = true_mod$ranef.coefs, ranef.ids = object$ranef.ids,
                         offset = object$offset, trial.size = object$trial.size, cutoffs = true_mod$cutoffs, powerparam = true_mod$powerparam))
               }
          if(new.ranefs) {
               message("New response-specific random intercepts used across simulated datasets (if the model involves them).")          
               out <- replicate(nsim, 
                    create.life(true.lv = true_mod$lv, lv.coefs = true_mod$lv.coefs, X = object$X, X.coefs = true_mod$X.coefs, traits = object$traits, 
                         traits.coefs = true_mod$traits.coefs, family = object$family, 
                         row.eff = object$row.eff, row.params = true_mod$row.params, row.ids = object$row.ids, 
                         ranef.params = true_mod$ranef.params, ranef.ids = object$ranef.ids,
                         offset = object$offset, trial.size = object$trial.size, cutoffs = true_mod$cutoffs, powerparam = true_mod$powerparam))
               }
          }
          
     if(new.lvs) {
          message("New latent variables generated for each simulated dataset (if the model involves them).")
          lv.control <- object$lv.control
          lv.control$lv.covparams <- true_mod$lv.covparams
          if(lv.control$type != "independent" & is.null(distmat))
               stop("For structured latent variables are to be generated, please supply distmat.")
          lv.control$distmat <- distmat

          
          if(!new.ranefs) {
               message("Same response-specific random intercepts used across simulated datasets (if the model involves them).")          
               out <- replicate(nsim, 
                    create.life(true.lv = NULL, lv.coefs = true_mod$lv.coefs, lv.control = lv.control, X = object$X, X.coefs = true_mod$X.coefs, 
                         traits = object$traits, traits.coefs = true_mod$traits.coefs, family = object$family, 
                         row.eff = object$row.eff, row.params = true_mod$row.params, row.ids = object$row.ids, 
                         true.ranef = true_mod$ranef.coefs, ranef.ids = object$ranef.ids,
                         offset = object$offset, trial.size = object$trial.size, cutoffs = true_mod$cutoffs, powerparam = true_mod$powerparam))
                    }
          if(new.ranefs) {
               message("New response-specific random intercepts used across simulated datasets (if the model involves them).")          
               out <- replicate(nsim, 
                    create.life(true.lv = NULL, lv.coefs = true_mod$lv.coefs, lv.control = lv.control, X = object$X, X.coefs = true_mod$X.coefs, 
                         traits = object$traits, traits.coefs = true_mod$traits.coefs, family = object$family, 
                         row.eff = object$row.eff, row.params = true_mod$row.params, row.ids = object$row.ids, 
                         ranef.params = true_mod$ranef.params, ranef.ids = object$ranef.ids,
                         offset = object$offset, trial.size = object$trial.size, cutoffs = true_mod$cutoffs, powerparam = true_mod$powerparam))
                    }
          }
            
     set.seed(NULL)
     return(out)
     }



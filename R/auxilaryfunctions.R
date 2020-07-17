##########
## Auxilary functions
##########
print.boral.startupInfo <- function() { 
     version <- packageVersion("boral")
     startup_note <- paste0("This is boral version ", version, ". If you recently updated boral, please check news(package = \"boral\") for the updates in the latest version.")
     packageStartupMessage(startup_note)
     }

.onAttach <- function(...) {
     print.boral.startupInfo()
     }

## Dunn-Smyth residuals
## Also create a confusion matrix for ordinal and multinomial data
ds.residuals <- function(object, est = "median") {  
     n <- object$n
     p <- object$p
     num.lv <- object$num.lv
     num.ord.levels <- object$num.ord.levels
     X <- object$X
     y <- object$y
     mus <- fitted.boral(object, est = est)

     if(any(object$family == "ordinal")) {
          message("One or more columns of y have ordinal responses. Constructing a single confusion matrix for these.")
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

          
     dsres_out <- matrix(NA,nrow = n, ncol = p)
     rownames(dsres_out) <- rownames(y)
     colnames(dsres_out) <- colnames(y)
     for(i in 1:n) { for(j in 1:p) {
          if(object$family[j] == "poisson") { 
               a <- ppois(as.vector(unlist(y[i,j]))-1, mus$out[i,j]); 
               b <- ppois(as.vector(unlist(y[i,j])), mus$out[i,j]); 		
               dsres_out[i,j] <- qnorm(runif(n = 1, min = a, max = b)) 
               }
          if(object$family[j] == "ztpoisson") { 
               a <- pztpois(as.vector(unlist(y[i,j]))-1, lambda = mus$out[i,j]); 
               b <- pztpois(as.vector(unlist(y[i,j])), lambda = mus$out[i,j]); 		
               dsres_out[i,j] <- qnorm(runif(n = 1, min = a, max = b)) 
               }
          if(object$family[j] == "negative.binomial") {
               if(est == "median") 
                    phis <- object$lv.coefs.median[,num.lv+2]+1e-5
               if(est == "mean") 
                    phis <- object$lv.coefs.mean[,num.lv+2]+1e-5
               a <- pnbinom(as.vector(unlist(y[i,j]))-1, mu=mus$out[i,j], size=1/phis[j]); 
               b <- pnbinom(as.vector(unlist(y[i,j])), mu=mus$out[i,j], size=1/phis[j])
               dsres_out[i,j] <- qnorm(runif(n = 1, min = a, max = b)) 
               }
          if(object$family[j] == "ztnegative.binomial") {
               if(est == "median") 
                    phis <- object$lv.coefs.median[,num.lv+2]+1e-5
               if(est == "mean") 
                    phis <- object$lv.coefs.mean[,num.lv+2]+1e-5
               a <- pztnbinom(as.vector(unlist(y[i,j]))-1, mu=mus$out[i,j], size=1/phis[j]); 
               b <- pztnbinom(as.vector(unlist(y[i,j])), mu=mus$out[i,j], size=1/phis[j])
               dsres_out[i,j] <- qnorm(runif(n = 1, min = a, max = b)) 
               }
          if(object$family[j] == "binomial") { 
               a <- pbinom(as.vector(unlist(y[i,j]))-1, object$trial.size[j], prob = mus$out[i,j]); 
               b <- pbinom(as.vector(unlist(y[i,j])), object$trial.size[j], prob = mus$out[i,j])
               dsres_out[i,j] <- qnorm(runif(n = 1, min = a, max = b)) 
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
               a <- plnorm(as.vector(unlist(y[i,j])), log(mus$out[i,j]), sdlog = (phis[j]))
               dsres_out[i,j] <- qnorm(a) 
               }
          if(object$family[j] == "tweedie") { 
               if(est == "median") 
                    phis <- object$lv.coefs.median[,num.lv+2]
               if(est == "mean") 
                    phis <- object$lv.coefs.mean[,num.lv+2]
               a <- b <- pTweedie(as.vector(unlist(y[i,j])), mu = mus$out[i,j], phi = phis[j], p = powerparam)
               if(as.vector(unlist(y[i,j])) == 0) 
                    a <- 0
               dsres_out[i,j] <- qnorm(runif(n = 1, min = a, max = b)) 
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
fitted.boral <- function(object, est = "median", linear.predictor = FALSE, ...) {
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
    
    
     if(linear.predictor)
          return(eta)
    
    
     index_multinom_cols <- which(object$family == "multinom")
     for(j in 1:p) {
          if(object$family[j] %in% c("binomial")) 
               fitted_out[,j] <- binomial(link = "probit")$linkinv(eta[,j])
          if(object$family[j] %in% c("beta")) 
               fitted_out[,j] <- exp(eta[,j])/(1+exp(eta[,j]))
          if(object$family[j] %in% c("poisson","lnormal","negative.binomial","tweedie","exponential","gamma")) 
               fitted_out[,j] <- exp(eta[,j])
          if(object$family[j] %in% c("ztpoisson")) 
               fitted_out[,j] <- ztpoisson()$linkinv(eta[,j])
          if(object$family[j] %in% c("ztnegative.binomial")) {
               if(est == "median") 
                    phis <- object$lv.coefs.median[,num.lv+2]+1e-5
               if(est == "mean") 
                    phis <- object$lv.coefs.mean[,num.lv+2]+1e-5
               fitted_out[,j] <- exp(eta[j,])/pnbinom(0, size = 1/phis[j], mu = exp(eta[j,]), lower.tail = FALSE, log.p = FALSE)
               }
          if(object$family[j] == "normal") 
               fitted_out[,j] <- eta[,j] 
     # 		if(object$family[j] == "multinom") {
     # 			if(est == "median") { if(!is.null(object$X.multinom.coefs.median)) eta2 <- eta[,j] + as.matrix(X)%*%object$X.multinom.coefs.median[which(index_multinom_cols == j),,] }
     # 			if(est == "mean") { if(!is.null(object$X.multinom.coefs.mean)) eta2 <- eta[,j] + as.matrix(X)%*%object$X.multinom.coefs.mean[which(index_multinom_cols == j),,] }
     # 			get_probs <- exp(eta2)/apply(exp(eta2),1,sum)	
     # 			fitted_multinom_probs[,j,] <- get_probs
     # 			}

          if(object$family[j] == "ordinal") {
               if(est == "median")
                    fitted_ordinal_probs[,j,] <- ordinal_conversion(n = n, lv = object$lv.median, lv.coefs.j = object$lv.coefs.median[j,], num.lv = num.lv, row.coefs = object$row.coefs, row.ids = object$row.ids, X = X, X.coefs.j = object$X.coefs.median[j,], cutoffs = object$cutoffs.median, est = "median")
               if(est == "mean")
                    fitted_ordinal_probs[,j,] <- ordinal_conversion(n = n, lv = object$lv.mean, lv.coefs.j = object$lv.coefs.mean[j,], num.lv = num.lv, row.coefs = object$row.coefs, row.ids = object$row.ids, X = X, X.coefs.j = object$X.coefs.mean[j,], cutoffs = object$cutoffs.mean, est = "mean")
               fitted_out[,j] <- apply(fitted_ordinal_probs[,j,],1,which.max) ## get max predicted probability
               }
          }	

     return(list(ordinal.probs = fitted_ordinal_probs, out = fitted_out))
     }

	
## Calculates DIC based on the conditional log-likelihood
get.dic <- function(jagsfit) { 
     jagsfit$BUGSoutput$DIC
     }
	

## Simple extraction of MCMC samples from fitted boral object
get.mcmcsamples <- function(object) {
     fit.mcmc <- object$jags.model$BUGSoutput
     if(is.null(fit.mcmc)) 
          stop("MCMC samples not found. Please use save.model = TRUE to save MCMC samples when using boral.")
     fit.mcmc <- mcmc(fit.mcmc$sims.matrix, start = 1, thin = object$mcmc.control$n.thin) ## Thanks to Guilliaume Blanchet for the original formatting!

     return(fit.mcmc)
     }

     

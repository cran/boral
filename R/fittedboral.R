## Fitted values
## For ordinal and multinomial data, returns a matrix of probabilities for each vector of rows
fitted.boral <- function(object, est = "median", include.ranef = TRUE, linear.predictor = FALSE, ...) {
     n <- object$n
     p <- object$p
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
          for(k in 1:ncol(object$row.ids)) 
               eta <- eta + matrix(object$row.coefs[[k]]$median[object$row.ids[,k]], nrow = n, ncol = p, byrow = FALSE)
          }
     if(!is.null(object$row.ids) && est == "mean") {
          for(k in 1:ncol(object$row.ids)) 
               eta <- eta + matrix(object$row.coefs[[k]]$mean[object$row.ids[,k]], nrow = n, ncol = p, byrow = FALSE) 
          }
    
     if(include.ranef) {
          if(!is.null(object$ranef.ids) && est == "median") {
               for(k0 in 1:ncol(object$ranef.ids)) 
                    eta <- eta + t(object$ranef.coefs.median[[k0]])[object$ranef.ids[,k0],] 
               }
          if(!is.null(object$ranef.ids) && est == "mean") {
               for(k0 in 1:ncol(object$ranef.ids)) 
                    eta <- eta + t(object$ranef.coefs.mean[[k0]])[object$ranef.ids[,k0],] 
               }
          }
          
          
     ##-------------------------
     ## Done!
     ##-------------------------
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

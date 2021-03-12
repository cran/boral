###############################
## Unseen functions
###############################
    
## Construct strings for the prior distributions used in makejagsboralmodel and makejagsboralnullmodel
construct_prior_strings <- function(x) {
     x$type[1] <- match.arg(x$type[1], choices = c("normal","cauchy","uniform"))
     x$type[2] <- match.arg(x$type[2], choices = c("normal","cauchy","uniform"))
     x$type[3] <- match.arg(x$type[3], choices = c("normal","cauchy","uniform"))
     x$type[4] <- match.arg(x$type[4], choices = c("halfnormal","halfcauchy","uniform"))

     
     if(x$type[1] == "normal") 
          prior.string1 <- paste0("dnorm(0,",1/x$hypparams[1],")")
     if(x$type[1] == "cauchy") 
          prior.string1 <- paste0("dt(0,",1/x$hypparams[1],",1)")
     if(x$type[1] == "uniform") 
          prior.string1 <- paste0("dunif(-",x$hypparams[1],",",x$hypparams[1],")")

     if(x$type[2] == "normal") {
          prior.string2 <- paste0("dnorm(0,", 1/x$hypparams[2], ")")
          prior.string22 <- paste0("dnorm(0,",1/x$hypparams[2],")I(0,)")
          }
     if(x$type[2] == "cauchy") {
          prior.string2 <- paste0("dt(0,",1/x$hypparams[2],",1)")
          prior.string22 <- paste0("dt(0,",1/x$hypparams[4],",1)I(0,)")
          }
     if(x$type[2] == "uniform") {
          prior.string2 <- paste0("dunif(-",x$hypparams[2],",",x$hypparams[2],")")
          prior.string22 <- paste0("dunif(0,",x$hypparams[4],")")
          }
     
     if(x$type[3] == "normal") 
          prior.string3 <- paste0("dnorm(0,", 1/x$hypparams[3], ")")
     if(x$type[3] == "cauchy") 
          prior.string3 <- paste0("dt(0,",1/x$hypparams[3],",1)")
     if(x$type[3] == "uniform") 
          prior.string3 <- paste0("dunif(-",x$hypparams[3],",",x$hypparams[3],")")

     if(x$type[4] == "uniform") 
          prior.string4 <- paste0("dunif(0,",x$hypparams[4],")")
     if(x$type[4] == "halfcauchy") 
          prior.string4 <- paste0("dt(0,",1/x$hypparams[4],",1)I(0,)")
     if(x$type[4] == "halfnormal") 
          prior.string4 <- paste0("dnorm(0,",1/x$hypparams[4],")I(0,)")
     #if(x$type[4] == "gamma") prior.string4 <- paste0("dgamma(",1/x$hypparams[4],",",1/x$hypparams[4],")")

     return(list(p1 = prior.string1, p2 = prior.string2, p22 = prior.string22, p3 = prior.string3, p4 = prior.string4))
     }


## Fill in empty components of prior.control and mcmc.control
fillin_prior_control <- function(x) {
     if(!("type" %in% names(x))) 
          x$type <- c("normal","normal","normal","uniform")
     if(!("hypparams" %in% names(x))) 
          x$hypparams <- c(10, 10, 10, 30)
     if(!("ssvs.index" %in% names(x))) 
          x$ssvs.index <- -1		
     if(!("ssvs.g" %in% names(x))) 
          x$ssvs.g <- 1e-6               
     if(!("ssvs.traitsindex" %in% names(x)))
          x$ssvs.traitsindex <- -1
          
     return(x)
     }

     
fillin_mcmc_control <- function(x) {
     if(!("n.burnin" %in% names(x))) 
          x$n.burnin <- 10000
     if(!("n.iteration" %in% names(x))) 
          x$n.iteration <- 40000
     if(!("n.thin" %in% names(x))) 
          x$n.thin <- 30
     if(!("seed" %in% names(x))) 
          x$seed <- NULL		
          
     return(x)
     }

     
## Given the lv, coefficients and cutoffs, return the multinomial probabilities for proportional odds regression for specific column of y
ordinal_conversion <- function(n, lv = NULL, lv.coefs.j = NULL, num.lv = NULL, 
     row.coefs = NULL, row.ids = NULL, X = NULL, X.coefs.j = NULL, offset.j = NULL, cutoffs, est = "median") {
	
     if(is.null(num.lv)) 
          num.lv <- 0
	
     etas <- matrix(0, nrow = n, ncol = length(cutoffs)) ## num.ord.levels - 1
     for(k in 1:length(cutoffs)) {
          etas[,k] <- rep(cutoffs[k],n) - lv.coefs.j[1]
          if(!is.null(lv)) 
               etas[,k] <- etas[,k] - as.matrix(lv)%*%lv.coefs.j[2:(num.lv+1)]
          if(!is.null(row.coefs)) { 
               if(est == "median") { 
                    for(k2 in 1:ncol(row.ids)) 
                         etas[,k] <- etas[,k] - row.coefs[[k2]]$median[row.ids[,k2]] 
                    } 
               if(est == "mean") { 
                    for(k2 in 1:ncol(row.ids)) 
                         etas[,k] <- etas[,k] - row.coefs[[k2]]$mean[row.ids[,k2]] 
                    }
               if(est == "ignore") { 
                    for(k2 in 1:ncol(row.ids)) 
                         etas[,k] <- etas[,k] - row.coefs[[k2]][row.ids[,k2]] 
                    }
               }
          if(!is.null(offset.j)) 
               etas[,k] <- etas[,k] - matrix(offset.j, ncol = 1)
          if(!is.null(X)) 
               etas[,k] <- etas[,k] - as.matrix(X)%*%X.coefs.j ## Don't forget the negative sign!
          }
          
     probs <- matrix(0,n,length(cutoffs)+1) ## num.ord.levels
     probs[,1] <- pnorm(etas[,1])
     for(k in 2:ncol(etas)) 
          probs[,k] <- pnorm(etas[,k]) - pnorm(etas[,k-1])
     probs[,length(cutoffs)+1] <- 1 - pnorm(etas[,length(cutoffs)])
     
     rm(etas)
     return(probs)
     }	


## Process Geweke's convergence diagnotics from a single chain MCMC fit	
process_geweke <- function(fit.mcmc, y, X = NULL, traits = NULL, family, num.lv, row.eff, row.ids, ranef.ids, 
     num.ord.levels = NULL, prior.control) { #type = "independent"
     
     p <- ncol(y)
     num.X <- 0
     if(!is.null(X)) 
          num.X <- ncol(X)
     num.traits <- 0
     if(!is.null(traits)) 
          num.traits <- ncol(traits)
        
     fit_geweke <- geweke.diag(fit.mcmc)[[1]] ## Takes first chain only
     out_gewekelist <- list(lv.coefs = matrix(fit_geweke[grep("lv.coefs", names(fit_geweke))], nrow = p))
     if(num.lv > 0) {
          out_gewekelist$lv.coefs <- matrix(out_gewekelist$lv.coefs[,-c(2:(num.lv+1))], nrow = p) ## Drop check on LV coefs
          fit_geweke <- fit_geweke[-grep("lvs",names(fit_geweke))] ## Drop check on lv
# 			if(type != "independent")
#                     out_gewekelist$lv.covparams <- fit_geweke[grep("lv.covparams",names(fit_geweke))] 
          }
     rownames(out_gewekelist$lv.coefs) <- colnames(y); 
     if(ncol(out_gewekelist$lv.coefs) > 1) {
        colnames(out_gewekelist$lv.coefs) <- c("beta0","Disperson") 
          } 
     else { 
          colnames(out_gewekelist$lv.coefs) <- c("beta0") 
          }

     if(row.eff != "none") {
          out_gewekelist$row.coefs <- vector("list", ncol(row.ids))
          names(out_gewekelist$row.coefs) <- colnames(row.ids)
          for(k in 1:ncol(row.ids))
               out_gewekelist$row.coefs[[k]] <- fit_geweke[grep(paste0("row.coefs.ID",k,"\\["), names(fit_geweke))]
          }
     if(row.eff == "random") {
          out_gewekelist$row.sigma <- vector("list", ncol(row.ids))
          names(out_gewekelist$row.sigma) <- colnames(row.ids)
          for(k in 1:ncol(row.ids))
               out_gewekelist$row.sigma[[k]] <- fit_geweke[grep(paste0("row.sigma.ID",k,"$"), names(fit_geweke))]
          }


     if(!is.null(ranef.ids)) {
          out_gewekelist$ranef.coefs <- vector("list", ncol(ranef.ids))
          names(out_gewekelist$ranef.coefs) <- colnames(ranef.ids)
          out_gewekelist$ranef.sigma <- matrix(NA, nrow = p, ncol = ncol(ranef.ids))
          rownames(out_gewekelist$ranef.sigma) <- colnames(y)
          colnames(out_gewekelist$ranef.sigma) <- colnames(ranef.ids)
          for(k0 in 1:ncol(ranef.ids)) {
               out_gewekelist$ranef.coefs[[k0]] <- matrix(fit_geweke[grep(paste0("ranef.coefs.ID",k0,"\\["), names(fit_geweke))], nrow = p)
               out_gewekelist$ranef.sigma[,k0] <- fit_geweke[grep(paste0("ranef.sigma.ID",k0,"\\["), names(fit_geweke))]
               }
          }

          
     if(num.X > 0) {
          out_gewekelist$X.coefs <- matrix(fit_geweke[grep("X.coefs", names(fit_geweke))], nrow = p)
          rownames(out_gewekelist$X.coefs) <- colnames(y); colnames(out_gewekelist$X.coefs) <- colnames(X)
          if(any(prior.control$ssvs.index > -1))
               out_gewekelist$X.coefs <- NULL ## Drop check on X coefs if SSVS is used
          }
            
     if(num.traits > 0) {
          out_gewekelist$traits.coefs <- cbind(fit_geweke[grep("traits.int", names(fit_geweke))],  
               matrix(fit_geweke[grep("traits.coefs", names(fit_geweke))], nrow = ncol(X)+1, ncol = ncol(traits)), 
               fit_geweke[grep("trait.sigma", names(fit_geweke))])
          rownames(out_gewekelist$traits.coefs) <- c("beta0",colnames(X)); 
          colnames(out_gewekelist$traits.coefs) <- c("kappa0",colnames(traits),"sigma")
          if(any(unlist(prior.control$ssvs.index) > -1))
               out_gewekelist$traits.coefs <- NULL ## Drop check on trait parameters if SSVS is used
          }

     if(any(family == "ordinal")) {
          out_gewekelist$cutoffs <- fit_geweke[grep("cutoffs", names(fit_geweke))]
          names(out_gewekelist$cutoffs) <- paste(1:(num.ord.levels - 1), "|", 2:num.ord.levels, sep = "") 

          if(sum(family == "ordinal") > 2)
               out_gewekelist$ordinal.sigma <- fit_geweke[grep("ordinal.sigma", names(fit_geweke))]
          }

     if(any(family == "tweedie"))
          out_gewekelist$powerparam <- fit_geweke[grep("powerparam", names(fit_geweke))]
    
    
     prop.outside <- table(2*pnorm(abs(unlist(out_gewekelist)), lower.tail = FALSE) < 0.05)/length(unlist(out_gewekelist))

     return(list(geweke.diag = out_gewekelist, prop.exceed = prop.outside))
     }


##-------------------
## Unused functions
##------------------
function() {
     ## Extract rhats from multiple chained MCMC fit	
     rhats <- function (x, asc = FALSE) {
          if(asc) { 
               x$BUGSoutput$summary[order(x$BUGSoutput$summary[, "Rhat"]), "Rhat", drop = FALSE] 
               }
          else { 
               x$BUGSoutput$summary[, "Rhat", drop = FALSE] 
               }
          }


# ## Process the rhats from multiple chained MCMC fit	
# process.rhats <- function(sims.matrix) {
#  		combined_fit_mcmc <- as.mcmc(sims.matrix)
# 		fit.rhats <- rhats(jagsfit, asc = FALSE)
#     		make.rhatslist <- list(lv.coefs = matrix(fit.rhats[grep("lv.coefs", rownames(fit.rhats))], nrow = p))
#     		#if(num.lv > 0) { fit.rhats <- fit.rhats[-grep("lvs",rownames(fit.rhats)),] } ## Drop check on lv
#     		if(num.lv > 0) { make.rhatslist$lv.coefs <- as.matrix(make.rhatslist$lv.coefs[,-c(2:(num.lv+1))]) } ## Drop check on LV coefs
#     		rownames(make.rhatslist$lv.coefs) <- colnames(y); 
#     		if(ncol(make.rhatslist$lv.coefs) > 1) { colnames(make.rhatslist$lv.coefs) <- c("beta0","Disperson") } 
# 		else { colnames(make.rhatslist$lv.coefs) <- c("beta0") }
# 
#     		if(row.eff != "none") {
#     			make.rhatslist$row.coefs <- fit.rhats[grep("row.coefs", rownames(fit.rhats))]		
#     			names(make.rhatslist$row.coefs) <- rownames(y) 
#     			}
#     		if(row.eff == "random") {
#     			make.rhatslist$row.sigma <- fit.rhats[grep("row.sigma.ID", rownames(fit.rhats))]
#    				names(make.rhatslist$row.sigma) <- c("Row random effects sigma") 
#     			}
# 
#    		if(num.X > 0) {
#     			make.rhatslist$X.coefs <- matrix(fit.rhats[grep("X.coefs", rownames(fit.rhats))], nrow = p)
#     			rownames(make.rhatslist$X.coefs) <- colnames(y); colnames(make.rhatslist$X.coefs) <- colnames(X) 
#     			}
# 		
# 		if(num.traits > 0) {
#     			make.rhatslist$traits.coefs <- cbind(fit.rhats[grep("trait.int", rownames(fit.rhats))], matrix(fit.rhats[grep("traits.coefs", rownames(fit.rhats))], nrow = ncol(X)+1), fit.rhats[grep("trait.sigma", rownames(fit.rhats))])
#     			rownames(make.rhatslist$traits.coefs) <- c("beta0",colnames(X)); colnames(make.rhatslist$traits.coefs) <- c("kappa0",colnames(traits),"sigma")
#     			}
#  
#    		if(any(family == "ordinal")) {
#    			make.rhatslist$cutoffs <- fit.rhats[grep("cutoffs", rownames(fit.rhats))]
#    			names(make.rhatslist$cutoffs) <- paste(1:(num.ord.levels - 1), "|", 2:num.ord.levels, sep = "") 
#  			if(sum(family == "ordinal") > 2) {
#  				make.rhatslist$ordinal.sigma <- fit.rhats[grep("ordinal.sigma", rownames(fit.rhats))]
#  				names(make.rhatslist$ordinal.sigma) <- "Rpecies-specific random intercept sigma for ordinal responses" 
#  				}
#    			}
# 
#    		if(any(family == "tweedie")) {
#    			make.rhatslist$powerparam <- fit.rhats[grep("powerparam", rownames(fit.rhats))]
#    			names(make.rhatslist$powerparam) <- "Common power parameter" 
#    			}
# 	
# 	return(make.rhatslist)
# 	}
	
	

}

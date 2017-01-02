###############################
## Unseen functions
##########################

## Given the lvs, coefficients and cutoffs, return the multinomial probabilities for proportional odds regression for specific element of y
ordinal.conversion.spp <- function(n, lv = NULL, lv.coefs.j = NULL, num.lv = NULL, row.coefs = NULL, row.ids = NULL, X = NULL, X.coefs.j = NULL, cutoffs, est = "median") {
	if(is.null(num.lv)) num.lv <- 0
	
	etas <- matrix(0, nrow = n, ncol = length(cutoffs)) ## num.ord.levels - 1
	for(k in 1:length(cutoffs)) {
		etas[,k] <- rep(cutoffs[k],n) - lv.coefs.j[1]
		if(!is.null(lv)) etas[,k] <- etas[,k] - as.matrix(lv)%*%lv.coefs.j[2:(num.lv+1)]
		if(!is.null(row.coefs)) { 
			if(est == "median") { for(k2 in 1:ncol(row.ids)) etas[,k] <- etas[,k] - row.coefs[[k2]]$median[row.ids[,k2]] } 
			if(est == "mean") { for(k2 in 1:ncol(row.ids)) etas[,k] <- etas[,k] - row.coefs[[k2]]$mean[row.ids[,k2]] }
			if(est == "ignore") { for(k2 in 1:ncol(row.ids)) etas[,k] <- etas[,k] - row.coefs[[k2]][row.ids[,k2]] }
			}
		if(!is.null(X)) etas[,k] <- etas[,k] - as.matrix(X)%*%X.coefs.j ## Don't forget the negative sign!
		}
	probs <- matrix(0,n,length(cutoffs)+1) ## num.ord.levels
	probs[,1] <- pnorm(etas[,1])
	for(k in 2:ncol(etas)) { probs[,k] <- pnorm(etas[,k]) - pnorm(etas[,k-1]) }		
	probs[,length(cutoffs)+1] <- 1 - pnorm(etas[,length(cutoffs)])

	rm(etas); return(probs)
	}	


## Extract rhats from multiple chained MCMC fit	
rhats <- function (x, asc = FALSE) {
	if(asc) { x$BUGSoutput$summary[order(x$BUGSoutput$summary[, "Rhat"]), "Rhat", drop = FALSE] }
	else { x$BUGSoutput$summary[, "Rhat", drop = FALSE] }
	}


# ## Process the rhats from multiple chained MCMC fit	
# process.rhats <- function(sims.matrix) {
#  		combined.fit.mcmc <- as.mcmc(sims.matrix)
# 		fit.rhats <- rhats(jagsfit, asc = FALSE)
#     		make.rhatslist <- list(lv.coefs = matrix(fit.rhats[grep("all.params", rownames(fit.rhats))], nrow = p))
#     		#if(num.lv > 0) { fit.rhats <- fit.rhats[-grep("lvs",rownames(fit.rhats)),] } ## Drop check on LVs
#     		if(num.lv > 0) { make.rhatslist$lv.coefs <- as.matrix(make.rhatslist$lv.coefs[,-c(2:(num.lv+1))]) } ## Drop check on LV coefs
#     		rownames(make.rhatslist$lv.coefs) <- colnames(y); 
#     		if(ncol(make.rhatslist$lv.coefs) > 1) { colnames(make.rhatslist$lv.coefs) <- c("beta0","Disperson") } 
# 		else { colnames(make.rhatslist$lv.coefs) <- c("beta0") }
# 
#     		if(row.eff != "none") {
#     			make.rhatslist$row.coefs <- fit.rhats[grep("row.params", rownames(fit.rhats))]		
#     			names(make.rhatslist$row.coefs) <- rownames(y) 
#     			}
#     		if(row.eff == "random") {
#     			make.rhatslist$row.sigma <- fit.rhats[grep("row.ranef.sigma", rownames(fit.rhats))]
#    				names(make.rhatslist$row.sigma) <- c("Row random effects sigma") 
#     			}
# 
#    		if(num.X > 0) {
#     			make.rhatslist$X.coefs <- matrix(fit.rhats[grep("X.params", rownames(fit.rhats))], nrow = p)
#     			rownames(make.rhatslist$X.coefs) <- colnames(y); colnames(make.rhatslist$X.coefs) <- colnames(X) 
#     			}
# 		
# 		if(num.traits > 0) {
#     			make.rhatslist$traits.coefs <- cbind(fit.rhats[grep("trait.int", rownames(fit.rhats))], matrix(fit.rhats[grep("traits.params", rownames(fit.rhats))], nrow = ncol(X)+1), fit.rhats[grep("sigma.trait", rownames(fit.rhats))])
#     			rownames(make.rhatslist$traits.coefs) <- c("beta0",colnames(X)); colnames(make.rhatslist$traits.coefs) <- c("kappa0",colnames(traits),"sigma")
#     			}
#  
#    		if(any(family == "ordinal")) {
#    			make.rhatslist$cutoffs <- fit.rhats[grep("alpha", rownames(fit.rhats))]
#    			names(make.rhatslist$cutoffs) <- paste(1:(num.ord.levels - 1), "|", 2:num.ord.levels, sep = "") 
#  			if(sum(family == "ordinal") > 2) {
#  				make.rhatslist$ordinal.sigma <- fit.rhats[grep("ordinal.ranef.sigma", rownames(fit.rhats))]
#  				names(make.rhatslist$ordinal.sigma) <- "Species-specific random intercept sigma for ordinal responses" 
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
	
	
## Sets up part of the JAGS script corresponding to family for responses; used in make.jagsboralmodel and make.jagsboralnullmodel. 
setup.resp.families.lv <- function(p, complete.family, num.lv, row.eff, row.ids, num.X, complete.trial.size, index.tweed.cols, index.ord.cols) {
	resp.family.script <- NULL
	
	for(j in 1:p) {			
		if(complete.family[j] != "multinom") {
			if(length(unique(complete.family)) == 1) {
				if(j == 1) {
					if(num.lv == 0) linpred.string <- paste("eta[i,j] <- 0", sep = "")
					if(num.lv > 0) linpred.string <- paste("eta[i,j] <- inprod(all.params[j,2:(num.lv+1)],lvs[i,])", sep = "")
					if(row.eff != "none") {
						for(k in 1:ncol(row.ids)) linpred.string <- paste(linpred.string, " + row.params.", colnames(row.ids)[k],"[row.ids[i,",k,"]]", sep ="")
						}
					if(num.X > 0) linpred.string <- paste(linpred.string, " + inprod(X.params[j,],X[i,])", sep ="")
					resp.family.script <- c(resp.family.script, paste("\t\t for(j in 1:p) { ", linpred.string, " }", sep = ""))
					}
				if(j > 1) { }
				}	
			if(length(unique(complete.family)) > 1) {		
				if(num.lv == 0) linpred.string <- paste("eta[i,",j, "] <- 0", sep = "")
				if(num.lv > 0) linpred.string <- paste("eta[i,",j, "] <- inprod(all.params[",j, ",2:(num.lv+1)],lvs[i,])", sep = "")
				if(row.eff != "none") {
					for(k in 1:ncol(row.ids)) linpred.string <- paste(linpred.string, " + row.params.", colnames(row.ids)[k],"[row.ids[i,",k,"]]", sep ="")
					}
				if(num.X > 0) linpred.string <- paste(linpred.string, " + inprod(X.params[", j, ",],X[i,])", sep ="")
				resp.family.script <- c(resp.family.script, paste("\t\t ", linpred.string, sep = ""))
				}				
			}			

			
		if(complete.family[j] == "negative.binomial") {
			if(length(unique(complete.family)) == 1) {
				if(j == 1) {
					#resp.family.script <- c(resp.family.script, paste("\t\t for(j in 1:p) { u[i,j] ~ dnorm(0, 1/all.params[",j, ",2]) }", sep = ""))
					resp.family.script <- c(resp.family.script, paste("\t\t for(j in 1:p) { u[i,j] ~ dgamma(1/all.params[j,num.lv+2], 1/all.params[j,num.lv+2]) }", sep = ""))
					resp.family.script <- c(resp.family.script, paste("\t\t for(j in 1:p) { y[i,j] ~ dpois(exp(all.params[j,1] + eta[i,j ])*(u[i,j])) } ## Parameterizing the NB as a multiplicative random effect models\n", sep = ""))
					}
				if(j > 1) { }
				}	
			if(length(unique(complete.family)) > 1) {		
				#resp.family.script <- c(resp.family.script, paste("\t\t u[i,",j, "] ~ dnorm(0, 1/all.params[",j, ",2])", sep = ""))
				resp.family.script <- c(resp.family.script, paste("\t\t u[i,",j, "] ~ dgamma(1/all.params[", j, ",num.lv+2], 1/all.params[",j, ",num.lv+2])", sep = ""))
				resp.family.script <- c(resp.family.script, paste("\t\t y[i,",j, "] ~ dpois(exp(all.params[", j, ",1] + eta[i,",j, "])*(u[i,", j, "])) ## Parameterizing the NB as a multiplicative random effect models, with size\n", sep = ""))
				}
			}
			
		if(complete.family[j] == "negative.binomial2") {
			if(length(unique(complete.family)) == 1) {
				if(j == 1) {
					resp.family.script <- c(resp.family.script, paste("\t\t for(j in 1:p) { u[i,j] <- 1/(1 + all.params[j,num.lv+2]*exp(all.params[j,1] + eta[i,j])) }", sep = ""))
					resp.family.script <- c(resp.family.script, paste("\t\t for(j in 1:p) { y[i,j] ~ dnegbin(u[i,j],1/all.params[j,num.lv+2]) } ## Parameterizing the NB as a function of prob and size \n", sep = ""))
					}
				if(j > 1) { }
				}
			if(length(unique(complete.family)) > 1) {		
				resp.family.script <- c(resp.family.script, paste("\t\t u[i,",j, "] <- 1/(1 + all.params[",j, ",num.lv+2]*exp(all.params[", j, ",1] + eta[i,",j, "]))", sep = ""))
				resp.family.script <- c(resp.family.script, paste("\t\t y[i,",j, "] ~ dnegbin(u[i,", j, "],1/all.params[",j, ",num.lv+2]) ## Parameterizing the NB as a function of prob and size \n", sep = ""))
				}
			}
			
		if(complete.family[j] == "normal") {
			if(length(unique(complete.family)) == 1) {
				if(j == 1) {
					resp.family.script <- c(resp.family.script, paste("\t\t for(j in 1:p) { y[i,j] ~ dnorm(all.params[j,1] + eta[i,j],pow(all.params[j,num.lv+2],-2)) } \n", sep = ""))
					}
				if(j > 1) { }
				}
			if(length(unique(complete.family)) > 1) {		
				resp.family.script <- c(resp.family.script, paste("\t\t y[i,", j, "] ~ dnorm(all.params[", j, ",1] + eta[i,", j, "],pow(all.params[", j, ",num.lv+2],-2)) \n", sep = ""))
				}
			}

		if(all(complete.family == "bernoulli")) { ## If all data are Bernoulli, then use step parameterization
			if(length(unique(complete.family)) == 1) {
				if(j == 1) {
					resp.family.script <- c(resp.family.script, paste("\t\t for(j in 1:p) { Z[i,j] ~ dnorm(all.params[j,1] + eta[i,j],(1-sum(all.params[",j,",2:(num.lv+1)]^2))) }", sep = ""))
					resp.family.script <- c(resp.family.script, paste("\t\t for(j in 1:p) { y[i,j] ~ dbern(step(Z[i,j])) }\n", sep = ""))
					}
				if(j > 1) { }
				}
			if(length(unique(complete.family)) > 1) {		
				resp.family.script <- c(resp.family.script, paste("\t\t Z[i,",j, "] ~ dnorm(all.params[", j, ",1] + eta[i,",j, "],(1-sum(all.params[",j,",2:(num.lv+1)]^2)))", sep = ""))
				resp.family.script <- c(resp.family.script, paste("\t\t y[i,",j, "] ~ dbern(step(Z[i,",j, "]))\n", sep = ""))
				}
			}
			
		if(complete.family[j] == "binomial") {
			resp.family.script <- c(resp.family.script, paste("\t\t y[i,",j, "] ~ dbin(phi(all.params[", j, ",1] + eta[i,",j, "]),",complete.trial.size[j],")\n", sep = ""))
			}
			
		if(complete.family[j] == "exponential") {
			if(length(unique(complete.family)) == 1) {
				if(j == 1) {
					resp.family.script <- c(resp.family.script, paste("\t\t for(j in 1:p) { y[i,j] ~ dexp(pow(exp(all.params[j,1] + eta[i,j]),-1)) }\n", sep = ""))
					}
				if(j > 1) { }
				}
			if(length(unique(complete.family)) > 1) {		
				resp.family.script <- c(resp.family.script, paste("\t\t y[i,", j, "] ~ dexp(pow(exp(all.params[", j, ",1] + eta[i,", j, "]),-1))\n", sep = ""))
				}
			}
			
		if(complete.family[j] == "gamma") {
			if(length(unique(complete.family)) == 1) {
				if(j == 1) {
					resp.family.script <- c(resp.family.script, paste("\t\t for(j in 1:p) { y[i,j] ~ dgamma(exp(all.params[j,1] + eta[i,j])*all.params[j,num.lv+2], all.params[j,num.lv+2]) } \n", sep = ""))
					}
				if(j > 1) { }
				}
			if(length(unique(complete.family)) > 1) {		
				resp.family.script <- c(resp.family.script, paste("\t\t y[i,", j, "] ~ dgamma(exp(all.params[", j, ",1] + eta[i,", j, "])*all.params[", j, ",num.lv+2], all.params[", j, ",num.lv+2])\n", sep = ""))
				}
			}
			
		if(complete.family[j] == "beta") {
			if(length(unique(complete.family)) == 1) {
				if(j == 1) {
					resp.family.script <- c(resp.family.script, paste("\t\t for(j in 1:p) { y[i,j] ~ dbeta(ilogit(all.params[j,1] + eta[i,j])*all.params[j,num.lv+2],(1-ilogit(all.params[j,1] + eta[i,j]))*all.params[j,num.lv+2]) }\n", sep = ""))
					}
				if(j > 1) { }
				}
			if(length(unique(complete.family)) > 1) {		
				resp.family.script <- c(resp.family.script, paste("\t\t y[i,", j, "] ~ dbeta(ilogit(all.params[", j, ",1] + eta[i,", j, "])*all.params[", j, ",num.lv+2],(1-ilogit(all.params[", j, ",1] + eta[i,", j, "]))*all.params[", j, ",num.lv+2])\n", sep = ""))
				}
			}
			
		if(complete.family[j] == "poisson") {
			if(length(unique(complete.family)) == 1) {
				if(j == 1) {
					resp.family.script <- c(resp.family.script, paste("\t\t for(j in 1:p) { y[i,j] ~ dpois(exp(all.params[j,1] + eta[i,j])) }\n", sep = ""))
					}
				if(j > 1) { }
				}
			if(length(unique(complete.family)) > 1) {		
				resp.family.script <- c(resp.family.script, paste("\t\t y[i,", j, "] ~ dpois(exp(all.params[", j, ",1] + eta[i,", j, "]))\n", sep = ""))
				}
			}	
			
		if(complete.family[j] == "lnormal") {
			if(length(unique(complete.family)) == 1) {
				if(j == 1) {
					resp.family.script <- c(resp.family.script, paste("\t\t for(j in 1:p) { y[i,j] ~ dlnorm(all.params[j,1] + eta[i,j],pow(all.params[j,num.lv+2],-2)) } \n", sep = ""))
					}
				if(j > 1) { }
				}
			if(length(unique(complete.family)) > 1) {		
				resp.family.script <- c(resp.family.script, paste("\t\t y[i,", j, "] ~ dlnorm(all.params[", j, ",1] + eta[i,", j, "],pow(all.params[", j, ",num.lv+2],-2)) \n", sep = ""))
				}
			}	
			
		if(complete.family[j] == "tweedie") {
			resp.family.script <- c(resp.family.script, paste("\t\t lambdanum[i,", which(index.tweed.cols == j), "] <- pow(exp(all.params[", j, ",1] + eta[i,", j, "]),2-powerparam)/(all.params[", j, ",num.lv+2]*(2-powerparam))", sep = ""))
			resp.family.script <- c(resp.family.script, paste("\t\t numfish[i,", which(index.tweed.cols == j), "] ~ dpois(lambdanum[i,", which(index.tweed.cols == j), "])", sep = ""))
			resp.family.script <- c(resp.family.script, paste("\t\t choose.shape[i,", which(index.tweed.cols == j), ",1] <- numfish[i,", which(index.tweed.cols == j), "]*(2-powerparam)/(powerparam-1)", sep = "")) ## If y > 0, then conditional on numfish, y is sum of independent gammas
			resp.family.script <- c(resp.family.script, paste("\t\t choose.rate[i,", which(index.tweed.cols == j), ",1] <- 1/(all.params[", j, ",num.lv+2]*(powerparam-1)*pow(exp(all.params[", j, ",1] + eta[i,", j, "]),powerparam-1))", sep = ""))
			resp.family.script <- c(resp.family.script, paste("\t\t choose.shape[i,", which(index.tweed.cols == j), ",2] <- 1", sep = "")) ## If y = 0, then Tweedie dist equals probability of Poisson equal 0
			resp.family.script <- c(resp.family.script, paste("\t\t choose.rate[i,", which(index.tweed.cols == j), ",2] <- exp(-lambdanum[i,", which(index.tweed.cols == j), "])", sep = ""))
			resp.family.script <- c(resp.family.script, paste("\t\t y[i,", j, "] ~ dgamma(choose.shape[i,", which(index.tweed.cols == j), ",1+equals(y[i,", which(index.tweed.cols == j), "],0)],choose.rate[i,", j, ",1+equals(y[i,", j, "],0)]) \n", sep = ""))
			}
			
		if(complete.family[j] == "ordinal") {
			if(length(index.ord.cols) == p) {
				if(j == 1) { 
					resp.family.script <- c(resp.family.script, paste("\t\t for(j in 1:p) { \n\t\t\t prob[i,j,1] <- phi(alpha[1]-eta[i,j]-all.params[j,1])", sep = ""))
					resp.family.script <- c(resp.family.script, paste("\t\t\t for(k in 2:(num.ord.levels-1)) { \n\t\t\t\t prob[i,j,k] <- phi(alpha[k]-eta[i,j]-all.params[j,1]) - phi(alpha[k-1]-eta[i,j]-all.params[j,1]) \n\t\t\t\t }", sep = ""))
					resp.family.script <- c(resp.family.script, paste("\t\t\t prob[i,j,num.ord.levels] <- 1-phi(alpha[num.ord.levels-1]-eta[i,j]-all.params[j,1])", sep = ""))
					resp.family.script <- c(resp.family.script, paste("\t\t\t y[i,j] ~ dcat(prob[i,j,]) \n\t\t\t } \n", sep = ""))
					}
				if(j > 1) { }
				}
			if(length(index.ord.cols) < p) {
				resp.family.script <- c(resp.family.script, paste("\t\t prob[i,", which(index.ord.cols == j), ",1] <- phi(alpha[1]-eta[i,", j, "]-all.params[", j, ",1])", sep = ""))
				resp.family.script <- c(resp.family.script, paste("\t\t for(k in 2:(num.ord.levels-1)) { prob[i,", which(index.ord.cols == j), ",k] <- phi(alpha[k]-eta[i,", j, "]-all.params[", j, ",1]) - phi(alpha[k-1]-eta[i,", j, "]-all.params[", j, ",1]) }", sep = ""))
				resp.family.script <- c(resp.family.script, paste("\t\t prob[i,", which(index.ord.cols == j), ",num.ord.levels] <- 1-phi(alpha[num.ord.levels-1]-eta[i,", j, "]-all.params[", j, ",1])", sep = ""))
				resp.family.script <- c(resp.family.script, paste("\t\t y[i,", j, "] ~ dcat(prob[i,", which(index.ord.cols == j), ",])\n", sep = ""))
				}
			}
			
		if(complete.family[j] == "multinom") { 
			stop("You shouldn't have gotten here!") ## Coefficients for LVs are constrained to be same for all levels! Otherwise identifiability constraints are hard!
#    		mod.general.lv <- c(mod.general.lv, paste("\t\t for(k in 1:num.multinom.levels[",j,"]) {",sep=""))
# 			if(num.X == 0 & row.eff) mod.general.lv <- c(mod.general.lv, paste("\t\t\t mu[i,",which(index.multinom.cols == j),",k] <- exp(row.params[i] + inprod(all.params[",j,",2:(num.lv+1)],lvs[i,]))",sep="")) 
# 			if(num.X > 0 & row.eff) mod.general.lv <- c(mod.general.lv, paste("\t\t\t mu[i,",which(index.multinom.cols == j),",k] <- exp(row.params[i] + inprod(all.params[",j,",2:(num.lv+1)],lvs[i,]) + inprod(X.multinom.params[",which(index.multinom.cols == j),",,k],X[i,]))",sep="")) 
# 			if(num.X == 0 & !row.eff) mod.general.lv <- c(mod.general.lv, paste("\t\t\t mu[i,",which(index.multinom.cols == j),",k] <- exp(inprod(all.params[",j,",2:(num.lv+1)],lvs[i,]))",sep="")) 
# 			if(num.X == 0 & !row.eff) mod.general.lv <- c(mod.general.lv, paste("\t\t\t mu[i,",which(index.multinom.cols == j),",k] <- exp(inprod(all.params[",j,",2:(num.lv+1)],lvs[i,]) + inprod(X.multinom.params[",which(index.multinom.cols == j),",,k],X[i,]))",sep="")) 			
# 			mod.general.lv <- c(mod.general.lv, paste("\t\t\t prob[i,",which(index.multinom.cols == j),",k] <- mu[i,",which(index.multinom.cols == j),",k]/sum(mu[i,",which(index.multinom.cols == j),",]) }",sep="")) 
# 			mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,",j,"] ~ dcat(prob[i,",which(index.multinom.cols == j),",]+0.001)\n",sep="")) 
			}		
		}
		
	return(resp.family.script)
	}
	
	
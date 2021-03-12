###############################
## Sets up part of the JAGS script corresponding to family for responses; used in make.jagsboralmodel and make.jagsboralnullmodel. 
## Unseen function
###############################

setup_respfamilies <- function(p, complete.family, num.lv, row.eff, row.ids, ranef.ids,
     offset, num.X, complete.trial.size, index.tweed.cols, index.ord.cols) {

     respfamily_script <- NULL
	
     for(j in 1:p) {
          if(complete.family[j] != "multinom") {
               if(length(unique(complete.family)) == 1) {
                    if(j == 1) {
                         if(num.lv == 0) 
                              linpred_string <- paste0("eta[i,j] <- 0")
                         if(num.lv > 0) 
                              linpred_string <- paste0("eta[i,j] <- inprod(lv.coefs[j,2:(num.lv+1)],lvs[i,])")
                         if(row.eff != "none") {
                              for(k in 1:ncol(row.ids)) 
                                   linpred_string <- paste0(linpred_string, " + row.coefs.ID",k,"[row.ids[i,",k,"]]")
                              }
                         if(!is.null(ranef.ids)) {
                              for(k0 in 1:ncol(ranef.ids)) 
                                   linpred_string <- paste0(linpred_string, " + ranef.coefs.ID",k0,"[j,ranef.ids[i,",k0,"]]")
                              }
                         if(num.X > 0) 
                              linpred_string <- paste0(linpred_string, " + inprod(X.coefs[j,],X[i,])")
                         if(!is.null(offset)) 
                              linpred_string <- paste0(linpred_string, " + offset[i,j]")
                         respfamily_script <- c(respfamily_script, paste0("\t\t for(j in 1:p) { ", linpred_string, " }"))
                         }
                    if(j > 1) { }
                    }	
               if(length(unique(complete.family)) > 1) {
                    if(num.lv == 0) 
                         linpred_string <- paste0("eta[i,",j, "] <- 0")
                    if(num.lv > 0) 
                         linpred_string <- paste0("eta[i,",j, "] <- inprod(lv.coefs[",j, ",2:(num.lv+1)],lvs[i,])")
                    if(row.eff != "none") 
                         {
                         for(k in 1:ncol(row.ids)) 
                              linpred_string <- paste0(linpred_string, " + row.coefs.ID",k,"[row.ids[i,",k,"]]")
                         }
                    if(!is.null(ranef.ids)) {
                         for(k0 in 1:ncol(ranef.ids)) 
                              linpred_string <- paste0(linpred_string, " + ranef.coefs.ID",k0,"[",j,",ranef.ids[i,",k0,"]]")
                         }
                    if(num.X > 0) 
                         linpred_string <- paste0(linpred_string, " + inprod(X.coefs[",j, ",],X[i,])")
                    if(!is.null(offset)) 
                         linpred_string <- paste0(linpred_string, " + offset[i,",j,"]")
                    respfamily_script <- c(respfamily_script, paste0("\t\t ", linpred_string))
                    }  
               }			
                
                
                
          if(complete.family[j] == "negative.binomial") {
               if(length(unique(complete.family)) == 1) {
                    if(j == 1) {
                         #respfamily_script <- c(respfamily_script, paste0("\t\t for(j in 1:p) { u[i,j] ~ dnorm(0, 1/lv.coefs[",j, ",2]) }"))
                         respfamily_script <- c(respfamily_script, paste0("\t\t for(j in 1:p) { u[i,j] ~ dgamma(1/lv.coefs[j,num.lv+2], 1/lv.coefs[j,num.lv+2]) }"))
                         respfamily_script <- c(respfamily_script, paste0("\t\t for(j in 1:p) { y[i,j] ~ dpois(exp(lv.coefs[j,1] + eta[i,j ])*(u[i,j])) } ## Parameterizing the NB as a multiplicative random effect models\n"))
                         }
                    if(j > 1) { }
                    }	
               if(length(unique(complete.family)) > 1) {
                    #respfamily_script <- c(respfamily_script, paste0("\t\t u[i,",j, "] ~ dnorm(0, 1/lv.coefs[",j, ",2])"))
                    respfamily_script <- c(respfamily_script, paste0("\t\t u[i,",j, "] ~ dgamma(1/lv.coefs[",j, ",num.lv+2], 1/lv.coefs[",j, ",num.lv+2])"))
                    respfamily_script <- c(respfamily_script, paste0("\t\t y[i,",j, "] ~ dpois(exp(lv.coefs[",j, ",1] + eta[i,",j, "])*(u[i,",j, "])) ## Parameterizing the NB as a multiplicative random effect models, with size\n"))
                    }
               }
                
		if(complete.family[j] == "ztnegative.binomial") {
			if(length(unique(complete.family)) == 1) {
				if(j == 1) {
					respfamily_script <- c(respfamily_script, paste0("\t\t for(j in 1:p) { u[i,j] <- 1/(1 + lv.coefs[j,num.lv+2]*exp(lv.coefs[j,1] + eta[i,j])) }"))
					respfamily_script <- c(respfamily_script, paste0("\t\t for(j in 1:p) { y[i,j] ~ dnegbin(u[i,j], 1/lv.coefs[j,num.lv+2]) T(1,) } \n"))
					}
				if(j > 1) { }
				}
			if(length(unique(complete.family)) > 1) {		
				respfamily_script <- c(respfamily_script, paste0("\t\t u[i,",j, "] <- 1/(1 + lv.coefs[",j, ",num.lv+2]*exp(lv.coefs[",j, ",1] + eta[i,",j, "]))"))
				respfamily_script <- c(respfamily_script, paste0("\t\t y[i,",j, "] ~ dnegbin(u[i,",j, "], 1/lv.coefs[",j, ",num.lv+2]) T(1,) \n"))
				}
			}
                
          if(complete.family[j] == "normal") {
               if(length(unique(complete.family)) == 1) {
                    if(j == 1) {
                         respfamily_script <- c(respfamily_script, paste0("\t\t for(j in 1:p) { y[i,j] ~ dnorm(lv.coefs[j,1] + eta[i,j],pow(lv.coefs[j,num.lv+2],-2)) } \n"))
                         }
                    if(j > 1) { }
                    }
               if(length(unique(complete.family)) > 1) {
                    respfamily_script <- c(respfamily_script, paste0("\t\t y[i,",j, "] ~ dnorm(lv.coefs[",j, ",1] + eta[i,",j, "],pow(lv.coefs[",j, ",num.lv+2],-2)) \n"))
                    }
               }

          if(complete.family[j] == "binomial") {
               respfamily_script <- c(respfamily_script, paste0("\t\t y[i,",j, "] ~ dbin(phi(lv.coefs[",j, ",1] + eta[i,",j, "]),",complete.trial.size[j],")\n"))
               }
                
          if(complete.family[j] == "exponential") {
               if(length(unique(complete.family)) == 1) {
                    if(j == 1) 
                         {
                         respfamily_script <- c(respfamily_script, paste0("\t\t for(j in 1:p) { y[i,j] ~ dexp(pow(exp(lv.coefs[j,1] + eta[i,j]),-1)) }\n"))
                         }
                    if(j > 1) { }
                    }
               if(length(unique(complete.family)) > 1) {		
                    respfamily_script <- c(respfamily_script, paste0("\t\t y[i,",j, "] ~ dexp(pow(exp(lv.coefs[",j, ",1] + eta[i,",j, "]),-1))\n"))
                    }
               }
                
          if(complete.family[j] == "gamma") {
               if(length(unique(complete.family)) == 1) {
                    if(j == 1) {
                         respfamily_script <- c(respfamily_script, paste0("\t\t for(j in 1:p) { y[i,j] ~ dgamma(exp(lv.coefs[j,1] + eta[i,j])*lv.coefs[j,num.lv+2], lv.coefs[j,num.lv+2]) } \n"))
                         }
                    if(j > 1) { }
                    }
               if(length(unique(complete.family)) > 1) {		
                    respfamily_script <- c(respfamily_script, paste0("\t\t y[i,",j, "] ~ dgamma(exp(lv.coefs[",j, ",1] + eta[i,",j, "])*lv.coefs[",j, ",num.lv+2], lv.coefs[",j, ",num.lv+2])\n"))
                    }
               }
                
          if(complete.family[j] == "beta") {
               if(length(unique(complete.family)) == 1) {
                    if(j == 1) 
                         {
                         respfamily_script <- c(respfamily_script, paste0("\t\t for(j in 1:p) { y[i,j] ~ dbeta(ilogit(lv.coefs[j,1] + eta[i,j])*lv.coefs[j,num.lv+2],(1-ilogit(lv.coefs[j,1] + eta[i,j]))*lv.coefs[j,num.lv+2]) }\n"))
                         }
                    if(j > 1) { }
                    }
               if(length(unique(complete.family)) > 1) {		
                    respfamily_script <- c(respfamily_script, paste0("\t\t y[i,",j, "] ~ dbeta(ilogit(lv.coefs[",j, ",1] + eta[i,",j, "])*lv.coefs[",j, ",num.lv+2],(1-ilogit(lv.coefs[",j, ",1] + eta[i,",j, "]))*lv.coefs[",j, ",num.lv+2])\n"))
                    }
               }
                
          if(complete.family[j] == "poisson") {
               if(length(unique(complete.family)) == 1) {
                    if(j == 1) {
                         respfamily_script <- c(respfamily_script, paste0("\t\t for(j in 1:p) { y[i,j] ~ dpois(exp(lv.coefs[j,1] + eta[i,j])) }\n"))
                         }
                    if(j > 1) { }
                    }
               if(length(unique(complete.family)) > 1) 
                    {		
                    respfamily_script <- c(respfamily_script, paste0("\t\t y[i,",j, "] ~ dpois(exp(lv.coefs[",j, ",1] + eta[i,",j, "]))\n"))
                    }
               }	
                
          if(complete.family[j] == "ztpoisson") {
               if(length(unique(complete.family)) == 1) {
                    if(j == 1) {
                         respfamily_script <- c(respfamily_script, paste0("\t\t for(j in 1:p) { y[i,j] ~ dpois(exp(lv.coefs[j,1] + eta[i,j])) T(1,) }\n"))
                         }
                    if(j > 1) { }
                    }
               if(length(unique(complete.family)) > 1) 
                    {		
                    respfamily_script <- c(respfamily_script, paste0("\t\t y[i,",j, "] ~ dpois(exp(lv.coefs[",j, ",1] + eta[i,",j, "])) T(1,) \n"))
                    }
               }	

          if(complete.family[j] == "lnormal") {
               if(length(unique(complete.family)) == 1) {
                    if(j == 1) {
                         respfamily_script <- c(respfamily_script, paste0("\t\t for(j in 1:p) { y[i,j] ~ dlnorm(lv.coefs[j,1] + eta[i,j],pow(lv.coefs[j,num.lv+2],-2)) } \n"))
                         }
                    if(j > 1) { }
                    }
               if(length(unique(complete.family)) > 1) {		
                    respfamily_script <- c(respfamily_script, paste0("\t\t y[i,",j, "] ~ dlnorm(lv.coefs[",j, ",1] + eta[i,",j, "],pow(lv.coefs[",j, ",num.lv+2],-2)) \n"))
                    }
               }  	
                
          if(complete.family[j] == "tweedie") {
               respfamily_script <- c(respfamily_script, paste0("\t\t lambdanum[i,", which(index.tweed.cols == j), "] <- pow(exp(lv.coefs[",j, ",1] + eta[i,",j, "]),2-powerparam)/(lv.coefs[",j, ",num.lv+2]*(2-powerparam))"))
               respfamily_script <- c(respfamily_script, paste0("\t\t numfish[i,", which(index.tweed.cols == j), "] ~ dpois(lambdanum[i,", which(index.tweed.cols == j), "])"))
               respfamily_script <- c(respfamily_script, paste0("\t\t choose.shape[i,", which(index.tweed.cols == j), ",1] <- numfish[i,", which(index.tweed.cols == j), "]*(2-powerparam)/(powerparam-1)")) ## If y > 0, then conditional on numfish, y is sum of independent gammas
               respfamily_script <- c(respfamily_script, paste0("\t\t choose.rate[i,", which(index.tweed.cols == j), ",1] <- 1/(lv.coefs[",j, ",num.lv+2]*(powerparam-1)*pow(exp(lv.coefs[",j, ",1] + eta[i,",j, "]),powerparam-1))"))
               respfamily_script <- c(respfamily_script, paste0("\t\t choose.shape[i,", which(index.tweed.cols == j), ",2] <- 1")) ## If y = 0, then Tweedie dist equals probability of Poisson equal 0
               respfamily_script <- c(respfamily_script, paste0("\t\t choose.rate[i,", which(index.tweed.cols == j), ",2] <- exp(-lambdanum[i,", which(index.tweed.cols == j), "])"))
               respfamily_script <- c(respfamily_script, paste0("\t\t y[i,",j, "] ~ dgamma(choose.shape[i,", which(index.tweed.cols == j), ",1+equals(y[i,", which(index.tweed.cols == j), "],0)],choose.rate[i,",j, ",1+equals(y[i,",j, "],0)]) \n"))
               }
                
          if(complete.family[j] == "ordinal") {
               if(length(index.ord.cols) == p) {
                    if(j == 1) { 
                         respfamily_script <- c(respfamily_script, paste0("\t\t for(j in 1:p) { \n\t\t\t prob[i,j,1] <- phi(cutoffs[1]-eta[i,j]-lv.coefs[j,1])"))
                         respfamily_script <- c(respfamily_script, paste0("\t\t\t for(k in 2:(num.ord.levels-1)) { \n\t\t\t\t prob[i,j,k] <- phi(cutoffs[k]-eta[i,j]-lv.coefs[j,1]) - phi(cutoffs[k-1]-eta[i,j]-lv.coefs[j,1]) \n\t\t\t\t }"))
                         respfamily_script <- c(respfamily_script, paste0("\t\t\t prob[i,j,num.ord.levels] <- 1-phi(cutoffs[num.ord.levels-1]-eta[i,j]-lv.coefs[j,1])"))
                         respfamily_script <- c(respfamily_script, paste0("\t\t\t y[i,j] ~ dcat(prob[i,j,]) \n\t\t\t } \n"))
                         }
                    if(j > 1) { }
                    }
               if(length(index.ord.cols) < p) {
                    respfamily_script <- c(respfamily_script, paste0("\t\t prob[i,", which(index.ord.cols == j), ",1] <- phi(cutoffs[1]-eta[i,",j, "]-lv.coefs[",j, ",1])"))
                    respfamily_script <- c(respfamily_script, paste0("\t\t for(k in 2:(num.ord.levels-1)) { prob[i,", which(index.ord.cols == j), ",k] <- phi(cutoffs[k]-eta[i,",j, "]-lv.coefs[",j, ",1]) - phi(cutoffs[k-1]-eta[i,",j, "]-lv.coefs[",j, ",1]) }"))
                    respfamily_script <- c(respfamily_script, paste0("\t\t prob[i,", which(index.ord.cols == j), ",num.ord.levels] <- 1-phi(cutoffs[num.ord.levels-1]-eta[i,",j, "]-lv.coefs[",j, ",1])"))
                    respfamily_script <- c(respfamily_script, paste0("\t\t y[i,",j, "] ~ dcat(prob[i,", which(index.ord.cols == j), ",])\n"))
                    }
               }
                
          if(complete.family[j] == "multinom") { 
               stop("You shouldn't have gotten here!") ## Coefficients for lv are constrained to be same for all levels! Otherwise identifiability constraints are hard!
#    		model_script <- c(model_script, paste0("\t\t for(k in 1:num.multinom.levels[",j,"]) {",sep=""))
# 			if(num.X == 0 & row.eff) model_script <- c(model_script, paste0("\t\t\t mu[i,",which(index_multinom_cols == j),",k] <- exp(row.coefs[i] + inprod(lv.coefs[",j,",2:(num.lv+1)],lvs[i,]))",sep="")) 
# 			if(num.X > 0 & row.eff) model_script <- c(model_script, paste0("\t\t\t mu[i,",which(index_multinom_cols == j),",k] <- exp(row.coefs[i] + inprod(lv.coefs[",j,",2:(num.lv+1)],lvs[i,]) + inprod(X.multinom.params[",which(index_multinom_cols == j),",,k],X[i,]))",sep="")) 
# 			if(num.X == 0 & !row.eff) model_script <- c(model_script, paste0("\t\t\t mu[i,",which(index_multinom_cols == j),",k] <- exp(inprod(lv.coefs[",j,",2:(num.lv+1)],lvs[i,]))",sep="")) 
# 			if(num.X > 0 & !row.eff) model_script <- c(model_script, paste0("\t\t\t mu[i,",which(index_multinom_cols == j),",k] <- exp(inprod(lv.coefs[",j,",2:(num.lv+1)],lvs[i,]) + inprod(X.multinom.params[",which(index_multinom_cols == j),",,k],X[i,]))",sep="")) 			
# 			model_script <- c(model_script, paste0("\t\t\t prob[i,",which(index_multinom_cols == j),",k] <- mu[i,",which(index_multinom_cols == j),",k]/sum(mu[i,",which(index_multinom_cols == j),",]) }",sep="")) 
# 			model_script <- c(model_script, paste0("\t\t y[i,",j,"] ~ dcat(prob[i,",which(index_multinom_cols == j),",]+0.001)\n",sep="")) 
               }		
          }
            
     return(respfamily_script)
     }
	

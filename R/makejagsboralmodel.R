make.jagsboralmodel <- function(family, num.X = 0, X.ind = NULL, num.traits = 0, which.traits = NULL,  
     lv.control = list(num.lv = 2, type = "independent"), 
     row.eff = "none", row.ids = NULL, ranef.ids = NULL, 
     offset = NULL, trial.size = 1, n, p, model.name = NULL, 
     prior.control = list(type = c("normal","normal","normal","uniform"), hypparams = c(10, 10, 10, 30), ssvs.index = -1, ssvs.g = 1e-6, ssvs.traitsindex = -1), num.lv = NULL) {
               
     check_which_traits(num.traits = num.traits, which.traits = which.traits, num.X = num.X, makejagsboralfile_messages = TRUE)
     if(is.null(which.traits)) { 
          which.traits <- vector("list",num.X+1)
          for(k in 1:(num.X+1)) 
               which.traits[[k]] <- 0 
          }           

     lv.control <- check_lv_control(num.lv = num.lv, lv.control = lv.control, need.distmat = FALSE)
     num.lv <- lv.control$num.lv
    
     complete_trial_size <- check_trial_size(family = family, trial.size = trial.size, create.complete.trial.size = TRUE, y = matrix(NA,nrow=1,ncol=p))

     complete_family <- check_family(family = family, y = matrix(1,nrow=1,ncol=p), traits = NULL) ## Done just to produce the complete_family vector
    
     if(row.eff != "none" && is.null(row.ids)) {
          row.ids <- matrix(1:n, ncol = 1)
          message("row.ids assumed to be matrix with one column and elements 1,2,...n i.e., a row-specific intercept.")
          }
     if(!is.null(row.ids)) {
          row.ids <- as.matrix(row.ids)
          if(nrow(row.ids) != n) 
               stop("Number of rows in the matrix row.ids should be equal to n.")
          if(is.null(colnames(row.ids))) 
               colnames(row.ids) <- paste0("ID", 1:ncol(row.ids))
          }

     
     ranef.ids <- check_ranef_ids(ranef.ids = ranef.ids, y = matrix(0,n,p))

     if(!is.null(offset)) { 
          if(!is.matrix(offset)) 
               stop("offset could be a matrix with the same dimensions as y.")
          }
            
               
     prior.control <- fillin_prior_control(x = prior.control)
     check_prior_control(prior.control = prior.control)
    
     if(length(prior.control$ssvs.index) == 1 & num.X > 0) 
          prior.control$ssvs.index <- rep(prior.control$ssvs.index, num.X)
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
        
     X.ind <- check_X_ind(X.ind = X.ind, p = p, num.X = num.X, prior.control = prior.control)
        
     index.ord.cols <- which(complete_family == "ordinal")
     index.tweed.cols <- which(complete_family == "tweedie")
    
    
     ##-------------------
     ## Checks done; starting writing JAGS script!
     ##-------------------

     
     model_script <- paste0("## JAGS model written for boral version ", packageDescription("boral")$Version, " on ", as.character(Sys.time()), " ##\n\n model {")
     model_script <- c(model_script, "\t ## Data Level ## \n\t for(i in 1:n) {")
    
     write.resp.script <- setup_respfamilies(p = p, complete.family = complete_family, num.lv = num.lv, row.eff = row.eff, row.ids = row.ids, ranef.ids = ranef.ids,
          offset = offset, num.X = num.X, complete.trial.size = complete_trial_size, index.tweed.cols = index.tweed.cols, index.ord.cols = index.ord.cols)
     model_script <- c(model_script, write.resp.script)	
     model_script <- c(model_script, paste0("\t\t }"))
     rm(write.resp.script)

     
     model_script <- c(model_script, paste0("\t ## Latent variables ##"))    
     if(lv.control$type == "independent")
          model_script <- c(model_script, paste0("\t for(i in 1:n) { for(k in 1:num.lv) { lvs[i,k] ~ dnorm(0,1) } } \n\n\t ## Process level and priors ##"))
     if(lv.control$type == "exponential")
          model_script <- c(model_script, paste0("\t for(k in 1:num.lv) { lvs[1:n,k] ~ dmnorm(zero.lvs,invSigma.lvs) } \n\t for(k1 in 1:n) { for(k2 in 1:n) { Sigma.lvs[k1,k2] <- exp(-distmat[k1,k2]/lv.covparams[1]) } } \n\t invSigma.lvs <- inverse(Sigma.lvs) \n\n\t ## Process level and priors ##"))
     if(lv.control$type == "squared.exponential")
          model_script <- c(model_script, paste0("\t for(k in 1:num.lv) { lvs[1:n,k] ~ dmnorm(zero.lvs,invSigma.lvs) } \n\t for(k1 in 1:n) { for(k2 in 1:n) { Sigma.lvs[k1,k2] <- exp(-pow(distmat[k1,k2]/lv.covparams[1],2)) } } \n\t invSigma.lvs <- inverse(Sigma.lvs) \n\n\t ## Process level and priors ##"))
     if(lv.control$type == "powered.exponential")
          model_script <- c(model_script, paste0("\t for(k in 1:num.lv) { lvs[1:n,k] ~ dmnorm(zero.lvs,invSigma.lvs) } \n\t for(k1 in 1:n) { for(k2 in 1:n) { Sigma.lvs[k1,k2] <- exp(-pow(distmat[k1,k2]/lv.covparams[1],lv.covparams[2])) } } \n\t invSigma.lvs <- inverse(Sigma.lvs) \n\n\t ## Process level and priors ##"))
     if(lv.control$type == "spherical")
          model_script <- c(model_script, paste0("\t for(k in 1:num.lv) { lvs[1:n,k] ~ dmnorm(zero.lvs,invSigma.lvs) } \n\t for(k1 in 1:n) { for(k2 in 1:n) { Sigma.lvs[k1,k2] <- step(lv.covparams[1] - distmat[k1,k2])*(1 - 1.5*distmat[k1,k2]/lv.covparams[1] + 0.5*pow(distmat[k1,k2]/lv.covparams[1],3)) } } \n\t invSigma.lvs <- inverse(Sigma.lvs) \n\n\t ## Process level and priors ##"))          
     #if(lv.control$type == "cauchy") ## DOES NOT WORK VERY WELL!
          #model_script <- c(model_script, paste0("\t for(k in 1:num.lv) { lvs[1:n,k] ~ dmnorm(zero.lvs,invSigma.lvs) } \n\t for(k1 in 1:n) { for(k2 in 1:n) { Sigma.lvs[k1,k2] <- pow(1+ pow(distmat[k1,k2]/lv.covparams[1],2), -lv.covparams[2]) } } \n\t invSigma.lvs <- inverse(Sigma.lvs) \n\n\t ## Process level and priors ##"))
     ## Matern not implemented to due complications/lack of direct availability of a BesselK function
    
     ## Build prior strings for all priors distributions
     prior.strings <- construct_prior_strings(x = prior.control)
    

     ## Code for response-specific intercept. Note this is set up different to how X variables are set up to save some coding space!
     ## No traits or traits included but not regressed against intercept
     if(num.traits == 0 || (num.traits > 0 & which.traits[[1]][1] == 0)) { 
          ## Not ordinal columns, then as per usual
          if(length(index.ord.cols) == 0) 
               model_script <- c(model_script, paste0("\t for(j in 1:p) { lv.coefs[j,1] ~ ", prior.strings$p1, " } ## Separate response intercepts")) 
          ## If 1 ordinal column, then intercept for this column equal 0
          if(length(index.ord.cols) == 1) {
               model_script <- c(model_script, paste0("\t lv.coefs[",index.ord.cols, ",1] <- 0 ## Single ordinal response intercept"))
               for(j in (1:p)[-index.ord.cols]) 
                    model_script <- c(model_script, paste0("\t lv.coefs[", j, ",1] ~ ", prior.strings$p1, "All other response intercepts"))
               }
          ## More than 1 ordinal column, then set up random intercept for this response
          if(length(index.ord.cols) > 1) {
               if(length(index.ord.cols) == p)
                    model_script <- c(model_script, paste0("\t for(j in 1:p) { lv.coefs[j,1] ~ dnorm(0,pow(ordinal.sigma,-2)) } ## Random intercept for all ordinal response"))
               else {
                    for(j in index.ord.cols)
                         model_script <- c(model_script, paste0("\t lv.coefs[",j, ",1] ~ dnorm(0,pow(ordinal.sigma,-2)) ## Random intercept for all ordinal response"))	
                    for(j in (1:p)[-index.ord.cols]) 
                         model_script <- c(model_script, paste0("\t lv.coefs[", j, ",1] ~ ", prior.strings$p1, "All other response intercepts"))
                    }
               model_script <- c(model_script, paste0("\t ordinal.sigma ~ ", prior.strings$p4))
               }
          if((num.traits > 0 & which.traits[[1]][1] == 0)) { 
               model_script <- c(model_script, paste0("\t traits.int[1] <- 0; for(l in 1:num.traits) { traits.coefs[1,l] <- 0 } \n\t trait.sigma[1] <- 0 ## Traits not used for intercept")) 
               }
          }
     ## Traits included in model and regressed against intercept
     if(num.traits > 0 & all(which.traits[[1]] > 0)) {
          ## If there are 0 or > 1 ordinal columns, then regress all intercepts against traits
          if(length(index.ord.cols) != 1) { 
               model_script <- c(model_script, paste0("\t for(j in 1:p) { lv.coefs[j,1] ~ dnorm(traits.int[1] + inprod(traits[j,],traits.coefs[1,1:num.traits]),pow(trait.sigma[1],-2)) } ## response intercepts regressed against traits"))
               }
          ## If there is 1 ordinal column, do not regress this intercept against trait	
          if(length(index.ord.cols) == 1) { 
               model_script <- c(model_script, paste0("\t lv.coefs[",index.ord.cols, ",1] <- 0 ## Ordinal response intercept"))
               for(j in (1:p)[-index.ord.cols]) 
                    model_script <- c(model_script, paste0("\t lv.coefs[", j, ",1] ~ dnorm(traits.int[1] + inprod(traits[",j,",],traits.coefs[1,1:num.traits]),pow(trait.sigma[1],-2)) ## All other intercepts"))
                    }
               model_script <- c(model_script, paste0("\t traits.int[1] ~ ", prior.strings$p1)) 
               for(l in which.traits[[1]]) {
                    if(prior.control$ssvs.traitsindex[[1]][which(which.traits[[1]] == l)] == -1)
                         model_script <- c(model_script, paste0("\t traits.coefs[",1,",",l,"] ~ ", prior.strings$p1, " ## Traits used for intercept")) 
                    if(prior.control$ssvs.traitsindex[[1]][which(which.traits[[1]] == l)] == 0) {
                         ssvs.prior.string <- paste0("dnorm(0,pow(", prior.control$hypparams[1], "*((1-ssvs.traitscoefs1",l,")*", prior.control$ssvs.g, " + ssvs.traitscoefs1", l, "),-1)); ssvs.traitscoefs1", l, " ~ dbern(0.5)")
                         model_script <- c(model_script, paste0("\t traits.coefs[",1,",",l,"] ~ ", ssvs.prior.string, " ## Traits used for intercept")) 
                         }
                    }
               if(length((1:num.traits)[-which.traits[[1]]]) > 0) {
                    for(l in (1:num.traits)[-which.traits[[1]]]) 
                         {
                         model_script <- c(model_script, paste0("\t traits.coefs[",1,",",l,"] <- 0 ## Traits not used for intercept")) 
                         }
                    }		
          model_script <- c(model_script, paste0("\t trait.sigma[1] ~ ", prior.strings$p4)) 
          }		 
            
            
     if(any(complete_family == "tweedie")) 
          model_script <- c(model_script, paste0("\t powerparam ~ dunif(1,2) ## Tweedie power parameter"))
     if(any(complete_family == "ordinal")) { 
          model_script <- c(model_script, paste0("\t for(k in 1:(num.ord.levels-1)) { cutoffs0[k] ~ ", prior.strings$p1, " }"))
          model_script <- c(model_script, paste0("\t cutoffs[1:(num.ord.levels-1)] <- sort(cutoffs0) ## Ordinal cutoffs"))
          }
            
            
     ## Priors on row effects
     if(row.eff == "fixed") {
          for(k in 1:ncol(row.ids)) 
               {
               model_script <- c(model_script, paste0("\n\t row.coefs.ID", k, "[1] <- 0"))
               model_script <- c(model_script, paste0("\n\t for(i in 2:n.ID[", k, "]) { row.coefs.ID", k, "[i] ~ ", prior.strings$p1, " } "))
               }
          }
     if(row.eff == "random") {
          for(k in 1:ncol(row.ids)) 
               {
               model_script <- c(model_script, paste0("\n\t for(i in 1:n.ID[", k, "]) { row.coefs.ID", k, "[i] ~ dnorm(0, pow(row.sigma.ID", k, ",-2)) } "))
               model_script <- c(model_script, paste0("\t row.sigma.ID", k, " ~ ", prior.strings$p4))
               }
          }

          
     ## Priors on response-specific random intercepts
     if(!is.null(ranef.ids)) {
          for(k0 in 1:ncol(ranef.ids)) {
               model_script <- c(model_script, paste0("\n\t for(j in 1:p) { for(i in 1:n.ranefID[", k0, "]) { ranef.coefs.ID", k0, "[j,i] ~ dnorm(0, pow(ranef.sigma.ID", k0, "[j],-2)) } }"))
               model_script <- c(model_script, paste0("\t for(j in 1:p) { ranef.sigma.ID", k0, "[j] ~ ",prior.strings$p4, " }"))
               }
          }     

          
     ## Priors on latent variables if required, controlled by prior.control$hypparams[2]
     if(lv.control$type %in% c("exponential","squared.exponential","spherical")) {
          model_script <- c(model_script, paste0("\t lv.covparams[1] ~ ", prior.strings$p22))
          }
     if(lv.control$type %in% c("powered.exponential")) {
          model_script <- c(model_script, paste0("\t lv.covparams[1] ~ ", prior.strings$p22))
          model_script <- c(model_script, paste0("\t lv.covparams[2] ~ dunif(0,2)"))
          }
    
     ## Priors on Latent variable coefficients, controlled by prior.control$hypparams[2]
     model_script <- c(model_script, paste0("\n\t for(i in 1:(num.lv-1)) { for(j in (i+2):(num.lv+1)) { lv.coefs[i,j] <- 0 } } ## Constraints to 0 on upper diagonal"))
     model_script <- c(model_script, paste0("\t for(i in 1:num.lv) { lv.coefs[i,i+1] ~ ", prior.strings$p22, " } ## Sign constraints on diagonal elements"))
     model_script <- c(model_script, paste0("\t for(i in 2:num.lv) { for(j in 2:i) { lv.coefs[i,j] ~ ", prior.strings$p2, " } } ## Free lower diagonals"))
     model_script <- c(model_script, paste0("\t for(i in (num.lv+1):p) { for(j in 2:(num.lv+1)) { lv.coefs[i,j] ~ ", prior.strings$p2, " } } ## All other elements"))

            
     ## Prior for X-coefficients, controlled by prior.control$hypparams[3]	
     if(num.X > 0) {
          model_script <- c(model_script, paste0("\n"))
          ## Traits not included in model
          if(num.traits == 0) { 
               for(i in 1:length(prior.control$ssvs.index)) {
                    if(prior.control$ssvs.index[i] == -1) {
                         if(!is.null(X.ind))
                         model_script <- c(model_script, paste0("\t for(j in 1:p) { X.coefs[j,", i, "] ~ ", prior.strings$p3, "I(-X.ind[j,",i,"],X.ind[j,",i,"]) } ")) 
                         if(is.null(X.ind))
                         model_script <- c(model_script, paste0("\t for(j in 1:p) { X.coefs[j,", i, "] ~ ", prior.strings$p3, " } ")) 
                         }
                    if(prior.control$ssvs.index[i] == 0) {
                         ssvs.prior.string <- paste0("dnorm(0,pow(", prior.control$hypparams[3], "*((1-ssvs.indX", i, "[j])*", prior.control$ssvs.g, " + ssvs.indX", i, "[j]),-1)); ssvs.indX", i, "[j] ~ dbern(0.5)")
                         model_script <- c(model_script, paste0("\t for(j in 1:p) { X.coefs[j,", i, "] ~ ", ssvs.prior.string, " }"))
                         }
                    if(prior.control$ssvs.index[i] > 0) {
                         ssvs.prior.string <- paste0("dnorm(0,pow(", prior.control$hypparams[3], "*((1-ssvs.gp", prior.control$ssvs.index[i], ")*", prior.control$ssvs.g, " + ssvs.gp", prior.control$ssvs.index[i], "),-1))")
                         model_script <- c(model_script, paste0("\t for(j in 1:p) { X.coefs[j,", i, "] ~ ", ssvs.prior.string, " } "))
                         }
                    } 			
               }
                    
          if(num.traits > 0) { for(i in 1:num.X) {
               ## Traits included but X coefs not regressed against them
               if(which.traits[[i+1]][1] == 0) { 
                    model_script <- c(model_script, paste0("\t for(j in 1:p) { X.coefs[j,", i, "] ~ ", prior.strings$p3, " } ## Coefficient not regressed against any traits")) 
                    model_script <- c(model_script, paste0("\t traits.int[",i+1,"] <- 0; trait.sigma[",i+1,"] <- 0; for(l in 1:num.traits) { traits.coefs[",i+1,",l] <- 0 } \n")) 
                    }
               ## Traits included and X coefs regressed against some of them
               if(all(which.traits[[i+1]] > 0)) { 
                    model_script <- c(model_script, paste0("\t for(j in 1:p) { X.coefs[j,", i, "] ~ dnorm(traits.int[",i+1,"] + inprod(traits[j,],traits.coefs[",i+1,",1:num.traits]),pow(trait.sigma[",i+1,"],-2)) } "))
                    model_script <- c(model_script, paste0("\t traits.int[",i+1,"] ~ ", prior.strings$p3)) 
                    for(l in which.traits[[i+1]]) {
                         if(prior.control$ssvs.traitsindex[[i+1]][which(which.traits[[i+1]] == l)] == -1)
                              model_script <- c(model_script, paste0("\t traits.coefs[",i+1,",",l,"] ~ ", prior.strings$p3, " ## Traits used for this X.coefs")) 
                         if(prior.control$ssvs.traitsindex[[i+1]][which(which.traits[[i+1]] == l)] == 0) {
                              ssvs.prior.string <- paste0("dnorm(0,pow(", prior.control$hypparams[3], "*((1-ssvs.traitscoefs",i+1,l,")*", prior.control$ssvs.g, " + ssvs.traitscoefs", i+1, l, "),-1)); ssvs.traitscoefs", i+1, l, " ~ dbern(0.5)")
                              model_script <- c(model_script, paste0("\t traits.coefs[",i+1,",",l,"] ~ ", ssvs.prior.string, " ## Traits used for this X.coefs")) 
                              }
                         }
                    if(length((1:num.traits)[-which.traits[[i+1]]]) > 0) {
                         for(l in (1:num.traits)[-which.traits[[i+1]]]) {
                              model_script <- c(model_script, paste0("\t traits.coefs[",i+1,",",l,"] <- 0 ## traits not used for this X.coefs")) 
                              }
                         }
                    model_script <- c(model_script, paste0("\t trait.sigma[",i+1,"] ~ ", prior.strings$p4, "\n"))
                    }
               } }
                            
          model_script <- c(model_script, paste0(""))
          if(any(prior.control$ssvs.index > 0)) {
               for(i in unique(prior.control$ssvs.index[prior.control$ssvs.index > 0])) 
                    model_script <- c(model_script, paste0("\t ssvs.gp", i, " ~ dbern(0.5)")) 
               }
          }	
# 	if(num.X > 0 & any(family == "multinom")) {
# 		model_script <- c(model_script, paste0("\t for(j in 1:",length(index_multinom_cols),") { for(i in 1:num.X) { X.multinom.params[j,i,1] <- 0 } } "))
# 		model_script <- c(model_script, paste0("\t for(j in 1:",length(index_multinom_cols),") { for(i in 1:num.X) { for(k in 2:num.multinom.levels) { X.multinom.params[j,i,k] ~ dnorm(0,",1/prior.control$hypparams[1],") } } } "))
# 		}
            
            
     ## Prior on dispersion parameters, controlled by prior.control$hypparams[4]
     if(!all(complete_family %in% c("poisson", "ztpoisson", "binomial", "ordinal", "multinom", "exponential"))) {
          model_script <- c(model_script, paste0("\t for(j in 1:p) { lv.coefs[j,num.lv+2] ~ ", prior.strings$p4, " } ## Dispersion parameters")) 
          }
                    
     model_script <- c(model_script, "\n\t }")

     
     if(!is.null(model.name)) { write(model_script, file = model.name) }
     if(is.null(model.name)) { write(model_script, file = "jagsboralmodel.txt") }
     }



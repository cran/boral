make.jagsboralmodel <- function(family, num.X = 0, num.traits = 0, which.traits = NULL, num.lv = 2, row.eff = "none", row.ids = NULL, offset = NULL, trial.size = 1, n, p, model.name = NULL, prior.control = list(type = c("normal","normal","normal","uniform"), hypparams = c(10, 10, 10, 30), ssvs.index = -1, ssvs.g = 1e-6)) {
	if(num.X == 0 & num.traits > 0) 
		stop("num.traits > 0 suggests traits are to be regressed against covariates X, so please set num.X > 0.") 
 	if(num.traits > 0 & is.null(which.traits)) 
		stop("If num.traits > 0, then please supply which.traits to inform what traits are regressed against which covariates.") 
 	if(!is.null(which.traits) & ((num.X+1) != length(which.traits))) 
		stop("which.traits should be a list with length 1+num.X.") 
 	if(!is.null(which.traits) & any(sapply(which.traits,length) > num.traits)) 
		stop("Each element in the list which.traits should have at most num.traits elements.") 
 	if(!is.null(which.traits) & any(prior.control$ssvs.index > -1)) 
		stop("Current version of boral only supports ssvs.index = -1 when traits are supplied...sorry!")
 	if(is.null(which.traits)) { 
		which.traits <- vector("list",num.X+1)
		for(k in 1:length(num.X+1)) 
			which.traits[[k]] <- 0 
		} 

	
	if(any(family == "binomial") & !(length(trial.size) %in% c(1, length(family)))) 
        stop("trial.size needs to be specified if any columns are binomially distributed; can either be a single element or a vector equal to p.")
	if(any(family == "binomial") & length(trial.size) == 1) {
		complete.trial.size <- rep(0, p)
		complete.trial.size[which(family == "binomial")] <- trial.size 
		}
	if(any(family == "binomial") & length(trial.size) == p) 
		complete.trial.size <- trial.size
	if(all(family != "binomial")) 
		complete.trial.size <- rep(0, p)

	if(length(family) == 1) 
		complete.family <- rep(family, p)
	if(length(family) == p) 
		complete.family <- family
	if(length(family) != p & length(family) != 1) 
		stop("Number of elements in family must either one or p")
	if(all(complete.family == "binomial") & all(complete.trial.size == 1)) 
		family <- rep("bernoulli",p)

	
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

		
	if(!is.null(offset)) { 
		if(!is.matrix(offset)) 
			stop("offset could be a matrix with the same dimensions as y.")
		}
		
     prior.control <- fillin.prior.control(x = prior.control)
		
	if(length(prior.control$hypparams) != 4 || length(prior.control$type) != 4) 
		stop("prior.control$type and prior.control$hypparams must be a vector of four elements. Please see boral help file as to what the elements correspond to.\n")
	if(!all(prior.control$type[-4] %in% c("normal","uniform","cauchy")))	
		stop("At least one of the first three elements of prior.control$type is not supported in current version of boral...sorry!")
	if(!(prior.control$type[4] %in% c("uniform","halfcauchy")))	
		stop("The fourth element of prior.control$type is not supported in current version of boral...sorry!")
	if(length(prior.control$ssvs.index) == 1 & num.X > 0) prior.control$ssvs.index <- rep(prior.control$ssvs.index, num.X)
	
	index.ord.cols <- which(complete.family == "ordinal")
	index.tweed.cols <- which(complete.family == "tweedie")
	
	
	## Checks done; starting writing JAGS script!

	
	mod.general.lv <- paste0("model {")
	mod.general.lv <- c(mod.general.lv, "\t ## Data Level ## \n\t for(i in 1:n) {")
	
	write.resp.script <- setup.resp.families.lv(p = p, complete.family = complete.family, num.lv = num.lv, row.eff = row.eff, row.ids = row.ids, offset = offset, num.X = num.X, complete.trial.size = complete.trial.size, index.tweed.cols = index.tweed.cols, index.ord.cols = index.ord.cols)
	mod.general.lv <- c(mod.general.lv, write.resp.script)	
	mod.general.lv <- c(mod.general.lv, paste0("\t\t } \n\n\t ## Latent variables ## \n\t for(i in 1:n) { for(k in 1:num.lv) { lvs[i,k] ~ dnorm(0,1) } } \n\n\t ## Process level and priors ##"))
	rm(write.resp.script)
    
	## Build prior strings for all priors distributions
     prior.strings <- construct.prior.strings(x = prior.control)

     
	## No traits or traits included but enviro coefs not regressed against them 
	if(num.traits == 0 || (num.traits > 0 & all(which.traits[[1]] == 0))) { 
		if(length(index.ord.cols) == 0) 
			mod.general.lv <- c(mod.general.lv, paste0("\t for(j in 1:p) { lv.coefs[j,1] ~ ", prior.strings$p1, " } ## Separate species intercepts")) 
		if(length(index.ord.cols) == 1) {
			mod.general.lv <- c(mod.general.lv, paste0("\t lv.coefs[",index.ord.cols, ",1] <- 0 ## Ordinal species intercept"))
			for(j in (1:p)[-index.ord.cols]) 
				mod.general.lv <- c(mod.general.lv, paste0("\t lv.coefs[", j, ",1] ~ ", prior.strings$p1))
			}
		if(length(index.ord.cols) > 1) {
			if(length(index.ord.cols) == p)
				mod.general.lv <- c(mod.general.lv, paste0("\t for(j in 1:p) { lv.coefs[j,1] ~ dnorm(0,pow(ordinal.sigma,-2)) } ## Ordinal species random intercept"))
			else {
				for(j in index.ord.cols)
					mod.general.lv <- c(mod.general.lv, paste0("\t lv.coefs[",j, ",1] ~ dnorm(0,pow(ordinal.sigma,-2)) ## Ordinal species random intercept"))	
				for(j in (1:p)[-index.ord.cols]) 
					mod.general.lv <- c(mod.general.lv, paste0("\t lv.coefs[", j, ",1] ~ ", prior.strings$p1))
				}				
			mod.general.lv <- c(mod.general.lv, paste0("\t ordinal.sigma ~ ", prior.strings$p4))
			}
		if((num.traits > 0 & all(which.traits[[1]] == 0))) { 
			mod.general.lv <- c(mod.general.lv, paste0("\t traits.int[1] <- 0; for(l in 1:num.traits) { traits.params[1,l] <- 0 } \n\t trait.sigma[1] <- 0")) 
			}
		}
	## Traits included and enviro coefs regressed against (some of) them
	if(num.traits > 0 & all(which.traits[[1]] > 0)) { 
		if(length(index.ord.cols) != 1) { ## If there are 0 or > 1 ordinal columns, then regress all intercepts against traits
			mod.general.lv <- c(mod.general.lv, paste0("\t for(j in 1:p) { lv.coefs[j,1] ~ dnorm(traits.int[1] + inprod(traits[j,],traits.params[1,1:num.traits]),pow(trait.sigma[1],-2)) } ## Species intercepts regressed against traits"))
			}
		if(length(index.ord.cols) == 1) { ## If there is 1 ordinal column, do not regress this intercept against trait
			mod.general.lv <- c(mod.general.lv, paste0("\t lv.coefs[",index.ord.cols, ",1] <- 0 ## Ordinal species intercept"))
			for(j in (1:p)[-index.ord.cols]) 
				mod.general.lv <- c(mod.general.lv, paste0("\t lv.coefs[", j, ",1] ~ dnorm(traits.int[1] + inprod(traits[",j,",],traits.params[1,1:num.traits]),pow(trait.sigma[1],-2))"))
			}
		for(l in which.traits[[1]]) {
			mod.general.lv <- c(mod.general.lv, paste0("\t traits.params[",1,",",l,"] ~ ", prior.strings$p1, " ## traits used for this spp intercept")) 
			}
		if(length((1:num.traits)[-which.traits[[1]]]) > 0) {
			for(l in (1:num.traits)[-which.traits[[1]]]) {
				mod.general.lv <- c(mod.general.lv, paste0("\t traits.params[",1,",",l,"] <- 0 ## traits not used for this X.coefs")) 
				}
			}
		
		mod.general.lv <- c(mod.general.lv, paste0("\t traits.int[1] ~ ", prior.strings$p1)) 
		mod.general.lv <- c(mod.general.lv, paste0("\t trait.sigma[1] ~ ", prior.strings$p4)) 
		}			 
		
		
	if(any(complete.family == "tweedie")) mod.general.lv <- c(mod.general.lv, paste0("\t powerparam ~ dunif(1,2) ## Tweedie power parameter"))
	if(any(complete.family == "ordinal")) { 
		mod.general.lv <- c(mod.general.lv, paste0("\t for(k in 1:(num.ord.levels-1)) { cutoffs0[k] ~ ", prior.strings$p1, " }"))
		mod.general.lv <- c(mod.general.lv, paste0("\t cutoffs[1:(num.ord.levels-1)] <- sort(cutoffs0) ## Ordinal cutoffs"))
		}
		
		
	## Priors on row effects
	if(row.eff == "fixed") {
		for(k in 1:ncol(row.ids)) 
			mod.general.lv <- c(mod.general.lv, paste0("\n\t for(i in 1:n.ID[", k, "]) { row.coefs.", colnames(row.ids)[k], "[i] ~ ", prior.strings$p1, " } "))
		}
	if(row.eff == "random") {
		for(k in 1:ncol(row.ids)) {
			mod.general.lv <- c(mod.general.lv, paste0("\n\t for(i in 1:n.ID[", k, "]) { row.coefs.", colnames(row.ids)[k], "[i] ~ dnorm(0, pow(row.sigma.", colnames(row.ids)[k], ",-2)) } "))
			mod.general.lv <- c(mod.general.lv, paste0("\t row.sigma.", colnames(row.ids)[k], " ~ ", prior.strings$p4))
			}
		#mod.general.lv <- c(mod.general.lv, paste0("\t row.ranef.mean ~ ", prior.string))
		}

		
	## Priors on Latent variable coefficients, controlled by prior.control$hypparams[2]
	mod.general.lv <- c(mod.general.lv, paste0("\n\t for(i in 1:(num.lv-1)) { for(j in (i+2):(num.lv+1)) { lv.coefs[i,j] <- 0 } } ## Constraints to 0 on upper diagonal"))
     mod.general.lv <- c(mod.general.lv, paste0("\t for(i in 1:num.lv) { lv.coefs[i,i+1] ~ ", prior.strings$p22, " } ## Sign constraints on diagonal elements"))
     mod.general.lv <- c(mod.general.lv, paste0("\t for(i in 2:num.lv) { for(j in 2:i) { lv.coefs[i,j] ~ ", prior.strings$p2, " } } ## Free lower diagonals"))
     mod.general.lv <- c(mod.general.lv, paste0("\t for(i in (num.lv+1):p) { for(j in 2:(num.lv+1)) { lv.coefs[i,j] ~ ", prior.strings$p2, " } } ## All other elements"))

		
	## Prior for X-coefficients, controlled by prior.control$hypparams[3]	
	if(num.X > 0) {
		mod.general.lv <- c(mod.general.lv, paste0("\n"))
	
		if(num.traits == 0) { ## Traits not included in model, so SSVS is permitted
			for(i in 1:length(prior.control$ssvs.index)) {
				if(prior.control$ssvs.index[i] == -1) {
					mod.general.lv <- c(mod.general.lv, paste0("\t for(j in 1:p) { X.coefs[j,", i, "] ~ ", prior.strings$p3, " } ")) 
					}
				if(prior.control$ssvs.index[i] == 0) {
					ssvs.prior.string <- paste0("dnorm(0,pow(", prior.control$hypparams[3], "*((1-ssvs.indcoefs", i, "[j])*", prior.control$ssvs.g, " + ssvs.indcoefs", i, "[j]),-1))")
					mod.general.lv <- c(mod.general.lv, paste0("\t for(j in 1:p) { X.coefs[j,", i, "] ~ ", ssvs.prior.string, "; ssvs.indcoefs", i, "[j] ~ dbern(0.5) }"))
					}
				if(prior.control$ssvs.index[i] > 0) {
					ssvs.prior.string <- paste0("dnorm(0,pow(", prior.control$hypparams[3], "*((1-ssvs.gpcoefs", prior.control$ssvs.index[i], ")*", prior.control$ssvs.g, " + ssvs.gpcoefs", prior.control$ssvs.index[i], "),-1))")
					mod.general.lv <- c(mod.general.lv, paste0("\t for(j in 1:p) { X.coefs[j,", i, "] ~ ", ssvs.prior.string, " } "))
					}
				} 				
			}
			
			
		if(num.traits > 0) { for(i in 1:num.X) {
			if(all(which.traits[[i+1]] == 0)) { ## Traits included but enviro coefs not regressed against them
				mod.general.lv <- c(mod.general.lv, paste0("\t for(j in 1:p) { X.coefs[j,", i, "] ~ ", prior.strings$p3, " } ## Coefficient not regressed against any traits")) 
				mod.general.lv <- c(mod.general.lv, paste0("\t traits.int[",i+1,"] <- 0; trait.sigma[",i+1,"] <- 0; for(l in 1:num.traits) { traits.params[",i+1,",l] <- 0 } \n")) 
				}
				
			if(all(which.traits[[i+1]] > 0)) { ## Traits included and environ coefs regressed against some of them
				mod.general.lv <- c(mod.general.lv, paste0("\t for(j in 1:p) { X.coefs[j,", i, "] ~ dnorm(traits.int[",i+1,"] + inprod(traits[j,],traits.params[",i+1,",1:num.traits]),pow(trait.sigma[",i+1,"],-2)) } "))
				for(l in which.traits[[i+1]]) {
					mod.general.lv <- c(mod.general.lv, paste0("\t traits.params[",i+1,",",l,"] ~ ", prior.strings$p3, " ## traits used for this X.coefs")) 
					}
				if(length((1:num.traits)[-which.traits[[i+1]]]) > 0) {
					for(l in (1:num.traits)[-which.traits[[i+1]]]) {
						mod.general.lv <- c(mod.general.lv, paste0("\t traits.params[",i+1,",",l,"] <- 0 ## traits not used for this X.coefs")) 
						}
					}
				mod.general.lv <- c(mod.general.lv, paste0("\t traits.int[",i+1,"] ~ ", prior.strings$p3)) 
				mod.general.lv <- c(mod.general.lv, paste0("\t trait.sigma[",i+1,"] ~ ", prior.strings$p4, "\n"))
				}
			} }
				
		mod.general.lv <- c(mod.general.lv, paste0(""))
		if(any(prior.control$ssvs.index > 0)) {
			for(i in unique(prior.control$ssvs.index[prior.control$ssvs.index > 0])) 
				mod.general.lv <- c(mod.general.lv, paste0("\t ssvs.gpcoefs", i, " ~ dbern(0.5)")) 
			}
		}	
# 	if(num.X > 0 & any(family == "multinom")) {
# 		mod.general.lv <- c(mod.general.lv, paste0("\t for(j in 1:",length(index.multinom.cols),") { for(i in 1:num.X) { X.multinom.params[j,i,1] <- 0 } } "))
# 		mod.general.lv <- c(mod.general.lv, paste0("\t for(j in 1:",length(index.multinom.cols),") { for(i in 1:num.X) { for(k in 2:num.multinom.levels) { X.multinom.params[j,i,k] ~ dnorm(0,",1/prior.control$hypparams[1],") } } } "))
# 		}
		
		
	## Prior on dispersion parameters, controlled by prior.control$hypparams[4]
	if(!all(complete.family %in% c("poisson", "binomial", "ordinal", "multinom", "exponential","bernoulli"))) {
	 	mod.general.lv <- c(mod.general.lv, paste0("\t for(j in 1:p) { lv.coefs[j,num.lv+2] ~ ", prior.strings$p4, " } ## Dispersion parameters")) 
	 	}
			
	mod.general.lv <- c(mod.general.lv, "\n\t }")

	
	if(!is.null(model.name)) { write(mod.general.lv, file = model.name) }
	if(is.null(model.name)) { write(mod.general.lv, file = "jagsboralmodel.txt") }
	}



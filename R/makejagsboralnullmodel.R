make.jagsboralnullmodel <- function (family, num.X = 0, X.ind = NULL, num.traits = 0, which.traits = NULL, 
     row.eff = "none", row.ids = NULL, offset = NULL, trial.size = 1, n, p, model.name = NULL, 
     prior.control = list(type = c("normal","normal","normal","uniform"), hypparams = c(10, 10, 10, 30), ssvs.index = -1, ssvs.g = 1e-6, ssvs.traitsindex = -1)) {

     if(num.X == 0 & num.traits > 0) 
        stop("num.traits > 0 suggests traits are to be regressed against covariates X, so please set num.X > 0") 
    if(num.traits > 0 & is.null(which.traits)) 
        stop("If num.traits > 0, then please supply which.traits to inform what traits are regressed against each covariate") 
    if(!is.null(which.traits) & ((num.X+1) != length(which.traits))) 
        stop("which.traits should be a list with length 1+num.X") 
    if(!is.null(which.traits) & any(sapply(which.traits,length) > num.traits)) 
        stop("Each element in the list which.traits should have at most num.traits elements") 
    if(is.null(which.traits)) { 
        which.traits <- vector("list",num.X+1)
        for(k in 1:(num.X+1)) 
            which.traits[[k]] <- 0 
        } 

    
    if(any(family == "binomial") & !(length(trial.size) %in% c(1, length(family)))) 
        stop("trial.size needs to be specified if any columns are binomially distributed; can either be a single element or a vector equal to p")
    if(any(family == "binomial") & length(trial.size) == 1) {
        complete_trial_size <- rep(0, p)
        complete_trial_size[which(family == "binomial")] <- trial.size 
        }
    if(any(family == "binomial") & length(trial.size) == p) 
        complete_trial_size <- trial.size
    if(all(family != "binomial")) 
        complete_trial_size <- rep(0, p)

    if(length(family) == 1) 
        complete_family <- rep(family, p)
    if(length(family) == p) 
        complete_family <- family
    if(length(family) != p & length(family) != 1) 
        stop("Number of elements in family must either one or p")
    if(all(complete_family == "binomial") & all(complete_trial_size == 1)) 
        family <- rep("bernoulli",p)

    
    if(row.eff != "none" && is.null(row.ids)) {
        row.ids <- matrix(1:n, ncol = 1)
        message("row.ids assumed to be matrix with one column and elements 1,2,...n i.e., a row-specific intercept")
        }
    if(!is.null(row.ids)) {
        row.ids <- as.matrix(row.ids)
        if(nrow(row.ids) != n) 
            stop("Number of rows in the matrix row.ids should be equal to n")
        if(is.null(colnames(row.ids))) 
            colnames(row.ids) <- paste0("ID", 1:ncol(row.ids))
        }
    
    
    if(!is.null(offset)) { 
        if(!is.matrix(offset)) 
            stop("offset could be a matrix with the same dimensions as y")
        }

            
    prior.control <- fillin.prior.control(x = prior.control)
    if(length(prior.control$hypparams) != 4 || length(prior.control$type) != 4) 
        stop("prior.control$type and prior.control$hypparams must be a vector of four elements. Please see boral help file as to what the elements correspond to.\n")
    if(!all(prior.control$type[-4] %in% c("normal","uniform","cauchy")))	
        stop("At least one of the first three elements of prior.control$type is not supported in current version of boral...sorry!")
    if(!(prior.control$type[4] %in% c("uniform","halfcauchy")))	
        stop("The fourth element of prior.control$type is not supported in current version of boral...sorry!")
    if(length(prior.control$ssvs.index) == 1 & num.X > 0) prior.control$ssvs.index <- rep(prior.control$ssvs.index, num.X)
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

    X.ind <- check.X.ind(X.ind = X.ind, p = p, num.X = num.X, prior.control = prior.control)

    index.ord.cols <- which(complete_family == "ordinal")
    index.tweed.cols <- which(complete_family == "tweedie")

    
    ## Checks done; starting writing JAGS script!
    
    
    model_script <- paste0("## JAGS model written for boral version ", packageDescription("boral")$Version, " on ", as.character(Sys.time()), " ##\n\n model {")
    model_script <- c(model_script, "\t ## Data Level ## \n\t for(i in 1:n) {")

    write.resp.script <- setup.resp.families.lv(p = p, complete.family = complete_family, num.lv = 0, row.eff = row.eff, row.ids = row.ids, offset = offset, num.X = num.X, complete.trial.size = complete_trial_size, index.tweed.cols = index.tweed.cols, index.ord.cols = index.ord.cols)
    model_script <- c(model_script, write.resp.script)
    model_script <- c(model_script, paste0("\t\t } \n\n\t ## Process level and priors ##"))
    
    
    ## Build prior strings for all priors distributions
    prior.strings <- construct.prior.strings(x = prior.control)

    
    ## Code for column-specific intercept. Note this is set up different to how X variables are set up to save some coding space!
    ## No traits or traits included but not regressed against intercept
    if(num.traits == 0 || (num.traits > 0 & which.traits[[1]][1] == 0)) { 
        ## Not ordinal columns, then as per usual
        if(length(index.ord.cols) == 0) 
            model_script <- c(model_script, paste0("\t for(j in 1:p) { lv.coefs[j,1] ~ ", prior.strings$p1, " } ## Separate species intercepts")) 
        ## If 1 ordinal column, then intercept for this column equal 0
        if(length(index.ord.cols) == 1) {
            model_script <- c(model_script, paste0("\t lv.coefs[",index.ord.cols, ",1] <- 0 ## Single ordinal species intercept"))
            for(j in (1:p)[-index.ord.cols]) 
                model_script <- c(model_script, paste0("\t lv.coefs[", j, ",1] ~ ", prior.strings$p1, "All other species intercepts"))
            }
        ## More than 1 ordinal column, then set up random intercept for this species
        if(length(index.ord.cols) > 1) {
            if(length(index.ord.cols) == p)
                model_script <- c(model_script, paste0("\t for(j in 1:p) { lv.coefs[j,1] ~ dnorm(0,pow(ordinal.sigma,-2)) } ## Random intercept for all ordinal species"))
                else {
                    for(j in index.ord.cols)
                        model_script <- c(model_script, paste0("\t lv.coefs[",j, ",1] ~ dnorm(0,pow(ordinal.sigma,-2)) ## Random intercept for all ordinal species"))	
                    for(j in (1:p)[-index.ord.cols]) 
                        model_script <- c(model_script, paste0("\t lv.coefs[", j, ",1] ~ ", prior.strings$p1, "All other species intercepts"))
                    }
                model_script <- c(model_script, paste0("\t ordinal.sigma ~ ", prior.strings$p4))
                }
            if((num.traits > 0 & which.traits[[1]][1] == 0)) { 
                model_script <- c(model_script, paste0("\t traits.int[1] <- 0; for(l in 1:num.traits) { traits.coefs[1,l] <- 0 } \n\t trait.sigma[1] <- 0 ## Traits not used for intercept")) 
                }
        }
    ## Traits included and intercept regressed against them
    if(num.traits > 0 & all(which.traits[[1]] > 0)) { 
        ## If there are 0 or > 1 ordinal columns, then regress all intercepts against traits
        if(length(index.ord.cols) != 1) { 
            model_script <- c(model_script, paste0("\t for(j in 1:p) { lv.coefs[j,1] ~ dnorm(traits.int[1] + inprod(traits[j,],traits.coefs[1,1:num.traits]),pow(trait.sigma[1],-2)) } ## Species intercepts regressed against traits"))
            }
        ## If there is 1 ordinal column, do not regress this intercept against trait	
        if(length(index.ord.cols) == 1) { 
            model_script <- c(model_script, paste0("\t lv.coefs[",index.ord.cols, ",1] <- 0 ## Ordinal species intercept"))
            for(j in (1:p)[-index.ord.cols]) 
                model_script <- c(model_script, paste0("\t lv.coefs[", j, ",1] ~ dnorm(traits.int[1] + inprod(traits[",j,",],traits.coefs[1,1:num.traits]),pow(trait.sigma[1],-2)) ## All other intercepts"))
            }
        for(l in which.traits[[1]]) {
            if(prior.control$ssvs.traitsindex[[1]][which.traits[[1]] == l] == -1)
                model_script <- c(model_script, paste0("\t traits.coefs[",1,",",l,"] ~ ", prior.strings$p1, " ## Traits used for intercept")) 
            if(prior.control$ssvs.traitsindex[[1]][which.traits[[1]] == l] == 0) {
                ssvs.prior.string <- paste0("dnorm(0,pow(", prior.control$hypparams[1], "*((1-ssvs.traitscoefs1",l,")*", prior.control$ssvs.g, " + ssvs.traitscoefs1", l, "),-1)); ssvs.traitscoefs1", l, " ~ dbern(0.5)")
                model_script <- c(model_script, paste0("\t traits.coefs[",1,",",l,"] ~ ", ssvs.prior.string, " ## Traits used for intercept")) 
                }
            }
        if(length((1:num.traits)[-which.traits[[1]]]) > 0) {
            for(l in (1:num.traits)[-which.traits[[1]]]) {
                model_script <- c(model_script, paste0("\t traits.coefs[",1,",",l,"] <- 0 ## Traits not used for intercept")) 
                }
            }		
        model_script <- c(model_script, paste0("\t traits.int[1] ~ ", prior.strings$p1)) 
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
            model_script <- c(model_script, paste0("\n\t for(i in 1:n.ID[", k, "]) { row.coefs.ID", k, "[i] ~ ", prior.strings$p1, " } "))
        }
    if(row.eff == "random") {
        for(k in 1:ncol(row.ids)) {
            model_script <- c(model_script, paste0("\n\t for(i in 1:n.ID[", k, "]) { row.coefs.ID", k, "[i] ~ dnorm(0, pow(row.sigma.ID", k, ",-2)) } "))
            model_script <- c(model_script, paste0("\t row.sigma.ID", k, " ~ ",prior.strings$p4))
            }
        #model_script <- c(model_script, paste0("\t row.ranef.mean ~ ", prior.strings$p1))
        }

            
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
            if(all(which.traits[[i+1]] == 0)) { 
                    model_script <- c(model_script, paste0("\t for(j in 1:p) { X.coefs[j,", i, "] ~ ", prior.strings$p3, " } ## Coefficient not regressed against any traits")) 
                    model_script <- c(model_script, paste0("\t traits.int[",i+1,"] <- 0; trait.sigma[",i+1,"] <- 0; for(l in 1:num.traits) { traits.coefs[",i+1,",l] <- 0 } \n")) 
                    }
            ## Traits included and X coefs regressed against some of them				
            if(all(which.traits[[i+1]] > 0)) { 
                model_script <- c(model_script, paste0("\t for(j in 1:p) { X.coefs[j,", i, "] ~ dnorm(traits.int[",i+1,"] + inprod(traits[j,],traits.coefs[",i+1,",1:num.traits]),pow(trait.sigma[",i+1,"],-2)) } "))
                for(l in which.traits[[i+1]]) {
                    if(prior.control$ssvs.traitsindex[[i+1]][which.traits[[i+1]] == l] == -1)
                        model_script <- c(model_script, paste0("\t traits.coefs[",i+1,",",l,"] ~ ", prior.strings$p3, " ## Traits used for this X.coefs")) 
                        if(prior.control$ssvs.traitsindex[[i+1]][which.traits[[i+1]] == l] == 0) {
                            ssvs.prior.string <- paste0("dnorm(0,pow(", prior.control$hypparams[3], "*((1-ssvs.traitscoefs",i+1,l,")*", prior.control$ssvs.g, " + ssvs.traitscoefs", i+1, l, "),-1)); ssvs.traitscoefs", i+1, l, " ~ dbern(0.5)")
                            model_script <- c(model_script, paste0("\t traits.coefs[",i+1,",",l,"] ~ ", ssvs.prior.string, " ## Traits used for this X.coefs")) 
                            }
                    }
                if(length((1:num.traits)[-which.traits[[i+1]]]) > 0) {
                    for(l in (1:num.traits)[-which.traits[[i+1]]]) {
                        model_script <- c(model_script, paste0("\t traits.coefs[",i+1,",",l,"] <- 0 ## traits not used for this X.coefs")) 
                        }
                    }
                model_script <- c(model_script, paste0("\t traits.int[",i+1,"] ~ ", prior.strings$p3)) 
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
    if(!all(complete_family %in% c("poisson", "binomial", "ordinal", "multinom", "exponential","bernoulli"))) {
        model_script <- c(model_script, paste0("\t for(j in 1:p) { lv.coefs[j,num.lv+2] ~ ", prior.strings$p4, " } ## Dispersion parameters")) 
        }
                        
    model_script <- c(model_script, "\n\t }")

    if(!is.null(model.name)) { write(model_script, file = model.name) }
    if(is.null(model.name)) { write(model_script, file = "jagsboralmodel.txt") }
    }

    
    

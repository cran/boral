##############
## Latent Variable Ordination and Regression using MCMC 
## Site effects fitted as fixed effects
## NB parameterized as V = mu + phi*mu^2
## Ordinal data handled as propotional odds regression, with same cutoff points for all spp. 
## Multinomial regression is available, but CURRENTLY NOT ALLOWED DUE TO COMPUTATION TIME
## This should make sense if the same ordinal scale is applied to all species
##############
# rm(list = ls())
# library(R2jags); 
# library(mvtnorm); 
# library(mvabund); 
# library(coda); 
# library(MASS)
# library(fishMod)
# source("auxilaryfunctions.R")
# #  
#data(spider)
#y <- spider$abun
#X <- spider$x; family = "negative.binomial"; num.lv = 2; site.eff = T; n.burnin = 4000; n.iteration = 24000; n.thin = 5; save.model = FALSE; calc.ics = TRUE; trial.size <- NULL; seed <- 123; hypparams = c(100,100); 

#library(ade4)
#data(dunedata)
#y = dunedata$veg+1 ## Shift levels up to start at 1
#X = model.matrix(~A1 + factor(use) - 1, data = dunedata$envir)
#family = rep("ordinal",30)
#num.lv = 2; site.eff = FALSE; n.burnin = 4000; n.iteration = 24000; n.thin = 4; save.model = TRUE; seed = 7; calc.ics = TRUE; trial.size = NULL; num.ord.levels <- 5; hypparams = c(100,100); 
#X <- matrix(rnorm(30*4),30,4)
#true.beta <- cbind(matrix(rnorm(length(family)*(ncol(X)+1)),length(family),ncol(X)+1),NA); 
#true.beta[nrow(true.beta),1] <- -sum(true.beta[-nrow(true.beta),1])
#true.ordinal <- seq(-0.5,0.5,length=num.ord.levels-1)
#y <- create.life(lv.coefs = true.beta[,c(1,ncol(true.beta))], X = X, X.coefs = true.beta[,-c(1,ncol(true.beta))], family = family, cutoffs = true.ordinal)

# n = 60; p <- 30
# X <- matrix(rnorm(n*2),n,2); beta <- cbind(matrix(rnorm(p*3),p,3),runif(p,0,5)); true.power <- 1.6
# mu <- exp(cbind(1,X)%*%t(beta[,1:3]))
# y <- matrix(NA,n,p)
# for(j in 1:ncol(y)) { y[,j] <- rTweedie(nrow(y), mu = mu[,j], phi = beta[j,4], p = true.power) }
# family = "tweedie"
# num.lv = 0; site.eff = FALSE; n.burnin = 4000; n.iteration = 24000; n.thin = 4; save.model = TRUE; seed = 1; calc.ics = TRUE; trial.size = NULL; hypparams = c(100,100); 
# 
# library(FD)
# data(tussock)
# y <- tussock$trait[,c("height","LDMC","leafN","leafsize","SLA","seedmass","clonality","resprouting","lifespan")]
# y$LDMC <- y$LDMC/1000 ## change to g/g
# y$leafsize <- y$leafsize/100 ## change to cm^2
# y[,"resprouting"] <- as.numeric(y[,"resprouting"])-1 ## 0 = no; 1 = yes
# levels(y[,"clonality"]) = c(3,2,1)
# y[,"clonality"] <- as.numeric(levels(y[,"clonality"]))[y[,"clonality"]]
# levels(y[,"lifespan"]) = c(0,0,1)
# y[,"lifespan"] <- as.numeric(levels(y[,"lifespan"]))[y[,"lifespan"]]
# y <- y[-which(is.na(rowSums(y))),]
# 
# family = c("lnormal","normal","normal","lnormal","lnormal","lnormal","ordinal","binomial","binomial")
# num.lv = 2; site.eff = FALSE; n.burnin = 4000; n.iteration = 24000; n.thin = 4; save.model = TRUE; seed = 123; calc.ics = TRUE; trial.size = 1; hypparams = c(100,100); X <- NULL

# n = 30; s <- 30; num.multinom.levels <- 4
# X <- matrix(rnorm(n*2),n,2)
# X.coefs <- rbind(matrix(0,s-1,2),c(1,2))
# X.multinom.coefs <- array(NA,dim=c(s-1,2,num.multinom.levels))
# for(k in 1:num.multinom.levels) { X.multinom.coefs[,,k] <- rnorm((s-1)*2) }
# site.coefs <- runif(n)
# lv.coefs <- cbind(matrix(runif(s,-3,-1),s,1),2)
# family = c(rep("multinom",s-1),"normal")
# num.lv = 2; site.eff = TRUE; n.burnin = 4000; n.iteration = 24000; n.thin = 4; save.model = TRUE; seed = 1; calc.ics = TRUE; trial.size = 1; hypparams = c(100,100); 
# y <- create.life(lv.coefs = lv.coefs, X = X, X.coefs = X.coefs, X.multinom.coefs = X.multinom.coefs, family = family, site.coefs = site.coefs)

boral <- function(y, ...) UseMethod("boral")

## Model is g(mu_{ij}) = site + theta0 + LV_i*theta_j + X_i*beta_j
boral.default <- function(y, X = NULL, family, trial.size = NULL, num.lv = 0, site.eff = FALSE, n.burnin = 4000, n.iteration = 24000, n.thin = 5, save.model = FALSE, seed = 123, calc.ics = TRUE, hypparams = c(100,100), ...) {
 	if(is.null(dim(y))) { cat("Converting y into a one column matrix.\n"); y <- matrix(y,ncol=1) }
 	if(!is.null(X) & is.null(dim(X))) { cat("Converting X into a one column matrix\n"); X <- matrix(X,ncol=1) }
 
 	if(num.lv > 3) warnings("We won't stop you, but please consider if you really want more than three latent variables in the model!")
 	if(length(family) != ncol(y) & length(family) != 1) { stop("Number of elements in family must either one or the # of columns in y") }
 	if(length(family) == 1) family <- rep(family, ncol(y))
 	if(!all(family %in% c("negative.binomial","poisson","binomial","normal","lnormal","tweedie","ordinal"))) 
 		stop("At least one of the elements in family is not compatible with current version of boral...sorry!") 
 	if(any(family == "binomial") & is.null(trial.size)) 
 		stop("trial.size needs to be specified if any columns are binomially distributed; can either be a single element or a vector equal to the # of columns in y. The latter will assume the specified trial size for all rows labelled binomial in the family argument.")
 	if(any(family == "ordinal")) {
 		if(sum(y[,family == "ordinal"] == 0) > 0) stop("For ordinal data, please shift minimum level to 1.") 
 		print("It is assumed all ordinal columns have the same number of levels -- please see help file as to the motivation behind this.")
 		print("boral may take a ridiculously long time to fit ordinal data models. Apolgies in advance!\n") }
 	if(any(family == "binomial")) {
 		if(!all(unlist(y[,family=="binomial"]) %in% c(0,1)))
 			stop("Any columns of y specified to be binomial must only contain 0/1 elements") }
 		
 		
 	if(!is.null(X)) { X.eff <- TRUE } else { X.eff <- FALSE } 
 	X.eff <- as.numeric(X.eff); 
 	if(!X.eff) { num.X <- 0 } else { num.X <- ncol(X) }
 
 	if(any(family == "binomial") & length(trial.size) == 1) { 
 		complete.trial.size <- rep(0,ncol(y)); complete.trial.size[which(family == "binomial")] <- trial.size }
 	if(any(family == "binomial") & length(trial.size) == ncol(y)) { complete.trial.size <- trial.size }
 	if(all(family != "binomial")) { complete.trial.size <- rep(0,ncol(y)) }
 
 	if(all(family != "ordinal")) { num.ord.levels <- 0; index.ord.cols <- NULL }
 	if(any(family == "ordinal")) { num.ord.levels <- max(y[,family == "ordinal"]); index.ord.cols <- which(family == "ordinal") }
 	
  	if(all(family != "multinom")) { num.multinom.levels <- 0; index.multinom.cols <- NULL }
  	if(any(family == "multinom")) { num.multinom.levels <- apply(y[,family == "multinom"],2,max); index.multinom.cols <- which(family == "multinom") }
 	
 	site.eff <- as.numeric(site.eff); 
 	n <- nrow(y); p <- ncol(y); 
 	n.chains = 1; ## Run one chain only to avoid arbitrary rotation problems
 	
 	if(num.lv > 0) make.jagsboralmodel(family, X.eff, site.eff, n, p, hypparams)
 	if(num.lv == 0) make.jagsboralnullmodel(family, X.eff, site.eff, n, p, hypparams)
 	
 	jags.data <- list("y","n","p","num.lv","num.X","complete.trial.size","num.ord.levels","num.multinom.levels")
 	if(X.eff) jags.data <-  c("X",jags.data)
 	if(any(family == "ordinal")) { ones <- matrix(1,n,p); jags.data <- c(jags.data,"ones") }
 
 	jags.params <- c("all.params")
 	if(num.lv > 0) jags.params <- c(jags.params, "lvs"); if(site.eff) jags.params <- c(jags.params,"site.params")
 	if(X.eff & any(family != "multinom")) jags.params <- c(jags.params,"X.params")
 	if(X.eff & any(family == "multinom")) jags.params <- c(jags.params,"X.multinom.params")
 	if(any(family == "tweedie")) jags.params <- c(jags.params,"powerparam")
 	if(any(family == "ordinal")) jags.params <- c(jags.params,"alpha")
 
 	jags.inits <- NULL
 	if(any(family %in% "tweedie")) { 
 		initial.list <- list("numfish" = matrix(1,n,p))
 		if(any(family %in% "ordinal")) initial.list$alpha0 <- seq(-1,1,length=num.ord.levels-1)
 		jags.inits <- function() { return(initial.list) } } 
 	set.seed(seed)
 
 	jagsfit <- suppressWarnings(jags(data=jags.data, inits=jags.inits, jags.params, model.file = "jagsboralmodel.txt", n.iter = n.iteration, n.burnin = n.burnin, n.thin = n.thin, n.chains = n.chains, DIC = T))
	## Format into big matrix; also deals with the conflict of as.mcmc converting to lists or the MCMC samples themselves
	fit.mcmcBase <- jagsfit$BUGSoutput
	fit.mcmc <- mcmc(fit.mcmcBase$sims.matrix, start = 1, thin = n.thin) ## Thanks to Guilliaume Blanchet for this!
# 	fit.mcmc <- as.mcmc(jagsfit)[[1]]
 	if(n.chains == 1) combined.fit.mcmc <- fit.mcmc
 	rm(fit.mcmc)
#	if(length(as.mcmc(jagsfit)) > 1) fit.mcmc <- as.mcmc(jagsfit) 
# 	if(n.chains > 1) {
# 		for(k in 1:n.chains) {
# 			if(k == 1) { combined.fit.mcmc <- fit.mcmc[[1]]; fit.mcmc[[1]] <- NA }
# 			if(k > 1) { combined.fit.mcmc <- rbind(combined.fit.mcmc, fit.mcmc[[k]]); fit.mcmc[[k]] <- NA }
# 			} }

 	## Flip dispersion parameters for NB returns phi_j, for normal and lognormal return sigma^2_j
 	sel.thetas <- grep("all.params",colnames(combined.fit.mcmc))
 	sel.thetas2 <- as.numeric(sel.thetas[(length(sel.thetas)-p+1):length(sel.thetas)])
 	combined.fit.mcmc[,sel.thetas2] <- 1/combined.fit.mcmc[,sel.thetas2]
 	if(any(family == "tweedie")) combined.fit.mcmc[,sel.thetas2[family == "tweedie"]] <- 1/combined.fit.mcmc[,sel.thetas2[family == "tweedie"]] ## Flip it back for tweedie
 	if(any(family %in% c("poisson","binomial","ordinal","multinom"))) combined.fit.mcmc[,sel.thetas2[family %in% c("poisson","binomial","ordinal","multinom")]] <- 0
 
  	## For any multinomial columns, set the corresponding rows in X.coefs to zero
  	if(any(family == "multinom") & X.eff == 1) {
  		for(k in index.multinom.cols) {
  			sel.multinom.col <- grep(paste("X.params\\[",k,",+",sep=""),colnames(combined.fit.mcmc))
  			combined.fit.mcmc[,sel.multinom.col] <- 0
  			}
  		}
 	
 	## Make output beautiful
 	if(is.null(colnames(y))) colnames(y) <- 1:ncol(y); if(is.null(rownames(y))) rownames(y) <- 1:nrow(y)
 	if(X.eff) { if(is.null(colnames(X))) colnames(X) <- 1:ncol(X); if(is.null(rownames(X))) rownames(X) <- 1:nrow(X) }
 
 	out.fit <- list(
 		lv.coefs.median = matrix(apply(combined.fit.mcmc[,grep("all.params",colnames(combined.fit.mcmc))],2,median),nrow=p),
 		lv.coefs.iqr = matrix(apply(combined.fit.mcmc[,grep("all.params",colnames(combined.fit.mcmc))],2,IQR),nrow=p),
 
 		lv.coefs.mean = matrix(apply(combined.fit.mcmc[,grep("all.params",colnames(combined.fit.mcmc))],2,mean),nrow=p), 
 		lv.coefs.sd = matrix(apply(combined.fit.mcmc[,grep("all.params",colnames(combined.fit.mcmc))],2,sd),nrow=p)
 		)
 	rownames(out.fit$lv.coefs.median) <- rownames(out.fit$lv.coefs.iqr) <- rownames(out.fit$lv.coefs.mean) <- rownames(out.fit$lv.coefs.sd) <- colnames(y)
 
 	if(num.lv > 0) {
 		out.fit$lv.median = matrix(apply(combined.fit.mcmc[,grep("lvs", colnames(combined.fit.mcmc))],2,median),nrow=n)
 		out.fit$lv.iqr = matrix(apply(combined.fit.mcmc[,grep("lvs", colnames(combined.fit.mcmc))],2,IQR),nrow=n)
		out.fit$lv.mean = matrix(apply(combined.fit.mcmc[,grep("lvs", colnames(combined.fit.mcmc))],2,mean),nrow=n)
		out.fit$lv.sd = matrix(apply(combined.fit.mcmc[,grep("lvs", colnames(combined.fit.mcmc))],2,sd),nrow=n) 
 
 		rownames(out.fit$lv.median) <- rownames(out.fit$lv.iqr) <- rownames(out.fit$lv.mean) <- rownames(out.fit$lv.sd) <- rownames(y)
 		colnames(out.fit$lv.median) <- colnames(out.fit$lv.iqr) <- colnames(out.fit$lv.mean) <- colnames(out.fit$lv.sd) <- paste("LV",1:num.lv,sep="")
 		colnames(out.fit$lv.coefs.median) <- colnames(out.fit$lv.coefs.iqr) <- colnames(out.fit$lv.coefs.mean) <- colnames(out.fit$lv.coefs.sd) <- c("theta0",paste("theta",1:num.lv,sep=""),"Dispersion")
 		}	
 	if(num.lv == 0) {	
 		colnames(out.fit$lv.coefs.median) <- colnames(out.fit$lv.coefs.iqr) <- colnames(out.fit$lv.coefs.mean) <- colnames(out.fit$lv.coefs.sd) <- c("theta0","Dispersion")
 		}
 		
 	if(site.eff) {
 		out.fit$site.coefs.median <- apply(combined.fit.mcmc[,grep("site.params", colnames(combined.fit.mcmc))],2,median)
 		out.fit$site.coefs.iqr <- apply(combined.fit.mcmc[,grep("site.params", colnames(combined.fit.mcmc))],2,IQR)
 		out.fit$site.coefs.mean <- apply(combined.fit.mcmc[,grep("site.params", colnames(combined.fit.mcmc))],2,mean)
 		out.fit$site.coefs.sd <- apply(combined.fit.mcmc[,grep("site.params", colnames(combined.fit.mcmc))],2,sd)
 
 		names(out.fit$site.coefs.median) <- names(out.fit$site.coefs.iqr) <- names(out.fit$site.coefs.mean) <- names(out.fit$site.coefs.sd) <- rownames(y) 
 		}
 
 	if(X.eff & any(family != "multinom")) {
 		out.fit$X.coefs.median <- matrix(apply(combined.fit.mcmc[,grep("X.params", colnames(combined.fit.mcmc))],2,median),nrow=p)
 		out.fit$X.coefs.iqr <- matrix(apply(combined.fit.mcmc[,grep("X.params", colnames(combined.fit.mcmc))],2,IQR),nrow=p)
 		out.fit$X.coefs.mean <- matrix(apply(combined.fit.mcmc[,grep("X.params", colnames(combined.fit.mcmc))],2,mean),nrow=p)
 		out.fit$X.coefs.sd <- matrix(apply(combined.fit.mcmc[,grep("X.params", colnames(combined.fit.mcmc))],2,sd),nrow=p)
 
 		rownames(out.fit$X.coefs.median) <- rownames(out.fit$X.coefs.iqr) <- rownames(out.fit$X.coefs.mean) <- rownames(out.fit$X.coefs.sd) <- colnames(y) 
 		colnames(out.fit$X.coefs.median) <- colnames(out.fit$X.coefs.iqr) <- colnames(out.fit$X.coefs.mean) <- colnames(out.fit$X.coefs.sd) <- colnames(X) 
 		}
 
  	if(X.eff & any(family == "multinom")) {
  		out.fit$X.multinom.coefs.median <- array(apply(combined.fit.mcmc[,grep("X.multinom.params", colnames(combined.fit.mcmc))],2,median),dim=c(length(index.multinom.cols),num.X,num.multinom.levels))
  		out.fit$X.multinom.coefs.iqr <- array(apply(combined.fit.mcmc[,grep("X.multinom.params", colnames(combined.fit.mcmc))],2,IQR),dim=c(length(index.multinom.cols),num.X,num.multinom.levels))
  		out.fit$X.multinom.coefs.mean <- array(apply(combined.fit.mcmc[,grep("X.multinom.params", colnames(combined.fit.mcmc))],2,mean),dim=c(length(index.multinom.cols),num.X,num.multinom.levels))
  		out.fit$X.multinom.coefs.sd <- array(apply(combined.fit.mcmc[,grep("X.multinom.params", colnames(combined.fit.mcmc))],2,sd),dim=c(length(index.multinom.cols),num.X,num.multinom.levels))
  
  		dimnames(out.fit$X.multinom.coefs.median) <- dimnames(out.fit$X.multinom.coefs.iqr) <- dimnames(out.fit$X.multinom.coefs.mean) <- dimnames(out.fit$X.multinom.coefs.sd) <- list("1" = index.multinom.cols, "2" = colnames(X), "level" = 1:num.multinom.levels)
  		}
 
 	if(any(family == "ordinal")) {
 		out.fit$cutoffs.median <- apply(combined.fit.mcmc[,grep("alpha", colnames(combined.fit.mcmc))],2,median)
 		out.fit$cutoffs.iqr <- apply(combined.fit.mcmc[,grep("alpha", colnames(combined.fit.mcmc))],2,IQR)
 		out.fit$cutoffs.mean <- apply(combined.fit.mcmc[,grep("alpha", colnames(combined.fit.mcmc))],2,mean)
 		out.fit$cutoffs.sd <- apply(combined.fit.mcmc[,grep("alpha", colnames(combined.fit.mcmc))],2,sd)
 		
 		names(out.fit$cutoffs.median) <- names(out.fit$cutoffs.iqr) <- names(out.fit$cutoffs.mean) <- names(out.fit$cutoffs.sd) <- paste(1:(num.ord.levels-1),"|",2:num.ord.levels,sep="")
 		}
 		
 	if(any(family == "tweedie")) {
 		out.fit$powerparam.median <- median(combined.fit.mcmc[,grep("powerparam", colnames(combined.fit.mcmc))])
 		out.fit$powerparam.iqr <- IQR(combined.fit.mcmc[,grep("powerparam", colnames(combined.fit.mcmc))])
 		out.fit$powerparam.mean <- mean(combined.fit.mcmc[,grep("powerparam", colnames(combined.fit.mcmc))])
 		out.fit$powerparam.sd <- sd(combined.fit.mcmc[,grep("powerparam", colnames(combined.fit.mcmc))])
 		
 		names(out.fit$powerparam.median) <- names(out.fit$powerparam.iqr) <- names(out.fit$powerparam.mean) <- names(out.fit$powerparam.sd) <- "Common power parameter"
 		}
 		
	get.hpds <- get.hpdintervals(y, X, combined.fit.mcmc, num.lv)
 	out.fit$hpdintervals <- get.hpds
 	
 	if(calc.ics) { 
 		cat("Calculating Information criteria\n")
 		get.ics <- get.measures(y, X, family, complete.trial.size, site.eff, num.lv, combined.fit.mcmc)
 		ics <- c(get.dic(jagsfit), get.ics$waic, get.ics$eaic, get.ics$ebic, get.ics$aic.median, get.ics$bic.median, get.ics$comp.lm, get.ics$num.params)
 		names(ics) <- c("Conditional DIC","WAIC","EAIC","EBIC","AIC at post. median", "BIC at post. median","Compound L-M at post. median","# of parameters")
		out.fit$ics <- ics 		
 		}
 
 	if(save.model) { out.fit$jags.model <- jagsfit }
 
 	out.fit$call <- match.call()
 	out.fit$n <- n; out.fit$p <- p; 
	#if(any(family == "tweedie")) y[,which(family == "tweedie")] <- y[,which(family == "tweedie")] - 0.001
 	if(!is.null(X)) out.fit$X <- X; if(is.null(X)) out.fit$X <- NULL
 	out.fit$y <- y; 
 	out.fit$family <- family; 
 	out.fit$num.lv <- num.lv; 
 	out.fit$num.X <- num.X; 
	out.fit$site.eff <- site.eff; 
 	out.fit$calc.ics <- calc.ics
 	out.fit$trial.size <- complete.trial.size
 	out.fit$hypparams <- hypparams
 	out.fit$num.ord.levels <- num.ord.levels
 	#out.fit$num.multinom.levels <- num.multinom.levels; 
 	
 	class(out.fit) <- "boral"
 	if(!save.model) { if(file.exists("jagsboralmodel.txt")) file.remove("jagsboralmodel.txt") }
 
 	return(out.fit) }	

################	
lvsplot <- function(x, jitter = FALSE, a = 1, newplot = TRUE, ...) {
 	if(x$num.lv > 2) stop("Manual plotting required for plotting beyond 2 latent variables")
 	if(x$num.lv == 0) stop("No latent variables to plot.")
 
 	n <- nrow(x$lv.med)
 	if(newplot == TRUE)
 		par(cex = a, cex.axis = a, cex.lab = a+0.5, mar = c(5,5,3,1), mfrow = c(1,2), cex.main = a+0.5, ...) 
 	if(newplot == FALSE)
 		par(cex = a, cex.axis = a, cex.lab = a+0.5, mar = c(5,5,3,1), cex.main = a+0.5, ...) 
 
 	if(x$num.lv == 1) {
 		plot(1:n, x$lv.med, xlab = "Row number", ylab = "Latent variable 1", main = "Plot of the latent variable posterior medians", cex = 1.2*a, type = "n", ...)
 		text(x = 1:n, y = x$lv.med, label = 1:n, cex = 1.2*a)
 		plot(1:n, x$lv.mean, xlab = "Row number", ylab = "Latent variable 1", main = "Plot of the latent variable posterior means", cex = 1.2*a, type = "n", ...)
 		text(x = 1:n, y = x$lv.mean, label = 1:n, cex = 1.2*a) 
 		}	
 	if(x$num.lv == 2) {
		plot(x$lv.med, xlab = "Latent variable 1", ylab = "Latent variable 2", main = "Plot of the latent variable posterior medians", cex = 1.2*a, type = "n", ...)
 		if(!jitter) text(x$lv.med, label = 1:n, cex = 1.2*a)
 		if(jitter) text(jitter(x$lv.med[,1]), jitter(x$lv.med[,2]), label = 1:n, cex = 1.2*a)
 		plot(x$lv.mean, xlab = "Latent variable 1", ylab = "Latent variable 2", main = "Plot of the latent variable posterior means", cex = 1.2*a, type = "n", ...)
 		if(!jitter) text(x$lv.mean, label = 1:n, cex = 1.2*a)
 		if(jitter) text(jitter(x$lv.mean[,1]), jitter(x$lv.mean[,2]), label = 1:n, cex = 1.2*a)
 		}	
 	}

print.boral <- function(x, ...) {
 	cat("Call:\n")
 	print(x$call)
 	cat("\n")
 	cat("Response matrix attributes\n \t# of rows:", x$n, "\n\t# of columns:", x$p, "\n") 
 	cat("Model attributes\n \tColumn families:", x$family, "\n\t# of latent variables:", x$num.lv, "\n\tRow effect included (0=no/1=yes)?", as.numeric(x$site.eff), "\n") 
 	if(any(x$family == "binomial")) cat("Trial sizes used (columns with binomial families):", x$trial.size,"\n")
 	if(any(x$family == "ordinal")) cat("Number of levels for ordinal data:", x$num.ord.levels,"\n")
 	#cat("Hyperparameters (variance in normal priors of coefficients, maximum in uniform prior for dispersion):", x$hypparams, "\n") 
 	if(x$num.X > 0) cat("Model matrix with", x$num.X, "covariates also fitted\n\n")
# 	cat("Output attributes\n")
# 	print(attributes(x))
 	}
	
print.summary.boral <- function(x, ...) {
 	cat("Call:\n")
 	print(x$call)
 	cat("\n")
 	
 	if(x$est == "median") { cat("Median point estimates\n\n LV coefficients (thetas) and dispersion parameter\n"); print(x$coefficients); cat("\n") }
 	if(x$est == "mean") { cat("Mean point estimates\n\n LV coefficients (thetas) and dispersion parameter\n"); print(x$coefficients); cat("\n") }		
 	
 	if(!is.null(x$site.coefficients)) { cat("Row coefficients\n"); print(x$site.coefficients); cat("\n") }
 	if(!is.null(x$X.coefficients)) { cat("X coefficients (betas)\n"); print(x$X.coefficients); cat("\n") }
 	if(!is.null(x$X.multinom.coefficients)) { cat("There are also coefficients corresponding to multinomial columns which have not been printed"); }
 	
 	if(any(x$family == "ordinal")) { cat("Proportional odds (Cumulative logit) regression intercepts\n"); print(x$cutoffs); cat("\n") }
 	if(any(x$family == "tweedie")) { cat("Tweedie power parameter \n"); print(x$powerparam); cat("\n") }
 
 	if(x$calc.ics) {
 		cat("DIC (pD = var(deviance)/2):", as.vector(unlist(x$ics[1])), "\n")
 		cat("WAIC:", as.vector(unlist(x$ics[2])), "\n")
 		cat("EAIC:", as.vector(unlist(x$ics[3])), "\n")
 		cat("EBIC:", as.vector(unlist(x$ics[4])), "\n")
 		cat("AIC at posterior median:", as.vector(unlist(x$ics[5])), "\n")
 		cat("BIC at posterior median:", as.vector(unlist(x$ics[6])), "\n")
 		cat("Compound Laplace-Metropolis estimator at posterior median:", as.vector(unlist(x$ics[7])), "\n")
		}
 	}	
		
summary.boral <- function(object, est = "median", ...) {
 	if(est == "median") {
 		gather.output <- list(call = object$call, coefficients = round(object$lv.coefs.med,3))
 		if(object$site.eff == TRUE) gather.output$site.coefficients = round(object$site.coefs.med,3)
 		if(object$num.X > 0) gather.output$X.coefficients = round(object$X.coefs.med,3)
 		if(any(object$family == "ordinal")) gather.output$cutoffs = round(object$cutoffs.median,3)
 		if(any(object$family == "tweedie")) gather.output$powerparam = round(object$powerparam.median,3)
 		if(!is.null(object$X.multinom.coefs.median)) gather.output$X.multinom.coefficients = round(object$X.multinom.coefs.median,3) }
 
 	if(est == "mean") {
 		gather.output <- list(call = object$call, coefficients = round(object$lv.coefs.mean,3))
 		if(object$site.eff == TRUE) gather.output$site.coefficients = round(object$site.coefs.mean,3)
 		if(object$num.X > 0) gather.output$X.coefficients = round(object$X.coefs.mean,3)
 		if(any(object$family == "ordinal")) gather.output$cutoffs = round(object$cutoffs.mean,3)
 		if(any(object$family == "tweedie")) gather.output$powerparam = round(object$powerparam.mean,3)
 		if(!is.null(object$X.multinom.coefs.mean)) gather.output$X.multinom.coefficients = round(object$X.multinom.coefs.mean,3) }
 
 	gather.output$est <- est
 	gather.output$calc.ics <- object$calc.ics
	gather.output$trial.size <- object$trial.size
 	gather.output$num.ord.levels <- object$num.ord.levels
 	
 	if(object$calc.ics) gather.output$ics = object$ics
 	class(gather.output) <- "summary.boral"
 	gather.output 
 	}
 			
plot.boral <- function(x, est = "median", jitter = FALSE, a = 1, ...) {
 	if(all(x$family %in% c("ordinal","multinom"))) 
 		stop("Residuals are not defined, and therefore residual analysis cannot be performed, if all columns of y are ordinal")
 	get.mus <- fitted.boral(x, est = est)$out
 	get.etas <- get.mus
 	get.ds.res <- ds.residuals(object = x, est = est)$residuals
 	for(j in 1:ncol(x$y)) {
 		if(x$family[j] == "binomial") get.etas[,j] <- log(get.mus[,j]/(1-get.mus[,j]))
 		if(x$family[j] %in% c("poisson","lnormal","negative.binomial","tweedie")) get.etas[,j] <- log(get.mus[,j]+1e-5)
 		if(x$family[j] == "normal") get.etas[,j] <- (get.mus[,j]) }
 
 	par(ask = T, cex = a, mar = c(5,5,1,1), cex.main = a) 
 	palette(rainbow(ncol(get.etas)))
 
 	matplot(get.etas, get.ds.res, ylab = "Dunn-Smyth Residuals", xlab = "Linear Predictors", type="n")
 	for(i in 1:ncol(get.etas)) { points(get.etas[,i], get.ds.res[,i], col=palette()[i]) }
 	abline(h=0, lty = 2, lwd = 2)
	
# 	matplot(get.mus, get.ds.res, ylab = "Dunn-Smyth Residuals", xlab = "Fitted Values", type="n")
# 	for(i in 1:ncol(get.mus)) { points(get.mus[,i], get.ds.res[,i], col=palette()[i]) }
# 	abline(h=0, lty = 2, lwd = 2)father

	matplot(get.ds.res, ylab = "Dunn-Smyth Residuals", xlab = "Row Index",type="n", xaxt = "n")
 	axis(side = 1, at = 1:nrow(x$y), labels = rownames(x$lv.mean), cex.axis = 1)
 	for (i in 1:ncol(get.mus)) { points(seq(1,nrow(x$y)),get.ds.res[,i], col=palette()[i]) }
 	abline(0,0,lty=2)
 
 	matplot(t(get.ds.res), ylab = "Dunn-Smyth Residuals", xlab = "Column Index", type="n", xaxt = "n")
 	axis(side = 1, at = 1:ncol(x$y), labels = rownames(x$coefs.mean), cex.axis = 1)
 	for(i in 1:ncol(get.mus)) { points(rep(i,nrow(get.etas)), get.ds.res[,i], col=palette()[i]) }
 	abline(h=0, lty = 2, lwd = 2)
 
	get.ds.res2 <- as.vector(unlist(get.ds.res))
 	qqnorm(get.ds.res2[is.finite(get.ds.res2)], main = "Normal Quantile Plot")
 	
 	palette("default")
 	}
 	
###############
## For some reason I cannot import mcmc function from coda 
###############
# fh.mcmc <- function (data = NA, start = 1, end = numeric(0), thin = 1) 
# {
#     if (is.matrix(data)) {
#         niter <- nrow(data)
#         nvar <- ncol(data)
#     }
#     else if (is.data.frame(data)) {
#         if (!all(sapply(data, is.numeric))) {
#             stop("Data frame contains non-numeric values")
#         }
#         data <- as.matrix(data)
#         niter <- nrow(data)
#         nvar <- ncol(data)
#     }
#     else {
#         niter <- length(data)
#         nvar <- 1
#     }
#     thin <- round(thin)
#     if (length(start) > 1) 
#         stop("Invalid start")
#     if (length(end) > 1) 
#         stop("Invalid end")
#     if (length(thin) != 1) 
#         stop("Invalid thin")
#     if (missing(end)) 
#         end <- start + (niter - 1) * thin
#     else if (missing(start)) 
#         start <- end - (niter - 1) * thin
#     nobs <- floor((end - start)/thin + 1)
#     if (niter < nobs) 
#         stop("Start, end and thin incompatible with data")
#     else {
#         end <- start + thin * (nobs - 1)
#         if (nobs < niter) 
#             data <- data[1:nobs, , drop = FALSE]
#     }
#     data
# }

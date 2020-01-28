coefsplot <- function(covname, x, labely = NULL, est = "median", ...) 
     {
     if(!is.null(labely)) if(!(length(labely) == nrow(x$X.coefs.median) || length(labely)==1)) 
          stop("If labely is not NULL, then it must be either of length one or a vector as long as the number of rows in x$X.coefs.median (number of species). Thanks!")
     if(!(covname %in% colnames(x$X.coefs.mean))) 
          stop("covname not found among the covariates in the boral object x")
     
     col.seq <- rep("black", length(x$hpdintervals$X.coefs[,covname,"lower"]))
     col.seq[x$hpdintervals$X.coefs[,covname,"lower"] < 0 & x$hpdintervals$X.coefs[,covname,"upper"] > 0] <- "grey"
     
     At.y <- rev(1:nrow(x$X.coefs.median)) ## So that the graphs plots in the same order as rownames of x$X.coefs.median

     if(est == "median")
          plot(x = x$X.coefs.median[,covname], y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = covname, xlim = c(min(x$hpdintervals$X.coefs[,covname,"lower"]), max(x$hpdintervals$X.coefs[,covname,"upper"])), pch = "x", ...)
     if(est == "mean")
          plot(x = x$X.coefs.mean[,covname], y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = covname, xlim = c(min(x$hpdintervals$X.coefs[,covname,"lower"]), max(x$hpdintervals$X.coefs[,covname,"upper"])), pch = "x", ...)

     segments(x0 = x$hpdintervals$X.coefs[,covname,"lower"], y0 = At.y, x1 = x$hpdintervals$X.coefs[,covname,"upper"], y1 = At.y, col = col.seq, ...)  
     abline(v=0, lty=3)
     
     if(is.null(labely)) 
          axis(2, at=At.y, labels = rownames(x$X.coefs.mean), las=1, ...) 
     if(!is.null(labely)) 
          {
          if(length(labely) == nrow(x$X.coefs.mean)) axis(2, at=At.y, labels=labely, las=1, ...) 
          if(length(labely) == 1) mtext(text = labely, side = 2, line = 3, las = 3, ...)
          } 

     if(exists("ssvs.indcoefs.mean", x)) 
          {
          message("Posterior probabilities of inclusion for ", covname, ":")
          print(round(x$ssvs.indcoefs.mean[,covname],3))
          message()
          }
     }

		

plot.boral <- function(x, est = "median", jitter = FALSE, ...) 
     {
     #if(all(x$family %in% c("ordinal","multinom"))) stop("Residuals are not defined, and therefore residual analysis cannot be performed, if all columns of y are ordinal")
     if(any(x$family %in% c("ordinal","multinom"))) 
          warning("For all columns of y that are ordinal, the first plot constructed is of Dunn-Smyth residuals against fitted values (i.e., the level with the highest predicted probability). Note this can make things very confusing to interpret if only SOME of the columns in the response matrix are ordinal")
     get_mus <- fitted.boral(x, est = est)$out
     get_etas <- get_mus
     get_ds_res <- ds.residuals(object = x, est = est)
     print(get_ds_res$agree.ordinal)
     get_ds_res <- get_ds_res$residuals
          
     for(j in 1:ncol(x$y)) 
               {
               if(x$family[j] %in% c("beta")) 
                    get_etas[,j] <- log((get_mus[,j]+1e-5)/(1-get_mus[,j]+1e-5))
               if(x$family[j] %in% c("binomial")) 
                    get_etas[,j] <- qnorm(get_mus[,j]+1e-5)
               if(x$family[j] %in% c("poisson","lnormal","negative.binomial","tweedie","exponential")) 
                    get_etas[,j] <- log(get_mus[,j]+1e-5)
               if(x$family[j] == "normal") 
                    get_etas[,j] <- (get_mus[,j]) 
               if(x$family[j] == "ordinal") { } ## Fitted values are the class with highest probability, which is already attained from fitted.boral
               }


     #.pardefault <- par(no.readonly = TRUE)	
     #par(ask = TRUE, cex = a, mar = c(5,5,2,1), cex.lab = 0.8*a, cex.main = a, las = 1, ...) 
     palette(rainbow(ncol(get_etas)))

     matplot(get_etas, get_ds_res, ylab = "Dunn-Smyth Residuals", xlab = "Linear Predictors", type="n", ...)
     for(i in 1:ncol(get_etas)) 
          points(get_etas[,i], get_ds_res[,i], col=palette()[i], ...)
     abline(h=0, lty = 2, lwd = 2)
     
     # 	matplot(get_mus, get_ds_res, ylab = "Dunn-Smyth Residuals", xlab = "Fitted Values", type="n")
     # 	for(i in 1:ncol(get_mus)) { points(get_mus[,i], get_ds_res[,i], col=palette()[i]) }
     # 	abline(h=0, lty = 2, lwd = 2)father

     matplot(get_ds_res, ylab = "Dunn-Smyth Residuals", xlab = "Row index",type="n", xaxt = "n", ...)
     axis(side = 1, at = 1:nrow(x$y), labels = rownames(x$lv.mean), ...)
     for (i in 1:ncol(get_mus)) 
          { 
          points(seq(1,nrow(x$y)),get_ds_res[,i], col=palette()[i], ...) 
          }
     abline(0,0,lty=2)

     matplot(t(get_ds_res), ylab = "Dunn-Smyth Residuals", xlab = "Column index", type="n", xaxt = "n", ...)
     axis(side = 1, at = 1:ncol(x$y), labels = rownames(x$coefs.mean), ...)
     for(i in 1:ncol(get_mus)) 
          { 
          points(rep(i,nrow(get_etas)), get_ds_res[,i], col=palette()[i], ...) 
          }
     abline(h=0, lty = 2, lwd = 2)

     get_ds_res2 <- as.vector(unlist(get_ds_res))
     qqnorm(get_ds_res2[is.finite(get_ds_res2)], main = "Normal Quantile Plot", ...)
     #qqline(y = get_ds_res2[is.finite(get_ds_res2)], ...)
     
     palette("default")
     #par(.pardefault) 	
     }

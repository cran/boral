coefsplot   <- function(covname, object, fourthcorner = FALSE, labely = NULL, est = "median", ...) {
     if(fourthcorner == TRUE) {
          covname <- labely <- NULL
          }     
     if(!exists("traits.coefs.median", object) & fourthcorner == TRUE) {
          warning("Fourth corner coefficients (traits.coefs) were not found in the object; reverting to fourthcorner = FALSE...")
          fourthcorner <- FALSE
          }

     if(fourthcorner == FALSE) {
          if(!is.null(labely)) if(!(length(labely) == nrow(object$X.coefs.median) || length(labely)==1)) 
               stop("If labely is not NULL, then it must be either of length one or a vector as long as the number of rows in object$X.coefs.median (number of response). Thanks!")
          if(!(covname %in% colnames(object$X.coefs.mean))) 
               stop("covname not found among the covariates in the boral object.")
          }
               
     ## Caterpillar plot of selected X coefficients
     if(fourthcorner == FALSE) {
          col.seq <- rep("black", length(object$hpdintervals$X.coefs[,covname,"lower"]))
          col.seq[object$hpdintervals$X.coefs[,covname,"lower"] < 0 & object$hpdintervals$X.coefs[,covname,"upper"] > 0] <- "grey"
          
          At.y <- rev(1:nrow(object$X.coefs.median)) ## So that the graphs plots in the same order as rownames of object$X.coefs.median

          if(est == "median")
               plot(x = object$X.coefs.median[,covname], y = At.y, yaxt = "n", ylab = "", xlab = "", col = col.seq, main = covname, 
                    xlim = c(min(object$hpdintervals$X.coefs[,covname,"lower"]), max(object$hpdintervals$X.coefs[,covname,"upper"])), pch = "x", ...)
          if(est == "mean")
               plot(x = object$X.coefs.mean[,covname], y = At.y, yaxt = "n", ylab = "", xlab = "", col = col.seq, main = covname, 
                    xlim = c(min(object$hpdintervals$X.coefs[,covname,"lower"]), max(object$hpdintervals$X.coefs[,covname,"upper"])), pch = "x", ...)

          segments(x0 = object$hpdintervals$X.coefs[,covname,"lower"], y0 = At.y, x1 = object$hpdintervals$X.coefs[,covname,"upper"], y1 = At.y, col = col.seq, ...)  
          abline(v=0, lty=3)
          
          if(is.null(labely)) 
               axis(2, at=At.y, labels = rownames(object$X.coefs.mean), las=1, ...) 
          if(!is.null(labely)) {
               if(length(labely) == nrow(object$X.coefs.mean)) 
                    axis(2, at = At.y, labels=labely, las=1, ...) 
               if(length(labely) == 1) 
                    mtext(text = labely, side = 2, line = 3, las = 3, ...)
               } 

          if(exists("ssvs.indcoefs.mean", object)) {
               message("Posterior probabilities of inclusion for ", covname, ":")
               print(round(object$ssvs.indcoefs.mean[,covname],3))
               message()
               }
          }


     ## Caterpillar plot of all fourth corner coefficients 
     if(fourthcorner == TRUE) {
          col.seq <- rep("black", length(object$hpdintervals$traits.coefs[,,"lower"]))
          col.seq[object$hpdintervals$traits.coefs[,,"lower"] < 0 & object$hpdintervals$traits.coefs[,,"upper"] > 0] <- "grey"
                    
          if(est == "median") {
               gettable <- as.data.frame.table(object$traits.coefs.median[,1:(object$num.traits+1)])
               getlower <- as.data.frame.table(object$hpdintervals$traits.coefs[,1:(object$num.traits+1),"lower"])
               getupper <- as.data.frame.table(object$hpdintervals$traits.coefs[,1:(object$num.traits+1),"upper"])
               }
          if(est == "mean") {
               gettable <- as.data.frame.table(object$traits.coefs.mean[,1:(object$num.traits+1)])
               getlower <- as.data.frame.table(object$hpdintervals$traits.coefs[,1:(object$num.traits+1),"lower"])
               getupper <- as.data.frame.table(object$hpdintervals$traits.coefs[,1:(object$num.traits+1),"upper"])
               }
               
          At.y <- rev(1:nrow(gettable))
          plot(x = gettable$Freq, y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = "", xlim = c(min(getlower$Freq), max(getupper$Freq)), pch = "x", 
               main = "Fourth-corner coefficient plot", ...)
          segments(x0 = getlower$Freq, y0 = At.y, x1 = getupper$Freq, y1 = At.y, col = col.seq, ...)  
          abline(v=0, lty=3)
          
          #axis(2, at = At.y, labels = apply(gettable[,1:2], 1, paste0, collapse=":"), las=1, ...) 
          axis(2, at = At.y, labels = gettable[,1], las=1, ...) 
          axis(4, at = At.y, labels = gettable[,2], las=1, ...) 
          

          if(exists("ssvs.traitscoefs.mean", object)) {
               message("Posterior probabilities of inclusion for fourth corner coefficients")
               print(round(object$ssvs.traitscoefs.mean,3))
               message()
               }
          }
     
     }
     
     
ranefsplot   <- function(sel.spp = NULL, object, ordered = FALSE, est = "median", ...) {
     if(!exists("ranef.ids", object)) {
          stop("Response-specific random intercepts were not found in the object.")
          }
     
     n <- nrow(object$y)
     p <- ncol(object$y)
     if(is.null(sel.spp)) { 
          sel.spp <- colnames(object$y)
          message("Random intercept plots for all response, one page at a time.") 
          }
     if(!is.null(sel.spp)) { 
          sel.spp <- as.vector(sel.spp)
          if(any(duplicated(sel.spp))) 
               stop("sel.spp should only contain unique values")
          
          if(is.numeric(sel.spp)) {
               if(any(sel.spp > p) || any(sel.spp < 1)) 
                    stop("If sel.spp is a numeric vector, then it should be an index of which columns of object$y are to be plotted.")
               
               final_spp_names <- colnames(object$y)[sel.spp]
               }
          if(is.character(sel.spp)) {
               if(!all(sel.spp %in% colnames(object$y))) 
                    stop("If sel.spp is a vector of response names, then it should only contain elements in colnames(object$y).")
               
               final_spp_names <- sel.spp
               }
          
          message("Only random intercept plots for response ", paste0(sel.spp, collapse = "/"), " will be constructed.") 
          }

          
     ## Caterpillar plot of (selected) predicted random intercepts
     how_many_ranefs <- ncol(object$ranef.ids)
     optimal_mfrow <- n2mfrow(how_many_ranefs)
    
     par(las = 1, ask = TRUE, mfrow = optimal_mfrow)
     for(k0 in 1:length(final_spp_names)) {
          sel_spp <- which(final_spp_names[k0] == colnames(object$y))
          
          for(k1 in 1:how_many_ranefs) {
               if(est == "median") 
                    cw_spp_estimates <- as.vector(object$ranef.coefs.median[[k1]][sel_spp,])
               if(est == "mean")
                    cw_spp_estimates <- as.vector(object$ranef.coefs.mean[[k1]][sel_spp,])
                    
               cw_hpd_intervals <- object$hpdintervals$ranef.coefs[[k1]][sel_spp,,]
               col_seq <- rep("black", nrow(cw_hpd_intervals))
               col_seq[cw_hpd_intervals[,"lower"] < 0 & cw_hpd_intervals[,"upper"] > 0] <- "grey"
               
               if(ordered) {
                    get_order <- order(cw_spp_estimates)
                    cw_spp_estimates <- cw_spp_estimates[get_order]
                    cw_hpd_intervals <- cw_hpd_intervals[get_order,]
                    }
               
               At_y <- rev(1:length(cw_spp_estimates)) 
               plot(x = cw_spp_estimates, y = At_y, yaxt = "n", ylab = "", xlab = "", col = col_seq, main = "", 
                    xlim = c(min(cw_hpd_intervals), max(cw_hpd_intervals)), pch = "x", ...)
               segments(x0 = cw_hpd_intervals[,"lower"], y0 = At_y, x1 = cw_hpd_intervals[,"upper"], y1 = At_y, col = col_seq, ...)  
               abline(v = 0, lty = 3)
               axis(2, at = At_y, labels = rownames(cw_hpd_intervals), las=1, ...) 
               mtext(colnames(object$ranef.ids)[k1], side = 3, line = 0.7)
               }
          
          mtext(final_spp_names[k0], side = 3, cex = 1.5, line = 2, font = 4)
          }
     
     }
     
     

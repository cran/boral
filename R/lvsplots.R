lvsplot <- function(x, jitter = FALSE, biplot = TRUE, ind.spp = NULL, alpha = 0.5, main = NULL, 
     est = "median", which.lvs = c(1,2), return.vals = FALSE, ...) {

    if(x$num.lv == 0) stop("No latent variables to plot")
    if(length(which.lvs) != 2) stop("which.lvs should be a vector of length 2, indicating which axes to plot. This argument is ignored if x$num.lv = 1")
    if(x$num.lv > 2 & any(which.lvs > x$num.lv)) stop("Fewer latent variables than chosen by which.lvs")

    n <- nrow(x$lv.median); p <- nrow(x$lv.coefs.median)
    if(!is.null(ind.spp)) { if(ind.spp > p) { ind.spp <- p } }
    if(biplot == TRUE & !is.null(ind.spp)) { 
        message("Only the first ", ind.spp, " ``most important'' latent variable coefficients included in biplot") 
        }
    if(biplot == TRUE & is.null(ind.spp)) { 
        ind.spp <- p; message("All latent variable coefficients included in biplot") 
        }

    if(x$num.lv == 1) {
        choose_lvs <- x$lv.median; 
        choose_lvs_coefs <- x$lv.coefs.median[,2]
        main <- "Plot of latent variable posterior medians"
        if(est == "mean") { 
            choose_lvs <- x$lv.mean
            choose_lvs_coefs <- x$lv.coefs.mean[,2] 
            }

        if(!biplot) {
            if(is.null(main) & est == "median") 
                main <- "Plot of latent variable posterior medians" 
            if(is.null(main) & est == "mean") 
                main <- "Plot of latent variable posterior means"
            plot(1:n, choose_lvs, xlab = "Row index", ylab = "Latent variable 1", main = main, type = "n", ...)
            if(!jitter) 
                text(x = 1:n, y = x$lv.median, label = 1:n, ...)
            if(jitter) 
                text(x = 1:n, y = jitter(x$lv.median), label = 1:n, ...)
            }

        if(biplot) {
            if(is.null(main) & est == "median") 
                main <- "Biplot of latent variable posterior medians"
            if(is.null(main) & est == "mean") 
                main <- "Biplot of latent variable posterior means" 
            plot(1:n, choose_lvs, xlab = "Row index", ylab = "Latent variable 1", main = main, type = "n", ...)
            if(!jitter) 
                text(x = 1:n, y = x$lv.median, label = 1:n, ...)
            if(jitter) 
                text(x = 1:n, y = jitter(x$lv.median), label = 1:n, ...)
            text(choose_lv_coefs, label = rownames(x$lv.coefs.mean), col = "red", ...)	
            }
        }


    if(x$num.lv > 1) {
        testcov <- tcrossprod(x$lv.median, x$lv.coefs.median[,2:(x$num.lv+1)])
        if(est == "mean") 
            testcov <- tcrossprod(x$lv.mean, x$lv.coefs.mean[,2:(x$num.lv+1)])

        do_svd <- svd(testcov, x$num.lv, x$num.lv)   		
        choose_lvs <- scale(do_svd$u*matrix(do_svd$d[1:x$num.lv]^alpha,nrow=x$n,ncol=x$num.lv,byrow=TRUE),center = TRUE, scale = FALSE)
        choose_lv_coefs <- scale(do_svd$v*matrix(do_svd$d[1:x$num.lv]^(1-alpha),nrow=x$p,ncol=x$num.lv,byrow=TRUE), center = TRUE, scale = FALSE)
        
        if(!biplot) {
            if(is.null(main) & est == "median") 
                main <- "Plot of latent variable posterior medians"
            if(is.null(main) & est == "mean") 
                main <- "Plot of latent variable posterior meas, axis"
            plot(choose_lvs, xlab = paste("Latent variable", which.lvs[1]), ylab = paste("Latent variable", which.lvs[2]), main = main, type = "n", ...)
            if(!jitter) 
                text(choose_lvs[,which.lvs], label = 1:n, ...)
            if(jitter) 
                text(jitter(choose_lvs[,which.lvs[1]]), jitter(choose_lvs[,which.lvs[2]]), label = 1:n, ...)
            }

        if(biplot) {
            if(is.null(main) & est == "median") 
                main <- "Biplot of latent variable posterior medians"
            if(is.null(main) & est == "mean") 
                main <- "Biplot of latent variable posterior means"
            largest.lnorms <- order(rowSums(choose_lv_coefs^2),decreasing=TRUE)[1:ind.spp] ## Open question as to whether largest.lnorms should be based on which.lvs only
                
            plot(rbind(choose_lvs,choose_lv_coefs)[,which.lvs], xlab = paste("Latent variable", which.lvs[1]), ylab = paste("Latent variable", which.lvs[2]), main = main, type = "n", xlim = 1.1*range(rbind(choose_lvs,choose_lv_coefs)[,which.lvs[1]]), ylim = 1.1*range(rbind(choose_lvs,choose_lv_coefs)[,which.lvs[2]]), ...)
            if(!jitter) 
                text(choose_lvs[,which.lvs], label = 1:n, ...)
            if(jitter) 
                text(jitter(choose_lvs[,which.lvs[1]]), jitter(choose_lvs[,which.lvs[2]]), label = 1:n, ...)
            text(choose_lv_coefs[largest.lnorms,which.lvs], label = rownames(x$lv.coefs.mean[largest.lnorms,]), col = "red", ...)	
            }
        }	

    out <- list(scaled.lvs = choose_lvs, scaled.lv.coefs = choose_lv_coefs)
    if(return.vals) return(out)
    }


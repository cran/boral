plot.boral <- function(x, est = "median", include.ranef = TRUE, jitter = FALSE, ...) {
     #if(all(x$family %in% c("ordinal","multinom"))) 
          #stop("Residuals are not defined, and therefore residual analysis cannot be performed, if all columns of y are ordinal")
     if(any(x$family %in% c("ordinal","multinom"))) 
          warning("For all columns of y that are ordinal, the first plot constructed is of Dunn-Smyth residuals against fitted values (i.e., the level with the highest predicted probability). Note this can make things very confusing to interpret if only SOME of the columns in the response matrix are ordinal")
     
     get_etas <- fitted.boral(x, est = est, include.ranef = include.ranef, linear.predictor = TRUE) 
     get_ds_res <- ds.residuals(object = x, est = est, include.ranef = include.ranef)
     print(get_ds_res$agree.ordinal)
     get_ds_res <- get_ds_res$residuals
          

     palette(rainbow(ncol(get_etas)))

     matplot(get_etas, get_ds_res, ylab = "Dunn-Smyth Residuals", xlab = "Linear Predictors", type="n", ...)
     for(i in 1:ncol(get_etas)) 
          points(get_etas[,i], get_ds_res[,i], col=palette()[i], ...)
     abline(h=0, lty = 2, lwd = 2)
     
     matplot(get_ds_res, ylab = "Dunn-Smyth Residuals", xlab = "Row index",type="n", xaxt = "n", ...)
     axis(side = 1, at = 1:nrow(x$y), labels = rownames(x$lv.mean), ...)
     for (i in 1:ncol(get_etas)) { 
          points(seq(1,nrow(x$y)),get_ds_res[,i], col=palette()[i], ...) 
          }
     abline(0,0,lty=2)

     matplot(t(get_ds_res), ylab = "Dunn-Smyth Residuals", xlab = "Column index", type="n", xaxt = "n", ...)
     axis(side = 1, at = 1:ncol(x$y), labels = rownames(x$coefs.mean), ...)
     for(i in 1:ncol(get_etas)) { 
          points(rep(i,nrow(get_etas)), get_ds_res[,i], col=palette()[i], ...) 
          }
     abline(h=0, lty = 2, lwd = 2)

     get_ds_res2 <- as.vector(unlist(get_ds_res))
     qqnorm(get_ds_res2[is.finite(get_ds_res2)], main = "Quantile-Quantile Plot", ...)
     #qqline(y = get_ds_res2[is.finite(get_ds_res2)], ...)
     
     palette("default")
     }

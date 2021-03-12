print.summary.boral <- function(x, ...) {
     message("Call:\n")
     print(x$call)
     message()
    
     if(x$est == "median")
          message("Median point estimates\n\n LV coefficients (thetas) and dispersion parameter (if applicable)"); 
     if(x$est == "mean")
          message("Mean point estimates\n\n LV coefficients (thetas) and dispersion parameter (if applicable)")
     print(x$coefficients)
     message() 
 	
     if(!is.null(x$lv.covparams)) { 
          message("LV covariance structure parameters\n")
          print(x$lv.covparams) 
          message() 
          }
     if(!is.null(x$X.coefficients)) { 
          message("X coefficients (betas)\n")
          print(x$X.coefficients)
          message() 
          }
     if(!is.null(x$X.multinom.coefficients)) 
          message("There are also coefficients corresponding to multinomial columns which have not been printed.")
     if(!is.null(x$traits.coefficients)) { 
          message("Trait coefficients")
          print(x$traits.coefficients)
          message() 
          }
     
     if(!is.null(x$row.coefficients)) { 
          message("Row coefficients\n")
          print(x$row.coefficients) 
          message() 
          }
     if(!is.null(x$ranef.coefficients)) { 
          message("There are also response-specific random intercepts which have not been printed.")
          message() 
          }

     if(any(x$family == "ordinal")) { 
          message("Proportional odds (Cumulative probit) cutoffs") 
          print(x$cutoffs)
          message() 
          }
     if(any(x$family == "tweedie")) { 
          message("Tweedie power parameter")
          print(x$powerparam) 
          message() 
          }
     }	
		
		
summary.boral <- function(object, est = "median", ...) {
     if(est == "median") {
          gather_output <- list(call = object$call, coefficients = round(object$lv.coefs.median,3))
          if(object$num.lv > 0)
               gather_output$lvs <- round(object$lv.median,3)
          if(object$lv.control$type != "independent")
               gather_output$lv.covparams <- round(object$lv.covparams.median,3)
          if(object$row.eff != "none") {
               for(k0 in 1:ncol(object$row.ids)) 
                    gather_output$row.coefficients[[k0]] = round(object$row.coefs[[k0]]$median,3)
               }
          if(!is.null(object$ranef.ids)) {
               for(k0 in 1:ncol(object$ranef.ids)) 
                    gather_output$ranef.coefficients[[k0]] = round(object$ranef.coefs.median[[k0]],3)
               }
          if(object$num.X > 0) 
               gather_output$X.coefficients <- round(object$X.coefs.median,3)
          if(object$num.traits > 0) 
               gather_output$traits.coefficients <- round(object$traits.coefs.median,3)
          if(any(object$family == "ordinal")) 
               gather_output$cutoffs <- round(object$cutoffs.median,3)
          if(any(object$family == "tweedie")) 
               gather_output$powerparam <- round(object$powerparam.median,3)
          if(!is.null(object$X.multinom.coefs.median)) 
               gather_output$X.multinom.coefficients <- round(object$X.multinom.coefs.median,3) 
          }

    if(est == "mean") {
          gather_output <- list(call = object$call, coefficients = round(object$lv.coefs.mean,3))
          if(object$num.lv > 0)
               gather_output$lvs <- round(object$lv.mean,3)
          if(object$lv.control$type != "independent")
               gather_output$lv.covparams <- round(object$lv.covparams.mean,3)
          if(object$row.eff != "none") {
               for(k0 in 1:ncol(object$row.ids)) 
                    gather_output$row.coefficients[[k0]] = round(object$row.coefs[[k0]]$mean,3)
               }
          if(!is.null(object$ranef.ids)) {
               for(k0 in 1:ncol(object$ranef.ids)) 
                    gather_output$ranef.coefficients[[k0]] = round(object$ranef.coefs.mean[[k0]],3)
               }
          if(object$num.X > 0) 
               gather_output$X.coefficients <- round(object$X.coefs.mean,3)
          if(object$num.traits > 0) 
               gather_output$traits.coefficients <- round(object$traits.coefs.mean,3)
          if(any(object$family == "ordinal")) 
               gather_output$cutoffs <- round(object$cutoffs.mean,3)
          if(any(object$family == "tweedie")) 
               gather_output$powerparam <- round(object$powerparam.mean,3)
          if(!is.null(object$X.multinom.coefs.mean)) 
               gather_output$X.multinom.coefficients <- round(object$X.multinom.coefs.mean,3) 
          }


     gather_output$est <- est
     gather_output$trial.size <- object$trial.size
     gather_output$num.ord.levels <- object$num.ord.levels
     gather_output$prior.control$ssvs.index <- object$prior.control$ssvs.index 


     if(any(object$prior.control$ssvs.index == 0)) 
          gather_output$ssvs.indcoefs.prob <- round(object$ssvs.indcoefs.mean,3)
     if(any(object$prior.control$ssvs.index > 0)) 
          gather_output$ssvs.gpcoefs.prob <- round(object$ssvs.gpcoefs.mean,3) 
     if(any(unlist(object$prior.control$ssvs.traitsindex) == 0)) 
          gather_output$ssvs.traitscoefs.prob <- round(object$ssvs.traitscoefs.mean,3) 

     class(gather_output) <- "summary.boral"
     gather_output 
     }
 			

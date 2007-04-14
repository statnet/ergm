ergm.initialfit<-function(theta0, MLestimate, Clist, Clist.miss, m, 
                          verbose=FALSE, ...) {
# Process input for call to ergm.mple or some other alternative fitting
# method.  If the user wishes only to obtain the fit from this method
# (MLestimate==FALSE), this fit is returned immediately upon return to
# ergm function and the function terminates.  Otherwise (MLestimate==TRUE),
# we also check to see whether the theta0 value is a valid vector of
# initial parameters.
  fitmethod <- match("MPLE", theta0)
  if (is.na(fitmethod) & MLestimate) { # theta0 should be a start vector
    theta0 <- as.vector(theta0)
    if (length(theta0) != Clist$nparam |
        length(theta0)!=length(m$coef.names)) {
      cat("theta0 is", theta0, "\n", "Clist$nparam is",Clist$nparam, "\n")
      stop(paste("Invalid starting parameter vector theta0;",
                 "unrecognized option or wrong number of arguments."))
    }
    fit <- list(coef=theta0) # anything else that should be in the list? 
  } else if (is.na(fitmethod) & !MLestimate) { # Error!
    stop(paste("Unrecognized fitting method", theta0,
               "used in conjuction with MLestimate=FALSE.\n"))
  } else {
    if (fitmethod==1) {  #  MPLE
      fit <- ergm.mple(Clist, Clist.miss, m, verbose=verbose, ...)
    }
  }
  fit
}

ergm.initialfit<-function(theta0, MLestimate, 
                          formula, nw, meanstats,
			  m, 
                          MPLEtype="glm", initial.loglik=NULL,
                          conddeg=NULL, MCMCparams=NULL, MHproposal=NULL,
                          force.MPLE=FALSE,
                          verbose=FALSE, ...) {
# Process input for call to ergm.mple or some other alternative fitting
# method.  If the user wishes only to obtain the fit from this method
# (MLestimate==FALSE), this fit is returned immediately upon return to
# ergm function and the function terminates.  Otherwise (MLestimate==TRUE),
# we also check to see whether the theta0 value is a valid vector of
# initial parameters.
  fitmethod <- match("MPLE", theta0)
  if(is.na(fitmethod) && any(m$offset)) { # theta0 has an offset
   force.MPLE <- TRUE
  }
#
  if(!is.null(conddeg)){
   formula.conddegmple <- ergm.update.formula(formula, ~ conddegmple + .)
   m.conddeg <- ergm.getmodel(formula.conddegmple, nw, drop=conddeg, initialfit=TRUE)
   fitmethod <- 1
   Clist <- ergm.Cprepare(nw, m.conddeg)
   Clist.miss <- ergm.design(nw, m.conddeg, verbose=verbose)
   Clist$meanstats=c(1,meanstats)
   conddeg <- list(m=m.conddeg, Clist=Clist, Clist.miss=Clist.miss)
  }
  #
  Clist <- ergm.Cprepare(nw, m)
  Clist.miss <- ergm.design(nw, m, verbose=verbose)
  Clist$meanstats=meanstats
  MCMCparams$Clist.miss=Clist.miss
#
# Check to see if we have "constraint
  if (is.na(fitmethod) && (MLestimate|any(m$offset))) { # theta0 should be a start vector
    
    theta0 <- as.vector(theta0)
    if (length(theta0) != Clist$nstats |
        length(theta0)!=length(m$coef.names)) {
      cat("theta0 is", theta0, "\n", "Clist$nstats is",Clist$nstats, "\n")
      stop(paste("Invalid starting parameter vector theta0;",
                 "unrecognized option or wrong number of parameters."))
    }

    if(force.MPLE){
      fit <- ergm.mple(Clist, Clist.miss, m, MPLEtype=MPLEtype,
                       theta0=theta0, conddeg=conddeg, 
		       MCMCparams=MCMCparams, MHproposal=MHproposal,
                       verbose=verbose, ...)
    }else{    
      if(!is.null(Clist.miss)){
        mle.lik <- -log(2)*(Clist$ndyads-network.naedgecount(nw))
      }else{
        mle.lik <- -log(2)*Clist$ndyads
      }
      if(!is.null(initial.loglik)){
        mle.lik <- mle.lik-initial.loglik
      }
      fit <- structure(list(coef=theta0, mle.lik=mle.lik),class="ergm")
    }
  } else if (is.na(fitmethod) & !MLestimate) { # Error!
    stop(paste("Unrecognized fitting method", theta0,
               "used in conjuction with MLestimate=FALSE.\n"))
  } else {
    if (fitmethod==1) {  #  MPLE
      fit <- ergm.mple(Clist, Clist.miss, m, MPLEtype=MPLEtype,
                       conddeg=conddeg, 
		       MCMCparams=MCMCparams, MHproposal=MHproposal,
		       verbose=verbose, ...)
    }
  }
  fit
}

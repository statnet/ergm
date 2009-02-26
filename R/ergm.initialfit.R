#  File ergm/R/ergm.initialfit.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
ergm.initialfit<-function(theta0, MLestimate, Clist, Clist.miss, m, 
                          MPLEtype="glm", initial.loglik=NULL,
                          force.MPLE=FALSE,
                          verbose=FALSE, ...) {
# Process input for call to ergm.mple or some other alternative fitting
# method.  If the user wishes only to obtain the fit from this method
# (MLestimate==FALSE), this fit is returned immediately upon return to
# ergm function and the function terminates.  Otherwise (MLestimate==TRUE),
# we also check to see whether the theta0 value is a valid vector of
# initial parameters.
  fitmethod <- match("MPLE", theta0)
  if (is.na(fitmethod) && MLestimate) { # theta0 should be a start vector
    
    theta0 <- as.vector(theta0)
    if (length(theta0) != Clist$nstats |
        length(theta0)!=length(m$coef.names)) {
      cat("theta0 is", theta0, "\n", "Clist$nstats is",Clist$nstats, "\n")
      stop(paste("Invalid starting parameter vector theta0;",
                 "unrecognized option or wrong number of parameters."))
    }

    if(force.MPLE){
      fit <- ergm.mple(Clist, Clist.miss, m, MPLEtype=MPLEtype,
                       theta0=theta0,
                       verbose=verbose, ...)
    }else{    
      if(!is.null(Clist.miss)){
        mle.lik <- -log(2)*(Clist$ndyads-Clist.miss$nedges)
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
                       verbose=verbose, ...)
    }
  }
  fit
}

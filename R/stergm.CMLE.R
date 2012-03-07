#  File ergm/R/stergm.CMLE.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
################################################################################
# The <stergm> function fits stergms from a specified formation and dissolution
# formula returning approximate MLE's based on MCMC estimation.
################################################################################

stergm.CMLE <- function(nw, formation, dissolution, times, offset.coef.form, offset.coef.diss,
                        eval.loglik,
                        estimate,
                        control,
                        verbose) {

  # Translate the "estimate" from the stergm() argument to the ergm() argument.
  estimate <- switch(estimate,
                     CMLE = "MLE",
                     CMPLE = "MPLE")
  if(estimate=="MPLE") warning("MPLE for stergm() is extremely imprecise at this time!")

  if(is.null(times)){
    if(inherits(nw, "network.list") || is.list(nw)){
      times  <- c(1,2)
      warning("Time points not specified for a list. Modeling transition from the first to the second network. This behavior may change in the future.")
    }else if(inherits(nw,"networkDynamic")){
      times  <- c(0,1)
      warning("Time points not specified for a networkDynamic. Modeling transition from time 0 to 1.")
    }
  }
  
  if(length(times)<2) stop("Time points whose transition is to be modeled was not specified.")
  if(length(times)>2) stop("Only two time points (one transition) are supported at this time.")

  if(inherits(nw, "network.list") || is.list(nw)){
    y0 <- nw[[times[1]]]
    y1 <- nw[[times[2]]]
  }else if(inherits(nw,"networkDynamic")){
    require(networkDynamic) # This is needed for the "%t%.network" function
    y0 <- nw %t% times[1]
    y1 <- nw %t% times[2]
  }

  if(!is.network(y0) || !is.network(y1)) stop("nw must be a networkDynamic, a network.list, or a list of networks.")
  
  
  # Construct the formation and dissolution networks; the
  # network.update cannot be used to copy attributes from y0 to y.form and
  # y.diss, since NA edges will be lost.
  
  y.form <- y0 | y1
  y.form <- nvattr.copy.network(y.form, y0)
  formation <- ergm.update.formula(formation, y.form~.)

  y.diss <- y0 & y1
  y.diss <- nvattr.copy.network(y.diss, y0)
  dissolution <- ergm.update.formula(dissolution, y.diss~.)

  # Apply initial values passed to control.stergm() the separate controls, if necessary.
  if(is.null(control$CMLE.control.form$init)) control$CMLE.control.form$init <- control$init.form
  if(is.null(control$CMLE.control.diss$init)) control$CMLE.control.diss$init <- control$init.diss
  
  # Now, call the ergm()s:
  cat("Fitting formation:\n")
  fit.form <- ergm(formation, constraints=~atleast(y0), offset.coef=offset.coef.form, eval.loglik=eval.loglik, estimate=estimate, control=control$CMLE.control.form, verbose=verbose)
  cat("Fitting dissolution:\n")
  fit.diss <- ergm(dissolution, constraints=~atmost(y0), offset.coef=offset.coef.diss, eval.loglik=eval.loglik, estimate=estimate, control=control$CMLE.control.diss, verbose=verbose)

  # Construct the output list. Conveniently, this is mainly a list consisting of two ergms.
  
  list(network=nw, times=times, formation=formation, dissolution=dissolution, formation.fit=fit.form, dissolution.fit=fit.diss, estimate=estimate)
}

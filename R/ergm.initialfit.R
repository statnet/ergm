#  File ergm/R/ergm.initialfit.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################

ergm.initialfit<-function(init, initial.is.final,
                          formula, nw, target.stats,
                          m, method = NULL,
                          MPLEtype="glm",
                          conddeg=NULL, control=NULL, MHproposal=NULL,
                          verbose=FALSE, ...) {
  method <- match.arg(method, c("MPLE","zeros"))
 
  # conddeg, whatever it does.
  if(!is.null(conddeg)){
   formula.conddegmple <- ergm.update.formula(formula, ~ conddegmple + .)
   m.conddeg <- ergm.getmodel(formula.conddegmple, nw, initialfit=TRUE)
   method <- "MPLE"
   Clist <- ergm.Cprepare(nw, m.conddeg)
   Clist.miss <- ergm.design(nw, m.conddeg, verbose=FALSE)
   m$target.stats=c(1,target.stats)
   conddeg <- list(m=m.conddeg, Clist=Clist, Clist.miss=Clist.miss)
  }
  
  Clist <- ergm.Cprepare(nw, m)
  Clist.miss <- ergm.design(nw, m, verbose=FALSE)
  m$target.stats<-target.stats
  control$Clist.miss<-Clist.miss

  # Respect init elements that are not offsets if it's only a starting value.
  if(!initial.is.final){ 
    m$etamap$offsettheta[!is.na(init)] <- TRUE
  }

  if(initial.is.final || any(is.na(init))){
    # If we aren't running MCMC after this or not all of init has been
    # supplied by the user, use MPLE.   
    # Also make sure that any initial values specified by the user are respected.
    fit <- switch(method,
                  MPLE = ergm.mple(Clist, Clist.miss, m, MPLEtype=MPLEtype,
                    init=init, conddeg=conddeg, 
                    control=control, MHproposal=MHproposal,
                    verbose=verbose, ...),
                  zeros = structure(list(coef=ifelse(is.na(init),0,init)),class="ergm"),
                  stop(paste("Invalid method specified for initial parameter calculation. Available methods are ",paste.and(formals()$method),".",sep=""))
                  )
  }else{
    # If this is just the initial value, *and* the user has supplied
    # all elements for init, just echo init.
    fit <- structure(list(coef=init),class="ergm")
  }
  fit
}

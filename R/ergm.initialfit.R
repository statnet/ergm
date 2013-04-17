#  File R/ergm.initialfit.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
####################################################################################
# The <ergm.initialfit> function fits an initial ergm object using either ML or MPL
# estimation.  If initial parameters are provided in 'init' and 'MLestimate' is 
# TRUE, the number of parameters in 'init' is checked for correctness.
# 
# --PARAMETERS--
#   init        :  either a vector whose first entry is "MPLE" or a vector
#                    of initial coefficients
#   MLestimate    :  whether a MLestimate should be used (T or F); 
#                       if TRUE, this may be overriden by 'force.MPLE'=TRUE  or
#                                init=c("MPLE", ...)
#                       if FALSE, 'init' must have "MPLE" as its 1st entry to 
#                                  avoid an error
#   formula       :  a formula of the form (nw ~ term(s)) 
#   nw            :  a network object, presumably that of 'formula'
#   target.stats     :  the mean statistics
#   m             :  the model as returned by <ergm.getmodel>
#   MPLEtype      :  the method for MPL estimation as either "glm", "penalized",
#                    or "logitreg"; this is ignored if ML estimation is used;
#                    default="glm" 
#   initial.loglik:  the initial log likelihood; default=NULL
#   conddeg       :  a formula for the conditional degree terms
#   control    :  a list of parameters for tuning the MCMC sampling;
#                    the only recognized component is 'samplesize'
#   MHproposal    :  an MHproposal object as returned by <getMHproposal>
#   force.MPLE    :  whether MPL estimation should be forced instead of ML 
#                    estimation (T or F); this is ignored if 'MLestimate'=FALSE
#                    or "MPLE" is an entry into 'init'; default=FALSE
#   verbose       :  whether the MPL estimation should be verbose (T or F); 
#                    default=FALSE
#   ...           :  addtional parameters that are used with MPL estimation;
#                    the only recognized parameeter is 'compressflag' which
#                    compresses the design matrix used by <ergm.mple>        
#
# --RETURNED--
#    an ergm object as one of the following lists
#     if MLE  -- a list with 2 components
#                  coef   : 'init'
#                  mle.lik:  the MLE likelihood
#    if MPLE -- the list returned by <ergm.mple>
#
######################################################################################

ergm.initialfit<-function(init, initial.is.final,
                          formula, nw, target.stats,
                          m, reference=~Bernoulli, method = NULL,
                          MPLEtype="glm",
                          conddeg=NULL, control=NULL, MHproposal=NULL, MHproposal.obs=NULL,
                          verbose=FALSE, ...) {
  method <- match.arg(method, ergm.init.methods(MHproposal$reference$name))
 
  # conddeg, whatever it does.
  if(method=="MPLE" && !is.null(conddeg)){
   formula.conddegmple <- ergm.update.formula(formula, . ~ conddegmple + .)
   m.conddeg <- ergm.getmodel(formula.conddegmple, nw, initialfit=TRUE)
   Clist <- ergm.Cprepare(nw, m.conddeg)
   Clist.miss <- ergm.design(nw, m.conddeg, verbose=FALSE)
   m$target.stats=c(1,target.stats)
   conddeg <- list(m=m.conddeg, Clist=Clist, Clist.miss=Clist.miss)
  }

  Clist <- ergm.Cprepare(nw, m)
  Clist.miss <- ergm.Cprepare(NVL(get.miss.dyads(MHproposal$arguments$constraints, MHproposal.obs$arguments$constraints), is.na(nw)), m)
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

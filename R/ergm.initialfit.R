####################################################################################
# The <ergm.initialfit> function fits an initial ergm object using either ML or MPL
# estimation.  If initial parameters are provided in 'theta0' and 'MLestimate' is 
# TRUE, the number of parameters in 'theta0' is checked for correctness.
# 
# --PARAMETERS--
#   theta0        :  either a vector whose first entry is "MPLE" or a vector
#                    of initial coefficients
#   MLestimate    :  whether a MLestimate should be used (T or F); 
#                       if TRUE, this may be overriden by 'force.MPLE'=TRUE  or
#                                theta0=c("MPLE", ...)
#                       if FALSE, 'theta0' must have "MPLE" as its 1st entry to 
#                                  avoid an error
#   formula       :  a formula of the form (nw ~ term(s)) 
#   nw            :  a network object, presumably that of 'formula'
#   meanstats     :  the mean statistics
#   m             :  the model as returned by <ergm.getmodel>
#   MPLEtype      :  the method for MPL estimation as either "glm", "penalized",
#                    or "logitreg"; this is ignored if ML estimation is used;
#                    default="glm" 
#   initial.loglik:  the initial log likelihood; default=NULL
#   conddeg       :  a formula for the conditional degree terms
#   MCMCparams    :  a list of parameters for tuning the MCMC sampling;
#                    the only recognized component is 'samplesize'
#   MHproposal    :  an MHproposal object as returned by <getMHproposal>
#   force.MPLE    :  whether MPL estimation should be forced instead of ML 
#                    estimation (T or F); this is ignored if 'MLestimate'=FALSE
#                    or "MPLE" is an entry into 'theta0'; default=FALSE
#   verbose       :  whether the MPL estimation should be verbose (T or F); 
#                    default=FALSE
#   ...           :  addtional parameters that are used with MPL estimation;
#                    the only recognized parameeter is 'compressflag' which
#                    compresses the design matrix used by <ergm.mple>        
#
# --RETURNED--
#    an ergm object as one of the following lists
#     if MLE  -- a list with 2 components
#                  coef   : 'theta0'
#                  mle.lik:  the MLE likelihood
#    if MPLE -- the list returned by <ergm.mple>
#
######################################################################################

ergm.initialfit<-function(theta0, initial.is.final,
                          formula, nw, meanstats,
                          m, method = "MPLE",
                          MPLEtype="glm",
                          conddeg=NULL, MCMCparams=NULL, MHproposal=NULL,
                          verbose=FALSE, ...) {
  method <- match.arg(method)
 
  # conddeg, whatever it does.
  if(!is.null(conddeg)){
   formula.conddegmple <- ergm.update.formula(formula, ~ conddegmple + .)
   m.conddeg <- ergm.getmodel(formula.conddegmple, nw, initialfit=TRUE)
   method <- "MPLE"
   Clist <- ergm.Cprepare(nw, m.conddeg)
   Clist.miss <- ergm.design(nw, m.conddeg, verbose=FALSE)
   Clist$meanstats=c(1,meanstats)
   conddeg <- list(m=m.conddeg, Clist=Clist, Clist.miss=Clist.miss)
  }
  
  Clist <- ergm.Cprepare(nw, m)
  Clist.miss <- ergm.design(nw, m, verbose=FALSE)
  Clist$meanstats=meanstats
  MCMCparams$Clist.miss=Clist.miss

  # Respect theta0 elements that are not offsets if it's only a starting value.
  if(!initial.is.final){ 
    m$etamap$offsettheta[!is.na(theta0)] <- TRUE
  }

  if(initial.is.final || any(is.na(theta0))){
    # If we aren't running MCMC after this or not all of theta0 has been
    # supplied by the user, use MPLE.   
    # Also make sure that any initial values specified by the user are respected.
    fit <- switch(method,
                  MPLE = ergm.mple(Clist, Clist.miss, m, MPLEtype=MPLEtype,
                    theta0=theta0, conddeg=conddeg, 
                    MCMCparams=MCMCparams, MHproposal=MHproposal,
                    verbose=verbose, ...),
                  error(paste("Invalid method specified for initial parameter calculation. Available methods are ",paste(formals()$method,collapse=","),".",sep=""))
                  )
  }else{
    # If this is just the initial value, *and* the user has supplied
    # all elements for theta0, just echo theta0.
    fit <- structure(list(coef=theta0),class="ergm")
  }
  fit
}

#  File R/ergm.initialfit.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
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
#   nw            :  a network object, presumably that of 'formula'
#   target.stats     :  the mean statistics
#   m             :  the model as returned by <ergm_model>
#   MPLEtype      :  the method for MPL estimation as either "glm", "penalized",
#                    or "logitreg"; this is ignored if ML estimation is used;
#                    default="glm" 
#   initial.loglik:  the initial log likelihood; default=NULL
#   control    :  a list of parameters for tuning the MCMC sampling;
#                    the only recognized component is 'samplesize'
#   proposal    :  an proposal object as returned by <getproposal>
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
ergm.initialfit<-function(init,
                          s, s.obs,
                          control=NULL,
                          verbose=FALSE, ...) {
  if(control$init.method!="skip" && any(is.NA(init))){
    # Respect init elements that are not offsets if it's only a starting value.
    s$model$etamap$offsettheta[!is.NA(init)] <- TRUE

    # Also make sure that any initial values specified by the user are respected.
    switch(control$init.method,
           MPLE = {
             if(control$MPLE.constraints.ignore) s$proposal$arguments$constraints <- s$proposal$arguments$constraints[".attributes"] # Drop all constraints except for .attributes .
             control$MPLE.samplesize <- control$init.MPLE.samplesize
             ergm.mple(s, s.obs,
                       init=init,
                       control=control,
                       verbose=verbose, ...)
           },
           zeros = structure(list(coefficients=.constrain_init(s$model, ifelse(is.NA(init),0,init)))),
           CD = ergm.CD.fixed(.constrain_init(s$model, ifelse(is.NA(init),0,init)),
                              s, s.obs, control, verbose, ...),
           stop("Invalid method specified for initial parameter calculation.")
    )
  }else{
    # If this is just the initial value, *and* the user has supplied
    # all elements for init, just echo init.
    structure(list(coefficients=init))
  }
}

.constrain_init <- function(m, init){
  if(is(m, "ergm_model")) m <- m$etamap
  init.no <- init[!m$offsettheta]
  maxtheta <- m$maxtheta[!m$offsettheta]
  init.no <- pmin(init.no, maxtheta - deInf(pmax(abs(maxtheta),1)*sqrt(.Machine$double.eps)))
  mintheta <- m$mintheta[!m$offsettheta]
  init.no <- pmax(init.no, mintheta + deInf(pmax(abs(mintheta),1)*sqrt(.Machine$double.eps)))
  init[!m$offsettheta] <- init.no
  init
}

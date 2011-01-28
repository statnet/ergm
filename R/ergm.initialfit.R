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
#   Clist         :  the list of parameters needed for ML or MPL estimation and 
#                    returned by <ergm.Cprepare>
#   Clist.miss    :  the list of parameters for the network of missing edges
#                    needed for ML or MPL estimation and returned by <ergm.design> 
#   m             :  the model as returned by <ergm.getmodel>
#   MPLEtype      :  the method for MPL estimation as either "glm", "penalized",
#                    or "logitreg"; this is ignored if ML estimation is used;
#                    default="glm" 
#   initial.loglik:  the initial log likelihood; default=NULL
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

ergm.initialfit<-function(theta0, MLestimate, Clist, Clist.miss, m, 
                          MPLEtype="glm", initial.loglik=NULL,
                          conddeg=NULL, MCMCparams=NULL, MHproposal=NULL,
                          force.MPLE=FALSE,
                          verbose=FALSE, ...) {
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
                       theta0=theta0, conddeg=conddeg, 
		       MCMCparams=MCMCparams, MHproposal=MHproposal,
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
                       conddeg=conddeg, 
		       MCMCparams=MCMCparams, MHproposal=MHproposal,
		       verbose=verbose, ...)
    }
  }
  fit
}

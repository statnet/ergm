## A function to compute and return the log-likelihood of an ERGM MLE.
logLik.ergm<-function(object, add=FALSE, force.reeval=FALSE, eval.loglik=add || force.reeval, control=control.logLik.ergm(), ...){
  check.control.class()
  
  control.transfer <- c("MCMC.burnin", "MCMC.interval", "MCMC.prop.weights",
"MCMC.prop.args", "MCMC.packagenames", "MCMC.init.maxedges", "MCMC.samplesize",
"obs.MCMC.burnin", "obs.MCMC.interval", "obs.MCMC.samplesize","warn.dyads")
  for(arg in control.transfer)
    if(is.null(control[[arg]]))
      control[arg] <- list(object$control[[arg]])

  control <- set.control.class("control.ergm.bridge")

  # "object" has an element control.
  loglik.control<-control

  
  out<-with(object,
            if(!force.reeval && !is.null(object$mle.lik)) mle.lik
            else{
              if(!eval.loglik) stop(nologLik.message(deparse(substitute(object))))
              ## If valued, compute a path sample from reference measure.
              if(!is.null(object$response)) ergm.bridge.0.llk(formula,response=response,reference=reference,constraints=constraints,coef=coef(object),control=loglik.control,llkonly=FALSE,...)
              else{
                ## If dyad-independent, just go from the deviance.
                if(is.dyad.independent(object) && is.null(object$sample) && is.dyad.ind.constraints(object$constrained, object$constrained.obs)) -object$glm$deviance/2 - -object$glm$null.deviance/2
                else
                  ## If dyad-dependent, bridge from a dyad-independent model.
                  ergm.bridge.dindstart.llk(formula,reference=reference,constraints=constraints,coef=coef(object),control=loglik.control,llkonly=FALSE,...)
              }
            }
            )
  
  if(is.numeric(out)){
    llk<-out
  }else{
    llk<-out$llk
    attr(llk,"br")<-out
  }
  class(llk)<-"logLik"
  attr(llk,"df")<-length(coef(object))
  attr(llk,"nobs")<-{
    # FIXME: We need a more general framework for handling constrained
    # and partially observed network "degrees of freedom".

    if(!is.dyad.ind.constraints(object$constrained, object$constrained.obs)
       && loglik.control$warn.dyads){
      warning("The number of observed dyads in this network is ill-defined due to complex constraints on the sample space.")
    }
    
    network.dyadcount(object$network) - network.edgecount(NVL(get.miss.dyads(object$constrained, object$constrained.obs),is.na(object$network)))
  }
  
  if(add){
    object$mle.lik<-llk
    object    
  } else llk
}

nologLik.message<-function(objname){
  paste("Log-likelihood was not estimated for this fit.\nTo get deviances, AIC, and/or BIC from fit `",objname,"` run \n  > ",objname,"<-logLik(",objname,", add=TRUE)\nto add it to the object or rerun this function with eval.loglik=TRUE.\n",sep="")
}

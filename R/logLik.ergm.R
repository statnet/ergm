#  File ergm/R/logLik.ergm.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
## A function to compute and return the log-likelihood of an ERGM MLE.
logLik.ergm<-function(object, add=FALSE, force.reeval=FALSE, eval.loglik=add || force.reeval, control=control.logLik.ergm(), ...){

  control.transfer <- c("MCMC.burnin", "MCMC.interval", "MCMC.prop.weights", "MCMC.prop.args", "MCMC.packagenames", "MCMC.init.maxedges","MCMC.samplesize")
  for(arg in control.transfer)
    if(is.null(control[[arg]]))
      control[arg] <- list(object$control[[arg]])


  # "object" has an element control.
  loglik.control<-control

  
  out<-with(object,
            if(!force.reeval && !is.null(object$mle.lik)) mle.lik
            else{
              if(!eval.loglik) stop(nologLik.message(deparse(substitute(object))))
              ## If valued, compute a path sample from reference measure.
                ## If dyad-independent, just go from the deviance.
                if(is.dyad.independent(object) && is.null(object$sample)) -object$glm$deviance/2
                else
                  ## If dyad-dependent, bridge from a dyad-independent model.
                  ergm.bridge.dindstart.llk(formula,constraints=constraints,coef=coef(object),control=loglik.control,llkonly=FALSE,...)
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
    #
    # Note that bd is always going to be in there.
    constr <- names(object$constrained)
    if(length(constr)==1){
      network.dyadcount(object$network)-network.naedgecount(object$network)
    }else if(length(constr)>=2){
      warning("The number of observed dyads in this network is ill-defined due to complex constraints on the sample space.")    
      network.dyadcount(object$network)-network.naedgecount(object$network)
    }else{
      switch(constr,
             # The following assumes that the atleast and/or the atmost network cannot have missing dyads.
             atleast = network.dyadcount(object$network)-network.naedgecount(object$network | object$constraint$atleast$nw) - network.edgecount(object$constraint$atleast$nw),
             atmost = network.edgecount(object$constraint$atmost$nw) - network.naedgecount(object$network & object$constraint$atleast$nw),
             {
               warning("The number of observed dyads in this network is ill-defined due to constraints on the sample space.")
               network.dyadcount(object$network)-network.naedgecount(object$network)
             })
    }
  }    
  
  if(add){
    object$mle.lik<-llk
    object    
  } else llk
}

nologLik.message<-function(objname){
  paste("Log-likelihood was not estimated for this fit.\nTo get deviances, AIC, and/or BIC from fit `",objname,"` run \n  > ",objname,"<-logLik(",objname,", add=TRUE)\nto add it to the object or rerun this function with eval.loglik=TRUE.\n",sep="")
}

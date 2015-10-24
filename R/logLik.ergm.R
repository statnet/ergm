## A function to compute and return the log-likelihood of an ERGM MLE.
logLik.ergm<-function(object, add=FALSE, force.reeval=FALSE, eval.loglik=add || force.reeval, control=control.logLik.ergm(), ...){

  if(!force.reeval && !is.null(object$mle.lik)) return(object$mle.lik)

  # Then, we need to recalculate...
  
  check.control.class()
  
  control.transfer <- c("MCMC.burnin", "MCMC.interval", "MCMC.prop.weights",
"MCMC.prop.args", "MCMC.packagenames", "MCMC.init.maxedges", "MCMC.samplesize",
"obs.MCMC.burnin", "obs.MCMC.interval", "obs.MCMC.samplesize","warn.dyads","MPLE.type","MPLE.max.dyad.types","parallel","parallel.type","parallel.version.check"
)
  for(arg in control.transfer)
    if(is.null(control[[arg]]))
      control[arg] <- list(object$control[[arg]])

  control <- set.control.class("control.ergm.bridge")

  # "object" has an element control.
  loglik.control<-control
  
  out<-with(object,
            {
              if(!eval.loglik) stop(nologLik.message(deparse(substitute(object))))
              
              ## If dyad-independent or MPLE, just go from the deviance.
              if(object$estimate=="MPLE"
                 || (is.dyad.independent(object)
                     && is.null(object$sample)
                     && is.null(object$response)))
			 if(control$MPLE.type=="penalized")
				 object$glm$loglik - object$glm.null$loglik else
                -object$glm$deviance/2 - -object$glm$null.deviance/2
              else
                ## If dyad-dependent but not valued and has a dyad-independent constraint, bridge from a dyad-independent model.
                if(is.dyad.independent(object$constrained, object$constrained.obs)
                   && is.null(object$response))
                  ergm.bridge.dindstart.llk(formula,reference=reference,constraints=constraints,coef=coef(object),control=loglik.control,llkonly=FALSE,...)
                else
                  ## If valued or has dyad-dependent constraint, compute a path sample from reference measure.
                  ergm.bridge.0.llk(formula,response=object$response,reference=reference,constraints=constraints,coef=coef(object),control=loglik.control,llkonly=FALSE,...)
            }
            )
  
  if(is.numeric(out)){
    llk<-out
  }else{
    llk<-out$llk
    attr(llk,"br")<-out
  }
  if(!inherits(llk,"logLik")){
    class(llk)<-"logLik"
    attr(llk,"df")<-length(coef(object))
    attr(llk,"nobs")<- if(is.null(object$null.lik)){ # If we can steal the number of observations from the null model...
      # FIXME: We need a more general framework for handling
      # constrained and partially observed network "degrees of
      # freedom". PROGRESS: We can handle dyad-independent ones fine,
      # now.
      
      if(!is.dyad.independent(object$constrained, object$constrained.obs)
         && loglik.control$warn.dyads){
        warning("The number of observed dyads in this network is ill-defined due to complex constraints on the sample space.")
      }
      
      network.dyadcount(object$network,FALSE) - network.edgecount(NVL(get.miss.dyads(object$constrained, object$constrained.obs),network.initialize(1)))
    }else attr(object$null.lik,"nobs")
  }


  if(!is.null(object$null.lik) && !is.na(object$null.lik)){ # If Null likelihood is defined, shift the MLE likelihood.
    llk[] <- c(llk + object$null.lik) # The brackets are important to preserve attr()s on llk, and the c() is important to strip the ones from the sum.
  }
  
  if(add){
    object$mle.lik<-llk
    object    
  } else llk
}

nologLik.message<-function(objname){
  paste("Log-likelihood was not estimated for this fit.\nTo get deviances, AIC, and/or BIC from fit `",objname,"` run \n  > ",objname,"<-logLik(",objname,", add=TRUE)\nto add it to the object or rerun this function with eval.loglik=TRUE.\n",sep="")
}

logLikNull <- function(object, ...) UseMethod("logLikNull")

logLikNull.ergm <- function(object, control=control.logLik.ergm(), ...){
  if(!is.null(object$null.lik)) object$null.lik

  nobs <- if(is.null(object$mle.lik)) network.dyadcount(object$network,FALSE) - network.edgecount(NVL(get.miss.dyads(object$constrained, object$constrained.obs),network.initialize(1))) else attr(object$mle.lik,"nobs")
  
  llk <-
    if(!is.null(object$response)){
      cat("Note: Null model likelihood calculation is not implemented for valued ERGMs at this time.\n")
      NA
    }else if(!is.dyad.independent(object$constrained, object$constrained.obs)){
      cat("Note: The constraint on the sample space is not dyad-independent. Null model likelihood is only implemented for dyad-independent constraints at this time. Number of observations is similarly ill-defined.\n")
      NA
    }else nobs * log(1/2)
  

  class(llk)<-"logLik"
  attr(llk,"df")<-0
  attr(llk,"nobs")<-nobs

  llk
}

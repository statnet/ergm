## A function to compute and return the log-likelihood of an ERGM MLE.
logLik.ergm<-function(object,nsteps=10,add=FALSE,eval.loglik=add,...){
  out<-with(object,
            if(!is.null(object$mle.lik)) mle.lik
            else{
              if(!eval.loglik) stop(nologLik.message(deparse(substitute(object))))
              ## If valued, compute a path sample from reference measure.
              if(!is.null(object$response)) ergm.bridge.0.llk(formula,response=response,reference=reference,constraints=constraints,theta=coef(object),nsteps=nsteps,llkonly=FALSE,...)
              else{
                ## If dyad-independent, just go from the deviance.
                if(is.dyad.independent(object)) -object$glm$deviance/2
                else
                  ## If dyad-dependent, bridge from a dyad-independent model.
                  ergm.bridge.dindstart.llk(formula,reference=reference,constraints=constraints,theta=coef(object),nsteps=nsteps,llkonly=FALSE,...)
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
  attr(llk,"nobs")<-network.dyadcount(object$network)-network.naedgecount(object$network)
  if(add){
    object$mle.lik<-llk
    object    
  } else llk
}

nologLik.message<-function(objname){
  paste("Log-likelihood was not estimated for this fit.\nTo get deviances, AIC, and/or BIC, run \n  > ",objname,"<-logLik.ergm(",objname,", add=TRUE)\nto add it to the object or rerun this function with eval.loglik=TRUE.\n",sep="")
}

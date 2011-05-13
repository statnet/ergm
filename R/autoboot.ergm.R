## An experimental routine to compute bootstrap-like distribution of
## parameter estimates for a fitted ergm.
##
## It takes a dyad-dependent ERGM object, `object`, a bootstrap
## resample size `S` that defaults to MCMCsamplesize of the ergm
## object, and an optional control.ergm list.  Note that the
## optimization step takes a fair amount of time (linear in `S` and
## otherwise depending on MCMCsamplesize and number of canonical
## statistics), so `S` should probably be smaller than the default for
## most uses.
##
## The return value is simply a data frame with coefficient names in
## the columns and replicates in the rows.
##
## The basic idea is that for each draw of an ERGM's sufficient
## statistics under the MLE (contained in ergm$sample), take 1 MCMC
## MLE optimization step from the actual MLE to evaluate what the MLE
## would have been had that value been the observed statistic.
## Provided these new MLEs are not too far from actual MLE, this
## should be equivalent to fitting an ERGM to each of a sample of
## networks drawn from the MLE distribution, a parametric bootstrap
## procedure.
##
## It should work perfectly well for curved models, but does not, at
## present, support partially observed ERGMs and may or may not work
## for models that have been drop-ed.

autoboot.ergm<-function(object, S=nrow(object$sample), control=control.ergm()){
  
  if(is.dyad.independent(object))
    stop(paste("ERGM fit `",deparse(substitute(object)),"` is dyad-independent. Standard errors returned by the object are exact.",sep=""))

  m<-ergm.getmodel(object$formula, ergm.getnetwork(object$formula), drop=FALSE, response=object$response)

  samp<-object$sample
  if(S>nrow(samp))stop("Bootstrap sample size cannot be greater than the ERGM sample size.")

  theta0<-coef(object)
  theta.boot<-apply(samp[sample(seq_len(nrow(samp)),S),],1,function(stat){
    v<-ergm.estimate(theta0=theta0,model=m,statsmatrix=sweep(samp,2,stat),statsmatrix.obs=NULL,
                     epsilon=control$epsilon,
                     nr.maxit=control$nr.maxit,
                     nr.reltol=control$nr.reltol,
                     calc.mcmc.se=FALSE, hessianflag=FALSE,
                     trustregion=+Inf, method=control$method,
                     metric=control$metric,
                     compress=control$compress, verbose=FALSE,
                     estimateonly=TRUE)
    v$coef
  })
  as.data.frame(t(theta.boot))
}

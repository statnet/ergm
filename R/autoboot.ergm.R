## An experimental routine to compute bootstrap-like distribution of
## parameter estimates for a fitted ergm.
##
## It takes an ERGM object estimated using MCMC MLE, `object`, a
## bootstrap resample size `R` that defaults to MCMCsamplesize of the
## ergm object, and an optional control.ergm list.  Note that the
## optimization step takes a fair amount of time (linear in `R` and
## otherwise depending on MCMCsamplesize and number of canonical
## statistics), so `R` should probably be smaller than the default for
## most uses.
##
## The return value is simply a data frame with coefficient names in
## the columns and replicates in the rows.
##
## The basic idea is that for each draw of an ERGM's sufficient
## statistics under the MLE (contained in the reweighted ergm$sample),
## take 1 MCMC MLE optimization step from the actual MLE to evaluate
## what the MLE would have been had that value been the observed
## statistic.  Provided these new MLEs are not too far from actual
## MLE, this should be equivalent to fitting an ERGM to each of a
## sample of networks drawn from the MLE distribution, a parametric
## bootstrap procedure.
##
## It should work perfectly well for curved models, but does not, at
## present, support partially observed ERGMs and may or may not work
## for models that have been drop-ed.

autoboot.ergm<-function(object, R, control=control.ergm()){

  if(is.null(object$sample)){
    if(is.dyad.independent(object))
      stop(paste("ERGM fit `",deparse(substitute(object)),"` is dyad-independent. Standard errors returned by the object are exact.",sep=""))
    else stop(paste("ERGM fit `",deparse(substitute(object)),"` is missing the sample.",sep=""))
  }

  m<-ergm.getmodel(object$formula, ergm.getnetwork(object$formula), drop=FALSE, response=object$response)

  samp<-object$sample

  ## "Compress" duplicate statistics:
  library(coda)
  samp<-compress.data.frame(as.data.frame(samp))
  samp.fr<-samp$frequencies
  samp<-as.matrix(samp$rows)
  S<-nrow(samp)
  
  ## Because the MCMC sample that came with the ergm is drawn not from
  ## the MLE, but from the iteration prior to the MLE, we need to
  ## reweight the realizations for the bootstrap sample:
  
  theta.samp<-object$MCMCtheta
  eta.samp<-ergm.eta(theta.samp, m$etamap)
  theta.mle<-coef(object)
  eta.mle<-ergm.eta(theta.mle, m$etamap)

  eta.diff<-eta.mle-eta.samp

  samp.w<-exp(samp %*% cbind(eta.diff))*samp.fr
  samp.w<-samp.w/sum(samp.w)

  ## Now, the expensive part of this procedure is the many
  ## optimization steps, so resample and compress:

  resamp<-sample.int(S,R,replace=TRUE,prob=samp.w)
  resamp.w<-tabulate(resamp, S)
  resamp.l<-which(resamp.w!=0)
    
  theta.boot<-apply(samp[resamp.l,],1,function(stat){
    v<-ergm.estimate(theta0=theta.samp,model=m,statsmatrix=sweep(samp,2,stat),statsmatrix.obs=NULL,
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

  resamp.ind<-rep(seq_along(resamp.l),resamp.w[resamp.w!=0])
  t(theta.boot)[resamp.ind,]
}

## Compress a data frame by eliminating duplicate rows while keeping
## track of their frequency.
compress.data.frame<-function(x){
  x<-sort(x)
  firsts<-which(!duplicated(x))
  freqs<-diff(c(firsts,nrow(x)+1))
  x<-x[firsts,]
  list(rows=x,frequencies=freqs)
}

## Sorts rows of a data frame in lexicographic order.
sort.data.frame<-function(x){
  x[do.call(order,sapply(seq_along(x),function(i)x[[i]],simplify=FALSE)),]
}

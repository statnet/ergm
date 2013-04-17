#  File R/autoboot.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
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

autoboot.ergm<-function(object, R, verbose=FALSE, control=object$control){
  check.control.class("ergm")
  if(is.null(object$sample) || any(is.na(object$sample))){
    if(is.dyad.independent(object))
      stop(paste("ERGM fit `",deparse(substitute(object)),"` is dyad-independent. Standard errors returned by the object are exact.",sep=""))
    else stop(paste("ERGM fit `",deparse(substitute(object)),"` is missing the sample.",sep=""))
  }

  m<-ergm.getmodel(object$formula, ergm.getnetwork(object$formula), response=object$response)

  samp<-object$sample

  if(verbose) cat("Compressing the statistics: ")
  
  ## "Compress" duplicate statistics:
  library(coda)
  samp<-compress.data.frame(as.data.frame(samp))
  samp.fr<-samp$frequencies
  samp<-as.matrix(samp$rows)
  S<-nrow(samp)

  if(verbose) cat(S,"/",nrow(object$sample)," distinct configurations.\n",sep="")
  
  ## Because the MCMC sample that came with the ergm is drawn not from
  ## the MLE, but from the iteration prior to the MLE, we need to
  ## reweight the realizations for the bootstrap sample:

  if(verbose) cat("Reweighting observations and sampling: ")
  
  theta.samp<-object$MCMCtheta
  eta.samp<-ergm.eta(theta.samp, m$etamap)
  theta.mle<-coef(object)
  eta.mle<-ergm.eta(theta.mle, m$etamap)

  eta.diff<-eta.mle-eta.samp

  if(verbose>1) cat("\nMCMC sample at [",theta.samp,"] for MLE at [",theta.mle,"], for a canonical difference of [",eta.diff,"].\n")
  
  samp.w<-exp(samp %*% cbind(eta.diff))*samp.fr
  samp.w<-samp.w/sum(samp.w)

  ## Now, the expensive part of this procedure is the many
  ## optimization steps, so resample and compress:

  resamp<-sample.int(S,R,replace=TRUE,prob=samp.w)
  resamp.w<-tabulate(resamp, S)
  resamp.l<-which(resamp.w!=0)

  if(verbose) cat(length(resamp.l), "distinct configurations.\n")

  if(verbose) cat("Running estimation:\n")
  
  theta.boot<-apply(samp[resamp.l,],1,function(stat){
    v<-ergm.estimate(init=theta.samp,model=m,statsmatrix=sweep(samp,2,stat),statsmatrix.obs=NULL,
                     epsilon=control$epsilon,
                     nr.maxit=control$MCMLE.NR.maxit,
                     nr.reltol=control$MCMLE.NR.reltol,
                     calc.mcmc.se=FALSE, hessianflag=FALSE,
                     trustregion=+Inf, method=control$MCMLE.method,
                     metric=control$MCMLE.metric,
                     compress=control$MCMC.compress, verbose=max(verbose-2,0),
                     estimateonly=TRUE)
    v$coef
  })

  resamp.ind<-rep(seq_along(resamp.l),resamp.w[resamp.w!=0])
  if(verbose) cat("Done.\n")  
  t(theta.boot)[resamp.ind,]
}

autoboot.se.ergm<-function(object, theta.boot=NULL, R, verbose=FALSE, control=object$control){
  check.control.class("ergm")
  if(is.null(theta.boot)) theta.boot<-autoboot.ergm(object, R, verbose=FALSE, control=control.ergm())

  sqrt(apply(sweep(theta.boot,2,coef(object))^2,2,mean))
}


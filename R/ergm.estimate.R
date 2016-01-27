#  File R/ergm.estimate.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################
##################################################################################
# The <ergm.estimate> function searches for and returns a maximizer of the
# log-likelihood function. This function is observation-process capable.
#
# --PARAMETERS--
#   init          : the vector of theta parameters that produced 'statsmatrix'
#   model           : the model, as returned by <ergm.getmodel>
#   statsmatrix     : the matrix of observed statistics that has already had the
#                     "observed statistics" vector subtracted out (i.e., the
#                     "observed stats" are assumed to be zero here)
#   statsmatrix.obs: the corresponding statsmatrix for the observation process;
#                     default=NULL
#   nr.maxit        : the maximum number of iterations to use within the <optim>
#                     rountine; default=1000
#   nr.reltol       : the relative tolerance to use within the <optim> rountine;
#                     default=sqrt(.Machine$double.eps)
#   metric          : the name of a metric to use, as either "Likelihood",
#                     "lognormal", "Median.Likelihood", or "EF.Likelihood";
#                     default="lognormal"
#   method          : the method to be used by the <optim> routine;
#                     default="Nelder-Mead"
#   compress        : whether the matrix of statistics should be compressed
#                     via <ergm.sufftoprob>; default=FALSE
#   calc.mcmc.se    : whether to calculate the standard errors induced by the
#                     MCMC algorithm; default=TRUE
#   hessainflag     : whether the Hessian matrix of the likelihood function
#                     should be computed; default=TRUE
#   verbose         : whether the progress of the estimation should be printed;
#                     default=FALSE
#   trace           : a non-negative interger specifying how much tracing
#                     information should be printed by the <optim> routine;
#                     default=6*'verbose'
#   trustregion     : the trust region parameter for the likelihood functions
#   dampening       : (logical) should likelihood dampening be used?
#  dampening.min.ess: effective sample size below which dampening is used
#   dampening.level : proportional distance from boundary of the convex hull
#                     move
#   estimateonly    : whether only the estimates (vs. the estimates and the
#                     standard errors) should be calculated; default=FALSE
#
#
# --IGNORED PARAMETERS--
#   epsilon         : ??; default=1e-10
#
#
# --RETURNED--
#   an ergm object as a list containing several components; for details
#   see the return list in the <ergm> function header (<ergm.estimate>=^)
#
###################################################################################         

ergm.estimate<-function(init, model, statsmatrix, statsmatrix.obs=NULL,
                        epsilon=1e-10, nr.maxit=1000, nr.reltol=sqrt(.Machine$double.eps),
                        metric="lognormal",
                        method="Nelder-Mead", compress=FALSE,
                        calc.mcmc.se=TRUE, hessianflag=TRUE,
                        verbose=FALSE, trace=6*verbose,
                        trustregion=20, 
                        dampening=FALSE,
                        dampening.min.ess=100,
                        dampening.level=0.1,
                        cov.type="normal",# cov.type="robust", 
                        estimateonly=FALSE, ...) {
  # If there is an observation process to deal with, statsmatrix.obs will not be NULL;
  # in this case, do some preprocessing.  Otherwise, skip ahead.
  obsprocess <- !is.null(statsmatrix.obs)
  if (obsprocess) {
    if (compress) { # See comment below for explanation of "compress"
      statsmatrix0.obs <- ergm.sufftoprob(statsmatrix.obs, compress=TRUE)
      probs.obs <- statsmatrix0.obs[,ncol(statsmatrix0.obs)]
      statsmatrix0.obs <- statsmatrix0.obs[,-ncol(statsmatrix0.obs), drop=FALSE]      
    }
    statsmatrix0.obs <- statsmatrix.obs
    probs.obs <- rep(1/nrow(statsmatrix0.obs),nrow(statsmatrix0.obs))
  } else { # these objects should be defined as NULL so they exist later on
    probs.obs <- NULL
    xsim.obs <- NULL
  }

  # Now check to see whether to compress the statsmatrix by searching for
  # nonunique rows.  After compression, rows should be unique and each row
  # has a 'prob' weight telling what proportion of the original rows match it.
  if(compress){
    if (verbose) { cat("Compressing the matrix of sampled sufficient statistcs.\n") }
    statsmatrix0 <- ergm.sufftoprob(statsmatrix,compress=TRUE)
    probs <- statsmatrix0[,ncol(statsmatrix0)]
    statsmatrix0 <- statsmatrix0[,-ncol(statsmatrix0), drop=FALSE]
  } else {
    statsmatrix0 <- statsmatrix
    probs <- rep(1/nrow(statsmatrix0),nrow(statsmatrix0))
  }

  # It is assumed that the statsmatrix0 matrix has already had the
  # "observed statistics" subtracted out.  Another way to say this is that
  # when ergm.estimate is called, the "observed statistics" should equal
  # zero when measured on the scale of the statsmatrix0 statistics.
  # Here, we recenter the statsmatrix0 matrix by subtracting some measure
  # of center (e.g., the column means).  Since this shifts the scale, the
  # value of xobs (playing the role of "observed statistics") must be
  # adjusted accordingly.
# av <- apply(sweep(statsmatrix0,1,probs,"*"), 2, sum)
  if(cov.type=="robust"){
   av <- apply(statsmatrix0,2,wtd.median,weight=probs)
   V=try(
      covMcd(statsmatrix0[,!model$etamap$offsetmap,drop=FALSE])$cov,
       silent=TRUE)
   if(inherits(V,"try-error")){
    V=cov(statsmatrix0[,!model$etamap$offsetmap,drop=FALSE])
   }
  }else{
   av <- apply(statsmatrix0, 2, weighted.mean, w=probs)
   V=cov(statsmatrix0[,!model$etamap$offsetmap,drop=FALSE])
  }
  xsim <- sweep(statsmatrix0, 2, av,"-")
  xobs <-  -av 
  # Do the same recentering for the statsmatrix0.obs matrix, if appropriate.
  # Note that xobs must be adjusted too.
  if(obsprocess) {
    if(cov.type=="robust"){
     V.obs=try(
        covMcd(statsmatrix0.obs[,!model$etamap$offsetmap,drop=FALSE])$cov,
               silent=TRUE)
     if(inherits(V.obs,"try-error")){
      V.obs=cov(statsmatrix0.obs[,!model$etamap$offsetmap,drop=FALSE])
      av.obs <- apply(statsmatrix0.obs,2,wtd.median,weight=probs.obs)
     }
    }else{
     V.obs=cov(statsmatrix0.obs[, !model$etamap$offsetmap, drop=FALSE])
     av.obs <- apply(statsmatrix0.obs, 2, weighted.mean, w=probs.obs)
    }
#   av.obs <- apply(sweep(statsmatrix0.obs,1,probs.obs,"*"), 2, sum)
    xsim.obs <- sweep(statsmatrix0.obs, 2, av.obs,"-")
    xobs <- av.obs-av
  }
  
  # Convert init (possibly "curved" parameters) to eta0 (canonical parameters)
  eta0 <- ergm.eta(init, model$etamap)

  # Choose appropriate loglikelihood, gradient, and Hessian functions
  # depending on metric chosen and also whether obsprocess==TRUE
  # Also, choose varweight multiplier for covariance term in loglikelihood
  # where 0.5 is the "true" value but this can be increased or decreased
  varweight <- 0.5
  if (verbose) { cat("Using", metric, "metric (see control.ergm function).\n") }
  if (obsprocess) {
    loglikelihoodfn <- switch(metric,
                              Likelihood=llik.fun.obs,
                              lognormal=llik.fun.obs,
                              logtaylor=llik.fun.obs,
                              Median.Likelihood=llik.fun.obs,
                              EF.Likelihood=llik.fun.obs,
                              llik.fun.obs.robust)
    gradientfn <- switch(metric,
                         Likelihood=llik.grad.obs,
                         lognormal=llik.grad.obs,
                         logtaylor=llik.grad.obs,
                         Median.Likelihood=llik.grad.obs,
                         EF.Likelihood=llik.grad.obs,
                         llik.grad.obs)
    Hessianfn <- switch(metric,
                        Likelihood=llik.hessian.obs,
                        lognormal=llik.hessian.obs,
                        logtaylor=llik.hessian.obs,
                        Median.Likelihood=llik.hessian.obs,
                        EF.Likelihood=llik.hessian.obs,
                        llik.hessian.obs)
  } else {
    loglikelihoodfn <- switch(metric,
                              Likelihood=llik.fun,
                              lognormal=llik.fun,
                              logtaylor=llik.fun.logtaylor,
                              Median.Likelihood=llik.fun.median,
                              EF.Likelihood=llik.fun.EF,
                              llik.fun2)
    gradientfn <- switch(metric,
                         Likelihood=llik.grad,
                         lognormal=llik.grad,
                         logtaylor=llik.grad,
                         Median.Likelihood=llik.grad,
                         EF.Likelihood=llik.grad,
                         llik.grad)
    Hessianfn <- switch(metric,
                        Likelihood=llik.hessian,
                        lognormal=llik.hessian,
                        logtaylor=llik.hessian,
                        Median.Likelihood=llik.hessian,
                        EF.Likelihood=llik.hessian,
                        llik.hessian)
  }
  
  # Now find maximizer of approximate loglikelihood ratio l(eta) - l(eta0).
  # First: If we're using the lognormal approximation, the maximizer is
  # closed-form.  We can't use the closed-form maximizer if we are
  # dealing with a curved exponential family.
  if (all(model$etamap$canonical!=0) && 
      (metric=="lognormal" || metric=="Likelihood")) {
    if (obsprocess) {
      if (verbose) { cat("Using log-normal approx with missing (no optim)\n") }
      Lout <- list(hessian = -(V-V.obs))
    } else {
      if (verbose) { cat("Using log-normal approx (no optim)\n") }
      Lout <- list(hessian = -V)
    }
    Lout$par <- try(eta0[!model$etamap$offsetmap] 
                    - solve(Lout$hessian, xobs[!model$etamap$offsetmap]),
                    silent=TRUE)
    # If there's an error, first try a robust matrix inverse.  This can often
    # happen if the matrix of simulated statistics does not ever change for one
    # or more statistics.
    if(inherits(Lout$par,"try-error")){
      Lout$par <- try(eta0[!model$etamap$offsetmap] 
                      - ginv(Lout$hessian) %*% 
                      xobs[!model$etamap$offsetmap],
                      silent=TRUE)
    }
    # If there's still an error, use the Matrix package to try to find an 
    # alternative Hessian approximant that has no zero eigenvalues.
    if(inherits(Lout$par,"try-error")){
      if (obsprocess) {
        Lout <- list(hessian = -(as.matrix(nearPD(V-V.obs)$mat)))
      }else{
        Lout <- list(hessian = -(as.matrix(nearPD(V)$mat)))
      }
      Lout$par <- eta0[!model$etamap$offsetmap] - solve(Lout$hessian, xobs[!model$etamap$offsetmap])
    }
    Lout$convergence <- 0 # maybe add some error-checking here to get other codes
    Lout$value <- 0.5*crossprod(xobs[!model$etamap$offsetmap],
                                Lout$par - eta0[!model$etamap$offsetmap])
    hessianflag <- TRUE # to make sure we don't recompute the Hessian later on
  } else {
    # "guess" will be the starting point for the optim search algorithm.
    # But only the non-offset values are relevant; the others will be
    # passed to the likelihood functions by way of the etamap$init element
    # of the model object.  NB:  This is a really ugly way to do this!  Change it?
    guess <- init[!model$etamap$offsettheta]
    model$etamap$init <- init
    
    loglikelihoodfn.trust<-function(trustregion=20, ...){
      value<-loglikelihoodfn(trustregion=trustregion, ...)
      grad<-gradientfn(trustregion=trustregion, ...)
      hess<-Hessianfn(...)
      hess[upper.tri(hess)]<-t(hess)[upper.tri(hess)]
#      print(value)
#      print(grad)
#      print(hess)
      list(value=value,gradient=as.vector(grad),hessian=hess)
    }

    if (verbose) { cat("Optimizing loglikelihood\n") }
    Lout <- try(trust(objfun=loglikelihoodfn.trust, parinit=guess,
                      rinit=1, 
                      rmax=100, 
                      parscale=rep(1,length(guess)), minimize=FALSE,
                      xobs=xobs,
                      xsim=xsim, probs=probs,
                      xsim.obs=xsim.obs, probs.obs=probs.obs,
                      varweight=varweight, trustregion=trustregion,
                      dampening=dampening,
                      dampening.min.ess=dampening.min.ess,
                      dampening.level=dampening.level,
                      eta0=eta0, etamap=model$etamap),
            silent=FALSE)
    Lout$par<-Lout$argument
#   if(Lout$value < trustregion-0.001){
#     current.scipen <- options()$scipen
#     options(scipen=3)
#     cat("the log-likelihood improved by",
#         format.pval(Lout$value,digits=4,eps=1e-4),"\n")
#     options(scipen=current.scipen)
#   }else{
#     cat("the log-likelihood did not improve.\n")
#   }
    if(inherits(Lout,"try-error") || Lout$value > max(199, trustregion) || Lout$value < -790) {
      if(!inherits(Lout,"try-error")) cat("Apparent likelihood improvement:", Lout$value, ".\n")
      cat("MLE could not be found. Trying Nelder-Mead...\n")
      Lout <- try(optim(par=guess, 
                        fn=llik.fun.median,
                        hessian=hessianflag,
                        method="Nelder-Mead",
                        control=list(trace=trace,fnscale=-1,maxit=100*nr.maxit,
                                     reltol=nr.reltol),
                        xobs=xobs, 
                        xsim=xsim, probs=probs, 
                        xsim.obs=xsim.obs, probs.obs=probs.obs,
                        varweight=varweight, trustregion=trustregion,
                        eta0=eta0, etamap=model$etamap),
              silent=FALSE)
      if(inherits(Lout,"try-error") || Lout$value > max(500, trustregion) ){
        cat(paste("No direct MLE exists!\n"))
      }
      if(Lout$convergence != 0 ){
        cat("Non-convergence after", nr.maxit, "iterations.\n")
      }
      cat("Nelder-Mead Log-likelihood ratio is ", Lout$value,"\n")
    }
  }

  theta <- init
  theta[!model$etamap$offsettheta] <- Lout$par
  names(theta) <- names(init)
  if (estimateonly) {
    # Output results as ergm-class object
    return(structure(list(coef=theta,
                          MCMCtheta=init,
                          samplesize=NROW(statsmatrix),
                          loglikelihood=Lout$value, 
                          failure=FALSE),
                        class="ergm"))
  } else {
    gradienttheta <- llik.grad(theta=Lout$par, xobs=xobs, xsim=xsim,
                          probs=probs, 
                          xsim.obs=xsim.obs, probs.obs=probs.obs,
                          varweight=varweight, eta0=eta0, etamap=model$etamap)
    gradient <- rep(NA, length=length(init))
    gradient[!model$etamap$offsettheta] <- gradienttheta
    #
    #  Calculate the auto-covariance of the MCMC suff. stats.
    #  and hence the MCMC s.e.
    #
    mc.se <- rep(NA, length=length(theta))
    mc.cov <- matrix(NA, length(theta), length(theta))
    covar <- NA
    if(!hessianflag){
      #  covar <- ginv(cov(xsim))
      #  Lout$hessian <- cov(xsim)
      Lout$hessian <- Hessianfn(theta=Lout$par, xobs=xobs, xsim=xsim,
                                probs=probs, 
                                xsim.obs=xsim.obs, probs.obs=probs.obs,
                                varweight=varweight,
                                eta0=eta0, etamap=model$etamap
                                )
    }
    
    covar <- matrix(NA, ncol=length(theta), nrow=length(theta))
    covar[!model$etamap$offsettheta,!model$etamap$offsettheta ] <- ginv(-Lout$hessian)
    dimnames(covar) <- list(names(theta),names(theta))
    He <- matrix(NA, ncol=length(theta), nrow=length(theta))
    He[!model$etamap$offsettheta,!model$etamap$offsettheta ] <- Lout$hessian
    dimnames(He) <- list(names(theta),names(theta))
    Lout$hessian <- He
    
    if(calc.mcmc.se){
      if (verbose) { cat("Starting MCMC s.e. computation.\n") }
      mc.cov <-
        if ((metric=="lognormal" || metric=="Likelihood")
            && length(model$etamap$curved)==0) {
          ergm.MCMCse.lognormal(theta=theta, init=init, 
                                statsmatrix=statsmatrix0, 
                                statsmatrix.obs=statsmatrix.obs,
                                H=V, H.obs=V.obs,
                                model=model)
        } else {
        ergm.MCMCse(theta=theta,init=init, 
                    statsmatrix=statsmatrix0,
                    statsmatrix.obs=statsmatrix.obs,
                    model=model)
      }
    }
    c0  <- loglikelihoodfn(theta=Lout$par, xobs=xobs,
                           xsim=xsim, probs=probs,
                           xsim.obs=xsim.obs, probs.obs=probs.obs,
                           varweight=0.5, eta0=eta0, etamap=model$etamap)
    c01 <- loglikelihoodfn(theta=Lout$par-Lout$par, xobs=xobs,
                           xsim=xsim, probs=probs,
                           xsim.obs=xsim.obs, probs.obs=probs.obs,
                           varweight=0.5, eta0=eta0, etamap=model$etamap)
    #
    # This is the log-likelihood calc from init=0
    #
    mcmcloglik <- -abs(c0 - c01)
    
    #   c1 <- theta1$loglikelihood
    #   c1  <- c01
    # loglikelihood <- mcmcloglik
    loglikelihood <- Lout$value
    
    #
    # Use the psuedo-likelihood as a base
    #
    iteration <- Lout$counts[1]
    #
    names(theta) <- names(init)
    
    ############################
    #
    #  Reconstruct the reduced form statsmatrix for curved 
    #  parametrizations
    #
    statsmatrix.all <- statsmatrix
    statsmatrix <- ergm.curved.statsmatrix(statsmatrix,theta,model$etamap)$sm
    statsmatrix0 <- statsmatrix
    if(all(dim(statsmatrix.all)==dim(statsmatrix))){
      statsmatrix.all <- NULL
    }

    # Output results as ergm-class object
    return(structure(list(coef=theta, sample=statsmatrix, sample.obs=statsmatrix.obs, 
                          iterations=iteration, #mcmcloglik=mcmcloglik,
                          MCMCtheta=init, 
                          loglikelihood=loglikelihood, gradient=gradient, hessian=Lout$hessian,
                          covar=covar, failure=FALSE,
                          mc.cov=mc.cov #, #acf=mcmcacf,
                          #fullsample=statsmatrix.all
                          ),
                        class="ergm"))
  }
}

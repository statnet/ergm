#  File R/ergm.estimate.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################
##################################################################################
# The <ergm.estimate> function searches for and returns a maximizer of the
# log-likelihood function. This function is observation-process capable.
#
# --PARAMETERS--
#   init          : the vector of theta parameters that produced 'statsmatrix'
#   model           : the model, as returned by <ergm_model>
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
#   calc.mcmc.se    : whether to calculate the standard errors induced by the
#                     MCMC algorithm; default=TRUE
#   hessainflag     : whether the Hessian matrix of the likelihood function
#                     should be computed; default=TRUE
#   verbose         : whether the progress of the estimation should be printed;
#                     default=FALSE
#   trace           : a non-negative interger specifying how much tracing
#                     information should be printed by the <optim> routine;
#                     default=6*'verbose'
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

ergm.estimate<-function(init, model, statsmatrices, statsmatrices.obs=NULL,
                        epsilon=1e-10, nr.maxit=1000, nr.reltol=sqrt(.Machine$double.eps),
                        metric="lognormal",
                        method="Nelder-Mead",
                        calc.mcmc.se=TRUE, hessianflag=TRUE,
                        verbose=FALSE, trace=6*verbose,
                        dampening=FALSE,
                        dampening.min.ess=100,
                        dampening.level=0.1,
                        steplen=1,
                        cov.type="normal",# cov.type="robust", 
                        estimateonly=FALSE, ...) {
  estimateonly <- estimateonly & !calc.mcmc.se
  # If there is an observation process to deal with, statsmatrices.obs
  # will not be NULL.
  obsprocess <- !is.null(statsmatrices.obs)

  # Construct an offsetless map and convert init (possibly "curved"
  # parameters) to eta0 (canonical parameters)
  etamap <- model$etamap
  etamap.no <- deoffset.etamap(etamap, init)
  eta0 <- ergm.eta(init[!etamap$offsettheta], etamap.no)


  statsmatrices.orig <- statsmatrices
  statsmatrices.orig.obs <- statsmatrices.obs

  statsmean <- colMeans.mcmc.list(statsmatrices.orig)
  if(!is.null(statsmatrices.orig.obs)){
    statsmatrices.obs <- lapply.mcmc.list(statsmatrices.orig.obs, .shift_scale_points, statsmean, steplen) # I.e., shrink each point of statsmatrix.obs towards the centroid of statsmatrix.
  }else{
    statsmatrices <- lapply.mcmc.list(statsmatrices.orig,sweep,2,(1-steplen)*statsmean,"-")
  }
  
  statsmatrix <- as.matrix(statsmatrices)
  if(obsprocess) statsmatrix.obs <- as.matrix(statsmatrices.obs)
  statsmatrix.orig <- as.matrix(statsmatrices.orig)
  if(obsprocess) statsmatrix.orig.obs <- as.matrix(statsmatrices.orig.obs)

    
  # Copy and compress the stats matrices after dropping the offset terms.
  xsim <- compress_rows(as.logwmatrix(statsmatrix[,!etamap$offsetmap, drop=FALSE]))
  lrowweights(xsim) <- -log_sum_exp(lrowweights(xsim)) # I.e., divide all weights by their sum.

  xsim.orig <- compress_rows(as.logwmatrix(statsmatrix.orig[,!etamap$offsetmap, drop=FALSE]))
  lrowweights(xsim.orig) <- -log_sum_exp(lrowweights(xsim.orig)) # I.e., divide all weights by their sum.

  if(obsprocess){
    xsim.obs <- compress_rows(as.logwmatrix(statsmatrix.obs[,!etamap$offsetmap, drop=FALSE]))
    lrowweights(xsim.obs) <- -log_sum_exp(lrowweights(xsim.obs)) 

    xsim.orig.obs <- compress_rows(as.logwmatrix(statsmatrix.orig.obs[,!etamap$offsetmap, drop=FALSE]))
    lrowweights(xsim.orig.obs) <- -log_sum_exp(lrowweights(xsim.orig.obs)) 
  }else{
    xsim.obs <- NULL
    xsim.orig.obs <- NULL
  }
  
  # It is assumed that the statsmatrix matrix has already had the
  # "observed statistics" subtracted out.  Another way to say this is
  # that when ergm.estimate is called, the "observed statistics"
  # should equal zero when measured on the scale of the statsmatrix
  # statistics.
  #' @importFrom robustbase covMcd
  if(cov.type=="robust"){
    tmp <- covMcd(decompress_rows(xsim, target.nrows=nrow(statsmatrix)))    
    av <- tmp$center
    V <- tmp$cov
  }else{
    av <- lweighted.mean(xsim, lrowweights(xsim))
    V <- lweighted.var(xsim, lrowweights(xsim))
  }
  
  xobs <- -av
  
  # Do the same recentering for the statsmatrix.obs matrix, if appropriate.
  # Note that xobs must be adjusted too.
  if(obsprocess) {
    if(cov.type=="robust"){
      tmp <- covMcd(decompress_rows(xsim.obs, target.nrows=nrow(statsmatrix.obs)))    
      av.obs <- tmp$center
      V.obs <- tmp$cov
    }else{
      av.obs <- lweighted.mean(xsim.obs, lrowweights(xsim.obs))
      V.obs <- lweighted.var(xsim.obs, lrowweights(xsim.obs))
    }

    xobs <- av.obs - av
  }
    
  # Choose appropriate loglikelihood, gradient, and Hessian functions
  # depending on metric chosen and also whether obsprocess==TRUE
  # Also, choose varweight multiplier for covariance term in loglikelihood
  # where 0.5 is the "true" value but this can be increased or decreased
  varweight <- 0.5
  if (verbose) { message("Using ", metric, " metric (see control.ergm function).") }
  if (obsprocess) {
    loglikelihoodfn <- switch(metric,
                              Likelihood=llik.fun.obs.lognormal,
                              lognormal=llik.fun.obs.lognormal,
                              logtaylor=llik.fun.obs.lognormal,
                              Median.Likelihood=llik.fun.obs.robust,
                              EF.Likelihood=llik.fun.obs.lognormal,
                              llik.fun.obs.IS)
    gradientfn <- switch(metric,
                         Likelihood=llik.grad.obs.IS,
                         lognormal=llik.grad.obs.IS,
                         logtaylor=llik.grad.obs.IS,
                         Median.Likelihood=llik.grad.obs.IS,
                         EF.Likelihood=llik.grad.obs.IS,
                         llik.grad.obs.IS)
    Hessianfn <- switch(metric,
                        Likelihood=llik.hessian.obs.IS,
                        lognormal=llik.hessian.obs.IS,
                        logtaylor=llik.hessian.obs.IS,
                        Median.Likelihood=llik.hessian.obs.IS,
                        EF.Likelihood=llik.hessian.obs.IS,
                        llik.hessian.obs.IS)
  } else {
    loglikelihoodfn <- switch(metric,
                              Likelihood=llik.fun.lognormal,
                              lognormal=llik.fun.lognormal,
                              logtaylor=llik.fun.logtaylor,
                              Median.Likelihood=llik.fun.median,
                              EF.Likelihood=llik.fun.EF,
                              llik.fun.IS)
    gradientfn <- switch(metric,
                         Likelihood=llik.grad.IS,
                         lognormal=llik.grad.IS,
                         logtaylor=llik.grad.IS,
                         Median.Likelihood=llik.grad.IS,
                         EF.Likelihood=llik.grad.IS,
                         llik.grad.IS)
    Hessianfn <- switch(metric,
                        Likelihood=llik.hessian.IS,
                        lognormal=llik.hessian.IS,
                        logtaylor=llik.hessian.IS,
                        Median.Likelihood=llik.hessian.IS,
                        EF.Likelihood=llik.hessian.IS,
                        llik.hessian.IS)
  }
  
  # Now find maximizer of approximate loglikelihood ratio l(eta) - l(eta0).
  # First: If we're using the lognormal approximation, the maximizer is
  # closed-form.  We can't use the closed-form maximizer if we are
  # dealing with a curved exponential family.
  if (!is.curved(model) &&
      (metric=="lognormal" || metric=="Likelihood") &&
      all(model$etamap$mintheta==-Inf) &&
      all(model$etamap$maxtheta==+Inf)) {
    if (obsprocess) {
      if (verbose) { message("Using log-normal approx with missing (no optim)") }
      # Here, setting posd.tol=0 ensures that the matrix is
      # nonnegative-definite: it is possible for some simulated
      # statistics not to change, but it is not possible for the
      # constrained sample to have a higher variance than the
      # unconstrained.
      Lout <- list(hessian = -as.matrix(snearPD(V-V.obs,posd.tol=0)$mat))
    } else {
      if (verbose) { message("Using log-normal approx (no optim)") }
      Lout <- list(hessian = -V)
    }
    Lout$par <- try(eta0 
                    + ssolve(-Lout$hessian, xobs),
                    silent=TRUE)
    # If there's an error, first try a robust matrix inverse.  This can often
    # happen if the matrix of simulated statistics does not ever change for one
    # or more statistics.
    if(inherits(Lout$par,"try-error")){
      Lout$par <- try(eta0 
                      + sginv(-Lout$hessian) %*%
                      xobs,
                      silent=TRUE)
    }
    # If there's still an error, use the Matrix package to try to find an 
    # alternative Hessian approximant that has no zero eigenvalues.
    if(inherits(Lout$par,"try-error")){
      if (obsprocess) {
        Lout <- list(hessian = -(as.matrix(snearPD(V-V.obs)$mat)))
      }else{
        Lout <- list(hessian = -(as.matrix(snearPD(V)$mat)))
      }
      Lout$par <- eta0 + ssolve(-Lout$hessian, xobs)
    }
    Lout$convergence <- 0 # maybe add some error-checking here to get other codes
    Lout$value <- 0.5*crossprod(xobs,
                                Lout$par - eta0)
    hessianflag <- TRUE # to make sure we don't recompute the Hessian later on
  } else {
    # "guess" will be the starting point for the optim search algorithm.
    guess <- init[!model$etamap$offsettheta]
    
    loglikelihoodfn.trust<-function(theta, ...){
      # Check for box constraint violation.
      if(anyNA(theta) ||
         any(theta < model$etamap$mintheta[!model$etamap$offsettheta]) ||
         any(theta > model$etamap$maxtheta[!model$etamap$offsettheta]))
        return(list(value=-Inf))
      
      value<-loglikelihoodfn(theta, ...)
      grad<-gradientfn(theta, ...)
      hess<-Hessianfn(theta, ...)
      hess[upper.tri(hess)]<-t(hess)[upper.tri(hess)]
#      message_print(value)
#      message_print(grad)
#      message_print(hess)
      list(value=value,gradient=as.vector(grad),hessian=hess)
    }

    if (verbose) { message("Optimizing loglikelihood") }
    #' @importFrom trust trust
    Lout <- try(trust(objfun=loglikelihoodfn.trust, parinit=guess,
                      rinit=1, 
                      rmax=100, 
                      parscale=rep(1,length(guess)), minimize=FALSE,
                      xsim=xsim,
                      xsim.obs=xsim.obs,
                      varweight=varweight,
                      dampening=dampening,
                      dampening.min.ess=dampening.min.ess,
                      dampening.level=dampening.level,
                      eta0=eta0, etamap=etamap.no),
                silent=FALSE)
    if(inherits(Lout,"try-error")) {
      message("MLE could not be found. Trying Nelder-Mead...")
      Lout <- try(optim(par=guess, 
                        fn=llik.fun.median,
                        hessian=hessianflag,
                        method="Nelder-Mead",
                        control=list(trace=trace,fnscale=-1,maxit=100*nr.maxit,
                                     reltol=nr.reltol),
                        xsim=xsim,
                        xsim.obs=xsim.obs,
                        varweight=varweight,
                        dampening=dampening,
                        dampening.min.ess=dampening.min.ess,
                        dampening.level=dampening.level,
                        eta0=eta0, etamap=etamap.no),
              silent=FALSE)
      if(inherits(Lout,"try-error")){
        message(paste("No direct MLE exists!"))
      }
      if(Lout$convergence != 0 ){
        message("Non-convergence after ", nr.maxit, " iterations.")
      }
      message("Nelder-Mead Log-likelihood ratio is ", Lout$value," ")
    } else Lout$par<-Lout$argument
  }

  theta <- init
  theta[!model$etamap$offsettheta] <- Lout$par
  names(theta) <- names(init)
  if (estimateonly) {
    # Output results as ergm-class object
    list(coefficients = setNames(theta, names(init)),
         MCMCtheta = init,
         samplesize = nrow(statsmatrix),
         loglikelihood = Lout$value,
         failure = FALSE)
  } else {
    gradienttheta <- llik.grad.IS(theta=Lout$par,
                        xsim=xsim,
                        xsim.obs=xsim.obs,
                        varweight=varweight,
                        eta0=eta0, etamap=etamap.no)
    gradient <- rep(NA, length=length(init))
    gradient[!model$etamap$offsettheta] <- gradienttheta
    #
    #  Calculate the auto-covariance of the MCMC suff. stats.
    #  and hence the MCMC s.e.
    #
    mc.se <- rep(NA, length=length(theta))
    mc.cov <- matrix(NA, length(theta), length(theta))
    covar <- NA
    if(!hessianflag || steplen!=1){
      Lout$hessian <- Hessianfn(theta=Lout$par,
                        xsim=xsim.orig,
                        xsim.obs=xsim.orig.obs,
                        varweight=varweight,
                        eta0=eta0, etamap=etamap.no
                        )
    }
    
    covar <- matrix(NA, ncol=length(theta), nrow=length(theta))
    covar[!model$etamap$offsettheta, !model$etamap$offsettheta] <- sginv(-Lout$hessian, tol=.Machine$double.eps^(3/4))
    dimnames(covar) <- list(names(theta),names(theta))
    He <- matrix(NA, ncol=length(theta), nrow=length(theta))
    He[!model$etamap$offsettheta, !model$etamap$offsettheta] <- Lout$hessian
    dimnames(He) <- list(names(theta),names(theta))
    Lout$hessian <- He
    
    if(calc.mcmc.se){
      if (verbose) message("Starting MCMC s.e. computation.")
      mcse.metric <-
        if ((metric == "lognormal" || metric == "Likelihood") && length(model$etamap$curved) == 0) "lognormal"
        else "IS"
      mc.cov <- ergm.MCMCse(model = model, theta = theta, init = init,
                            statsmatrices = statsmatrices,
                            statsmatrices.obs = statsmatrices.obs,
                            H = V, H.obs = V.obs,
                            metric = mcse.metric)
    }

    # Output results as ergm-class object
    list(
      coefficients = setNames(theta, names(init)),
      sample = statsmatrices, sample.obs = statsmatrices.obs,
      iterations = Lout$counts[1], MCMCtheta = init,
      loglikelihood = Lout$value, gradient = gradient, hessian = Lout$hessian,
      covar = covar, failure = FALSE, mc.cov = mc.cov
    )
  }
}

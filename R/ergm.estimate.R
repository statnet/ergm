##################################################################################
# The <ergm.estimate> function searches for and returns a maximizer of the
# log-likelihood function. This function is missing-data capable.
#
# --PARAMETERS--
#   theta0          : the vector of theta parameters that produced 'statsmatrix'
#   model           : the model, as returned by <ergm.getmodel>
#   statsmatrix     : the matrix of observed statistics that has already had the
#                     "observed statistics" vector subtracted out (i.e., the
#                     "observed stats" are assumed to be zero here)
#   statsmatrix.miss: the corresponding statsmatrix for the missing edges;
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

ergm.estimate<-function(theta0, model, statsmatrix, statsmatrix.miss=NULL,
                        epsilon=1e-10, nr.maxit=1000, nr.reltol=sqrt(.Machine$double.eps),
                        metric="lognormal",
                        method="Nelder-Mead", compress=FALSE,
                        calc.mcmc.se=TRUE, hessianflag=TRUE,
                        verbose=FALSE, trace=6*verbose,
                        trustregion=20, 
                        estimateonly=FALSE, ...) {
  # If there are missing data to deal with, statsmatrix.miss will not be NULL;
  # in this case, do some preprocessing.  Otherwise, skip ahead.
  missingflag <- !is.null(statsmatrix.miss)
  if (missingflag) {
    if (compress) { # See comment below for explanation of "compress"
      statsmatrix0.miss <- ergm.sufftoprob(statsmatrix.miss, compress=TRUE)
      probs.miss <- statsmatrix0.miss[,ncol(statsmatrix0.miss)]
      statsmatrix0.miss <- statsmatrix0.miss[,-ncol(statsmatrix0.miss), drop=FALSE]      
    }
    statsmatrix0.miss <- statsmatrix.miss
    probs.miss <- rep(1/nrow(statsmatrix0.miss),nrow(statsmatrix0.miss))
  } else { # these objects should be defined as NULL so they exist later on
    probs.miss <- NULL
    xsim.miss <- NULL
  }

  # Now check to see whether to compress the statsmatrix by searching for
  # nonunique rows.  After compression, rows should be unique and each row
  # has a 'prob' weight telling what proportion of the original rows match it.
  if(compress){
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
   av <- apply(statsmatrix0, 2, weighted.mean, w=probs)
   V=cov(statsmatrix0[,!model$etamap$offsetmap,drop=FALSE])
  xsim <- sweep(statsmatrix0, 2, av,"-")
  xobs <-  -av 
  # Do the same recentering for the statsmatrix0.miss matrix, if appropriate.
  # Note that xobs must be adjusted too.
  if(missingflag) {
     V.miss=cov(statsmatrix0.miss[, !model$etamap$offsetmap, drop=FALSE])
     av.miss <- apply(statsmatrix0.miss, 2, weighted.mean, w=probs.miss)
#   av.miss <- apply(sweep(statsmatrix0.miss,1,probs.miss,"*"), 2, sum)
    xsim.miss <- sweep(statsmatrix0.miss, 2, av.miss,"-")
    xobs <- av.miss-av
  }
  
  # Convert theta0 (possibly "curved" parameters) to eta0 (canonical parameters)
  eta0 <- ergm.eta(theta0, model$etamap)

  # Choose appropriate loglikelihood, gradient, and Hessian functions
  # depending on metric chosen and also whether missingflag==TRUE
  # Also, choose varweight multiplier for covariance term in loglikelihood
  # where 0.5 is the "true" value but this can be increased or decreased
  varweight <- 0.5
  if (verbose) { cat("Using", metric, "metric (see control.ergm function).\n") }
  if (missingflag) {
    loglikelihoodfn <- switch(metric,
                              Likelihood=llik.fun.miss,
                              lognormal=llik.fun.miss,
                              Median.Likelihood=llik.fun.miss,
                              EF.Likelihood=llik.fun.miss,
                              llik.fun.miss.robust)
    gradientfn <- switch(metric,
                         Likelihood=llik.grad.miss,
                         lognormal=llik.grad.miss,
                         Median.Likelihood=llik.grad.miss,
                         EF.Likelihood=llik.grad.miss,
                         llik.grad.miss)
    Hessianfn <- switch(metric,
                        Likelihood=llik.hessian.miss,
                        lognormal=llik.hessian.miss,
                        Median.Likelihood=llik.hessian.miss,
                        EF.Likelihood=llik.hessian.miss,
                        llik.hessian.miss)
  } else {
    loglikelihoodfn <- switch(metric,
                              Likelihood=llik.fun,
                              lognormal=llik.fun,
                              Median.Likelihood=llik.fun.median,
                              EF.Likelihood=llik.fun.EF,
                              llik.fun2)
    gradientfn <- switch(metric,
                         Likelihood=llik.grad,
                         lognormal=llik.grad,
                         Median.Likelihood=llik.grad,
                         EF.Likelihood=llik.grad,
                         llik.grad2)
    Hessianfn <- switch(metric,
                        Likelihood=llik.hessian,
                        lognormal=llik.hessian,
                        Median.Likelihood=llik.hessian,
                        EF.Likelihood=llik.hessian,
                        llik.hessian2)
  }
  
  # Now find maximizer of approximate loglikelihood ratio l(eta) - l(eta0).
  # First: If we're using the lognormal approximation, the maximizer is
  # closed-form.  We can't use the closed-form maximizer if we are
  # dealing with a curved exponential family.
  if (all(model$etamap$canonical==1) && 
      (metric=="lognormal" || metric=="Likelihood")) {
    if (missingflag) {
      if (verbose) { cat("Using log-normal approx with missing (no optim)\n") }
      Lout <- list(hessian = -(V-V.miss))
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
                      - robust.inverse(Lout$hessian) %*% 
                      xobs[!model$etamap$offsetmap],
                      silent=TRUE)
    }
    # If there's still an error, use the Matrix package to try to find an 
    # alternative Hessian approximant that has no zero eigenvalues.
    if(inherits(Lout$par,"try-error") && !inherits(try(library(Matrix)), "try-error")){
      if (missingflag) {
        Lout <- list(hessian = -(as.matrix(nearPD(V-V.miss)$mat)))
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
    # passed to the likelihood functions by way of the etamap$theta0 element
    # of the model object.  NB:  This is a really ugly way to do this!  Change it?
    guess <- theta0[!model$etamap$offsettheta]
    model$etamap$theta0 <- theta0
    
    if (verbose) { cat("Optimizing loglikelihood\n") }
    Lout <- try(optim(par=guess,
                      fn=loglikelihoodfn,   gr=gradientfn,
                      hessian=hessianflag,
                      method=method,
                      control=list(trace=trace, fnscale=-1,
                                   maxit=nr.maxit,reltol=nr.reltol),
                      xobs=xobs,
                      xsim=xsim, probs=probs,
                      xsim.miss=xsim.miss, probs.miss=probs.miss,
                      varweight=varweight, trustregion=trustregion,
                      eta0=eta0, etamap=model$etamap),
            silent=FALSE)
#   if(Lout$value < trustregion-0.001){
#     current.scipen <- options()$scipen
#     options(scipen=3)
#     cat("the log-likelihood improved by",
#         format.pval(Lout$value,digits=4,eps=1e-4),"\n")
#     options(scipen=current.scipen)
#   }else{
#     cat("the log-likelihood did not improve.\n")
#   }
    if(inherits(Lout,"try-error") || Lout$value > 199 || Lout$value < -790) {
      cat("MLE could not be found. Trying Nelder-Mead...\n")
      Lout <- try(optim(par=guess, 
                        fn=loglikelihoodfn,
                        hessian=hessianflag,
                        method="Nelder-Mead",
                        control=list(trace=trace,fnscale=-1,maxit=100*nr.maxit,
                                     reltol=0.01),
                        xobs=xobs, 
                        xsim=xsim, probs=probs, 
                        xsim.miss=xsim.miss, probs.miss=probs.miss,
                        varweight=varweight, trustregion=trustregion,
                        eta0=eta0, etamap=model$etamap),
              silent=FALSE)
      if(inherits(Lout,"try-error") || Lout$value > 500 ){
        cat(paste("No direct MLE exists!\n"))
      }
      if(Lout$convergence != 0 ){
        cat("Non-convergence after", nr.maxit, "iterations.\n")
      }
      cat("Nelder-Mead Log-likelihood ratio is ", Lout$value,"\n")
    }
  }

  theta <- theta0
  theta[!model$etamap$offsettheta] <- Lout$par
  names(theta) <- names(theta0)
  if (estimateonly) {
    # Output results as ergm-class object
    return(structure(list(coef=theta,
                          MCMCtheta=theta0,
                          samplesize=NROW(statsmatrix),
                          loglikelihood=Lout$value, 
                          failure=FALSE),
                        class="ergm"))
  } else {
    gradienttheta <- llik.grad(theta=Lout$par, xobs=xobs, xsim=xsim,
                          probs=probs, 
                          xsim.miss=xsim.miss, probs.miss=probs.miss,
                          varweight=varweight, eta0=eta0, etamap=model$etamap)
    gradient <- rep(NA, length=length(theta0))
    gradient[!model$etamap$offsettheta] <- gradienttheta
    #
    #  Calculate the auto-covariance of the MCMC suff. stats.
    #  and hence the MCMC s.e.
    #
    mc.se <- rep(NA, length=length(theta))
    covar <- NA
    if(!hessianflag){
     #  covar <- robust.inverse(cov(xsim))
     #  Lout$hessian <- cov(xsim)
     Lout$hessian <- Hessianfn(theta=theta, xobs=xobs, xsim=xsim,
                               probs=probs, 
                               xsim.miss=xsim.miss, probs.miss=probs.miss,
                               varweight=varweight,
                               eta0=eta0, etamap=model$etamap
                               )
     covar <- matrix(NA, ncol=length(theta), nrow=length(theta))
     covar[!model$etamap$offsettheta,!model$etamap$offsettheta ] <- robust.inverse(-Lout$hessian)
     dimnames(covar) <- list(names(theta),names(theta))
    }else{
     covar <- matrix(NA, ncol=length(theta), nrow=length(theta))
     covar[!model$etamap$offsettheta,!model$etamap$offsettheta ] <- robust.inverse(-Lout$hessian)
     dimnames(covar) <- list(names(theta),names(theta))
     He <- matrix(NA, ncol=length(theta), nrow=length(theta))
     He[!model$etamap$offsettheta,!model$etamap$offsettheta ] <- Lout$hessian
     dimnames(He) <- list(names(theta),names(theta))
     Lout$hessian <- He
    }
    if(calc.mcmc.se){
      if (verbose) { cat("Starting MCMC s.e. computation.\n") }
      if ((metric=="lognormal" || metric=="Likelihood")
          && length(model$etamap$curved)==0) {
        mc.se <- ergm.MCMCse.lognormal(theta=theta, theta0=theta0, 
                                       statsmatrix=statsmatrix0, 
                                       statsmatrix.miss=statsmatrix.miss,
                                       H=V, H.miss=V.miss,
                                       model=model)
      } else {
        MCMCse <- ergm.MCMCse(theta=theta,theta0=theta0, 
                             statsmatrix=statsmatrix0,
                             statsmatrix.miss=statsmatrix.miss,
                             model=model)
        mc.se <- MCMCse$mc.se
#       The next line forces the s.e. in summary.ergm to combine
#       the hessian of the likelihood plus the MCMC s.e.
#       covar <- MCMCse$mc.cov+covar
#       If the above line is commented out only the hessian of the likelihood 
#       is used.
      }
    }
    c0  <- loglikelihoodfn(theta=Lout$par, xobs=xobs,
                           xsim=xsim, probs=probs,
                           xsim.miss=xsim.miss, probs.miss=probs.miss,
                           varweight=0.5, eta0=eta0, etamap=model$etamap)
    c01 <- loglikelihoodfn(theta=Lout$par-Lout$par, xobs=xobs,
                           xsim=xsim, probs=probs,
                           xsim.miss=xsim.miss, probs.miss=probs.miss,
                           varweight=0.5, eta0=eta0, etamap=model$etamap)
    #
    # This is the log-likelihood calc from theta0=0
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
    names(theta) <- names(theta0)
    
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

    # Commented out for now -- this should probably not be added to the
    # ergm object because it is not always needed.
    #if (verbose) cat("Starting MCMC s.e. ACF computation.\n")
    #if(calc.mcmc.se){
    #  mcmcacf <- ergm.MCMCacf(statsmatrix0)
    #}else{
    #  mcmcacf <- covar-covar
    #}
    #if (verbose) cat("Ending MCMC s.e. ACF computation.\n")

    # Output results as ergm-class object
    return(structure(list(coef=theta, sample=statsmatrix, sample.miss=statsmatrix.miss, 
                          iterations=iteration, #mcmcloglik=mcmcloglik,
                          MCMCtheta=theta0, 
                          loglikelihood=loglikelihood, gradient=gradient,
                          covar=covar, samplesize=NROW(statsmatrix), failure=FALSE,
                          mc.se=mc.se#, #acf=mcmcacf,
                          #fullsample=statsmatrix.all
                          ),
                        class="ergm"))
  }
}

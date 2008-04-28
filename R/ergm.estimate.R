ergm.estimate<-function(theta0, model, statsmatrix, statsmatrix.miss=NULL,
                        epsilon=1e-10, nr.maxit=100, nr.reltol=sqrt(.Machine$double.eps),
                        metric="Likelihood",
                        method="Nelder-Mead", compress=FALSE,
                        calc.mcmc.se=TRUE, hessian=TRUE,
                        verbose=FALSE, trace=6*verbose,
                        trustregion=20, 
                        estimateonly=FALSE, ...) {
  samplesize <- dim(statsmatrix)[1]
  if(compress){
    statsmatrix0 <- ergm.sufftoprob(statsmatrix,compress=TRUE)
    probs <- statsmatrix0[,ncol(statsmatrix0)]
    statsmatrix0 <- statsmatrix0[,-ncol(statsmatrix0),drop=FALSE]
  }else{
    statsmatrix0 <- statsmatrix
    probs <- rep(1/nrow(statsmatrix0),nrow(statsmatrix0))
  }
  av <- apply(sweep(statsmatrix0,1,probs,"*"), 2, sum)
  xsim <- sweep(statsmatrix0, 2, av,"-")
  xobs <- - av
  
#
# Set up the initial estimate
#
  guess <- theta0[!model$etamap$offsettheta]
  if (verbose) cat("Converting theta0 to eta0\n")
  eta0 <- ergm.eta(theta0, model$etamap) #unsure about this
  model$etamap$theta0 <- theta0
#
# Log-Likelihood and gradient functions
#
  if(metric=="Likelihood"){
#   Default method
#   check degeneracy removed from here?
    penalty <- 0.5
  }else{
#
#   Simple convergence without log-normal modification to
#   log-likelihood computation
#
    llik.fun <- llik.fun2
    llik.grad <- llik.fun2
    llik.hessian <- llik.hessian2
    penalty <- 0.5
  }
  if (verbose) cat("Optimizing loglikelihood\n")
  Lout <- try(optim(par=guess, 
                    fn=llik.fun, #  gr=llik.grad,
                    hessian=hessian,
                    method=method,
                    control=list(trace=trace,fnscale=-1,maxit=nr.maxit,reltol=nr.reltol),
                    xobs=xobs,
                    xsim=xsim, probs=probs,
                    penalty=0.5, trustregion=trustregion,
                    eta0=eta0, etamap=model$etamap))
# if(verbose){cat("Log-likelihood ratio is", Lout$value,"\n")}
  if(Lout$value < trustregion-0.001){
   current.scipen <- options()$scipen
   options(scipen=3)
   cat("the log-likelihood improved by",
       format.pval(Lout$value,digits=4,eps=1e-4),"\n")
   options(scipen=current.scipen)
  }else{
   cat("the log-likelihood did not improve.\n")
  }
  if(inherits(Lout,"try-error") || Lout$value > 199 ||
     Lout$value < -790) {
    cat("MLE could not be found. Trying Nelder-Mead...\n")
    Lout <- try(optim(par=guess, 
                      fn=llik.fun,
                      hessian=hessian,
                      method="Nelder-Mead",
                      control=list(trace=trace,fnscale=-1,maxit=100*nr.maxit,
#                       reltol=nr.reltol,
                        reltol=0.01),
                      xobs=xobs, 
                      xsim=xsim, probs=probs, 
                      penalty=0.5, trustregion=trustregion,
                      eta0=eta0, etamap=model$etamap))
    if(inherits(Lout,"try-error") || Lout$value > 500 ){
      cat(paste("No direct MLE exists!\n"))
    }
    if(Lout$convergence != 0 ){
      cat("Non-convergence after", nr.maxit, "iterations.\n")
    }
    cat("Nelder-Mead Log-likelihood ratio is ", Lout$value,"\n")
  }
  theta <- theta0
  theta[!model$etamap$offsettheta] <- Lout$par
  names(theta) <- names(theta0)
  if (estimateonly) {
    # Output results as ergm-class object
    return(structure(list(coef=theta, 
                          MCMCtheta=theta0, 
                          samplesize=samplesize, 
                          failure=FALSE),
                        class="ergm"))
  } else {
    gradient <- llik.grad(theta=Lout$par, xobs=xobs, xsim=xsim,
                          probs=probs, 
                          penalty=0.5, eta0=eta0, etamap=model$etamap)
    gradient[model$etamap$offsettheta] <- 0
    #
    #  Calculate the auto-covariance of the MCMC suff. stats.
    #  and hence the MCMC s.e.
    #
    mc.se <- rep(NA, length=length(theta))
    covar <- NA
    if(!hessian){
      #  covar <- robust.inverse(cov(xsim))
      #  Lout$hessian <- cov(xsim)
      Lout$hessian <- llik.hessian(theta=theta, xobs=xobs, xsim=xsim,
                                   probs=probs, 
                                   penalty=0.5,
                                   eta0=eta0, etamap=model$etamap
                                   )
      Lout$hessian[,model$etamap$offsettheta] <- 0
      Lout$hessian[model$etamap$offsettheta,] <- 0
    }
    if(calc.mcmc.se){
      if (verbose) cat("Starting MCMC s.e. computation.\n")
        mcmcse <- ergm.MCMCse(theta, theta0, statsmatrix0,
                              statsmatrix.miss,
                              model=model)
        mc.se <- mcmcse$mc.se
        #   covar <- robust.inverse(-mcmcse$hessian)
        H <- mcmcse$hessian
        covar <- robust.inverse(-H)
        if(all(!is.na(diag(covar))) && all(diag(covar)<0)){covar <- -covar}
        mc.se[model$etamap$offsettheta] <- NA
    }
    if(inherits(covar,"try-error") | is.na(covar[1])){
      covar <- robust.inverse(-Lout$hessian)
    }
    covar[,model$etamap$offsettheta ] <- NA
    covar[ model$etamap$offsettheta,] <- NA
    c0  <- llik.fun(theta=Lout$par, xobs=xobs,
                    xsim=xsim, probs=probs,
                    penalty=0.5, eta0=eta0, etamap=model$etamap)
    #   VIP: Note added penalty for more skewness in the 
    #        values computed relative to 0
    #
    c01 <- llik.fun(theta=Lout$par-Lout$par, xobs=xobs,
                    xsim=xsim, probs=probs,
                    penalty=0.67, eta0=eta0, etamap=model$etamap)
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
#    if (verbose) cat("Starting MCMC s.e. ACF computation.\n")
#    if(calc.mcmc.se){
#      mcmcacf <- ergm.MCMCacf(statsmatrix0)
#    }else{                 
#      mcmcacf <- covar-covar
#    }
#    if (verbose) cat("Ending MCMC s.e. ACF computation.\n")
      
      # Output results as ergm-class object
    return(structure(list(coef=theta, sample=statsmatrix, 
                          iterations=iteration, #mcmcloglik=mcmcloglik,
                          MCMCtheta=theta0, 
                          loglikelihood=loglikelihood, gradient=gradient,
                          covar=covar, samplesize=samplesize, failure=FALSE,
                          mc.se=mc.se#, #acf=mcmcacf,
                          #fullsample=statsmatrix.all
                          ),
                        class="ergm"))
  }
}

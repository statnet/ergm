ergm.estimate<-function(theta0, model, xobs=NULL, statsmatrix, statsmatrix.miss=NULL,    #####Added xobs=NULL
                        epsilon=1e-10, nr.maxit=1000, nr.reltol=sqrt(.Machine$double.eps),
                        metric="Likelihood",
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
    # Center the statsmatrix0.miss matrix by subtracting the column mean vector:
    av.miss <- apply(sweep(statsmatrix0.miss,1,probs.miss,"*"), 2, sum)
    xsim.miss <- sweep(statsmatrix0.miss, 2, av.miss,"-")
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
  
  # Center the statsmatrix0 matrix by subtracting the column means.  Note that
  # this means that xobs, if not NULL, should already have these means subtracted
  # out!  THIS IS A BAD IDEA!!  It means that the statsmatrix and xobs are not
  # on the same scale.
  # To do:  Change this so that either centering
  # is not done or the centering is applied to the xobs vector as well in a 
  # sensible fashion.  Maybe the latter would be better for reasons of numerical
  # stability? 
  av <- apply(sweep(statsmatrix0,1,probs,"*"), 2, sum)
  xsim <- sweep(statsmatrix0, 2, av,"-")
  
  # fix xobs.  Current algorithm is stupid:  If xobs==NULL, then set xobs=-av
  # and otherwise, do nothing.  Should be:  If xobs==NULL, set it to zero.  Then,
  # subtract av and add av.miss.
  if (is.null(xobs)) {
    xobs <-  -av + ifelse(missingflag, av.miss, 0)
  }
  
  # Convert theta0 to eta0
  eta0 <- ergm.eta(theta0, model$etamap)
  
  # "guess" will be the starting point for the optim search algorithm.
  # But only the non-offset values are relevant; the others will be
  # passed to the likelihood functions by way of the etamap$theta0 element
  # of the model object.  NB:  This is a really ugly way to do this!  Change it?
  guess <- theta0[!model$etamap$offsettheta]
  model$etamap$theta0 <- theta0

  #
# Log-Likelihood and gradient functions
#

  # Choose appropriate loglikelihood, gradient, and Hessian functions
  # Also, choose varweight multiplier for covariance term in loglikelihood
  # where 0.5 is the "true" value but this can be increased or decreased
  if(metric=="Likelihood" && !missingflag) {
    loglikelihoodfn <- llik.fun
    gradientfn <- llik.grad
    Hessianfn <- llik.hessian
    varweight <- 0.5
  } else if (metric == "Likelihood" && missingflag) {
    loglikelihoodfn <- llik.fun.miss
    gradientfn <- llik.grad.miss
    Hessianfn <- llik.hessian.miss
    varweight <- 0.5
  } else {
    loglikelihoodfn <- llik.fun2
    gradientfn <- llik.grad2
    Hessianfn <- llik.hessian2
    varweight <- 0.5
  }

  if (verbose) cat("Optimizing loglikelihood\n")
  Lout <- try(optim(par=guess,                  
                    fn=loglikelihoodfn, #  gr=gradientfn,
                    hessian=hessianflag,
                    method=method,
                    control=list(trace=trace, fnscale=-1,
                                 maxit=nr.maxit,reltol=nr.reltol),
                    xobs=xobs,
                    xsim=xsim, probs=probs,
                    xsim.miss = ifelse(missingflag, xsim.miss, NULL),
                    probs.miss = ifelse(missingflag, xsim.miss, NULL),
                    varweight=varweight, trustregion=trustregion,
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
                      hessian=hessianflag,
                      method="Nelder-Mead",
                      control=list(trace=trace,fnscale=-1,maxit=100*nr.maxit,
                        reltol=0.01),
                      xobs=xobs, 
                      xsim=xsim, probs=probs, 
                      xsim.miss = ifelse(missingflag, xsim.miss, NULL),
                      probs.miss = ifelse(missingflag, xsim.miss, NULL),
                      varweight = varweight, trustregion=trustregion,
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
                          samplesize=NROW(statsmatrix),
                          failure=FALSE),
                        class="ergm"))
  } else {
    gradient <- llik.grad(theta=Lout$par, xobs=xobs, xsim=xsim,
                          probs=probs, 
                          xsim.miss = ifelse(missingflag, xsim.miss, NULL),
                          probs.miss = ifelse(missingflag, xsim.miss, NULL),
                          varweight=varweight, eta0=eta0, etamap=model$etamap)
    gradient[model$etamap$offsettheta] <- 0
    #
    #  Calculate the auto-covariance of the MCMC suff. stats.
    #  and hence the MCMC s.e.
    #
    mc.se <- rep(NA, length=length(theta))
    covar <- NA
    if(!hessianflag){
      #  covar <- robust.inverse(cov(xsim))
      #  Lout$hessian <- cov(xsim)
      Lout$hessian <- llik.hessian(theta=theta, xobs=xobs, xsim=xsim,
                                   probs=probs, 
                                   xsim.miss = ifelse(missingflag, xsim.miss, NULL),
                                   probs.miss = ifelse(missingflag, xsim.miss, NULL),
                                   varweight=varweight,
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
                    xsim.miss = ifelse(missingflag, xsim.miss, NULL),
                    probs.miss = ifelse(missingflag, xsim.miss, NULL),
                    varweight=0.5, eta0=eta0, etamap=model$etamap)
    #   VIP: Note added varweight for more skewness in the 
    #        values computed relative to 0
    #
    c01 <- llik.fun(theta=Lout$par-Lout$par, xobs=xobs,
                    xsim=xsim, probs=probs,
                    xsim.miss = ifelse(missingflag, xsim.miss, NULL),
                    probs.miss = ifelse(missingflag, xsim.miss, NULL),
                    varweight=0.67, eta0=eta0, etamap=model$etamap)
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



# pre-processing:  
#  (0) find sample size of statsmatrix;
#  (1) set up centered version of the statsmatrix; (choices:  compress?  missing?)
#  (2) obtain first guess at optimal
#  (3) convert theta0 to eta0
#  (4) modify value of xobs like statsmatrix has been modified
#  (5) change model$etamap$theta0 to theta0 (??)
#  (6) Determine appropriate set of likelihood and gradient functions
#  (7) call optim function
#  (8) check for errors in optim call (retry if errors found)
#  (9) report change in loglikelihood; report on non-convergence
#  (10) modify value of theta to optimal; add names
#  (11) return if estimateonly==TRUE
#  (12) calculate gradient vector
#  (13) If optim did not already return hessian matrix approximation, calculate hessian
#  (14) For both gradient and Hessian, zero out entries corresponding to model$etamap$offsettheta
#  (15) if calc.mcmc.se==TRUE, calculate MCMC S.E.
#  (16) set covar equal to the inverse of the Hessian
#  (17) for entries corresponding to model$etamap$offsettheta, set covar entries to NA
#  (18) calculate mcmcloglik approximation (currently using very ad-hoc method)



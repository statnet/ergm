ergm.estimate.only<-function(theta0, model, 
                        statsmatrix, statsmatrix.miss,
                        epsilon=1e-10, nr.maxit=100,
                        metric="Likelihood",
                        method="Nelder-Mead", compress=FALSE,
                        calc.mcmc.se=TRUE, hessian=TRUE,
                        verbose=FALSE, trace=6*verbose,
                        trustregion=20, ...) {
  samplesize <- dim(statsmatrix)[1]
  if(compress){
    statsmatrix0 <- ergm.sufftoprob(statsmatrix,compress=TRUE)
    probs <- statsmatrix0[,ncol(statsmatrix0)]
    statsmatrix0 <- statsmatrix0[,-ncol(statsmatrix0)]
    if(!is.null(statsmatrix.miss)){
      statsmatrix0.miss <- ergm.sufftoprob(statsmatrix.miss,compress=TRUE)
      probs.miss <- statsmatrix0.miss[,ncol(statsmatrix0.miss)]
      statsmatrix0.miss <- statsmatrix0.miss[,-ncol(statsmatrix0.miss)]
    }else{
      statsmatrix0.miss <- NULL
      probs.miss <- NULL
    }
  }else{
    statsmatrix0 <- statsmatrix
    probs <- rep(1/nrow(statsmatrix0),nrow(statsmatrix0))
    if(!is.null(statsmatrix.miss)){
      statsmatrix0.miss <- statsmatrix.miss
      probs.miss <- rep(1/nrow(statsmatrix0.miss),nrow(statsmatrix0.miss))
    }else{
      statsmatrix0.miss <- NULL
      probs.miss <- NULL
    }
  }
  av <- apply(sweep(statsmatrix0,1,probs,"*"), 2, sum)
  xsim <- sweep(statsmatrix0, 2, av,"-")
  if(!is.null(statsmatrix.miss)){
    av.miss <- apply(sweep(statsmatrix0.miss,1,probs.miss,"*"), 2, sum)
    xsim.miss <- sweep(statsmatrix0.miss, 2, av.miss,"-")
    xobs <- av.miss-av
  }else{
    xsim.miss <- NULL
    probs.miss <- NULL
    xobs <- - av
  }
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
   if(is.null(statsmatrix.miss)){
#   Default method
#   check degeneracy removed from here?
    penalty <- 0.5
   }else{
    llik.fun <- llik.fun.miss
    llik.grad <- llik.fun.miss
    llik.hessian <- llik.hessian.miss
    penalty <- 0.5
   }
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
                    control=list(trace=trace,fnscale=-1,maxit=nr.maxit),
                    xobs=xobs,
                    xsim=xsim, probs=probs,
                    xsim.miss=xsim.miss, probs.miss=probs.miss,
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
                        reltol=0.01),
                      xobs=xobs, 
                      xsim=xsim, probs=probs, 
                      xsim.miss=xsim.miss, probs.miss=probs.miss,
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

# Output results as ergm-class object
  structure(list(coef=theta, 
                 MCMCtheta=theta0, 
                 samplesize=samplesize, failure=FALSE),
            class="ergm") 
}

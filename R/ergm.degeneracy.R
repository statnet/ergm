ergm.degeneracy <- function(object, 
                          control=ergm.control(),
                          verbose=FALSE) {
  
  if(!is.ergm(object)){
    stop("A ergm object argument must be given.")
  }
  if(!is.null(object$mplefit$glm) & is.matrix(object$sample)){
   # So a MCMC fit
   if(object$loglikelihood>control$trustregion-0.0001){
    object$degeneracy <- control$trustregion
   }else{
    changeobs <- (-2*object$mplefit$glm$y+1)*model.matrix(object$mplefit$glm)
    object$degeneracy.type <- apply(changeobs, 1, ergm.compute.degeneracy,
         object$MCMCtheta, object$etamap, object$sample)
    wgts <- object$mplefit$glm$prior.weights
    names(wgts) <- "num.dyads"
    object$degeneracy.type <- t(rbind(object$degeneracy.type,wgts))
    object$degeneracy <- max(object$degeneracy.type[,1],na.rm=TRUE)
   }
  }else{
   # So a non-MCMC fit
   if(!is.null(object$glm)){
   # So the MPLE was fit
    # This is the change in log-likelihood for logistic regression
    object$degeneracy.type <- abs(model.matrix(object$glm) %*% object$glm$coef)
    wgts <- object$glm$prior.weights
    object$degeneracy.type <- cbind(object$degeneracy.type,wgts)
    colnames(object$degeneracy.type) <- c("delta.log.lik","num.dyads")
#   dr <- residuals(object$glm,type="deviance")^2 / object$glm$prior.weights
#   object$degeneracy.type <- rep(dr/2, fit$glm$prior.weights)
    object$degeneracy <- max(object$degeneracy.type[,1],na.rm=TRUE)
   }else{
    object$degeneracy <- control$trustregion
    object$degeneracy.type <- NULL
   }
  }
  object
}

ergm.compute.degeneracy<-function(xobs, theta0, etamap, statsmatrix,
                        epsilon=1e-10, nr.maxit=100,
                        verbose=FALSE, trace=6*verbose,
                        hessian=FALSE,
                        trustregion=20, ...) {
  samplesize <- dim(statsmatrix)[1]
  statsmatrix0 <- statsmatrix
  probs <- rep(1/nrow(statsmatrix0),nrow(statsmatrix0))
  statsmatrix0.miss <- NULL
  probs.miss <- NULL
  av <- apply(sweep(statsmatrix0,1,probs,"*"), 2, sum)
  xsim <- sweep(statsmatrix0, 2, av,"-")
  xsim.miss <- NULL
  probs.miss <- NULL
# xobs0 <- summary(model$formula)
  xobs <- -xobs - av
#
# Set up the initial estimate
#
  guess <- theta0
  if (verbose) cat("Converting theta0 to eta0\n")
  eta0 <- ergm.eta(theta0, etamap) #unsure about this
  etamap$theta0 <- theta0
#
# Log-Likelihood and gradient functions
#
  penalty <- 0.5
  if (verbose) cat("Optimizing loglikelihood\n")
  Lout <- try(optim(par=guess, 
                    fn=llik.fun, #  gr=llik.grad,
                    hessian=FALSE,
                    method="BFGS",
                    control=list(trace=trace,fnscale=-1,maxit=nr.maxit),
                    xobs=xobs,
                    xsim=xsim, probs=probs,
                    xsim.miss=xsim.miss, probs.miss=probs.miss,
                    penalty=0.5, trustregion=trustregion,
                    eta0=eta0, etamap=etamap))
  if(verbose){cat("Log-likelihood ratio is", Lout$value,"\n")}
  if(inherits(Lout,"try-error") || Lout$value > 199 ||
     Lout$value < -790) {
    cat("MLE could not be found. Degenerate!\n")
    cat("Nelder-Mead Log-likelihood ratio is ", Lout$value,"\n")
#   return(list(coef=theta0, 
#      loglikelihood=Lout$value))
  }
  theta <- Lout$par
  names(theta) <- names(theta0)
# c0  <- llik.fun(theta=Lout$par, xobs=xobs,
#                 xsim=xsim, probs=probs,
#                 xsim.miss=xsim.miss, probs.miss=probs.miss,
#                 penalty=0.5, eta0=eta0, etamap=etamap)
  loglikelihood <- Lout$value
  names(loglikelihood) <- "loglikelihood"

# loglikelihood
#
# Returns the change in the log-likelihood ratio if the dyad is toggled
# and the change in the coefficient
# list(coef=theta, 
#      loglikelihood=loglikelihood)
  c(loglikelihood, theta)
}
